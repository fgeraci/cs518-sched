/*
 *  linux/drivers/block/ll_rw_blk.c
 *
 * Copyright (C) 1991, 1992 Linus Torvalds
 * Copyright (C) 1994,      Karl Keyte: Added support for disk statistics
 * Elevator latency, (C) 2000  Andrea Arcangeli <andrea@suse.de> SuSE
 * Queue request tables / lock, selectable elevator, Jens Axboe <axboe@suse.de>
 * kernel-doc documentation started by NeilBrown <neilb@cse.unsw.edu.au> -  July2000
 * bio rewrite, highmem i/o, etc, Jens Axboe <axboe@suse.de> - may 2001
 */

/*
 * This handles all read/write requests to block devices
 */
#include <linux/sched.h>
#include <linux/kernel.h>
#include <linux/kernel_stat.h>
#include <linux/errno.h>
#include <linux/string.h>
#include <linux/config.h>
#include <linux/locks.h>
#include <linux/mm.h>
#include <linux/swap.h>
#include <linux/init.h>
#include <linux/smp_lock.h>
#include <linux/bootmem.h>
#include <linux/completion.h>
#include <linux/compiler.h>
#include <scsi/scsi.h>

#include <asm/system.h>
#include <asm/io.h>
#include <linux/blk.h>
#include <linux/highmem.h>
#include <linux/slab.h>
#include <linux/module.h>

/*
 * MAC Floppy IWM hooks
 */

#ifdef CONFIG_MAC_FLOPPY_IWM
extern int mac_floppy_init(void);
#endif

/*
 * For the allocated request tables
 */
static kmem_cache_t *request_cachep;

/*
 * The "disk" task queue is used to start the actual requests
 * after a plug
 */
DECLARE_TASK_QUEUE(tq_disk);

/* This specifies how many sectors to read ahead on the disk. */

int read_ahead[MAX_BLKDEV];

/* blk_dev_struct is:
 *	request_queue
 *	*queue
 */
struct blk_dev_struct blk_dev[MAX_BLKDEV]; /* initialized by blk_dev_init() */

/*
 * blk_size contains the size of all block-devices in units of 1024 byte
 * sectors:
 *
 * blk_size[MAJOR][MINOR]
 *
 * if (!blk_size[MAJOR]) then no minor size checking is done.
 */
int * blk_size[MAX_BLKDEV];

/*
 * blksize_size contains the size of all block-devices:
 *
 * blksize_size[MAJOR][MINOR]
 *
 * if (!blksize_size[MAJOR]) then 1024 bytes is assumed.
 */
int * blksize_size[MAX_BLKDEV];

/*
 * The following tunes the read-ahead algorithm in mm/filemap.c
 */
int * max_readahead[MAX_BLKDEV];

/*
 * How many reqeusts do we allocate per queue,
 * and how many do we "batch" on freeing them?
 */
int queue_nr_requests, batch_requests;
unsigned long blk_max_low_pfn, blk_max_pfn;
int blk_nohighio = 0;

/**
 * blk_get_queue: - return the queue that matches the given device
 * @dev:    device
 *
 * Description:
 *     Given a specific device, return the queue that will hold I/O
 *     for it. This is either a &struct blk_dev_struct lookup and a
 *     call to the ->queue() function defined, or the default queue
 *     stored in the same location.
 *
 **/
inline request_queue_t *blk_get_queue(kdev_t dev)
{
	struct blk_dev_struct *bdev = blk_dev + major(dev);

	if (bdev->queue)
		return bdev->queue(dev);
	else
		return &blk_dev[major(dev)].request_queue;
}

void blk_queue_prep_rq(request_queue_t *q, prep_rq_fn *pfn)
{
	q->prep_rq_fn = pfn;
}

/**
 * blk_queue_make_request - define an alternate make_request function for a device
 * @q:  the request queue for the device to be affected
 * @mfn: the alternate make_request function
 *
 * Description:
 *    The normal way for &struct bios to be passed to a device
 *    driver is for them to be collected into requests on a request
 *    queue, and then to allow the device driver to select requests
 *    off that queue when it is ready.  This works well for many block
 *    devices. However some block devices (typically virtual devices
 *    such as md or lvm) do not benefit from the processing on the
 *    request queue, and are served best by having the requests passed
 *    directly to them.  This can be achieved by providing a function
 *    to blk_queue_make_request().
 *
 * Caveat:
 *    The driver that does this *must* be able to deal appropriately
 *    with buffers in "highmemory". This can be accomplished by either calling
 *    bio_kmap() to get a temporary kernel mapping, or by calling
 *    blk_queue_bounce() to create a buffer in normal memory.
 **/
void blk_queue_make_request(request_queue_t * q, make_request_fn * mfn)
{
	/*
	 * set defaults
	 */
	q->max_phys_segments = MAX_PHYS_SEGMENTS;
	q->max_hw_segments = MAX_HW_SEGMENTS;
	q->make_request_fn = mfn;
	blk_queue_max_sectors(q, MAX_SECTORS);
	blk_queue_hardsect_size(q, 512);

	/*
	 * by default assume old behaviour and bounce for any highmem page
	 */
	blk_queue_bounce_limit(q, BLK_BOUNCE_HIGH);

	init_waitqueue_head(&q->queue_wait);
}

/**
 * blk_queue_bounce_limit - set bounce buffer limit for queue
 * @q:  the request queue for the device
 * @dma_addr:   bus address limit
 *
 * Description:
 *    Different hardware can have different requirements as to what pages
 *    it can do I/O directly to. A low level driver can call
 *    blk_queue_bounce_limit to have lower memory pages allocated as bounce
 *    buffers for doing I/O to pages residing above @page. By default
 *    the block layer sets this to the highest numbered "low" memory page.
 **/
void blk_queue_bounce_limit(request_queue_t *q, u64 dma_addr)
{
	unsigned long bounce_pfn = dma_addr >> PAGE_SHIFT;
	unsigned long mb = dma_addr >> 20;
	static request_queue_t *last_q;

	/*
	 * set appropriate bounce gfp mask -- unfortunately we don't have a
	 * full 4GB zone, so we have to resort to low memory for any bounces.
	 * ISA has its own < 16MB zone.
	 */
	if (dma_addr == BLK_BOUNCE_ISA) {
		init_emergency_isa_pool();
		q->bounce_gfp = GFP_NOIO | GFP_DMA;
	} else
		q->bounce_gfp = GFP_NOHIGHIO;

	/*
	 * keep this for debugging for now...
	 */
	if (dma_addr != BLK_BOUNCE_HIGH && q != last_q) {
		printk("blk: queue %p, ", q);
		if (dma_addr == BLK_BOUNCE_ANY)
			printk("no I/O memory limit\n");
		else
			printk("I/O limit %luMb (mask 0x%Lx)\n", mb, (long long) dma_addr);
	}

	q->bounce_pfn = bounce_pfn;
	last_q = q;
}


/**
 * blk_queue_max_sectors - set max sectors for a request for this queue
 * @q:  the request queue for the device
 * @max_sectors:  max sectors in the usual 512b unit
 *
 * Description:
 *    Enables a low level driver to set an upper limit on the size of
 *    received requests.
 **/
void blk_queue_max_sectors(request_queue_t *q, unsigned short max_sectors)
{
	q->max_sectors = max_sectors;
}

/**
 * blk_queue_max_phys_segments - set max phys segments for a request for this queue
 * @q:  the request queue for the device
 * @max_segments:  max number of segments
 *
 * Description:
 *    Enables a low level driver to set an upper limit on the number of
 *    physical data segments in a request.  This would be the largest sized
 *    scatter list the driver could handle.
 **/
void blk_queue_max_phys_segments(request_queue_t *q, unsigned short max_segments)
{
	q->max_phys_segments = max_segments;
}

/**
 * blk_queue_max_hw_segments - set max hw segments for a request for this queue
 * @q:  the request queue for the device
 * @max_segments:  max number of segments
 *
 * Description:
 *    Enables a low level driver to set an upper limit on the number of
 *    hw data segments in a request.  This would be the largest number of
 *    address/length pairs the host adapter can actually give as once
 *    to the device.
 **/
void blk_queue_max_hw_segments(request_queue_t *q, unsigned short max_segments)
{
	q->max_hw_segments = max_segments;
}

/**
 * blk_queue_max_segment_size - set max segment size for blk_rq_map_sg
 * @q:  the request queue for the device
 * @max_size:  max size of segment in bytes
 *
 * Description:
 *    Enables a low level driver to set an upper limit on the size of a
 *    coalesced segment
 **/
void blk_queue_max_segment_size(request_queue_t *q, unsigned int max_size)
{
	q->max_segment_size = max_size;
}

/**
 * blk_queue_hardsect_size - set hardware sector size for the queue
 * @q:  the request queue for the device
 * @size:  the hardware sector size, in bytes
 *
 * Description:
 *   This should typically be set to the lowest possible sector size
 *   that the hardware can operate on (possible without reverting to
 *   even internal read-modify-write operations). Usually the default
 *   of 512 covers most hardware.
 **/
void blk_queue_hardsect_size(request_queue_t *q, unsigned short size)
{
	q->hardsect_size = size;
}

/**
 * blk_queue_segment_boundary - set boundary rules for segment merging
 * @q:  the request queue for the device
 * @mask:  the memory boundary mask
 **/
void blk_queue_segment_boundary(request_queue_t *q, unsigned long mask)
{
	q->seg_boundary_mask = mask;
}

void blk_queue_assign_lock(request_queue_t *q, spinlock_t *lock)
{
	spin_lock_init(lock);
	q->queue_lock = lock;
}

static char *rq_flags[] = { "REQ_RW", "REQ_RW_AHEAD", "REQ_BARRIER",
			   "REQ_CMD", "REQ_NOMERGE", "REQ_STARTED",
			   "REQ_DONTPREP", "REQ_DRIVE_CMD", "REQ_DRIVE_TASK",
			   "REQ_PC", "REQ_BLOCK_PC", "REQ_SENSE",
			   "REQ_SPECIAL" };

void blk_dump_rq_flags(struct request *rq, char *msg)
{
	int bit;

	printk("%s: dev %02x:%02x: ", msg, major(rq->rq_dev), minor(rq->rq_dev));
	bit = 0;
	do {
		if (rq->flags & (1 << bit))
			printk("%s ", rq_flags[bit]);
		bit++;
	} while (bit < __REQ_NR_BITS);

	if (rq->flags & REQ_CMD)
		printk("sector %lu, nr/cnr %lu/%u\n", rq->sector,
						       rq->nr_sectors,
						       rq->current_nr_sectors);

	printk("\n");
}

/*
 * standard prep_rq_fn that builds 10 byte cmds
 */
int ll_10byte_cmd_build(request_queue_t *q, struct request *rq)
{
	int hard_sect = get_hardsect_size(rq->rq_dev);
	sector_t block = rq->hard_sector / (hard_sect >> 9);
	unsigned long blocks = rq->hard_nr_sectors / (hard_sect >> 9);

	if (!(rq->flags & REQ_CMD))
		return 0;

	memset(rq->cmd, 0, sizeof(rq->cmd));

	if (rq_data_dir(rq) == READ)
		rq->cmd[0] = READ_10;
	else 
		rq->cmd[0] = WRITE_10;

	/*
	 * fill in lba
	 */
	rq->cmd[2] = (block >> 24) & 0xff;
	rq->cmd[3] = (block >> 16) & 0xff;
	rq->cmd[4] = (block >>  8) & 0xff;
	rq->cmd[5] = block & 0xff;

	/*
	 * and transfer length
	 */
	rq->cmd[7] = (blocks >> 8) & 0xff;
	rq->cmd[8] = blocks & 0xff;

	return 0;
}

void blk_recount_segments(request_queue_t *q, struct bio *bio)
{
	struct bio_vec *bv, *bvprv = NULL;
	int i, nr_phys_segs, nr_hw_segs, seg_size, cluster;

	if (unlikely(!bio->bi_io_vec))
		return;

	cluster = q->queue_flags & (1 << QUEUE_FLAG_CLUSTER);
	seg_size = nr_phys_segs = nr_hw_segs = 0;
	bio_for_each_segment(bv, bio, i) {
		if (bvprv && cluster) {
			int phys, seg;

			if (seg_size + bv->bv_len > q->max_segment_size) {
				nr_phys_segs++;
				goto new_segment;
			}

			phys = BIOVEC_PHYS_MERGEABLE(bvprv, bv);
			seg = BIOVEC_SEG_BOUNDARY(q, bvprv, bv);
			if (!phys || !seg)
				nr_phys_segs++;
			if (!seg)
				goto new_segment;

			if (!BIOVEC_VIRT_MERGEABLE(bvprv, bv))
				goto new_segment;

			seg_size += bv->bv_len;
			bvprv = bv;
			continue;
		} else {
			nr_phys_segs++;
		}
new_segment:
		nr_hw_segs++;
		bvprv = bv;
		seg_size = bv->bv_len;
	}

	bio->bi_phys_segments = nr_phys_segs;
	bio->bi_hw_segments = nr_hw_segs;
	bio->bi_flags |= (1 << BIO_SEG_VALID);
}


inline int blk_phys_contig_segment(request_queue_t *q, struct bio *bio,
				   struct bio *nxt)
{
	if (!(q->queue_flags & (1 << QUEUE_FLAG_CLUSTER)))
		return 0;

	if (!BIOVEC_PHYS_MERGEABLE(__BVEC_END(bio), __BVEC_START(nxt)))
		return 0;
	if (bio->bi_size + nxt->bi_size > q->max_segment_size)
		return 0;

	/*
	 * bio and nxt are contigous in memory, check if the queue allows
	 * these two to be merged into one
	 */
	if (BIO_SEG_BOUNDARY(q, bio, nxt))
		return 1;

	return 0;
}

inline int blk_hw_contig_segment(request_queue_t *q, struct bio *bio,
				 struct bio *nxt)
{
	if (!(q->queue_flags & (1 << QUEUE_FLAG_CLUSTER)))
		return 0;

	if (!BIOVEC_VIRT_MERGEABLE(__BVEC_END(bio), __BVEC_START(nxt)))
		return 0;
	if (bio->bi_size + nxt->bi_size > q->max_segment_size)
		return 0;

	/*
	 * bio and nxt are contigous in memory, check if the queue allows
	 * these two to be merged into one
	 */
	if (BIO_SEG_BOUNDARY(q, bio, nxt))
		return 1;

	return 0;
}

/*
 * map a request to scatterlist, return number of sg entries setup. Caller
 * must make sure sg can hold rq->nr_phys_segments entries
 */
int blk_rq_map_sg(request_queue_t *q, struct request *rq, struct scatterlist *sg)
{
	struct bio_vec *bvec, *bvprv;
	struct bio *bio;
	int nsegs, i, cluster;

	nsegs = 0;
	cluster = q->queue_flags & (1 << QUEUE_FLAG_CLUSTER);

	/*
	 * for each bio in rq
	 */
	bvprv = NULL;
	rq_for_each_bio(bio, rq) {
		/*
		 * for each segment in bio
		 */
		bio_for_each_segment(bvec, bio, i) {
			int nbytes = bvec->bv_len;

			if (bvprv && cluster) {
				if (sg[nsegs - 1].length + nbytes > q->max_segment_size)
					goto new_segment;

				if (!BIOVEC_PHYS_MERGEABLE(bvprv, bvec))
					goto new_segment;
				if (!BIOVEC_SEG_BOUNDARY(q, bvprv, bvec))
					goto new_segment;

				sg[nsegs - 1].length += nbytes;
			} else {
new_segment:
				memset(&sg[nsegs],0,sizeof(struct scatterlist));
				sg[nsegs].page = bvec->bv_page;
				sg[nsegs].length = nbytes;
				sg[nsegs].offset = bvec->bv_offset;

				nsegs++;
			}
			bvprv = bvec;
		} /* segments in bio */
	} /* bios in rq */

	return nsegs;
}

/*
 * the standard queue merge functions, can be overridden with device
 * specific ones if so desired
 */

static inline int ll_new_mergeable(request_queue_t *q,
				   struct request *req,
				   struct bio *bio)
{
	int nr_phys_segs = bio_phys_segments(q, bio);

	if (req->nr_phys_segments + nr_phys_segs > q->max_phys_segments) {
		req->flags |= REQ_NOMERGE;
		q->last_merge = NULL;
		return 0;
	}

	/*
	 * A hw segment is just getting larger, bump just the phys
	 * counter.
	 */
	req->nr_phys_segments += nr_phys_segs;
	return 1;
}

static inline int ll_new_hw_segment(request_queue_t *q,
				    struct request *req,
				    struct bio *bio)
{
	int nr_hw_segs = bio_hw_segments(q, bio);
	int nr_phys_segs = bio_phys_segments(q, bio);

	if (req->nr_hw_segments + nr_hw_segs > q->max_hw_segments
	    || req->nr_phys_segments + nr_phys_segs > q->max_phys_segments) {
		req->flags |= REQ_NOMERGE;
		q->last_merge = NULL;
		return 0;
	}

	/*
	 * This will form the start of a new hw segment.  Bump both
	 * counters.
	 */
	req->nr_hw_segments += nr_hw_segs;
	req->nr_phys_segments += nr_phys_segs;
	return 1;
}

static int ll_back_merge_fn(request_queue_t *q, struct request *req, 
			    struct bio *bio)
{
	if (req->nr_sectors + bio_sectors(bio) > q->max_sectors) {
		req->flags |= REQ_NOMERGE;
		q->last_merge = NULL;
		return 0;
	}

	if (BIOVEC_VIRT_MERGEABLE(__BVEC_END(req->biotail), __BVEC_START(bio)))
		return ll_new_mergeable(q, req, bio);

	return ll_new_hw_segment(q, req, bio);
}

static int ll_front_merge_fn(request_queue_t *q, struct request *req, 
			     struct bio *bio)
{
	if (req->nr_sectors + bio_sectors(bio) > q->max_sectors) {
		req->flags |= REQ_NOMERGE;
		q->last_merge = NULL;
		return 0;
	}

	if (BIOVEC_VIRT_MERGEABLE(__BVEC_END(bio), __BVEC_START(req->bio)))
		return ll_new_mergeable(q, req, bio);

	return ll_new_hw_segment(q, req, bio);
}

static int ll_merge_requests_fn(request_queue_t *q, struct request *req,
				struct request *next)
{
	int total_phys_segments = req->nr_phys_segments +next->nr_phys_segments;
	int total_hw_segments = req->nr_hw_segments + next->nr_hw_segments;

	/*
	 * First check if the either of the requests are re-queued
	 * requests.  Can't merge them if they are.
	 */
	if (req->special || next->special)
		return 0;

	/*
	 * Will it become to large?
	 */
	if ((req->nr_sectors + next->nr_sectors) > q->max_sectors)
		return 0;

	total_phys_segments = req->nr_phys_segments + next->nr_phys_segments;
	if (blk_phys_contig_segment(q, req->biotail, next->bio))
		total_phys_segments--;

	if (total_phys_segments > q->max_phys_segments)
		return 0;

	total_hw_segments = req->nr_hw_segments + next->nr_hw_segments;
	if (blk_hw_contig_segment(q, req->biotail, next->bio))
		total_hw_segments--;
    
	if (total_hw_segments > q->max_hw_segments)
		return 0;

	/* Merge is OK... */
	req->nr_phys_segments = total_phys_segments;
	req->nr_hw_segments = total_hw_segments;
	return 1;
}

/*
 * "plug" the device if there are no outstanding requests: this will
 * force the transfer to start only after we have put all the requests
 * on the list.
 *
 * This is called with interrupts off and no requests on the queue.
 * (and with the request spinlock acquired)
 */
void blk_plug_device(request_queue_t *q)
{
	/*
	 * common case
	 */
	if (!elv_queue_empty(q))
		return;

	if (!test_and_set_bit(QUEUE_FLAG_PLUGGED, &q->queue_flags))
		queue_task(&q->plug_tq, &tq_disk);
}

/*
 * remove the plug and let it rip..
 */
static inline void __generic_unplug_device(request_queue_t *q)
{
	/*
	 * not plugged
	 */
	if (!test_and_clear_bit(QUEUE_FLAG_PLUGGED, &q->queue_flags))
		return;

	/*
	 * was plugged, fire request_fn if queue has stuff to do
	 */
	if (!elv_queue_empty(q))
		q->request_fn(q);
}

/**
 * generic_unplug_device - fire a request queue
 * @q:    The &request_queue_t in question
 *
 * Description:
 *   Linux uses plugging to build bigger requests queues before letting
 *   the device have at them. If a queue is plugged, the I/O scheduler
 *   is still adding and merging requests on the queue. Once the queue
 *   gets unplugged (either by manually calling this function, or by
 *   running the tq_disk task queue), the request_fn defined for the
 *   queue is invoked and transfers started.
 **/
void generic_unplug_device(void *data)
{
	request_queue_t *q = (request_queue_t *) data;
	unsigned long flags;

	spin_lock_irqsave(q->queue_lock, flags);
	__generic_unplug_device(q);
	spin_unlock_irqrestore(q->queue_lock, flags);
}

static int __blk_cleanup_queue(struct request_list *list)
{
	struct list_head *head = &list->free;
	struct request *rq;
	int i = 0;

	while (!list_empty(head)) {
		rq = list_entry(head->next, struct request, queuelist);
		list_del(&rq->queuelist);
		kmem_cache_free(request_cachep, rq);
		i++;
	}

	if (i != list->count)
		printk("request list leak!\n");

	list->count = 0;
	return i;
}

/**
 * blk_cleanup_queue: - release a &request_queue_t when it is no longer needed
 * @q:    the request queue to be released
 *
 * Description:
 *     blk_cleanup_queue is the pair to blk_init_queue().  It should
 *     be called when a request queue is being released; typically
 *     when a block device is being de-registered.  Currently, its
 *     primary task it to free all the &struct request structures that
 *     were allocated to the queue.
 * Caveat: 
 *     Hopefully the low level driver will have finished any
 *     outstanding requests first...
 **/
void blk_cleanup_queue(request_queue_t * q)
{
	int count = queue_nr_requests;

	count -= __blk_cleanup_queue(&q->rq[READ]);
	count -= __blk_cleanup_queue(&q->rq[WRITE]);

	if (count)
		printk("blk_cleanup_queue: leaked requests (%d)\n", count);

	elevator_exit(q, &q->elevator);

	memset(q, 0, sizeof(*q));
}

static int blk_init_free_list(request_queue_t *q)
{
	struct request_list *rl;
	struct request *rq;
	int i;

	INIT_LIST_HEAD(&q->rq[READ].free);
	INIT_LIST_HEAD(&q->rq[WRITE].free);
	q->rq[READ].count = 0;
	q->rq[WRITE].count = 0;

	/*
	 * Divide requests in half between read and write
	 */
	rl = &q->rq[READ];
	for (i = 0; i < queue_nr_requests; i++) {
		rq = kmem_cache_alloc(request_cachep, SLAB_KERNEL);
		if (!rq)
			goto nomem;

		/*
		 * half way through, switch to WRITE list
		 */
		if (i == queue_nr_requests / 2)
			rl = &q->rq[WRITE];

		memset(rq, 0, sizeof(struct request));
		rq->rq_status = RQ_INACTIVE;
		list_add(&rq->queuelist, &rl->free);
		rl->count++;
	}

	init_waitqueue_head(&q->rq[READ].wait);
	init_waitqueue_head(&q->rq[WRITE].wait);
	return 0;
nomem:
	blk_cleanup_queue(q);
	return 1;
}

static int __make_request(request_queue_t *, struct bio *);

/**
 * blk_init_queue  - prepare a request queue for use with a block device
 * @q:    The &request_queue_t to be initialised
 * @rfn:  The function to be called to process requests that have been
 *        placed on the queue.
 *
 * Description:
 *    If a block device wishes to use the standard request handling procedures,
 *    which sorts requests and coalesces adjacent requests, then it must
 *    call blk_init_queue().  The function @rfn will be called when there
 *    are requests on the queue that need to be processed.  If the device
 *    supports plugging, then @rfn may not be called immediately when requests
 *    are available on the queue, but may be called at some time later instead.
 *    Plugged queues are generally unplugged when a buffer belonging to one
 *    of the requests on the queue is needed, or due to memory pressure.
 *
 *    @rfn is not required, or even expected, to remove all requests off the
 *    queue, but only as many as it can handle at a time.  If it does leave
 *    requests on the queue, it is responsible for arranging that the requests
 *    get dealt with eventually.
 *
 *    The queue spin lock must be held while manipulating the requests on the
 *    request queue.
 *
 * Note:
 *    blk_init_queue() must be paired with a blk_cleanup_queue() call
 *    when the block device is deactivated (such as at module unload).
 **/
int blk_init_queue(request_queue_t *q, request_fn_proc *rfn, spinlock_t *lock)
{
	int ret;

	if (blk_init_free_list(q))
		return -ENOMEM;

	if ((ret = elevator_init(q, &q->elevator, elevator_linus))) {
		blk_cleanup_queue(q);
		return ret;
	}

	q->request_fn		= rfn;
	q->back_merge_fn       	= ll_back_merge_fn;
	q->front_merge_fn      	= ll_front_merge_fn;
	q->merge_requests_fn	= ll_merge_requests_fn;
	q->prep_rq_fn		= NULL;
	q->plug_tq.sync		= 0;
	q->plug_tq.routine	= &generic_unplug_device;
	q->plug_tq.data		= q;
	q->queue_flags		= (1 << QUEUE_FLAG_CLUSTER);
	q->queue_lock		= lock;
	
	blk_queue_segment_boundary(q, 0xffffffff);

	blk_queue_make_request(q, __make_request);
	blk_queue_max_segment_size(q, MAX_SEGMENT_SIZE);

	blk_queue_max_hw_segments(q, MAX_HW_SEGMENTS);
	blk_queue_max_phys_segments(q, MAX_PHYS_SEGMENTS);
	return 0;
}

#define blkdev_free_rq(list) list_entry((list)->next, struct request, queuelist)
/*
 * Get a free request. queue lock must be held and interrupts
 * disabled on the way in.
 */
static inline struct request *get_request(request_queue_t *q, int rw)
{
	struct request *rq = NULL;
	struct request_list *rl = q->rq + rw;

	if (!list_empty(&rl->free)) {
		rq = blkdev_free_rq(&rl->free);
		list_del(&rq->queuelist);
		rl->count--;
		rq->flags = 0;
		rq->rq_status = RQ_ACTIVE;
		rq->special = NULL;
		rq->q = q;
		rq->rl = rl;
	}

	return rq;
}

/*
 * No available requests for this queue, unplug the device.
 */
static struct request *get_request_wait(request_queue_t *q, int rw)
{
	DECLARE_WAITQUEUE(wait, current);
	struct request_list *rl = &q->rq[rw];
	struct request *rq;

	spin_lock_prefetch(q->queue_lock);

	generic_unplug_device(q);
	add_wait_queue(&rl->wait, &wait);
	do {
		set_current_state(TASK_UNINTERRUPTIBLE);
		if (rl->count < batch_requests)
			schedule();
		spin_lock_irq(q->queue_lock);
		rq = get_request(q, rw);
		spin_unlock_irq(q->queue_lock);
	} while (rq == NULL);
	remove_wait_queue(&rl->wait, &wait);
	current->state = TASK_RUNNING;
	return rq;
}

struct request *blk_get_request(request_queue_t *q, int rw, int gfp_mask)
{
	struct request *rq;

	BUG_ON(rw != READ && rw != WRITE);

	spin_lock_irq(q->queue_lock);
	rq = get_request(q, rw);
	spin_unlock_irq(q->queue_lock);

	if (!rq && (gfp_mask & __GFP_WAIT))
		rq = get_request_wait(q, rw);

	if (rq) {
		rq->flags = 0;
		rq->buffer = NULL;
		rq->bio = rq->biotail = NULL;
		rq->waiting = NULL;
	}
	return rq;
}

void blk_put_request(struct request *rq)
{
	blkdev_release_request(rq);
}

/* RO fail safe mechanism */

static long ro_bits[MAX_BLKDEV][8];

int is_read_only(kdev_t dev)
{
	int minor,major;

	major = major(dev);
	minor = minor(dev);
	if (major < 0 || major >= MAX_BLKDEV) return 0;
	return ro_bits[major][minor >> 5] & (1 << (minor & 31));
}

void set_device_ro(kdev_t dev,int flag)
{
	int minor,major;

	major = major(dev);
	minor = minor(dev);
	if (major < 0 || major >= MAX_BLKDEV) return;
	if (flag) ro_bits[major][minor >> 5] |= 1 << (minor & 31);
	else ro_bits[major][minor >> 5] &= ~(1 << (minor & 31));
}

void drive_stat_acct(struct request *rq, int nr_sectors, int new_io)
{
	unsigned int major = major(rq->rq_dev);
	int rw = rq_data_dir(rq);
	unsigned int index;

	index = disk_index(rq->rq_dev);
	if ((index >= DK_MAX_DISK) || (major >= DK_MAX_MAJOR))
		return;

	kstat.dk_drive[major][index] += new_io;
	if (rw == READ) {
		kstat.dk_drive_rio[major][index] += new_io;
		kstat.dk_drive_rblk[major][index] += nr_sectors;
	} else if (rw == WRITE) {
		kstat.dk_drive_wio[major][index] += new_io;
		kstat.dk_drive_wblk[major][index] += nr_sectors;
	} else
		printk(KERN_ERR "drive_stat_acct: cmd not R/W?\n");
}

/*
 * add-request adds a request to the linked list.
 * queue lock is held and interrupts disabled, as we muck with the
 * request queue list.
 */
static inline void add_request(request_queue_t * q, struct request * req,
			       struct list_head *insert_here)
{
	drive_stat_acct(req, req->nr_sectors, 1);

	/*
	 * elevator indicated where it wants this request to be
	 * inserted at elevator_merge time
	 */
	__elv_add_request(q, req, insert_here);
}

/*
 * Must be called with queue lock held and interrupts disabled
 */
void blkdev_release_request(struct request *req)
{
	struct request_list *rl = req->rl;
	request_queue_t *q = req->q;

	req->rq_status = RQ_INACTIVE;
	req->q = NULL;
	req->rl = NULL;

	if (q) {
		if (q->last_merge == &req->queuelist)
			q->last_merge = NULL;
	}

	/*
	 * Request may not have originated from ll_rw_blk. if not,
	 * it didn't come out of our reserved rq pools
	 */
	if (rl) {
		list_add(&req->queuelist, &rl->free);

		if (++rl->count >= batch_requests &&waitqueue_active(&rl->wait))
			wake_up(&rl->wait);
	}
}

/*
 * Has to be called with the request spinlock acquired
 */
static void attempt_merge(request_queue_t *q, struct request *req,
			  struct request *next)
{
	if (!rq_mergeable(req) || !rq_mergeable(next))
		return;

	/*
	 * not contigious
	 */
	if (req->sector + req->nr_sectors != next->sector)
		return;

	if (rq_data_dir(req) != rq_data_dir(next)
	    || !kdev_same(req->rq_dev, next->rq_dev)
	    || req->nr_sectors + next->nr_sectors > q->max_sectors
	    || next->waiting || next->special)
		return;

	/*
	 * If we are allowed to merge, then append bio list
	 * from next to rq and release next. merge_requests_fn
	 * will have updated segment counts, update sector
	 * counts here.
	 */
	if (q->merge_requests_fn(q, req, next)) {
		elv_merge_requests(q, req, next);

		blkdev_dequeue_request(next);

		req->biotail->bi_next = next->bio;
		req->biotail = next->biotail;

		req->nr_sectors = req->hard_nr_sectors += next->hard_nr_sectors;

		blkdev_release_request(next);
	}
}

static inline void attempt_back_merge(request_queue_t *q, struct request *rq)
{
	struct list_head *next = rq->queuelist.next;

	if (next != &q->queue_head)
		attempt_merge(q, rq, list_entry_rq(next));
}

static inline void attempt_front_merge(request_queue_t *q, struct request *rq)
{
	struct list_head *prev = rq->queuelist.prev;

	if (prev != &q->queue_head)
		attempt_merge(q, list_entry_rq(prev), rq);
}

/**
 * blk_attempt_remerge  - attempt to remerge active head with next request
 * @q:    The &request_queue_t belonging to the device
 * @rq:   The head request (usually)
 *
 * Description:
 *    For head-active devices, the queue can easily be unplugged so quickly
 *    that proper merging is not done on the front request. This may hurt
 *    performance greatly for some devices. The block layer cannot safely
 *    do merging on that first request for these queues, but the driver can
 *    call this function and make it happen any way. Only the driver knows
 *    when it is safe to do so.
 **/
void blk_attempt_remerge(request_queue_t *q, struct request *rq)
{
	unsigned long flags;

	spin_lock_irqsave(q->queue_lock, flags);
	attempt_back_merge(q, rq);
	spin_unlock_irqrestore(q->queue_lock, flags);
}

static int __make_request(request_queue_t *q, struct bio *bio)
{
	struct request *req, *freereq = NULL;
	int el_ret, rw, nr_sectors, cur_nr_sectors, barrier;
	struct list_head *insert_here;
	sector_t sector;

	sector = bio->bi_sector;
	nr_sectors = bio_sectors(bio);
	cur_nr_sectors = bio_iovec(bio)->bv_len >> 9;
	rw = bio_data_dir(bio);

	/*
	 * low level driver can indicate that it wants pages above a
	 * certain limit bounced to low memory (ie for highmem, or even
	 * ISA dma in theory)
	 */
	blk_queue_bounce(q, &bio);

	spin_lock_prefetch(q->queue_lock);

	barrier = test_bit(BIO_RW_BARRIER, &bio->bi_rw);

	spin_lock_irq(q->queue_lock);
again:
	req = NULL;
	insert_here = q->queue_head.prev;

	if (blk_queue_empty(q) || barrier) {
		blk_plug_device(q);
		goto get_rq;
	}

	el_ret = elv_merge(q, &req, bio);
	switch (el_ret) {
		case ELEVATOR_BACK_MERGE:
			BUG_ON(!rq_mergeable(req));
			if (!q->back_merge_fn(q, req, bio)) {
				insert_here = &req->queuelist;
				break;
			}

			elv_merge_cleanup(q, req, nr_sectors);

			req->biotail->bi_next = bio;
			req->biotail = bio;
			req->nr_sectors = req->hard_nr_sectors += nr_sectors;
			drive_stat_acct(req, nr_sectors, 0);
			attempt_back_merge(q, req);
			goto out;

		case ELEVATOR_FRONT_MERGE:
			BUG_ON(!rq_mergeable(req));
			if (!q->front_merge_fn(q, req, bio)) {
				insert_here = req->queuelist.prev;
				break;
			}

			elv_merge_cleanup(q, req, nr_sectors);

			bio->bi_next = req->bio;
			req->bio = bio;
			/*
			 * may not be valid. if the low level driver said
			 * it didn't need a bounce buffer then it better
			 * not touch req->buffer either...
			 */
			req->buffer = bio_data(bio);
			req->current_nr_sectors = cur_nr_sectors;
			req->hard_cur_sectors = cur_nr_sectors;
			req->sector = req->hard_sector = sector;
			req->nr_sectors = req->hard_nr_sectors += nr_sectors;
			drive_stat_acct(req, nr_sectors, 0);
			attempt_front_merge(q, req);
			goto out;

		/*
		 * elevator says don't/can't merge. get new request
		 */
		case ELEVATOR_NO_MERGE:
			/*
			 * use elevator hints as to where to insert the
			 * request. if no hints, just add it to the back
			 * of the queue
			 */
			if (req)
				insert_here = &req->queuelist;
			break;

		default:
			printk("elevator returned crap (%d)\n", el_ret);
			BUG();
	}

	/*
	 * Grab a free request from the freelist - if that is empty, check
	 * if we are doing read ahead and abort instead of blocking for
	 * a free slot.
	 */
get_rq:
	if (freereq) {
		req = freereq;
		freereq = NULL;
	} else if ((req = get_request(q, rw)) == NULL) {
		spin_unlock_irq(q->queue_lock);

		/*
		 * READA bit set
		 */
		if (bio->bi_rw & (1 << BIO_RW_AHEAD)) {
			set_bit(BIO_RW_BLOCK, &bio->bi_flags);
			goto end_io;
		}

		freereq = get_request_wait(q, rw);
		spin_lock_irq(q->queue_lock);
		goto again;
	}

	/*
	 * first three bits are identical in rq->flags and bio->bi_rw,
	 * see bio.h and blkdev.h
	 */
	req->flags = (bio->bi_rw & 7) | REQ_CMD;

	/*
	 * REQ_BARRIER implies no merging, but lets make it explicit
	 */
	if (barrier)
		req->flags |= (REQ_BARRIER | REQ_NOMERGE);

	req->errors = 0;
	req->hard_sector = req->sector = sector;
	req->hard_nr_sectors = req->nr_sectors = nr_sectors;
	req->current_nr_sectors = req->hard_cur_sectors = cur_nr_sectors;
	req->nr_phys_segments = bio_phys_segments(q, bio);
	req->nr_hw_segments = bio_hw_segments(q, bio);
	req->buffer = bio_data(bio);	/* see ->buffer comment above */
	req->waiting = NULL;
	req->bio = req->biotail = bio;
	req->rq_dev = bio->bi_dev;
	add_request(q, req, insert_here);
out:
	if (freereq)
		blkdev_release_request(freereq);
	spin_unlock_irq(q->queue_lock);
	return 0;

end_io:
	bio->bi_end_io(bio, nr_sectors);
	return 0;
}


/*
 * If bio->bi_dev is a partition, remap the location
 */
static inline void blk_partition_remap(struct bio *bio)
{
	int major, minor, drive, minor0;
	struct gendisk *g;
	kdev_t dev0;

	major = major(bio->bi_dev);
	if ((g = get_gendisk(bio->bi_dev))) {
		minor = minor(bio->bi_dev);
		drive = (minor >> g->minor_shift);
		minor0 = (drive << g->minor_shift); /* whole disk device */
		/* that is, minor0 = (minor & ~((1<<g->minor_shift)-1)); */
		dev0 = mk_kdev(major, minor0);
		if (!kdev_same(dev0, bio->bi_dev)) {
			bio->bi_dev = dev0;
			bio->bi_sector += g->part[minor].start_sect;
		}
		/* lots of checks are possible */
	}
}

/**
 * generic_make_request: hand a buffer to it's device driver for I/O
 * @bio:  The bio describing the location in memory and on the device.
 *
 * generic_make_request() is used to make I/O requests of block
 * devices. It is passed a &struct bio, which describes the I/O that needs
 * to be done.
 *
 * generic_make_request() does not return any status.  The
 * success/failure status of the request, along with notification of
 * completion, is delivered asynchronously through the bio->bi_end_io
 * function described (one day) else where.
 *
 * The caller of generic_make_request must make sure that bi_io_vec
 * are set to describe the memory buffer, and that bi_dev and bi_sector are
 * set to describe the device address, and the
 * bi_end_io and optionally bi_private are set to describe how
 * completion notification should be signaled.
 *
 * generic_make_request and the drivers it calls may use bi_next if this
 * bio happens to be merged with someone else, and may change bi_dev and
 * bi_sector for remaps as it sees fit.  So the values of these fields
 * should NOT be depended on after the call to generic_make_request.
 *
 * */
void generic_make_request(struct bio *bio)
{
	int major = major(bio->bi_dev);
	int minor = minor(bio->bi_dev);
	request_queue_t *q;
	sector_t minorsize = 0;
	int ret, nr_sectors = bio_sectors(bio);

	/* Test device or partition size, when known. */
	if (blk_size[major])
		minorsize = blk_size[major][minor];
	if (minorsize) {
		unsigned long maxsector = (minorsize << 1) + 1;
		unsigned long sector = bio->bi_sector;

		if (maxsector < nr_sectors || maxsector - nr_sectors < sector) {
			if (blk_size[major][minor]) {
				
				/* This may well happen - the kernel calls
				 * bread() without checking the size of the
				 * device, e.g., when mounting a device. */
				printk(KERN_INFO
				       "attempt to access beyond end of device\n");
				printk(KERN_INFO "%s: rw=%ld, want=%ld, limit=%Lu\n",
				       kdevname(bio->bi_dev), bio->bi_rw,
				       (sector + nr_sectors)>>1,
				       (long long) blk_size[major][minor]);
			}
			set_bit(BIO_EOF, &bio->bi_flags);
			goto end_io;
		}
	}

	/*
	 * Resolve the mapping until finished. (drivers are
	 * still free to implement/resolve their own stacking
	 * by explicitly returning 0)
	 *
	 * NOTE: we don't repeat the blk_size check for each new device.
	 * Stacking drivers are expected to know what they are doing.
	 */
	do {
		q = blk_get_queue(bio->bi_dev);
		if (!q) {
			printk(KERN_ERR
			       "generic_make_request: Trying to access nonexistent block-device %s (%Lu)\n",
			       kdevname(bio->bi_dev), (long long) bio->bi_sector);
end_io:
			bio->bi_end_io(bio, nr_sectors);
			break;
		}

		BUG_ON(bio_sectors(bio) > q->max_sectors);

		/*
		 * If this device has partitions, remap block n
		 * of partition p to block n+start(p) of the disk.
		 */
		blk_partition_remap(bio);

		ret = q->make_request_fn(q, bio);
		blk_put_queue(q);

	} while (ret);
}

/*
 * our default bio end_io callback handler for a buffer_head mapping.
 */
static int end_bio_bh_io_sync(struct bio *bio, int nr_sectors)
{
	struct buffer_head *bh = bio->bi_private;

	BIO_BUG_ON(nr_sectors != (bh->b_size >> 9));

	bh->b_end_io(bh, test_bit(BIO_UPTODATE, &bio->bi_flags));
	bio_put(bio);
	return 0;
}

/**
 * submit_bio: submit a bio to the block device layer for I/O
 * @rw: whether to %READ or %WRITE, or maybe to %READA (read ahead)
 * @bio: The &struct bio which describes the I/O
 *
 * submit_bio() is very similar in purpose to generic_make_request(), and
 * uses that function to do most of the work. Both are fairly rough
 * interfaces, @bio must be presetup and ready for I/O.
 *
 */
int submit_bio(int rw, struct bio *bio)
{
	int count = bio_sectors(bio);

	/*
	 * do some validity checks...
	 */
	BUG_ON(!bio->bi_end_io);

	BIO_BUG_ON(!bio->bi_size);
	BIO_BUG_ON(!bio->bi_io_vec);

	bio->bi_rw = rw;

	if (rw & WRITE)
		kstat.pgpgout += count;
	else
		kstat.pgpgin += count;

	generic_make_request(bio);
	return 1;
}

/**
 * submit_bh: submit a buffer_head to the block device layer for I/O
 * @rw: whether to %READ or %WRITE, or maybe to %READA (read ahead)
 * @bh: The &struct buffer_head which describes the I/O
 *
 **/
int submit_bh(int rw, struct buffer_head * bh)
{
	struct bio *bio;

	BUG_ON(!test_bit(BH_Lock, &bh->b_state));
	BUG_ON(!buffer_mapped(bh));
	BUG_ON(!bh->b_end_io);

	set_bit(BH_Req, &bh->b_state);

	/*
	 * from here on down, it's all bio -- do the initial mapping,
	 * submit_bio -> generic_make_request may further map this bio around
	 */
	bio = bio_alloc(GFP_NOIO, 1);

	bio->bi_sector = bh->b_blocknr * (bh->b_size >> 9);
	bio->bi_dev = bh->b_dev;
	bio->bi_io_vec[0].bv_page = bh->b_page;
	bio->bi_io_vec[0].bv_len = bh->b_size;
	bio->bi_io_vec[0].bv_offset = bh_offset(bh);

	bio->bi_vcnt = 1;
	bio->bi_idx = 0;
	bio->bi_size = bh->b_size;

	bio->bi_end_io = end_bio_bh_io_sync;
	bio->bi_private = bh;

	return submit_bio(rw, bio);
}

/**
 * ll_rw_block: low-level access to block devices
 * @rw: whether to %READ or %WRITE or maybe %READA (readahead)
 * @nr: number of &struct buffer_heads in the array
 * @bhs: array of pointers to &struct buffer_head
 *
 * ll_rw_block() takes an array of pointers to &struct buffer_heads,
 * and requests an I/O operation on them, either a %READ or a %WRITE.
 * The third %READA option is described in the documentation for
 * generic_make_request() which ll_rw_block() calls.
 *
 * This function provides extra functionality that is not in
 * generic_make_request() that is relevant to buffers in the buffer
 * cache or page cache.  In particular it drops any buffer that it
 * cannot get a lock on (with the BH_Lock state bit), any buffer that
 * appears to be clean when doing a write request, and any buffer that
 * appears to be up-to-date when doing read request.  Further it marks
 * as clean buffers that are processed for writing (the buffer cache
 * wont assume that they are actually clean until the buffer gets
 * unlocked).
 *
 * ll_rw_block sets b_end_io to simple completion handler that marks
 * the buffer up-to-date (if approriate), unlocks the buffer and wakes
 * any waiters.  As client that needs a more interesting completion
 * routine should call submit_bh() (or generic_make_request())
 * directly.
 *
 * Caveat:
 *  All of the buffers must be for the same device, and must also be
 *  a multiple of the current approved size for the device.
 *
 **/
void ll_rw_block(int rw, int nr, struct buffer_head * bhs[])
{
	unsigned int major;
	int correct_size;
	int i;

	if (!nr)
		return;

	major = major(bhs[0]->b_dev);

	/* Determine correct block size for this device. */
	correct_size = get_hardsect_size(bhs[0]->b_dev);

	/* Verify requested block sizes. */
	for (i = 0; i < nr; i++) {
		struct buffer_head *bh = bhs[i];
		if (bh->b_size & (correct_size - 1)) {
			printk(KERN_NOTICE "ll_rw_block: device %s: "
			       "only %d-char blocks implemented (%u)\n",
			       kdevname(bhs[0]->b_dev),
			       correct_size, bh->b_size);
			goto sorry;
		}
	}

	if ((rw & WRITE) && is_read_only(bhs[0]->b_dev)) {
		printk(KERN_NOTICE "Can't write to read-only device %s\n",
		       kdevname(bhs[0]->b_dev));
		goto sorry;
	}

	for (i = 0; i < nr; i++) {
		struct buffer_head *bh = bhs[i];

		/* Only one thread can actually submit the I/O. */
		if (test_and_set_bit(BH_Lock, &bh->b_state))
			continue;

		/* We have the buffer lock */
		atomic_inc(&bh->b_count);
		bh->b_end_io = end_buffer_io_sync;

		switch(rw) {
		case WRITE:
			if (!atomic_set_buffer_clean(bh))
				/* Hmmph! Nothing to write */
				goto end_io;
			__mark_buffer_clean(bh);
			break;

		case READA:
		case READ:
			if (buffer_uptodate(bh))
				/* Hmmph! Already have it */
				goto end_io;
			break;
		default:
			BUG();
	end_io:
			bh->b_end_io(bh, test_bit(BH_Uptodate, &bh->b_state));
			continue;
		}

		submit_bh(rw, bh);
	}
	return;

sorry:
	/* Make sure we don't get infinite dirty retries.. */
	for (i = 0; i < nr; i++)
		mark_buffer_clean(bhs[i]);
}

#ifdef CONFIG_STRAM_SWAP
extern int stram_device_init (void);
#endif

inline void blk_recalc_rq_segments(struct request *rq)
{
	struct bio *bio;
	int nr_phys_segs, nr_hw_segs;

	rq->buffer = bio_data(rq->bio);

	nr_phys_segs = nr_hw_segs = 0;
	rq_for_each_bio(bio, rq) {
		/* Force bio hw/phys segs to be recalculated. */
		bio->bi_flags &= ~(1 << BIO_SEG_VALID);

		nr_phys_segs += bio_phys_segments(rq->q, bio);
		nr_hw_segs += bio_hw_segments(rq->q, bio);
	}

	rq->nr_phys_segments = nr_phys_segs;
	rq->nr_hw_segments = nr_hw_segs;
}

inline void blk_recalc_rq_sectors(struct request *rq, int nsect)
{
	if (rq->flags & REQ_CMD) {
		rq->hard_sector += nsect;
		rq->hard_nr_sectors -= nsect;
		rq->sector = rq->hard_sector;
		rq->nr_sectors = rq->hard_nr_sectors;

		rq->current_nr_sectors = bio_iovec(rq->bio)->bv_len >> 9;
		rq->hard_cur_sectors = rq->current_nr_sectors;

		/*
		 * if total number of sectors is less than the first segment
		 * size, something has gone terribly wrong
		 */
		if (rq->nr_sectors < rq->current_nr_sectors) {
			printk("blk: request botched\n");
			rq->nr_sectors = rq->current_nr_sectors;
		}
	}
}

/**
 * end_that_request_first - end I/O on one buffer.
 * @req:      the request being processed
 * @uptodate: 0 for I/O error
 * @nr_sectors: number of sectors to end I/O on
 *
 * Description:
 *     Ends I/O on the first buffer attached to @req, and sets it up
 *     for the next buffer_head (if any) in the cluster.
 *     
 * Return:
 *     0 - we are done with this request, call end_that_request_last()
 *     1 - still buffers pending for this request
 **/

int end_that_request_first(struct request *req, int uptodate, int nr_sectors)
{
	int nsect, total_nsect;
	struct bio *bio;

	req->errors = 0;
	if (!uptodate)
		printk("end_request: I/O error, dev %s, sector %lu\n",
			kdevname(req->rq_dev), req->sector);

	total_nsect = 0;
	while ((bio = req->bio)) {
		nsect = bio_iovec(bio)->bv_len >> 9;

		BIO_BUG_ON(bio_iovec(bio)->bv_len > bio->bi_size);

		/*
		 * not a complete bvec done
		 */
		if (unlikely(nsect > nr_sectors)) {
			int residual = (nsect - nr_sectors) << 9;

			bio->bi_size -= residual;
			bio_iovec(bio)->bv_offset += residual;
			bio_iovec(bio)->bv_len -= residual;
			blk_recalc_rq_sectors(req, nr_sectors);
			blk_recalc_rq_segments(req);
			return 1;
		}

		/*
		 * account transfer
		 */
		bio->bi_size -= bio_iovec(bio)->bv_len;
		bio->bi_idx++;

		nr_sectors -= nsect;
		total_nsect += nsect;

		if (!bio->bi_size) {
			req->bio = bio->bi_next;

			if (unlikely(bio_endio(bio, uptodate, total_nsect)))
				BUG();

			total_nsect = 0;
		}

		if ((bio = req->bio)) {
			blk_recalc_rq_sectors(req, nsect);

			/*
			 * end more in this run, or just return 'not-done'
			 */
			if (unlikely(nr_sectors <= 0)) {
				blk_recalc_rq_segments(req);
				return 1;
			}
		}
	}

	return 0;
}

void end_that_request_last(struct request *req)
{
	if (req->waiting)
		complete(req->waiting);

	blkdev_release_request(req);
}

#define MB(kb)	((kb) << 10)

int __init blk_dev_init(void)
{
	struct blk_dev_struct *dev;
	int total_ram;

	request_cachep = kmem_cache_create("blkdev_requests",
					   sizeof(struct request),
					   0, SLAB_HWCACHE_ALIGN, NULL, NULL);

	if (!request_cachep)
		panic("Can't create request pool slab cache\n");

	for (dev = blk_dev + MAX_BLKDEV; dev-- != blk_dev;)
		dev->queue = NULL;

	memset(ro_bits,0,sizeof(ro_bits));
	memset(max_readahead, 0, sizeof(max_readahead));

	total_ram = nr_free_pages() << (PAGE_SHIFT - 10);

	/*
	 * Free request slots per queue.
	 * (Half for reads, half for writes)
	 */
	queue_nr_requests = 64;
	if (total_ram > MB(32))
		queue_nr_requests = 256;

	/*
	 * Batch frees according to queue length
	 */
	if ((batch_requests = queue_nr_requests / 4) > 32)
		batch_requests = 32;
	printk("block: %d slots per queue, batch=%d\n", queue_nr_requests, batch_requests);

	blk_max_low_pfn = max_low_pfn;
	blk_max_pfn = max_pfn;

#if defined(CONFIG_IDE) && defined(CONFIG_BLK_DEV_IDE)
	ide_init();		/* this MUST precede hd_init */
#endif
#if defined(CONFIG_IDE) && defined(CONFIG_BLK_DEV_HD)
	hd_init();
#endif
#if defined(__i386__)	/* Do we even need this? */
	outb_p(0xc, 0x3f2);
#endif

	return 0;
};

EXPORT_SYMBOL(end_that_request_first);
EXPORT_SYMBOL(end_that_request_last);
EXPORT_SYMBOL(blk_init_queue);
EXPORT_SYMBOL(blk_get_queue);
EXPORT_SYMBOL(blk_cleanup_queue);
EXPORT_SYMBOL(blk_queue_make_request);
EXPORT_SYMBOL(blk_queue_bounce_limit);
EXPORT_SYMBOL(generic_make_request);
EXPORT_SYMBOL(blkdev_release_request);
EXPORT_SYMBOL(generic_unplug_device);
EXPORT_SYMBOL(blk_attempt_remerge);
EXPORT_SYMBOL(blk_max_low_pfn);
EXPORT_SYMBOL(blk_queue_max_sectors);
EXPORT_SYMBOL(blk_queue_max_phys_segments);
EXPORT_SYMBOL(blk_queue_max_hw_segments);
EXPORT_SYMBOL(blk_queue_max_segment_size);
EXPORT_SYMBOL(blk_queue_hardsect_size);
EXPORT_SYMBOL(blk_queue_segment_boundary);
EXPORT_SYMBOL(blk_rq_map_sg);
EXPORT_SYMBOL(blk_nohighio);
EXPORT_SYMBOL(blk_dump_rq_flags);
EXPORT_SYMBOL(submit_bio);
EXPORT_SYMBOL(blk_queue_assign_lock);
EXPORT_SYMBOL(blk_phys_contig_segment);
EXPORT_SYMBOL(blk_hw_contig_segment);

EXPORT_SYMBOL(ll_10byte_cmd_build);
EXPORT_SYMBOL(blk_queue_prep_rq);
