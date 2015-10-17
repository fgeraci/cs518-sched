
/* /include/linux/sched.h */

extern void scheduler_tick(int user_tick, int system);
typedef struct task_struct task_t; // opaque of task_struct
typedef struct prio_array prio_array_t;

/* /kernel/sched.c */

/*
	System timer that the kernel periodically
	calls and MIGHT (check logic) marks processes as 
	needing rescheduling, if interrupted. This will
	trigger schedule()
	
	Caller: 
		timer.c 	- update_process_times
		sched.c		- sched_fork
	Implement
		sched.c
	Ref
		sched.h
*/

/* Queue struct and opaque */

typedef struct runqueue runqueue_t;

struct runqueue {
	spinlock_t lock;
    unsigned long 	nr_running, nr_switches, expired_timestamp,
					nr_uninterruptible, timestamp_last_tick;
					task_t *curr, *idle;
    struct mm_struct *prev_mm;
    prio_array_t *active, *expired, arrays[2];
    int best_expired_prio, prev_cpu_load[NR_CPUS];
 #ifdef CONFIG_NUMA
    atomic_t *node_nr_running;
    int prev_node_load[MAX_NUMNODES];
 #endif
    task_t *migration_thread;
    struct list_head migration_queue;
	atomic_t nr_iowait;
};

void schedule_tick(int user_ticks, int sys_ticks) {
	/* Logic here */
}

struct prio_array {
	/* - TODO - check data members exist
	int nr_active;
	unsigned long bitmap[BITMAP_SIZE];
	struct list_head queue[MAX_PRIO];
	*/
};

// 602. schedule(void) {
		// data members
		long *switch_count;
		task_t *prev, *next;		// CS518 - opaquing task_struct
		runqueue_t *rq;				// pointer to current queue
		prio_array *array;			// priority levels
		struct list_head *queue;	// defined OK - list.h
		unsigned long run_time;		// total runtime allowed
		unsigned long long now; 	// just now double long - 64 bit unsigned integer
		int idx;
		// end data members
// }

/* >>>>> EXTERNAL TO sched.c / sched.h - TIMER / FORK ... <<<<< */

/* /kernel/timer.c */

	// update_process_times
		scheduler_tick(user_tick, system);
	//