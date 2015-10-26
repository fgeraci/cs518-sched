#ifndef __LINUX_SMPLOCK_H
#define __LINUX_SMPLOCK_H

#include <linux/config.h>

#ifndef CONFIG_SMP

#define lock_kernel()				do { } while(0)
#define unlock_kernel()				do { } while(0)
#define release_kernel_lock(task, cpu)		do { } while(0)
#define reacquire_kernel_lock(task)		do { } while(0)
#define kernel_locked() 1
/*
* Re-acquire the kernel lock
*/
static inline void reacquire_kernel_lock(struct task_struct *task)
{
         if (unlikely(task->lock_depth >= 0))
                get_kernel_lock();
}

#else

#include <asm/smplock.h>

#endif /* CONFIG_SMP */

#endif
