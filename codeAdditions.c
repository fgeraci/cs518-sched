
/* /include/linux/sched.h */

extern void scheduler_tick(int user_tick, int system);

/* /kernel/sched.c */

/*
	System timer that the kernel periodically
	calls and marks processes as needing rescheduling.
	
	Caller: 
		timer.c 	- update_process_times
		sched.c		- sched_fork
	Implement
		sched.c
	Ref
		sched.h
*/
void schedule_tick(int user_ticks, int sys_ticks) {
	/* Logic here */
}


/* /kernel/timer.c */

	// update_process_times
		scheduler_tick(user_tick, system);
	//