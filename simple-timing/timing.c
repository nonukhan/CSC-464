

#include "timing.h"

// global variables to hold times
struct timeval *global_start_t = NULL, *global_stop_t = NULL;
int global_curr_t = 0, global_max_t = 0, global_time_file;

void init_timing(int max, char *time_file) {

	if( (global_time_file = open(time_file, O_RDWR | O_CREAT)) == -1) {
		fprintf(stderr, "failed to open timing file.\ntiming initialization Failed\n");
		return;
	}

	global_max_t = max;
	global_curr_t = 0;
	global_start_t = (struct timeval *) malloc(max*sizeof(struct timeval));
	global_stop_t = (struct timeval *) malloc(max*sizeof(struct timeval));
}


// gets the current time and places it in global_start_t 
void get_start_time() {	
	if (global_start_t) {
		gettimeofday(&global_start_t[global_curr_t], NULL);
	}
}


// gets the current time and places it in global_stop_t 
void get_finish_time() {
	if (global_stop_t) {
		gettimeofday(&global_stop_t[global_curr_t], NULL);
		global_curr_t++;
	}
}


// prints both global_start_t and global_stop_t to stdout. The line is preceded by T_GR_HOOK
void complete_timing() {
	int i;
	if (global_start_t && global_stop_t) {
		for (i=0; i<global_curr_t; i++) {
			fprintf(global_time_file, "%s %ld.%06ld, %ld.%06ld\n",
					T_GR_HOOK,
					(long int)global_start_t[i].tv_sec, 
					(long int)global_start_t[i].tv_usec, 
					(long int)global_stop_t[i].tv_sec, 
					(long int)global_stop_t[i].tv_usec);
		}
		free(global_start_t);
		free(global_stop_t);
	} else {
		fprintf(stderr, "Timing not initialized!\n");
	}
}
