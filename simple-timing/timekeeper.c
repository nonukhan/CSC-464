/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *	Stephen Tredger, March 2012
 *
 *	timekeeper.c
 *		Used to gather timestamps for profiling of various functions.
 *		timekeepers are structures that hold up to a set of MAX_TIMES start/stop
 *		timestamps.
 */

#include "timekeeper.h"


// wrapper for malloc to avoid clutter of error checking when mallocing
void * err_malloc(size_t len) {
	void *buff;
	
	if( !(buff = malloc(len)) ) {
		fprintf(stderr, "malloc(): failed to allocate %u bytes\n", (unsigned int)len);
		exit(T_EXIT_FAIL);
	}
	return buff;
}


// wrapper for calloc so we can call like malloc and don't have to have clutter
//	of error checking when mallocing elsewhere
void * clr_malloc(size_t len) {
	void *buff;
	
	if( !(buff = calloc(len, sizeof(char))) ) {
		fprintf(stderr, "malloc(): failed to allocate %u bytes\n", (unsigned int)len);
		exit(T_EXIT_FAIL);
	}
	return buff;
}


// allocate memory for and initializes a timekeeper stuct
timekeeper * create_timekeeper(char *filename, char *hook) {
	timekeeper *tk;
	size_t name_size;
	
	tk = (timekeeper *) clr_malloc(sizeof(timekeeper));	
	
	// if we have a filename place it in timefile, otherwise leave it NULL
	if (filename) {
		if( (name_size = strlen(filename)) > 0) {
			tk->timefile = (char *) clr_malloc(name_size+1);
			memcpy(tk->timefile, filename, name_size);
			tk->timefile[name_size] = '\0'; // make sure we have a NULL terminator
		}
	}
	
	// if we have a grep hook place it in
	if (hook) {
		if (strlen(hook) < GR_HOOK_LEN) {
			memcpy(tk->gr_hook, hook, strlen(hook));
			tk->gr_hook[strlen(hook)] = '\0';		
		} else {
			memcpy(tk->gr_hook, hook, GR_HOOK_LEN);
			tk->gr_hook[GR_HOOK_LEN] = '\0';
		}
	}
	
	return tk;
}


// frees timekeeper struct
void destroy_timekeeper(timekeeper *tk) {
	
	// nulling after a free may not be thread safe??????
	free(tk->timefile);
	tk->timefile = NULL;
	free(tk);
	tk = NULL;
}


// gets a timestamp and places into the start timeval array in tk. 
//	Returns if the number of start/stop timestamps are not equal
//	or we already have MAX_TIMES timestamps
void get_start_time(timekeeper *tk) {
	
	assert(tk);
	if (tk->start_t_num != tk->stop_t_num) {
		fprintf(stderr, "get_start_time(): previous timer not stopped!\n");
		return;
	}
	if (tk->start_t_num >= MAX_TIMES) {
		fprintf(stderr, "get_start_time(): timekeeper full!\n");
		return;
	}
	
	// get timestamp and increment start_t_num
	gettimeofday(&tk->start[tk->start_t_num++], NULL);
}


// gets a timestamp and places into the stop timeval array in tk. 
//	Returns if the number of start timestamps isn't exactly
//	one greater than the number of stop timestamps
void get_stop_time(timekeeper *tk) {
	
	assert(tk);
	if (tk->start_t_num != tk->stop_t_num+1) {
		fprintf(stderr, "get_stop_time(): timer not started!\n");
		return;
	}
	assert(tk->stop_t_num < MAX_TIMES);
	
	// get timestamp and increment stop_t_num
	gettimeofday(&tk->stop[tk->stop_t_num++], NULL);
}


// print timestamps (in secs) to timefile. If timefile is NULL or
//	can't be opened prints to stdout instead
//	Returns if the number of timestamps are unequal.
void print_timekeeper(timekeeper *tk) {
	int i;
	FILE *outfile;
	
	assert(tk);
	if (tk->start_t_num != tk->stop_t_num) {
		fprintf(stderr, "print_timekeeper(): unequal number of start/stop times!\n");
		return;
	}
	
	if (tk->timefile) {
		if( !(outfile = fopen(tk->timefile, "a+"))) {
			perror("print_timekeeper(): failed to open timing file.\nprinting to stdout\n");
			outfile = stdout;
		}
	} else {
		outfile = stdout;
	}
	
	for (i=0; i<tk->start_t_num; i++) {
		// print "<gr_hook> <start_t> <stop_t>" timefile
		fprintf(outfile, "%s %ld.%06ld, %ld.%06ld\n",
				tk->gr_hook,
				(long int)tk->start[i].tv_sec, 
				(long int)tk->start[i].tv_usec, 
				(long int)tk->stop[i].tv_sec, 
				(long int)tk->stop[i].tv_usec);
	}
}
