
#include <sys/time.h>
#include <stdio.h>

#define T_GR_HOOK "--->>"


void init_timing(int max, char *time_file);


// gets the current time and stores it for later use by the print_times() function.
//	Succesive calls to this function will overwrite the previous time!
void get_start_time();


// gets the current time and stores it for later use by the print_times() function.
//	Succesive calls to this function will overwrite the previous time!
void get_stop_time();


// prints times (sec) recorded by the get_start_time() and get_stop_time() 
//	functions to stdout. The line is preceded by T_GR_HOOK
void complete_timing();
