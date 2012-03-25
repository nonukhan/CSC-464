/* * * * * * * *
 * Stephen Tredger, March 2012
 *
 * baseline.c
 *	simple test to get overhead of collecting timestamps
 */


#include "timekeeper.h"

int main() {
	int i;
	timekeeper *tk;
	
	tk = create_timekeeper(NULL, NULL);
		
	for (i=0; i<MAX_TIMES; i++) {
		get_start_time(tk);
		get_stop_time(tk);
	}
	
	print_timekeeper(tk);
	destroy_timekeeper(tk);
	return 0;
}
