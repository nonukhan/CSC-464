/* * * * * * * *
 * Stephen Tredger, March 2012
 *
 * baseline.c
 *	simple test to get overhead of collecting timestamps
 */


#define MAXTIMES 1000

#include "timekeeper.h"

int main() {
	int i;
	timekeeper *tk;
	
	tk = create_timekeeper(MAXTIMES, NULL, NULL);
		
	for (i=0; i<MAXTIMES; i++) {
		get_start_time(tk);
		get_stop_time(tk);
	}
	
	print_timekeeper(tk);
	destroy_timekeeper(tk);
	return 0;
}
