void MPIInitialize(int *argc, char **argv);
boolean Host(int IAM);
boolean Node(int IAM);
long MPINodeWorkAmt(long datasets);
long MPINodeStartPosn(long datasets);
void MPICopyWholeFile(const char *xname, char *xnamestr, char *xprog);
void MPIOpenFile(FILE **fileptr, char *name, char *namestr, char *mode,
		 char *prog, char *name2);
void MPIOpenStdinFile(FILE *fileptr, char *name, char *namestr, char *mode,
		      char *prog, char *name2);
void MPICopyPartialFile(char *name, char *namestr, char *prog, long datasets);
