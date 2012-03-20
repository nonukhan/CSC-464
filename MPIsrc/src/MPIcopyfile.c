#include <mpi.h>
#include <string.h>
#include <unistd.h>  /* needed for unlink */
#include "phylip.h"


#define BSZ 256

void MPIInitialize(int *argc, char **argv)
{
  int size, rc;
  rc = MPI_Init(argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
      printf("\nMPI VERSION REQUIRES AT LEAST 2 PROCESSORS.\n");
      MPI_Finalize();
      exit(0); 
      }
  }

boolean Host(int IAM)
  {
  if ( IAM == 0)                                             
    return true;
  else return false;
  }

boolean Node(int IAM)
  {
  if ( IAM != 0)                                             
    return true;
  else return false;
  }

long MPINodeWorkAmt(long datasets)
   {
   int nodework, SIZE, IAM;
   MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
   MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
   nodework=datasets/SIZE;
   if ((datasets % SIZE) > IAM) nodework++;
   return (nodework);
   }

/****************************************************************/
/* Returns the unique node starting position for the input file */
/****************************************************************/
long MPINodeStartPosn(long datasets)
   {
   int nodework, startp, SIZE, IAM, i;
   MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
   MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
   if (Host(IAM)) startp=1;
   else {
      startp=1;
      for(i=0;i<IAM;i++) {
        nodework=datasets/SIZE;
        if ((datasets % SIZE) > i) nodework++;
        startp=startp+nodework;
        }
      }
   return (startp);
   }
   

/**************************************************************************/
/* This routine copies a file that exists on the host processor to all    */
/* processors. THE COPIED FILE IS THE SAME NAME AS THE ORIGINAL           */
/* FILENAME (*name) EXCEPT THAT HAS A .NODENUM ADDED. E.G.                */
/* If the original file is named FILENAME then                            */
/*    On Node 0 copied filename = FILENAME.0                              */
/*    On Node 1 copied filename = FILENAME.1                              */
/*    On Node 2 copied filename = FILENAME.2  etc.                        */
/**************************************************************************/
void MPICopyWholeFile(const char *name, char *namestr, char *prog)
  {
    int i;
    char *j;
    FILE *origfile;
    FILE *copiedfile;
    static MPI_Status status;
    static int IAM, SIZE;
    char namewnode[BSZ];
    char namewnode2[BSZ];
    char name2[BSZ];

    static char buffer[BSZ];
    for(i=0;i<BSZ;i++) buffer[i]=' ';

    MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);

    strncpy(name2,name,BSZ);

    sprintf(namewnode,"%s.%d",name2,IAM);
    strncpy(namewnode2,namewnode,BSZ);
    copiedfile = fopen(namewnode,"w");
    /*openfile(&copiedfile,namewnode,namestr,"w",prog,namewnode2);*/
    if (Host(IAM)) {
      origfile = fopen(name,"r");
      /* openfile(&origfile,name,namestr,"r",prog,name2); */
      /*
      printf("IAM HOST\n"); 
      fflush(stdout);
      */ 
      for(;;) {
        j = fgets(buffer, BSZ, origfile);
        /*printf("j=%s\n",j);*/
        if ( j == NULL) {
	  for (i=1; i<SIZE; i++) {
	    buffer[0]='\0';
	    MPI_Send(buffer,BSZ,MPI_CHAR,i,0,MPI_COMM_WORLD);
	    }
          /*FClose(origfile);*/
          fclose(origfile);
          break;
	  }
        else {
          fputs(buffer,copiedfile);
          for (i=1; i<SIZE; i++) {
	    MPI_Send(buffer,BSZ,MPI_CHAR,i,0,MPI_COMM_WORLD);
            /*printf("HOST SENT BUFFER=%s\n",buffer); fflush(stdout); */
	    }
          }
        }
      }
    if (Node(IAM)) {
      /*
      printf("IAM a node\n"); 
      fflush(stdout);
      */ 
      for(;;) {
	 /*printf("Node is about to receive"); fflush(stdout); */
         MPI_Recv(buffer,BSZ,MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
         /*
         printf("BUFFER=\n"); fflush(stdout);
         printf("%s",buffer); fflush(stdout);
         */
         if (buffer[0] == '\0') break;
         fputs(buffer, copiedfile);
         }
      }
    fclose(copiedfile);
    /*FClose(copiedfile); */
  }
/**************************************************************************/
/* This routine copies the processor specific filenames into a combined   */
/* output filename on the host processor.                                 */
/* THE COPIED FILE IS THE SAME NAME AS THE PROCESSOR SPECIFIC             */
/* FILENAME (*name) EXCEPT THAT IT DOES NOT HAVE A .NODENUM ADDED. E.G.   */
/* The copied  file is named FILENAME then                                */
/*    On Node 0 filename = FILENAME.0                                     */
/*    On Node 1 filename = FILENAME.1                                     */
/*    On Node 2 filename = FILENAME.2  etc.                               */
/**************************************************************************/
void MPICollectOutputFiles(const char *name, char *namestr, char *prog)
  {
    int i;
    char *j;
    FILE *origfile;
    FILE *copiedfile;
    static MPI_Status status;
    int IAM, SIZE;
    char namewnode[BSZ];
    char namewnode2[BSZ];
    char name2[BSZ];

    char buffer[BSZ];

    MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);

    strncpy(name2,name,BSZ);

    sprintf(namewnode,"%s.%d",name2,IAM);
    strncpy(namewnode2,namewnode,BSZ);
    copiedfile = fopen(namewnode,"r");
    /*openfile(&copiedfile,namewnode,namestr,"r",prog,namewnode2);*/
    if (Node(IAM)) {
      /*
      printf("IAM NODE\n"); 
      fflush(stdout);
      */ 
      /* This makes me wait until the host is ready for my data */
      MPI_Recv(buffer,BSZ,MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
      for(;;) {
        j = NULL;  /* This seems to clear up seg fault on J */
        j = fgets(buffer, BSZ, copiedfile);
        /*printf("j=%s\n",j);*/
        if ( j == NULL) {
            buffer[0]='\0';
	    MPI_Send(buffer,BSZ,MPI_CHAR,0,0,MPI_COMM_WORLD);
            break;
	  }
        else {
	    MPI_Send(buffer,BSZ,MPI_CHAR,0,0,MPI_COMM_WORLD);
	    }
        }
      /*FClose(copiedfile);*/
      fclose(copiedfile);
      }
    if (Host(IAM)) {
      origfile = fopen(name,"w");
      /* openfile(&origfile,name,namestr,"w",prog,name2); */
      /*
      printf("IAM Host\n"); 
      fflush(stdout);
      */ 
      /***********************************/
      /* First, copy what is on the HOST */
      /***********************************/
      for(;;) {
        j = NULL;  /* This seems to clear up seg fault on J */
        j = fgets(buffer, BSZ, copiedfile);
        if ( j == NULL) break;
        fputs(buffer, origfile);
        }
      /***********************************/
      /* Then, get stuff from the nodes  */
      /***********************************/
      for(i=1;i<SIZE;i++) {
	 /* Tell the node that I want it to start sending its data */
	 MPI_Send(buffer,BSZ,MPI_CHAR,i,0,MPI_COMM_WORLD);
         for(;;) {
           MPI_Recv(buffer,BSZ,MPI_CHAR,i,0,MPI_COMM_WORLD, &status);
           if (buffer[0] == '\0') break;
           fputs(buffer, origfile);
           }
      }
    fclose(origfile);
    /*FClose(origfile); */
    /*FClose(copiedfile);*/
    fclose(copiedfile);
    }
  }

/**************************************************************************/
/* This routine opens processor specific filenames. Processor specific    */
/* filenames begin with (*name) with the .PROCESSORNUM appended to the    */
/* end of the file. For example, if (*name) is named FILENAME then        */
/*    On Node 0 filename = FILENAME.0 is opened                           */
/*    On Node 1 filename = FILENAME.1 is opened                           */
/*    On Node 2 filename = FILENAME.2 is opened etc.                      */
/**************************************************************************/
void MPIOpenFile(FILE **fileptr, char *name, char *namestr, char *mode,
    char *prog, char *name2)
  {
    char buffer[BSZ];
    char namewnode[BSZ];
    int i;
    char *j;
    int IAM;
    MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
    sprintf(namewnode,"%s.%d",name,IAM);
    sprintf(name2,"%s.%d",name,IAM);
    openfile(fileptr,namewnode,namestr,mode,prog,name2);
  }

/**************************************************************************/
/* This routine opens the processor specific standard input file          */
/* Uses REOPEN to the stdin unit.  Processor specific                     */
/* filenames begin with (*name) with the .PROCESSORNUM appended to the    */
/* end of the file. For example, if (*name) is named FILENAME then        */
/*    On Node 0 filename = FILENAME.0 is opened                           */
/*    On Node 1 filename = FILENAME.1 is opened                           */
/*    On Node 2 filename = FILENAME.2 is opened etc.                      */
/**************************************************************************/
void MPIOpenStdinFile(FILE *fileptr, char *name, char *namestr, char *mode,
    char *prog, char *name2)
  {
    char buffer[BSZ];
    char namewnode[BSZ];
    int i;
    char *j;
    int IAM;
   
    MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
    sprintf(namewnode,"%s.%d",name,IAM);
    sprintf(name2,"%s.%d",name,IAM);
    freopen(namewnode,mode,fileptr);
    /*printf("%d opened up file name %s (%s) as fileptr %ld\n",IAM,namewnode,name2,fileptr); */
  }

void MPICleanupFiles(char *name)
{
   char namewnode[BSZ];
   int i;
   int IAM;

   MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
   sprintf(namewnode,"%s.%d",name,IAM);
   unlink(namewnode);
}


/**************************************************************************/
/* This routine copies a file that exists on the host processor and       */
/* divides that file across all processors. THE COPIED FILE IS THE SAME   */
/* NAME AS THE ORIGINAL FILENAME (*name) EXCEPT THAT HAS A .NODENUM ADDED */
/* E.G. If the original file is named FILENAME then                       */
/*    On Node 0 copied partial file is named = FILENAME.0                 */
/*    On Node 1 copied partial file is named = FILENAME.1                 */
/*    On Node 2 copied partial file is named = FILENAME.2  etc.           */
/**************************************************************************/
void MPICopyPartialFile(char *name, char *namestr, char *prog, long datasets)
  {
    char buffer[BSZ];
    char endbuffer[BSZ];
    char namewnode[BSZ];
    int i;
    int setnum;
    char *j;
    FILE *origfile;
    FILE *copiedfile;
    static MPI_Status status;
    int IAM, SIZE;
    long nodework; 
    int distnum,linestoread;
    long numspecies, alignwidth;
    MPI_Comm_rank(MPI_COMM_WORLD, &IAM);
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
    if (datasets < (long) SIZE) {
      printf("\nMPI VERSION REQUIRES NUMBER OF DATASETS GREATER THAN\n");
      printf("\nNUMBER OF PROCESSORS.\n");
      MPI_Finalize();
      exit(0);
    }
    sprintf(namewnode,"%s.%d",name,IAM);
    openfile(&copiedfile,namewnode,namestr,"w",prog,namewnode);
    if (Host(IAM)) {
      openfile(&origfile,name,namestr,"r",prog,name);
      nodework=datasets/ (long) SIZE;
      for (setnum=0; setnum<SIZE; setnum++) {
	if ((datasets % SIZE) > IAM) nodework++;
        /********************************************/
        /* Distribute the file here. This processor */
        /* gets "nodework" datasets                 */
        /********************************************/
        for (distnum=1;  distnum < nodework; distnum++) {
	    fscanf(infile, "%ld%ld", &numspecies, &alignwidth);    
            linestoread=((alignwidth / 60) * numspecies) + 
                         (alignwidth % 60) * numspecies;
	    for (i = 0; i < linestoread; i++) {
               j=fgets(buffer, BSZ, origfile);
               if ( j == NULL) printf("FATAL ERROR\n");
               if (setnum == 0) {
                  fputs(buffer,copiedfile);
	          }
               else {
	          MPI_Send(buffer,BSZ,MPI_CHAR,setnum,0,MPI_COMM_WORLD);
	          }
	       }
 	    }
        if (setnum != 0) {
	   MPI_Send("\0",BSZ,MPI_CHAR,setnum,0,MPI_COMM_WORLD);
	   }
        }
      FClose(origfile);
      }
    if (Node(IAM)) {
      for(;;) {
         MPI_Recv(buffer,BSZ,MPI_CHAR,0,0,MPI_COMM_WORLD, &status);
         if (buffer[0] == '\0') break;
         fputs(buffer, copiedfile);
         }
      }
    FClose(copiedfile);
  }


