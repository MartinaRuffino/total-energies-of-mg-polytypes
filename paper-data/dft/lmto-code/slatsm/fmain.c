/* Main and support programs to link unix shell with fortran programs. */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI
#include "mpi.h"
#endif

/* --- Parameters that define Interface between fortran and C function calls ---
   The following must be set as compiler switches

   CMAIN           is the name of the entry point function, e.g. CMAIN=main

   ... the following describe how fortran function names are mangled
   FC_UNDERSCORE   if 0, means fortran appends no underscore to a function
                   if 1, means fortran appends an underscore to a function
                   if 2, means fortran appends two underscores to function names that already contain underscores

   FC_UPPERCASE    if 1, means fortran converts a function name to upper case


   ... The following are used for extracting command-line arguments from fortran
       There are two function calls fmain.c supplies for the command-line argument:
       nargc()        returns the number of arguments
       gtargc(iarg,s) returns particular argument iarg in string s.


   NOPASSEDARGS=#  if #=0 then argc, argv are passed to program entry point CMAIN
                   Functions nargc and gtargc just return data from supplied argc, argv.

                   if #=1 then argc, argv are not passed to program entry point CMAIN
                   and in an initialization step fmain.c uses the function calls to
                   extract the information and store it locally.
                   fmain then functions in the same way as NOPASSEDARGS = 0.

   ... the following three tokens are used when NOPASSEDARGS=1
   NARGFCALL=fn    name of function call that returns # arguments
                   It need not be supplied; however, 
                   nargc() always returns 0.

   ADD_TO_NARGFCALL=strn  (optionally used in conjunction with NARGFCALL=fn)
                   number-of-arguments = fn + strn
                   Designed for implementations that return something different from # args, e.g.
                   # args = iargc() + 1

   ARGFCALL=fn     name of function call that returns one entry in argv
                   It need not be supplied; however, gtargc is then not defined.

*/

/* macro FC_FUNC converts a standard name to a fortranized version */
/* FC_FUNC must be bypassed for sun systems */
#if FC_UPPERCASE == 1
#  if FC_UNDERSCORE == 1
#    define FC_FUNC(x,X) X ## _
#  else
#    if FC_UNDERSCORE == 2
#      define FC_FUNC(x,X) X ## _  /*two underscores are added ONLY to names that already contain an underscore! */
#    else
#      define FC_FUNC(x,X) X
#    endif
#  endif
#else
#  if FC_UNDERSCORE == 1
#    define FC_FUNC(x,X) x ## _
#  else
#    if FC_UNDERSCORE == 2
#      define FC_FUNC(x,X) x ## _
#    else
#      define FC_FUNC(x,X) x
#    endif
#  endif
#endif

void FC_FUNC(cexit,CEXIT)();

static int ret_val=0, aargc=0; static char **pargv=NULL;

/* --- Entry point --- */
#if NOPASSEDARGS == 0
int CMAIN(argc,argv)
int argc; char *argv[];
#else
int CMAIN()
#endif

{

  /* Used when a copy of command-line arguments is to be made */
#define MAXARGS 300
#define MAX_ARG_LEN 256
  char *argv_copy[MAXARGS];

#if (FC_UNDERSCORE == 0 && FC_UPPERCASE == 1)
  void FMAIN();
#endif
#if (FC_UNDERSCORE == 1 && FC_UPPERCASE == 1)
  void FMAIN_();
#endif
#if (FC_UNDERSCORE == 2 && FC_UPPERCASE == 1)
  void FMAIN_();
#endif
#if (FC_UNDERSCORE == 0 && FC_UPPERCASE == 0)
  void fmain();
#endif
#if (FC_UNDERSCORE == 1 && FC_UPPERCASE == 0)
  void fmain_();
#endif
#if (FC_UNDERSCORE == 2 && FC_UPPERCASE == 0)
  void fmain_();
#endif

/* --------------- fortran-to-C linkage ---------------- */
#if NOPASSEDARGS == 0
  aargc = argc;
  pargv = argv;
#else
#ifdef NARGFCALL
  int NARGFCALL();  size_t bytes; int len, i,j; char pargi[MAX_ARG_LEN];
#ifdef ARGFCALL
#define HAVE_GTARGC
  void ARGFCALL();
#endif

  aargc = NARGFCALL();
#ifdef ADD_TO_NARGFCALL
  aargc += ADD_TO_NARGFCALL;
#endif

#ifdef ARGFCALL
  for (i = 0; i < aargc; i++) {

/*  Fortran call, using fortran conventions :
    CALL ARGFCALL(i,pargi)
*/
    len = MAX_ARG_LEN-1;
    ARGFCALL(&i,pargi,len);
    j = len;
    while (--j > 0 && pargi[j] == ' ') { len--; } pargi[len] = 0;

    len = 1 + strlen(pargi);  bytes = len * sizeof(char);
    if ( !(argv_copy[i] = malloc(bytes)) ) {
      printf("fmain.c unable to allocate char array\n");
      exit(-1);
    }
    strcpy(argv_copy[i], pargi);

/*      printf("argument %d is %s\n",i,pargi); */

  }
  pargv = argv_copy;
#endif   /* ARGFCALL */
#endif   /* NARGFCALL */
#endif   /* NOPASSEDARGS */

/*    printf("aargc assigned to %d\n",aargc); */

#ifdef MPI
  {
  int argc_mpi = aargc; char **pargv_mpi = pargv;
  char **argv_mpicopy;

/*    printf("argc is %d\n",aargc);   printf("argc_mpi is %d\n",argc_mpi); */

  MPI_Init(&argc_mpi,&pargv_mpi);

/*    printf("argc_mpi is after MPI_Init %d\n",argc_mpi); */

 {
  size_t bytes;
  int len, i, procid, master;

  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
  master = 0;
/*MPI_Bcast(&aargc,1,MPI_INT,master,MPI_COMM_WORLD); */
  MPI_Bcast(&argc_mpi,1,MPI_INT,master,MPI_COMM_WORLD);

  argv_mpicopy = malloc(argc_mpi * sizeof(argv_mpicopy));

  for (i = 0; i < argc_mpi; i++) {
    if (procid == master) {
      len = 1 + strlen(pargv[i]);
    }
    MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
    bytes = len * sizeof(char);
    if ( !(argv_mpicopy[i] = malloc(bytes)) ) {
      printf("Process %d: unable to allocate %d bytes\n",procid,(unsigned int)bytes);
      exit(-1);
    }
    if (procid == master) {
       strcpy(argv_mpicopy[i], pargv[i]);
    }
    MPI_Bcast(argv_mpicopy[i],len,MPI_CHAR,master,MPI_COMM_WORLD);
  }
 }
  aargc = argc_mpi;
  pargv = argv_mpicopy;
/*   printf("aargc is now %d:\n",aargc); */
  }
#endif

#ifdef AIX
  save_me();  /* so we won't get killed when page space is low */
#endif

/* --- Pass control to routine fmain --- */

/* this is function call fmain_(), appropriately mangled for fortran conventions*/
/*  fmain_(); */
  FC_FUNC(fmain,FMAIN)(); /*FC_FUNC does not work for Sun systems*/

 /* call cexit for cleanup */
 {
   int i1=0,i2=1;
   FC_FUNC(cexit,CEXIT)(&i1,&i2);
 }

 /* This is normal exit */
  exit (ret_val);
}

/* --- function cexit: if *ps is nonzero, exit with retval pv --- */
/*cexit_(pv,ps)*/

void FC_FUNC(cexit,CEXIT)(pv,ps)
int *pv,*ps;
{

  ret_val = *pv;
#ifdef CRAY2
  exit (ret_val);
#else
  if (*ps != 0) {

/*  printf("exiting with retval = %d ...\n",ret_val); */

    /* run MPI_Finalize if hasn't been already */
#ifdef MPI
    int flag, MPI_Finalized();
    MPI_Finalized(&flag);
    if ( flag == 0) MPI_Finalize();
#endif

    exit (ret_val);
  }
#endif
}

int FC_FUNC(nargc,NARGC)()
{
  return(aargc);
}

/* A fortran-callable 'getarg'.  It requires that command-line
   arguments are passed to the CMAIN entry point (NOPASSEDARGS == 0)
*/

#if NOPASSEDARGS == 0
#define HAVE_GTARGC
#endif
#ifdef HAVE_GTARGC
void FC_FUNC(gtargc,GTARGC)(iarg,ps,len)
int *iarg; char *ps; short len;
{
  int i,maxlen; char *pps;

/* to handle fortran bug ... */
  len = (len < 0) ? -len : len;

  if (*iarg > aargc)
    { puts("getarg: request for nonexistent command line arg");
      exit(-1);
    }

/*copy string to ps, filling with blanks if passed string longer ...*/
  maxlen = strlen(*(pargv + *iarg));

/*   printf("maxlen is %d\n",maxlen); */
/*   maxlen = strlen(pargv[*iarg]); */
/*   printf("maxlen is %d\n",maxlen); */

  maxlen = (maxlen < len) ? maxlen : len;
  for (i = -1, pps=ps ; ++i<maxlen ;) *pps++ = *(*(pargv + *iarg) + i);
  while (i < len) {*pps++ = ' '; i++;}

}
#endif

/* void fmain() { void ftime(); printf("hello, world"); ftime();} */

#if DEC
void s_abort()
{
  exit(-1);
}
#endif
