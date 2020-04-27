/*  Conditional compilation of FORTRAN programs

C#ifdef, C#ifndef, C#else, C#elseif, C#endif

This preprocessor, ccomp, is intended to provide a simplified, FORTRAN
compatible version of C conditional compilation.  FORTRAN statements
beginning with C# are preprocessor directives; the ones implemented
now are C#ifdef, C#ifndef, C#else, C#elseif, C#endif.  There is also
C#define to define a name, and C#include to include another file,
as discussed below.

At any point, the preprocessor is in an UNCOMMENT (i.e. true) or in a
COMMENT (i.e. false) state.  When ccomp encounters a new C#ifdef,
C#ifndef, C#else, C#elseif, C#endif, it re-evaluates its state, which
it does in the usual way, as follows.

Directives C#ifdef, C#ifndef, C#elseif are followed by an expression,
which for now consists only of a name, eg C#ifdef CRAY.  Names are set
by a prior C#define, eg C#define CRAY.  The expression is true when
the name exists, it is false otherwise.

Supposing that the current state is true, or UNCOMMENT.
C#ifdef expr, with expr true, retains the UNCOMMENT state.
              If expr evaluates to false it switches state to COMMENT.
C#ifndef expr is like C#ifdef, but with the opposite result of expr
              NB: C#ifdef/C#endif blocks may be nested;
              C#ifdef and C#ifndef increment the nesting level
C#elseif expr switches the state to COMMENT
C#else        switches the state to COMMENT also
C#endif       decrements the nesting level, and keeps the UNCOMMENT state

Supposing that the current state is false, or COMMENT.
C#ifdef expr  keeps the COMMENT state
C#ifndef expr also keeps the COMMENT state
C#elseif expr,if expr is true, switches the state to whatever it was
              in the prior nesting level.
              If expr is false, sets the state to COMMENT
C#else        switches the state to its value in the prior
              nesting level
C#endif       decrements the nesting level, sets the state to its value
              in the prior nesting level

How ccomp determines whether to modify code:

Whether the lines between a C#ifdef/C#else/C#endif blocks need to be
commented out or uncommented depends on whether they have been
commented out in a previous pass.  This information is contained in a
'C' following the directive, eg C#ifdefC, C#ifndefC, C#elseC,
C#elseifC, C#endifC.  The preprocessor will set this, it is best
advised to create any new blocks uncommented and let the preprocessor
do the commenting.

When in a false state, ccomp adds a comment to any non-directive line
(ones beginning with something other than C#); that wasn't already
commented (it inserts a comment character in front of the line and
lengthens the line length by one character).  When in a true state, it
uncomments any line that was previouslly commented by removing the
leading comment.

C#include fname includes file fnam in the output.  Both C#ifdef and
C#include are ignored when inside an C#ifdef/C#endif block that
evaluates to COMMENT.  (Setting the -i switch tells ccomp includes
lines from a C#include directive, regardless of its state.)

As with C, ccomp distinguishes case in the names.  Output is to
standard output, unless a destination filename is supplied.

There is a primitive facility to make logical expressions using the AND
(&) and OR (|) operators, such C#ifdef john & bill, or C#ifdef john |
bill, is provided.  Precedence of operators is strictly left to right,
so that john | bill & mike is equivalent to (john | bill) & mike,
whereas john & bill | mike is equivalent to (john & bill) | mike

History:

11 Jul 1996 (version 3) C#define is ignored when inside an
                        C#ifdef/C#endif block that is false.
 1 Jul 1996 (version 2) added C#include file-name  to includes file-name
                        unless conditional compilation is false.
*/

#include <stdio.h>

#define COMMENT   1
#define UNCOMMENT 0
/*#define IBM_VM 0 */

/* switches: */
int Options;
#define FILTER 2
#define PRINTALL 4
#define INCLALL 8

/* comment character */
char cchr='C', cchr2='c';

/* Rules of operation:
  A new conditional is encountered by a C#ifdef (or C#elseif or C#else)
  Every time a new conditional is encountered, the following is determined:

  Conditional compilation is in one of two states:
    COMMENT:   should be commented
    UNCOMMENT: should be uncommented
  Previous condition is in one of two states:
    COMMENT:   was commented
    UNCOMMENT: was uncommented

  The previous condition defaults to state UNCOMMENT, unless explicitly
  specified by appending a C to the conditional, eg: C#ifdef goes to C#ifdefC,
  etc.  This will be done by ccomp if it comments something through a
  pass, and is not normally done manually.
 */

struct nest_table
{
  int num;
  char truth[10];
  char truthe[10];
};

#define LEAVE 0
#define ADDC 1
#define REMOVEC 3

char *fexts();

char s[256];
short qchange;
int lineno;
char dlist[256],ulist[256];
char *pdlist,*pulist;
short ndef=0,nudef=0;

#define INIT_PTR_LIST {pdlist=dlist; pulist=ulist;}

void movmem();

main(argc,argv)
short argc; char *argv[];
{
  char srcnam[80],dstnam[80],olddir[80],*lpdlist=NULL;
  FILE *Source; FILE *Dest;
  short iarg=0;
  void ccomp(),Error();

  if (argc == 1) Error();

/*--- set initial options ---*/
  Options = PRINTALL;
  Dest = stdout;

  INIT_PTR_LIST;

/*--- get source file name and switches --- */
  do {if (iarg >= argc) Error();} while (sw(argv[++iarg]));
  if (Options & FILTER) Source = stdin; else strcpy(srcnam,argv[iarg++]);

#ifdef MSDOS
  if (!fexts(srcnam)) strcat(srcnam,".for");
#endif
#ifdef unix
  if (!fexts(srcnam)) strcat(srcnam,".f");
#endif

/* get the destination file name, if there is one */
#ifdef unix
  if (iarg != argc) {
    getcwd(olddir,80);
    strcpy(dstnam,argv[iarg++]);
    if (iarg != argc) Error();
    if (!chdir(dstnam)) {  /* if the destination name a directory ... */
      chdir(olddir);
      strcat(dstnam,"\\");
      strcat(dstnam,srcnam);
    }
    printf("ccomp: file %s to %s\n",srcnam,dstnam);
    Dest = fopen(dstnam,"w");
  }
#endif
#ifdef IBM_VM
  strcpy(dstnam,srcnam);
  if (!fexts(srcnam)) strcat(dstnam," fortfix");
  if (!fexts(srcnam)) strcat(srcnam," fortran");
  printf("ccomp: file %s to %s\n",srcnam,dstnam);
  Dest = fopen(dstnam,"w");
#endif
#if unix | IBM_VM
#else
  if (iarg != argc) {
    strcpy(dstnam,argv[iarg++]);
    if (iarg != argc) Error();
    printf("ccomp: file %s to %s\n",srcnam,dstnam);
    Dest = fopen(dstnam,"w");
  }
#endif

/*--- open the source file ---*/
  if (!(Options & FILTER)) {
    if ((Source=fopen(srcnam,"r")) == NULL) Error();
  }

  ccomp(Source,Dest,lpdlist);
  return (0);
}

#define max(a,b)   ((a) > (b) ? (a) : (b))

void ccomp(Source,Dest,lpdlist)
FILE *Source,*Dest; char *lpdlist;
{
  struct nest_table nest;
  int i,nmov,change_state;
  char cntrl[40],name[40],*ps,*psC,*psO;
  void errchk();

  s[0] = cchr;
  nest.truthe[0] = nest.truth[0] = UNCOMMENT;
  change_state = LEAVE;
  qchange = lineno = nest.num = 0;

/*Prints out all #define to stdout and sets up lpdlist */
  if (lpdlist == NULL) {
    for (i=ndef, lpdlist=dlist; i--;) {
      fprintf(Dest,"%c#define %s\n",cchr,lpdlist);
      lpdlist += strlen(lpdlist)+1;
    }
  }

  while(ps=s+1,fgets(ps,255,Source)) {
    FILE *new_Source; int already_wrote = 0;
    lineno++;
    if ((*ps==cchr2 || *ps==cchr) && *(ps+1)=='#') {
      psC = ps; name[0] = '\0'; sscanf(ps+2,"%s %s",cntrl,name);
      already_wrote = 0;

/* ...include directive ... */
      if (!strcmp(cntrl,"include") &&
          (nest.truth[nest.num] == UNCOMMENT || Options & INCLALL)) {
        if ((new_Source=fopen(name,"r")) == NULL) {
          char sout[200];
          sprintf(sout,"ccomp failed to open file %s\n",name);
          errchk(sout);
        }
        ccomp(new_Source, Dest, lpdlist); already_wrote = 1;
      }

/* ...define directive ... */
      if (!strcmp(cntrl,"define")) {
        if (among(ulist,nudef,name)) continue;
        if (nest.truth[nest.num] == UNCOMMENT)
	  { strcat(lpdlist,name); lpdlist += strlen(lpdlist)+1; ndef++; }
      }

/* ...ifndef directive ... */
      else if (!strncmp(cntrl,"ifndef",6)) {
        psC = ps+6;
/*      Push this condition onto stack */
        nest.truthe[nest.num+1] = nest.truth[nest.num];
        nest.num++;
        nest.truth[nest.num] =
          ((nest.truthe[nest.num] == COMMENT) || match(psC)) ?
          COMMENT: UNCOMMENT;
      }

/* ...ifdef directive ... */
      else if (!strncmp(cntrl,"ifdef",5)) {
        psC = ps+5;
/*      Push this condition onto stack */
        nest.truthe[nest.num+1] = nest.truth[nest.num];
        nest.num++;
        nest.truth[nest.num] =
          ((nest.truthe[nest.num] == COMMENT) || !match(psC)) ?
          COMMENT: UNCOMMENT;
      }

/* ...elseif directive ... */
      else if (!strncmp(cntrl,"elseif",6)) {
        psC = ps+6;
        if (!nest.num) errchk("ccomp: #elseif found before #ifdef");
        if (nest.truth[nest.num] == UNCOMMENT) nest.truthe[nest.num] = COMMENT;
        nest.truth[nest.num] =
          ((nest.truthe[nest.num] == COMMENT) || !match(psC)) ?
          COMMENT: UNCOMMENT;
      }

/* ...else directive ... */
      else if (!strncmp(cntrl,"else",4)) {
        psC = ps+4;
        if (!nest.num) errchk("ccomp: #else found before #ifdef");
        if (nest.truth[nest.num] == UNCOMMENT) nest.truthe[nest.num] = COMMENT;
        nest.truth[nest.num] =
          ((nest.truthe[nest.num] == COMMENT)) ?
          COMMENT: UNCOMMENT;
      }

/* ...endif directive ... */
      else if (!strncmp(cntrl,"endif",5)) {
        psC = ps+5;
        if (!nest.num--) errchk("ccomp: #endif found before #ifdef");
      }

/*    Check for a change of state */
      if (ps != psC) {
        psC += 2;
        change_state = LEAVE;
        if (((*psC == cchr) ? COMMENT: UNCOMMENT) ^ nest.truth[nest.num]) {
          psO = ps;
          nmov = psC - ps;
/* shift right, eg C#ifdefC to _C#ifdef;  or left, eg _C#ifdef becomes C#ifdefC */
          if (nest.truth[nest.num] == COMMENT) {
            change_state = ADDC;
            ps -= 1;
          }
          else {
            change_state = REMOVEC;
            ps += 1; psC += 1;
          }
          movmem(psO,ps,nmov);
          if (nest.truth[nest.num] == COMMENT) *(psC-1) = cchr;
        }
      }
    }
    else {
      if (!nest.num && change_state != LEAVE)
        errchk("ccomp: not within conditional compilation but changing state");
      switch(change_state) {
#ifdef DEBUG
        case ADDC: ps = "add comment\r\n"; break;
        case REMOVEC: ps = "remove comment\r\n"; break;
        case LEAVE: ps = "leave as is\r\n";
#else
        case ADDC: ps -= 1; *ps = cchr; break;
        case REMOVEC: if (*ps == cchr2 || *ps == cchr) {ps += 1; break;}
                      else errchk("ccomp: Line begins without comment");
#endif
      }
    }
    if (! already_wrote &&
        (Options & PRINTALL || nest.truth[nest.num] == UNCOMMENT))
      fputs(ps,Dest);
  }
}

#include <ctype.h>

match(psC)
char *psC;
{
  char nam[64];
  int state[20],op[20],i,j,res;
  void next_word();

  next_word(&psC,nam);
  i = -1;
  do {
    i++;
    next_word(&psC,nam);
    state[i] = (nam[0] == '!') ? ! among(dlist,ndef,nam+1): among(dlist,ndef,nam);
    next_word(&psC,nam);

  } while ((op[i]=1,!strcmp(nam,"|")) || (op[i]=2,!strcmp(nam,"&")));

  for (j=0,res=state[0]; j<i; j++) {
    if (op[j] == 1) res |= state[j+1];
    else if (op[j] == 2) res &= state[j+1];
  }

/*  printf("match %d %d\n",i,res);*/

  return (res);
}

void next_word(psC,nam)
char **psC,*nam;
{
  while (isspace(**psC)) (*psC)++;
  while (!isspace(**psC)) {if (!(*nam++ = **psC)) break; (*psC)++;}
  *nam = '\0';
}
among(list,nlist,name)
char *list,*name;
{
  while (nlist-- > 0) {
    if (!strcmp(name,list)) return(1);
    list += strlen(list)+1;
  }
  return(0);
}

void errchk (string)
  char *string;
{
  if (string) {
    fprintf(stderr,"%d %s",lineno,string);
    exit (-1);
  }
}

sw(s)
char *s;
{
  void Error();
  if (*s++ != '-') return(0);
  switch(*s++) {
    case 'f': Options |= FILTER; break;
    case 'd': sscanf(s,"%s",pdlist); pdlist += strlen(pdlist)+1; ndef++;
              break;
    case 'u': sscanf(s,"%s",pulist); pulist += strlen(pulist)+1; nudef++;
              break;
    case 'c': sscanf(s,"%c",&cchr); cchr2 = cchr;
              break;
    case 'o': Options ^= PRINTALL;
              break;
    case 'i': Options ^= INCLALL;
              break;
    case 's': exit(0);
              break;
    default: Error();
  }
  return(1);
}

void Error()
{
  printf("usage: ccomp [-cdou] source-file [dest-file]\n%s%s%s%s%s%s"
        ,"              -dname defines name\n"
        ,"              -uname undefs name\n"
        ,"              -cchar changes comment from 'C' to 'char'\n"
        ,"              -o     prints out only uncommented lines\n"
        ,"              -s     exits silently with value 0\n"
        ,"              -i     always incorporate lines from C#include file\n"
        );
  printf("              version 3\n");
  exit(1);
}

void movmem(ps,pd,len)
char *ps,*pd; int len;
{
  if (ps > pd) {while(len--) *pd++ = *ps++;}
  else {ps += len; pd += len; while(len--) *--pd = *--ps;}
}

#ifdef TEST
main()
{
  char *ps,*pd;
  ps = "abcdefghi";
  pd = ps+2;
  movmem(ps,pd,3); printf("ps and pd are %s %s\n",ps,pd);
  pd = "abcdefghi";
  ps = pd+2;
  movmem(ps,pd,3); printf("ps and pd are %s %s\n",ps,pd);
}
#endif

/* fexts finds an extension to a file name.
   returns NULL if no extension was found, a pointer to extension if found
 */

char *fexts(fname)
register char *fname;
{
  register short extlen;
  register char *pnam;

#ifdef MSDOS
#define EXTCHR '.'
#define MAXEXT 4
#endif
#ifdef unix
#define EXTCHR '.'
#define MAXEXT 10
#endif
#ifdef IBM_VM
#define EXTCHR ' '
#define MAXEXT 9
#endif

  if ((pnam=fname+(extlen=strlen(fname))) == fname) return(NULL);
  extlen = (extlen > MAXEXT) ? MAXEXT: extlen;
  while(extlen--) if (*--pnam == EXTCHR) return(pnam+1);
  return(NULL);
}
