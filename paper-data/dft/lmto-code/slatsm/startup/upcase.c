/* upcase: filter that copies stdin to stdout, raising the case */

#include <stdio.h>
#include <ctype.h>

main(argc,argv)
short argc; char *argv[];
{
  char fname[80];
  FILE *Source;
  short iarg;

/*--- get source file name and switches --- */
    makupcase(stdin);
}

#define MAXLEN 4096
makupcase(Source)
FILE *Source;
{
  char s[MAXLEN],*ps;

  while(fgets(s,MAXLEN,Source)) {
    ps = s-1; while (*++ps) if (islower(*ps)) *ps = toupper(*ps);
  }
  fputs(s,stdout);
}
