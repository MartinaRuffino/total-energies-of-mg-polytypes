/* locase: filter that copies stdin to stdout, lowering the case */

#include <stdio.h>
#include <ctype.h>

main(argc,argv)
short argc; char *argv[];
{
  char fname[80];
  FILE *Source;
  short iarg;

/*--- get source file name and switches --- */
    maklocase(stdin);
}

#define MAXLEN 4096
maklocase(Source)
FILE *Source;
{
  char s[MAXLEN],*ps;

  while(fgets(s,MAXLEN,Source)) {
    ps = s-1; while (*++ps) if (isupper(*ps)) *ps = tolower(*ps);
  }
  fputs(s,stdout);
}
