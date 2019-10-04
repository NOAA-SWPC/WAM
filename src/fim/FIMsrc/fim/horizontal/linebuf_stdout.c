#if ( defined FORTRANUNDERSCORE )
#define linebuf_stdout linebuf_stdout_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define linebuf_stdout linebuf_stdout__
#endif

#include <stdio.h>
void linebuf_stdout ()
{
  setlinebuf (stdout);
  printf ("linebuf_stdout: output will be line buffered from now on\n");
}
