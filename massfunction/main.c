#include <stdio.h>
#include <stdlib.h>
#include "../cinclude/io_tree.h"
#ifdef MPI
#include <mpi.h>
#include "../cinclude/mympi.h"
#endif

#define STRBUFFER 1024
long long *count;

FILE *check_fopen(const char *filename, const char *mode) {
  FILE *pFile;
  pFile = fopen (filename, mode);
  if (pFile == NULL) {
    fnprintf(stderr, STRBUFFER, "ERROR: Cannot open file %s\nStop at %s l:%d of %s\n",filename,__func__,__LINE__,__FILE__);
    exit(1);
  } else
    return pFile;
}

int main(int argc, char **argv) {
  int i;
  FILE *fp;
  int masscolumn = 20;
  int numbins = 50;
  float max_v = 1.e12;
  float min_v = 1.e7;
  char buffer[STRBUFFER];
  char format[STRBUFFER];
  /* ASCII */
  strcpy(format,"");
  for(i=0;i<masscolumn;i++)
    strcat(buffer,"%g ");
  
  fp = check_fopen(filename,"r");
  while((fgets(buffer,STRBUFFER,fp)) != NULL) 
    if((buffer[0] != "#") && (buffer[0] != "\n")) {
      sscanf(buffer,)
    }
  
  
  return 0;
}
