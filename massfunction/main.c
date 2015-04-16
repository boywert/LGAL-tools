#include <stdio.h>
#include <stdlib.h>
#include "../cinclude/io_tree.h"
#include "../cinclude/mympi.h"
#ifdef MPI
#include <mpi.h>
#endif

int main(int argc, char **argv) {
  int i;
  struct mpi_vars mympi;
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mympi.ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &mympi.NTask);
#else
  mympi.NTask = 1;
  mympi.ThisTask = 0;
#endif //MPI
  if(!mympi.ThisTask) {
  }
  return 0;
}
