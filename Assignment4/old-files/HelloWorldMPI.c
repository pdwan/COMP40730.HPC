 #include <stdio.h>
 #include <stdlib.h>
 #include <mpi.h>

 int main( int argc, char** argv )
 {
    int myId, numProc, len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Get_processor_name(name, &len);
    printf("Id: %d, Name: %s, NumProcs: %d\n", myId, name, numProc);
    MPI_Finalize();
    return 0;
 }
