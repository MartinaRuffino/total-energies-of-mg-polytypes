#include <cublas.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DIRECT

void m2fl(int m, int n, double* a, int lda, const char *fln) {
   FILE *fl = fopen(fln, "w");
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
         fprintf(fl, " %12.4f", a[i*lda + j]);
      }
      fprintf(fl, "\n");
   }
   fclose(fl);
}

extern "C"
void pcudgemm(char transa, char transb, int m, int n, int k, double alpha,
            double *a, int lda, double *b, int ldb, double beta, double *c, int ldc,
            MPI_Comm comm) {
//    clock_t t[8];
//    double dt;

//    t[0] = clock();

   int ip, np;
   MPI_Comm_rank(comm, &ip);
   MPI_Comm_size(comm, &np);

//    printf("pcudgemm%d/%d: %c %c %d %d %d\n", ip, np, transa, transb, m, n, k);


//     cublasStatus st;

   size_t s = sizeof(double);


//    cublasSetMatrix(lc, na, s, c, lc, dc, lc);
   int ix = ip*(n/np);
   int nn;

   if (ip == np-1) {
      nn = n - (np-1)*(n/np);
   } else {
      nn = n/np;
   }

//    t[1] = clock();

   cublasDgemm (transa, transb, m, nn, k, alpha, a, lda, &b[ix*ldc], ldb, beta, &c[ix*ldc], ldc);
//    cudaDeviceSynchronize();
//    t[2] = clock();
#ifndef DIRECT
   double *hc = (double *) malloc(m*n*sizeof(double));
//    double *hc = (double *) calloc(m*n, sizeof(double));
   cublasGetMatrix(m, nn, s, &c[ix*ldc], ldc, &hc[ix*m], m);
//    cudaDeviceSynchronize();
//    m2fl(m,n,hc,m,(ip==1?"hc1-1":"hc1-0"));

   int *recvcount = (int*) malloc(np*sizeof(int));
   int *displs    = (int*) malloc(np*sizeof(int));

   for (int i = 0; i < np; ++i) recvcount[i] = (n/np)*m; // do not replace with nn because all elements will be wrong on the (np-1)th process
   recvcount[np-1] = (n - (np-1)*(n/np))*m;

   for (int i = 0; i < np; ++i) displs[i] = i*(n/np)*m;

//    t[3] = clock();

   MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hc, recvcount, displs, MPI_DOUBLE, comm);

//    t[4] = clock();
//    m2fl(m,n,hc,m,(ip==1?"hc2-1":"hc2-0"));

   free(recvcount);
   free(displs);

//    t[5] = clock();

   cublasSetMatrix(m, n, s, hc, m, c, ldc);
//    cudaDeviceSynchronize();
//    t[6] = clock();

   free(hc);

#else

   int *recvcount = (int*) malloc(np*sizeof(int));
   int *displs    = (int*) malloc(np*sizeof(int));

   for (int i = 0; i < np; ++i) recvcount[i] = (n/np)*m; // do not replace with nn because all elements will be wrong on the (np-1)th process
   recvcount[np-1] = (n - (np-1)*(n/np))*m;

   for (int i = 0; i < np; ++i) displs[i] = i*(n/np)*ldc;

//    t[3] = clock();
//    printf("DIRECT MPI_Allgatherv\n");
   MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, c, recvcount, displs, MPI_DOUBLE, comm);

//    t[4] = clock();
//    m2fl(m,n,hc,m,(ip==1?"hc2-1":"hc2-0"));

   free(recvcount);
   free(displs);

//    MPI_Status stat;
//    if (ip == 0) {
//       MPI_Send(&c[0*(n/np)*ldc], (n/np)*m, MPI_DOUBLE, 1, 0, comm);
//       MPI_Recv(&c[1*(n/np)*ldc], (n - (np-1)*(n/np))*m, MPI_DOUBLE, 1, 1, comm, &stat);
//    } else if (ip == 1) {
//       MPI_Recv(&c[0*(n/np)*ldc], (n/np)*m, MPI_DOUBLE, 0, 0, comm, &stat);
//       MPI_Send(&c[1*(n/np)*ldc], (n - (np-1)*(n/np))*m, MPI_DOUBLE, 1,  1, comm);
//    }

#endif
//    t[7] = clock();

//    exit(0);
//     dt = t[1] - t[0]; printf("pcudgemm%d prep: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[2] - t[1]; printf("pcudgemm%d calc: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[3] - t[2]; printf("pcudgemm%d pcom: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[4] - t[3]; printf("pcudgemm%d comm: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[5] - t[4]; printf("pcudgemm%d frec: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[6] - t[5]; printf("pcudgemm%d setm: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[7] - t[6]; printf("pcudgemm%d freh: %12.8f\n", ip, dt/CLOCKS_PER_SEC);
//     dt = t[7] - t[0]; printf("pcudgemm%d totl: %12.8f\n", ip, dt/CLOCKS_PER_SEC);

}