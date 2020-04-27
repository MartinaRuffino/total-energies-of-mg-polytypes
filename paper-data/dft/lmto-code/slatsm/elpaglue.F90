subroutine pdsyel2(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_real_2stage: Solves the real eigenvalue problem with a 2 stage approach
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------

!#ifdef! ELPA
!   use elpa1
!   use elpa2
!#endif

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
   real*8, intent(inout) :: a(lda,*), ev(na), q(ldq,*)

!#ifdef! ELPA
!   call solve_evp_real_2stage(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
!#else
   call rx('ELPA not available. Put -dELPA in PPFLAGS, recompile and link with the lib')
!#endif

end subroutine pdsyel2



subroutine pzheel2(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_complex_2stage: Solves the complex eigenvalue problem with a 2 stage approach
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------

!#ifdef! ELPA
!   use elpa1
!   use elpa2
!#endif

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
   complex*16, intent(inout) :: a(lda,*), q(ldq,*)
   real*8, intent(inout) :: ev(na)

!#ifdef! ELPA
!   call solve_evp_complex_2stage(na, nev, a, lda, ev, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
!#else
   call rx('ELPA not available. Put -dELPA in PPFLAGS, recompile and link with the lib')
!#endif

end subroutine pzheel2
