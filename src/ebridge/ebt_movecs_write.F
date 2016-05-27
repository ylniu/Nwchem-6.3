      logical function movecs_write(rtdb, basis, filename, 
     $     scftype, title,
     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "bq.fh"
c
c     Temporary routine
c
      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
      integer basis             ! [input] Basis handle(-1 if not accessible)
      character*(*) filename    ! [input] File to write to
      character*(*) scftype     ! [input] Type of SCF (dft, scf, mcscf)
      character*(*) title       ! [input] Title of job that created vectors
      integer nbf               ! [input] No. of functions in basis
      integer nsets             ! [input] No. of sets of vectors
      integer nmo(nsets)        ! [input] No. of vectors in each set
      integer ldocc             ! [input]
      integer ldevals           ! [input]
      double precision occ(ldocc, nsets) ! [input] Occupation numbers 
      double precision evals(ldevals, nsets) ! [input] Eigenvalues
      integer g_vecs(nsets)     ! Global array with eigen-vectors
c
      integer unitno            ! Unit no. for writing
      parameter (unitno = 67)   ! These need to be managed !!!
      integer bq_handle         ! Handle of the BQ fragment
      integer lentit
      integer lenbas
      integer l_vecs, k_vecs
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      character*26 date
      character*32 geomsum, basissum, bqsum, key
      character*20 scftype20    ! DECL MUST match movecs_converged routine
      character*128 basis_name, trans_name
      double precision energy, enrep
      integer inntsize
c
      if (oscfps) call pstat_on(ps_vecio)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
c     Generate info to go into convergence and restart line
c
      if (basis.ne.-1 .and. rtdb.ne.-1) then
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
         if (.not. bas_geom(basis, geom))
     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_write: bad geometry',0, GEOM_ERR)
         if (.not. bas_name(basis, basis_name, trans_name))
     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
         if (trans_name .ne. ' ') basis_name = trans_name
         if (.not. bq_get_active(bq_handle)) then
              bqsum = ' '
         else
           if (.not. bq_namespace(bq_handle,bqsum)) then
              bqsum = ' '
           endif
         endif
c         
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        date = ' '
      else
         basissum = ' '
         geomsum  = ' '
         bqsum    = ' '
         date     = ' '
      endif
      scftype20 = scftype
c
c
      ok = 0
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it is given a value of zero.
c
      write(key,'(a,'':energy'')') scftype(1:inp_strlen(scftype))
      if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy)) then
        energy=0.0d0
      endif
c
c get nuclear repulsion energy
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $  call errquit('unable to get nuclear rep energy',geom, GEOM_ERR)
c
      if (ga_nodeid() .eq. 0) then
         if (nsets .gt. 2) then
            write(6,*) ' movecs_write: nsets > 2 ', nsets
            goto 10
         endif
         open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c        Information about convergence
c
         write(unitno, err=1001) basissum, geomsum, bqsum, scftype20,
     $                           date
c
c        Check that read routines are both consistent with this
c
         write(unitno, err=1001) scftype20
         lentit = max(1,inp_strlen(title)) ! 0 length record confuses f2c
         write(unitno, err=1001) lentit
         write(unitno, err=1001) title(1:lentit)
         lenbas = max(1,inp_strlen(basis_name))
         write(unitno, err=1001) lenbas
         write(unitno, err=1001) basis_name(1:lenbas)
         write(unitno, err=1001) nsets
         write(unitno, err=1001) nbf
         write(unitno, err=1001) (nmo(i),i=1,nsets)
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_write',l_vecs,k_vecs))
     $        call errquit('movecs_write: ma failed', nbf, MA_ERR)
         do iset = 1, nsets
            write(unitno, err=1001) (occ(j,iset),j=1,nbf)
            write(unitno, err=1001) (evals(j,iset),j=1,nbf)
            do i = 1, nmo(iset)
               call ga_get(g_vecs(iset), 1, nbf, i, i, dbl_mb(k_vecs),1)
               call swrite(unitno, dbl_mb(k_vecs), nbf)
*               write(unitno, err=1001) (dbl_mb(k_vecs+j), j=0,nbf-1)
            enddo
         enddo
         if (.not. ma_pop_stack(l_vecs))
     $        call errquit('movecs_write: ma pop failed', l_vecs,
     &       MA_ERR)
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it is given a value of zero.
c
         write(unitno, err=1001) energy, enrep
c         
         close(unitno,err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
c
      movecs_write = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote molecular orbitals to ',a/)
         call util_flush(luout)
      endif
      call ga_sync()
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_write: failing writing to ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
