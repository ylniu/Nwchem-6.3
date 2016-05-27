      subroutine ebt_tran_dip2(geom, rho, nbf, basis, coc, mtpval)
      implicit none
      !-------------------------------------------------------------------------
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
      !-------------------------------------------------------------------------
      integer geom, basis
      double precision coc(3)
      double precision mtpval(3)
      double precision rho(nbf, nbf)
      integer nbf
      integer nshell
      integer ishell, jshell, ijshell, mpole, ncomp
      integer maxbuf, maxscr, mtpint, nbf_max
      integer l_dens, k_dens
      integer l_scr, k_scr
      integer l_buf, k_buf
      integer ilo, ihi, jlo, jhi, jdim, idim
      integer me, nproc
      integer g_rho
      integer i, j
      !-------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,nbf,nbf,'MO_A',0,0,g_rho))
     $ call errquit('task_et g_rho create failed', 0, GA_ERR)
      !-------------------------------------------------------------------------
      do i=1, nbf
        call ga_put(g_rho, 1, nbf, i, i, rho(1,i), 1)
      end do
      !-------------------------------------------------------------------------
      me = ga_nodeid()
      nproc = ga_nnodes()
      mpole = 1
      call int_init_1emtpl(maxbuf,maxscr,basis,mpole)
c
      if (.not. bas_geom(basis, geom)) call errquit
     $   ('hnd_elfcon: bad basis', 555, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $   ('hnd_elfcon: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &   ('hnd_elfcon: bas_nbf_cn_max failed',555, BASIS_ERR)
c
      !-------------------------------------------------------------------------
		! For k_dens, nbf_max(=3, for f4) is enough.
		! But nbf(=36, for f4) is also OK.
		!
      if (.not.ma_push_get(mt_dbl,nbf*nbf,'dens',l_dens,k_dens))
     &    call errquit('hnd_elfcon: ma 1 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxscr,'scratch',l_scr,k_scr))
     &    call errquit('hnd_elfcon: ma 2 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxbuf,'int buf',l_buf,k_buf))
     &    call errquit('hnd_elfcon: ma 3 failed',911,MA_ERR)
c
c     Zero mtpval result array
c
      ncomp  = (mpole+1)*(mpole+2)/2
      mtpval = 0.D0
c      call dcopy(ncomp,0.0d0,0,mtpval,1)
c
c     ----- calculate multipole integral component(s) at all points -----
c
      ijshell = 0
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)
         idim = ihi - ilo + 1

         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi)) call
     &            errquit('hnd_elfcon: bas_cn2bfr',basis,BASIS_ERR)
               jdim = jhi - jlo + 1
c
c     Get the density patch, make the integrals and contract
c
               call ga_get(g_rho, ilo, ihi, jlo, jhi,
     &              dbl_mb(k_dens), idim)
c
               call dcopy(maxbuf,0.0d0,0,dbl_mb(k_buf),1)
      !-------------------------------------------------------------------------
               call int_mpolel(basis,ishell,basis,jshell,mpole,coc,
     &              maxscr,dbl_mb(k_scr),maxbuf,dbl_mb(k_buf),mtpint)
               call ebt_multi_reduce22(dbl_mb(k_buf),dbl_mb(k_dens),
     &              mtpval,idim,jdim,ncomp, ilo, jlo)
      !-------------------------------------------------------------------------
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c     Collect components from all the nodes for all points
c
      call ga_sync()
      call ga_dgop(msg_efgs_col,mtpval,ncomp,'+')
c
c     Clean up MA data blocks
c
      if (.not.ma_pop_stack(l_buf)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_buf failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_scr)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_scr failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_dens)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_dens failed',911,MA_ERR)
      !-------------------------------------------------------------------------
      if (.not.ga_destroy(g_rho))
     &  call errquit('could not destroy g_rho array', g_rho, GA_ERR)
      !-------------------------------------------------------------------------
      return
      end subroutine
c
c ------------------------------------------------------------------------------
c
      subroutine ebt_multi_reduce22
     & (mblock,block,rblock,idim,jdim,nblock,ilo,jlo)
c
      implicit none
      integer idim,jdim, nblock
      double precision mblock(nblock,jdim,idim), block(idim,jdim)
      double precision rblock(nblock)
c
      integer iblock,i,j, ilo, jlo
c
      do i = 1, idim
         do j = 1, jdim
            do iblock = 1, nblock
               rblock(iblock) = rblock(iblock)
     &           - block(i,j) * mblock(iblock,j,i) * 2.D0
            enddo
         enddo
      enddo

      return
      end
