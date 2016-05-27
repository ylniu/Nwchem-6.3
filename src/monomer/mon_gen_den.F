      logical function gen_den(nset, nbf, nocc, nv, nmo,
     $  g_vars, g_zai, g_mo, g_cimo)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
c
		integer nset, nbf, nocc(nset), nmo(nset)
		integer nv(nset)
		integer g_vars, g_mo, g_zai, g_cimo, unitno
		integer gamo, gimo, g_tmp
		integer nci, iroot, nroots, inci
		integer l_vecs, k_vecs
		integer iii, jjj, kkk
		integer ncoeff, idx_nc
		integer i, j, k, iset, me
c ------------------------------------------------------------------
		gen_den = .true.
		me = ga_nodeid()
!-------------------------------------------------------------------------------
		call ga_sync()
		call ga_get(g_vars, 1, 1, 1, 1, nroots, 1)
		call ga_get(g_vars, 4, 4, 1, 1, nci   , 1)
		iset = 1
!-------------------------------------------------------------------------------
		if (.not. ma_push_get(mt_dbl,nbf*nbf,'mon_movecs_read',
     $  l_vecs,k_vecs))
     $  call errquit('ma_push_get: k_vecs failed', nbf*nbf, MA_ERR)
!-------------------------------------------------------------------------------
      j=nci
		if (.not. ma_push_get(mt_dbl,j,'cis',ncoeff,idx_nc))
     $	call errquit('cis_write: cis failed', 1, MA_ERR)
!-------------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,nmo(iset),nv(iset)  ,'AMO',0,0,gamo))
     $        call errquit('task_mon ga_create failed', gamo, GA_ERR)
!-------------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,nmo(iset),nocc(iset),'IMO',0,0,gimo))
     $        call errquit('task_mon ga_create failed', gimo, GA_ERR)
!-------------------------------------------------------------------------------
		call ga_sync()
		do i = 1, nmo(iset)
			k = k_vecs + (i-1) * nmo(iset)
			call ga_get(g_mo, 1, nbf, i, i, dbl_mb(k), 1)

			if (me .eq. 0 ) then
				if (i<=nocc(iset)) then
					call ga_put(gimo, 1, nbf, i, i, dbl_mb(k), 1)
				else
					j=i-nocc(iset)
					call ga_put(gamo, 1, nbf, j, j, dbl_mb(k), 1)
				end if
			end if
			call ga_sync()
		enddo
!-------------------------------------------------------------------------------
		call mon_print(gamo)
		call mon_print(gimo)
!-------------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,nmo(iset),nocc(iset),'tmp1',0,0,g_tmp))
     $        call errquit('task_mon ga_create failed', g_tmp, GA_ERR)
      call ga_zero(g_tmp)
!-------------------------------------------------------------------------------
      call ga_dgemm('N','N',nmo(iset),nocc(iset),nv(iset),
     & 1.d0,gamo,g_zai,0.d0,g_tmp)
      if(ga_nodeid()==0) then
      write(*,'("g_tmp")')
      end if
      call mon_print(g_tmp)
!-------------------------------------------------------------------------------
c    ! niuyingli
c    g_zai : Z_{ai}
c    gamo  : C_{a \mu}
c    gimo  : C_{i \nv}
c    g_cimo: C_{a \mu}^T . Z_{ai} . C_{i \nv}
      call ga_dgemm('N','T',nmo(iset),nmo(iset),nocc(iset),
     & 1.d0,g_tmp,gimo,0.d0,g_cimo)
      call ga_sync()
!-------------------------------------------------------------------------------
      if (.not.ga_destroy(gamo))
     & call errquit('could not destroy gamo array',gamo, GA_ERR)
!-------------------------------------------------------------------------------
      if (.not.ga_destroy(gimo))
     & call errquit('could not destroy gimo array',gimo, GA_ERR)
!-------------------------------------------------------------------------------
      if (.not.ga_destroy(g_tmp))
     & call errquit('could not destroy g_tmp array',g_tmp, GA_ERR)
!-------------------------------------------------------------------------------
      return
      end
