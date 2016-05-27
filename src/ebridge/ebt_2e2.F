c     $Id: ebt_2e.F 19708 2010-10-29 18:04:21Z d3y133 $
		subroutine ebt_2e2(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &	g_den, den, Jab)
      implicit none
#include "errquit.fh"
c
c     calculates the two-electron contribution to Vba, Jab
c
c     nbf, basis, geom                                          [input]
c     max2e, mem2 (for int_2e4c)                                [input]
c     g_den (handles for alpha and beta densities)              [input]
c     Jab (two electron contrib to Vab)                         [output]
c
#include "rtdb.fh"
#include "schwarz.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "msgids.fh"   
#include "stdio.fh"   
c
      integer i,j,k,l, i1, j1, k1, l1, maxl1
      integer rtdb, me
      integer g_den(4,4)
      integer nbf, geom, basis,max2e, mem2
      integer l_gJ, k_gJ, l_scr, k_scr, l_den, k_den
      integer ish, jsh, ksh, lsh, nsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer nproc, max_sh_bf, blen
      integer maxlsh, shellnum, ng, numi, numj
      double precision Jabijkl
      integer nxtask, next
      external nxtask
      double precision smax, sij, skl
      double precision norm, G(10000), dij
      double precision tol2e
      double precision dabsmax, adenmax, bdenmax
		real(8) :: Jab(4,4,4,4)
		real(8) :: den(nbf,nbf,4,4)
      external  dabsmax
c     -------------------------------------------     
c     find the maximum value of the alpha density
c     -------------------------------------------     
c
		me      = ga_nodeid()
		Jab     = 0.D0
		bdenmax = 0.D0
		!-----------------------------------------------------------------------------
		call ga_distribution(g_den(1,1),me,ilo,ihi,jlo,jhi)
		!-----------------------------------------------------------------------------
		do i=1, 4
			do j=1, i
  				numi=ihi-ilo+1
   			numj=jhi-jlo+1
				if (.not.ma_push_get(MT_DBL,numi*numj,'tdens',l_den,k_den))
     $		call errquit('ebt_2e:ma_push_get failed for tdens',k_den,MA_ERR)
				call dfill(numi*numj, 0.0d0, dbl_mb(k_den), 1)
  				call ga_get(g_den(i,j), ilo, ihi, jlo, jhi, dbl_mb(k_den),numi)
				adenmax = dabsmax(numi*numj,dbl_mb(k_den))
				call ga_sync()
				call ga_dgop(msg_et_adens,adenmax,1,'max')
				bdenmax = max(adenmax, bdenmax)
			end do
		end do
		adenmax = bdenmax
c
c     -------------------------------------------------
c     initialize establish 2e and scr arrays
c     -------------------------------------------------
c     
      if (.not.ma_push_get(MT_DBL, max2e, '2e buffer', l_gJ, k_gJ))
     $   call errquit('ebt_2e:ma_push_get failed for 2e J buffer',k_gJ,
     &     MA_ERR)
c     
      if (.not.ma_push_get(MT_DBL, mem2, '2e scr', l_scr, k_scr))
     $   call errquit('ebt_2e:ma_push_get failed for 2e scratch',k_scr,
     &     MA_ERR)
c     
c     ------------------
c     get basis set info
c     ------------------
c     
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('ebt_2e: bas_numcont failed', basis, BASIS_ERR)
c     
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('ebt_2e: could not get max_sh_bf',0, BASIS_ERR)

      blen = min(nbf,36,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
		!-----------------------------------------------------------------
		! set and print screening parameters
		!-----------------------------------------------------------------
		call schwarz_init(geom, basis)
		smax = schwarz_max() 
		dij  = adenmax*adenmax
		!-----------------------------------------------------------------
 1212 format (1x,a30,1x,a1,1pE14.4)
		if (me==0) then
			write(LuOut,*)
			write(LuOut,*)'Start two-electron contribution '
			write(LuOut,1212)'Largest two-electron integral ',':',smax
			write(LuOut,1212)'Integral screening threshold  ',':',tol2e
			write(LuOut,*)
		endif
		!-----------------------------------------------------------------
		! begin loop over shells
		!-----------------------
		call ga_sync()
		shellnum = 0
		nproc    = ga_nnodes()
		next     = nxtask(nproc,1)
		!-----------------------------------------------------------------
		do ish = nsh, 1, -1
		!-----------------------------------------------------------------
			if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $		call errquit('ebt_2e: bas_cn2bfr', ish, BASIS_ERR)
			!--------------------------------------------------------------
			do jsh = ish, 1, -1
			!--------------------------------------------------------------
				if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $			call errquit('ebt_2e: bas_cn2bfr', jsh, BASIS_ERR)
				sij = schwarz_shell(ish,jsh)
				if (sij*smax*dij .gt. tol2e) then
					!--------------------------------------------------------
					do ksh = ish, 1, -1
					!--------------------------------------------------------
						if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $					call errquit('ebt_2e: bas_cn2bfr', ksh, BASIS_ERR)
						maxlsh = ksh
						if (ksh==ish) maxlsh = jsh
						!-----------------------------------------------------
						do lsh = maxlsh, 1, -1
						!-----------------------------------------------------
							if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $					call errquit('ebt_2e: bas_cn2bfr', lsh, BASIS_ERR)
							!--------------------------------------------------
 							if (shellnum==next) then
								!-----------------------------------------------
								skl = schwarz_shell(ksh,lsh)      
								if (sij*skl*smax > tol2e) then
									!--------------------------------------------
 									call int_2e4c(basis, ish, jsh, basis, 
     $								ksh, lsh, mem2, dbl_mb(k_scr), max2e, 
     $								dbl_mb(k_gJ))
									!
									!--------------------------------------------
									!
									do i1=4 , 1, -1
									do j1=i1, 1, -1
									do k1=i1, 1, -1
										maxl1=k1
										if (k1==i1) maxl1=j1
									do l1=maxl1, 1, -1
			!--------------------------------------------------------------
			call ebt_dens_prod2(nbf,ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     $		ish,jsh,ksh,lsh,i1,j1,k1,l1,norm, G, ng, den)
			!--------------------------------------------------------------
			if (sij*skl*smax*norm > tol2e) then
				!-----------------------------------------------------------
				call ebt_mult(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &			ish,jsh,ksh,lsh,G,ng,dbl_mb(k_gJ),Jabijkl)
				!-----------------------------------------------------------
				Jab(i1,j1,k1,l1) = Jab(i1,j1,k1,l1) + Jabijkl
				!-----------------------------------------------------------
			end if
			!--------------------------------------------------------------
									end do ! l1
									end do ! k1
									end do ! j1
									end do ! i1
									!
									!--------------------------------------------
									!
  								endif  !norm screening
								next = nxtask(nproc,1)
								!-----------------------------------------------
   						endif !shellnum
							shellnum = shellnum + 1
						enddo ! lsh
					enddo ! ksh
				endif  !sij x smax
			enddo ! jsh
		enddo ! ish
		!-----------------------------------------------------------------
		next = nxtask(-nproc,1)
		call ga_sync()
		!-----------------------------------------------------------------
		i=4**4
		call ga_dgop(msg_et_2e,Jab,i,'+')
		!-----------------------------------------------------------------
		! destroy arrays
		!-----------------------------------------------------------------
		call schwarz_tidy()
		!-----------------------------------------------------------------
		return
		end
