       subroutine ebt_mult2(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &   ish,jsh,ksh,lsh,G,ng,gJ,JJJ, nbf,den)
       implicit none
c
c  this routine multiplies the integrals and density product, G,
c  and returns the contribution to the 2e part of Vab.
c
c  variables are defined in ebt_dens_mult.F
c
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer ng, nbf
       integer i,j,k,l, i1, j1, k1, l1
       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl, maxl1
       double precision gJ(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
       double precision G(ng)
       double precision JJJ(4,4,4,4)
		real(8) :: dij, dkl, dijkl, delta
		real(8) :: den(nbf,nbf,4,4)
c
c  initialize JJJ and ng
c  =========================
c
       JJJ = 0.D0
       ng  = 0
c
c  loop over bf indices
c  ===================
c
		do i = ilo, ihi
			!----------------------------------------------------------------------------
			maxj = jhi
			if(jsh==ish) maxj = i
			!----------------------------------------------------------------------------
   		do j = jlo, maxj
				!-------------------------------------------------------------------------
				if (i==j) then
					dij = 1.0d0
				else
					dij = 2.0d0
				end if
				!-------------------------------------------------------------------------
				maxk = khi
    			if(ksh==ish .and. jsh==lsh) maxk = i
				!-------------------------------------------------------------------------
    			do k = klo, maxk
					!----------------------------------------------------------------------
					maxl = lhi
					if(lsh.eq.ksh) maxl = k
					if(ish.eq.ksh.AND.jsh.eq.lsh.AND.k.eq.i) maxl=j
					!----------------------------------------------------------------------
					do l = llo, maxl
						!-------------------------------------------------------------------
						if (k==l) then
							dkl = 1.0d0
						else
							dkl = 2.0d0
						end if
						!-------------------------------------------------------------------
						ng = ng + 1
						!-------------------------------------------------------------------
						write(*,'("ttt", 2i4, 3(2x, 4i6))') ga_nodeid(), ng,
     &					ilo,ihi,jlo, jhi, klo, khi, llo, lhi, i, j, k, l
						!-------------------------------------------------------------------
						dkl = 0.0d0
						if (k.eq.l) dkl = 1.0d0
						dijkl = 0.0d0
						if (i.eq.k.AND.j.eq.l) dijkl = 1.0d0
						delta = (1.0d0+dij)*(1.0d0+dkl)
						!ng = ng + 1
						!GG(ng)  = 2.0d0*((GGijkl + GGilkj)/delta)/(1.0d0+dijkl)
						!-------------------------------------------------------------------
						do i1=4 , 1, -1
							do j1=i1, 1, -1
								do k1=i1, 1, -1
									maxl1=k1
									if (k1==i1) maxl1=j1
									do l1=maxl1, 1, -1
										JJJ(i1,k1,j1,l1)=JJJ(i1,k1,j1,l1)
     &									+ gJ(l,k,j,i)*den(l,k,k1,i1)*den(j,i,l1,j1)
								!-------------------------------------------------------------
									end do ! l1
								end do ! k1
							end do ! j1
						end do ! i1
    				enddo ! l
   			enddo ! k
			enddo ! j
		enddo ! i
c
 1111 format (2(F20.16),4x,2(F20.16))
       return
       end


c $Id: ebt_JJJ_mult.F 21176 2011-10-10 06:35:49Z d3y133 $
