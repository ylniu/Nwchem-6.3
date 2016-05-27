		subroutine ebt_dens_prod(nbf,ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &	ish,jsh,ksh,lsh, norm,GG,ng,blen,
     &	pa_ij,pa_kl,pa_ji,pa_lk,pb_ij,pb_kl,pb_ji,pb_lk) 
		implicit none
c
c  this routine multiplies the densities and forms the array GG
c  (a product of densities) and calculates the maximum value, which
c  is used as a screening parameter.
c
c      ilo,ihi,jlo,jhi,klo,khi,llo,lhi - bf index ranges [input]
c      ish,jsh,ksh,lsh - shell indixes                   [input]
c      pa,pb - alpha and beta densities                  [input]
c      norm                                              [output]
c      GG - densit products                              [output]
c      ng - number of elements in GG                     [output]
c 
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer nbf,ng
       integer i,j,k,l, ii,jj,kk,ll
       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl, blen
       double precision Gijkl,Gjikl,Gijlk,Gjilk
       double precision GGijkl
       double precision dij,dkl,dijkl
       double precision delta
       double precision GG(10000), dabsmax, norm
       double precision pa_ij(blen,blen),pa_kl(blen,blen)
       double precision pa_ji(blen,blen),pa_lk(blen,blen)
       double precision pb_ij(blen,blen),pb_kl(blen,blen) 
       double precision pb_ji(blen,blen),pb_lk(blen,blen)
       external dabsmax
c
c
c  initialize ng
c  =============
c
       ng = 0
c
c  loop over bf indices
c  ====================
c
      do i = ilo, ihi
c
          maxj = jhi
          if(jsh.eq.ish) maxj = i
          do j = jlo, maxj
c
            maxk = khi
            if(ksh.eq.ish.AND.jsh.eq.lsh) maxk = i
            do k = klo, maxk
c
              maxl = lhi
              if(lsh.eq.ksh) maxl = k
              if(ish.eq.ksh.AND.jsh.eq.lsh.AND.k.eq.i) maxl=j
              do l = llo, maxl
c
              ii = i - ilo+1
              jj = j - jlo+1
              kk = k - klo+1
              ll = l - llo+1
c
                Gijkl =   pa_ij(ii,jj)*pb_kl(kk,ll)
     &                  + pb_ij(ii,jj)*pa_kl(kk,ll)
c
                Gjikl =   pa_ji(jj,ii)*pb_kl(kk,ll)
     &                  + pb_ji(jj,ii)*pa_kl(kk,ll)
c
                Gijlk =   pa_ij(ii,jj)*pb_lk(ll,kk)
     &                  + pb_ij(ii,jj)*pa_lk(ll,kk)
c
                Gjilk =   pa_ji(jj,ii)*pb_lk(ll,kk)
     &                  + pb_ji(jj,ii)*pa_lk(ll,kk)
c
                if (i==j) then
                   dij = 2.0d0
                else
                   dij = 1.0d0
                end if
c
                if (k==l) then
                   dkl = 2.0d0
                else
                   dkl = 1.0d0
                end if
c
                if (i==k .and. j==l) then
                   dijkl = 2.0d0
                else
                   dijkl = 1.0d0
                end if
c
                delta = dij*dkl/dijkl
c
                GGijkl = Gijkl + Gjikl + Gijlk + Gjilk
c
                ng = ng + 1
                GG(ng) = GGijkl / delta
c
               enddo
            enddo
         enddo
      enddo
c
      norm = dabsmax(ng,GG)
c
      return
      end


c $Id: et_dens_prod.F 21176 2011-10-10 06:35:49Z d3y133 $
