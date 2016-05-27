!        subroutine eet_dens_prod(Dnbf,Anbf,
!      &  ilo,ihi,jlo,jhi,
!      &  klo,khi,llo,lhi,
!      &  ish,jsh,ksh,lsh,
!      &  norm,GG,ng,Dblen,Ablen,
!      &  pa_ij,pa_ji,pb_lk,pb_kl)
       subroutine eet_dens_prod(Dnbf,Anbf,
     &  ilo,ihi,jlo,jhi,
     &  klo,khi,llo,lhi,
     &  ish,jsh,ksh,lsh,
     &  norm,GG,ng,
     &  numi,numj,numk,numl,
     &  pa_ij,pa_ji,pb_lk,pb_kl)
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
       integer Dnbf, Anbf,ng
       integer i,j,k,l, ii,jj,kk,ll
       integer numi,numj,numk,numl
       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl, Dblen, Ablen
       double precision Gijkl,Gjikl,Gijlk,Gjilk
       double precision GGijkl
       double precision dij,dkl,dijkl
       double precision delta
       double precision GG(10000), dabsmax, norm
       double precision pa_ij(numi,numj),pa_ji(numj,numi)
       double precision pb_lk(numl,numk),pb_kl(numk,numl)
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
            do k = klo, khi
c
              maxl = lhi
              if(lsh.eq.ksh) maxl = k
              do l = llo, maxl
c
              ii = i - ilo+1
              jj = j - jlo+1
              kk = k - klo+1
              ll = l - llo+1
c
              ng = ng + 1
              if (i==j .and. k==l) then
                Gijkl  = pa_ij(ii,jj)*pb_kl(kk,ll)
                GG(ng) = Gijkl
              else if (i==j .and. k/=l) then
                Gijkl  = pa_ij(ii,jj)*pb_kl(kk,ll)
                Gijlk  = pa_ij(ii,jj)*pb_lk(ll,kk)
                GG(ng) = Gijkl + Gijlk
              else if (i/=j .and. k==l) then
                Gijkl  = pa_ij(ii,jj)*pb_kl(kk,ll)
                Gjikl  = pa_ji(jj,ii)*pb_kl(kk,ll)
                GG(ng) = Gijkl + Gjikl
              else
                Gijkl  = pa_ij(ii,jj)*pb_kl(kk,ll)
                Gjikl  = pa_ji(jj,ii)*pb_kl(kk,ll)
                Gijlk  = pa_ij(ii,jj)*pb_lk(ll,kk)
                Gjilk  = pa_ji(jj,ii)*pb_lk(ll,kk)
                GG(ng) = Gijkl + Gjikl + Gijlk + Gjilk
              end if
c
               enddo
            enddo
         enddo
      enddo
      norm = dabsmax(ng,GG)
      return
      end


c $Id: eet_dens_prod.F 21176 2011-10-10 06:35:49Z d3y133 $
