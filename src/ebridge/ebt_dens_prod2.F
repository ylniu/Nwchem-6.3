		subroutine ebt_dens_prod2(nbf,ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &	ish,jsh,ksh,lsh,i1,j1,k1,l1,norm,GG,ng,den)
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
       integer i , j , k , l
       integer i1, j1, k1, l1
       integer ilo, ihi, jlo, jhi
       integer klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl, blen
       double precision Gijkl,Gjikl,Gijlk,Gjilk
       double precision GGijkl
       double precision dij,dkl,dijkl
       double precision delta
       double precision GG(1000), dabsmax, norm
       double precision den(nbf, nbf, 4, 4)
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
                Gijkl =   den(i,j,i1,j1)*den(k,l,k1,l1)
     &                  + den(i,j,k1,l1)*den(k,l,i1,j1)
c
                Gjikl =   den(j,i,i1,j1)*den(k,l,k1,l1)
     &                  + den(j,i,k1,l1)*den(k,l,i1,j1)
c
                Gijlk =   den(i,j,i1,j1)*den(l,k,k1,l1)
     &                  + den(i,j,k1,l1)*den(l,k,i1,j1)
c
                Gjilk =   den(j,i,i1,j1)*den(l,k,k1,l1)
     &                  + den(j,i,k1,l1)*den(l,k,i1,j1)
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
               enddo
            enddo
         enddo
      enddo
c
      norm = dabsmax(ng,GG)
c
      return
      end

c $Id: ebt_dens_prod.F 21176 2011-10-10 06:35:49Z d3y133 $
