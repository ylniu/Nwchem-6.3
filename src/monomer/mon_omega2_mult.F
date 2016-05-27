       subroutine mon_omega2_mult(
     & ilo,ihi,jlo,jhi,
     & klo,khi,llo,lhi,
     & ish,jsh,ksh,lsh,
     & G,ng,gJ,gK,omega2,omega_ex)
       implicit none
c
c  this routine multiplies the integrals and density product, G,
c  and returns the contribution to the 2e part of Vab.
c
c  variables are defined in mon_dens_mult.F
c
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
c
       integer ng
       integer i,j,k,l
       integer ilo, ihi, jlo, jhi
       integer klo, khi, llo, lhi
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl
       double precision gJ(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
       double precision gK(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
       double precision G(ng)
       double precision omega2, omega_ex
       double precision rtmp
c
c  initialize omega2 and ng
c  =========================
c
       omega2   = 0.d0
       omega_ex = 0.d0
       
       ng = 0
c
c  loop over bf indices
c  ===================
c
      do i = ilo, ihi
c
        maxj = jhi
        if(jsh.eq.ish) maxj = i
        do j = jlo, maxj
          do k = klo, khi
            maxl = lhi
            if(lsh.eq.ksh) maxl = k
            do l = llo, maxl
              ng = ng + 1
              if (abs(gJ(l,k,j,i)) > 1.D-40 .and.
     $            abs(G(ng))       > 1.D-40 ) then
                omega2   = omega2   + gJ(l,k,j,i)*G(ng)
!                omega_ex = omega_ex + gK(l,k,j,i)*G(ng)
                if (util_print("2e4c",print_debug)) then
                  rtmp = i*1.D6 + j*1.D4 + k*1.D2 + l
                  write(LuOut,1111) rtmp,i,j,k,l, gJ(l,k,j,i),
     $              G(ng), gJ(l,k,j,i)*G(ng)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
c
 1111 format (f11.1,4i4, 4x, 3(es20.12))
       return
       end


c $Id: mon_omega2_mult.F 21176 2011-10-10 06:35:49Z d3y133 $
