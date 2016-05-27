c     $Id: eet_2e.F 19708 2010-10-29 18:04:21Z d3y133 $
      subroutine eet_2e(rtdb,Dnbf, Anbf, Dbas, Abas, Dgeom, Ageom,
     & max2e,mem2,tol2e,g_pa,g_pb,omega2, Vex)
      implicit none
#include "errquit.fh"
c
c     calculates the two-electron contribution to Vba, omega2
c
c     nbf, basis, geom                                          [input]
c     max2e, mem2 (for int_2e4c)                                [input]
c     g_pa, g_pb (handles for alpha and beta densities)         [input]
c     omega2 (two electron contrib to Vab)                      [output]
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
      integer lpa_ij, kpa_ij, lpa_ji, kpa_ji
      integer lpb_kl, kpb_kl, lpb_lk, kpb_lk
      integer rtdb
      integer g_pa,g_pb
      integer max2e, mem2, me
      integer Dnbf, Anbf, Dbas, Abas, Dgeom, Ageom
      integer l_gJ, k_gJ, l_scr, k_scr, l_den, k_den
      integer l_gK, k_gK, l_scK, k_scK
      integer ish, jsh, ksh, lsh, nsh, Dnsh, Ansh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer nproc, Dblen, Ablen, Dmax_sh_bf, Amax_sh_bf
      integer maxlsh, shellnum, ng
      integer idim, jdim, kdim, ldim
      double precision omega2, omega2_val, Vex, omega_ex
      integer nxtask, next
      integer i, j, k, l
      external nxtask
      double precision sij, skl, Dsmax, Asmax
      double precision norm, G(10000), dij
      double precision tol2e
      double precision dabsmax, adenmax,bdenmax
      external  dabsmax
c     -------------------------------------------     
c     find the maximum value of the alpha density
c     -------------------------------------------     
c
      me = ga_nodeid()
      call ga_distribution(g_pa,me,ilo,ihi,jlo,jhi)
      idim=ihi-ilo+1
      jdim=jhi-jlo+1
c
      if (.not.ma_push_get(MT_DBL,idim*jdim,'temp dens',l_den,k_den))
     $  call errquit('eet_2e:ma_push_get failed for temp dens',k_den,
     &       MA_ERR)
      call dfill(idim*jdim, 0.0d0, dbl_mb(k_den), 1)
c
      call ga_get(g_pa, ilo, ihi, jlo, jhi, dbl_mb(k_den),1)
      adenmax = dabsmax(idim*jdim,dbl_mb(k_den))
c
      call ga_sync()
c
c     find the maximum of CIMO_D
      call ga_dgop(msg_eet_adens,adenmax,1,'max')
c
c     -------------------------------------------     
c     find the maximum value of the beta density
c     -------------------------------------------     
c
      call ga_distribution(g_pb,ga_nodeid(),ilo,ihi,jlo,jhi)   
      idim=ihi-ilo+1
      jdim=jhi-jlo+1
c
      call dfill(idim*jdim, 0.0d0, dbl_mb(k_den), 1)
c
      call ga_get(g_pb, ilo, ihi, jlo, jhi, dbl_mb(k_den),1)
      bdenmax = dabsmax(idim*jdim,dbl_mb(k_den))
c
      if (.not. ma_pop_stack(l_den))
     &  call errquit('eet_2e: pop l_den', l_den, MA_ERR)
c
      call ga_sync()
c
c     find the maximum of CIMO_A
      call ga_dgop(msg_eet_bdens,bdenmax,1,'max')
c
c     -------------------------------------------------
c     initialize omega2 and establish 2e and scr arrays
c     -------------------------------------------------
c
      omega2 = 0.d0
      Vex    = 0.d0
c
      if (.not.ma_push_get(MT_DBL, max2e, '2e buffer', l_gJ, k_gJ))
     $   call errquit('eet_2e:ma_push_get failed for 2e J buffer',k_gJ,
     &       MA_ERR)
c
      if (.not.ma_push_get(MT_DBL, mem2, '2e scr', l_scr, k_scr))
     $   call errquit('eet_2e:ma_push_get failed for 2e scratch',k_scr,
     &       MA_ERR)
c
      if (.not.ma_push_get(MT_DBL, max2e, '2e buffer', l_gK, k_gK))
     $   call errquit('eet_2e:ma_push_get failed for 2e J buffer',k_gK,
     &       MA_ERR)
c
      if (.not.ma_push_get(MT_DBL, mem2, '2e scr', l_scK, k_scK))
     $   call errquit('eet_2e:ma_push_get failed for 2e scratch',k_scK,
     &       MA_ERR)
c
c     ------------------
c     get basis set info
c     ------------------
c
      if (.not. bas_numcont(Dbas, Dnsh)) call errquit
     $     ('eet_2e: Dbas_numcont failed', Dbas, BASIS_ERR)
c
      if (.not. bas_nbf_cn_max(Dbas,Dmax_sh_bf))
     $     call errquit('eet_2e: could not get Dmax_sh_bf',0, BASIS_ERR)

      Dblen = min(Dnbf,36,6*Dmax_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
c
c     ------------------
c     get basis set info
c     ------------------
c
      if (.not. bas_numcont(Abas, Ansh)) call errquit
     $     ('eet_2e: Abas_numcont failed', Abas, BASIS_ERR)
c
      if (.not. bas_nbf_cn_max(Abas,Amax_sh_bf))
     $     call errquit('eet_2e: could not get Amax_sh_bf',0, BASIS_ERR)

      Ablen = min(Anbf,36,6*Amax_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
c
c     ----------------------------------
c     set and print screening parameters
c     ----------------------------------
c
      call schwarz_init(Dgeom, Dbas)
      !Dsmax= schwarz_max()
      call schwarz_init(Ageom, Abas)
      !Asmax= schwarz_max()
      dij = max(adenmax*adenmax,bdenmax*bdenmax)
c
 1212 format (1x,a30,1x,a1,1pE14.4)
      if (ga_nodeid().eq.0) then
      write(LuOut,*)
      write(LuOut,*)'Start two-electron contribution,  H2(RP)...'
      !write(LuOut,1212)'Density screening tolerance   ',':',dij
      !write(LuOut,1212)'Largest two-electron integral ',':',Dsmax
      !write(LuOut,1212)'Largest two-electron integral ',':',Asmax
      write(LuOut,1212)'Integral screening threshold  ',':',tol2e
      write(LuOut,*)
      endif
c
c     ---------------------------------
c     allocate temporary density arrays
c     ---------------------------------
c
      if (.not. ma_push_get(mt_dbl,Dblen*Dblen,'pa_ij', lpa_ij, kpa_ij))
     $     call errquit('eet_2e:could not allocate lpa_ij',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,Dblen*Dblen,'pa_ji', lpa_ji, kpa_ji))
     $     call errquit('eet_2e:could not allocate lpa_ji',0, MA_ERR)
c
      if (.not. ma_push_get(mt_dbl,Ablen*Ablen,'pb_kl', lpb_kl, kpb_kl))
     $     call errquit('eet_2e:could not allocate lpb_kl',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl,Ablen*Ablen,'pb_lk', lpb_lk, kpb_lk))
     $     call errquit('eet_2e:could not allocate lpb_lk',0, MA_ERR)
c
c     ----------------------     
c     begin loop over shells
c     ----------------------     
c     
      call ga_sync()
c ------------------------------------------------------------------------------
      if (util_print("2e4c",print_debug) .and. me==0) then
        write(LuOut,1515)
        write(LuOut,1616) "i(D)", "j(D)", "k(A)", "l(A)",
     $    "Int(2e4c)", "CIMO(D)*CIMO(A)", "Int*CIMO"
      end if
 1515 format (" Two electron integral output")
 1616 format (4a6, 3(2x,a20))
c ------------------------------------------------------------------------------
      shellnum = 0
      nproc = ga_nnodes()
      next = nxtask(nproc,1)
c
      do ish = Dnsh, 1, -1
        if (.not. bas_cn2bfr(Dbas, ish, ilo, ihi))
     $    call errquit('eet_2e: bas_cn2bfr', ish, BASIS_ERR)
c
        do jsh = ish, 1, -1
          if (.not. bas_cn2bfr(Dbas, jsh, jlo, jhi))
     $      call errquit('eet_2e: bas_cn2bfr', jsh, BASIS_ERR)
c
c         ???
          sij = schwarz_shell(ish,jsh)
          if (sij .gt. tol2e) then
c
            do ksh = Ansh, 1, -1
              if (.not. bas_cn2bfr(Abas, ksh, klo, khi))
     $          call errquit('eet_2e: bas_cn2bfr', ksh,BASIS_ERR)
c
              do lsh = ksh, 1, -1
                if (.not. bas_cn2bfr(Abas, lsh, llo, lhi))
     $            call errquit('eet_2e: bas_cn2bfr', lsh,BASIS_ERR)
c
                if (shellnum.eq.next) then
c
                  skl = schwarz_shell(ksh,lsh)
c
                  if (sij*skl.gt.tol2e) then
                    idim=ihi-ilo+1
                    jdim=jhi-jlo+1
                    kdim=khi-klo+1
                    ldim=lhi-llo+1
c
                    call ga_get(g_pa,ilo,ihi,jlo,jhi,
     $                dbl_mb(kpa_ij),idim)
                    call ga_get(g_pa,jlo,jhi,ilo,ihi,
     $                dbl_mb(kpa_ji),jdim)
                    call ga_get(g_pb,klo,khi,llo,lhi,
     $                dbl_mb(kpb_kl),kdim)
                    call ga_get(g_pb,llo,lhi,klo,khi,
     $                dbl_mb(kpb_lk),ldim)
c
                    call eet_dens_prod(Dnbf,Anbf,
     $                ilo,ihi,jlo,jhi,
     $                klo,khi,llo,lhi,
     $                ish,jsh,ksh,lsh,
     $                norm,G,ng,
     $                idim,jdim,kdim,ldim,
     $                dbl_mb(kpa_ij),dbl_mb(kpa_ji),
     $                dbl_mb(kpb_lk),dbl_mb(kpb_kl))
c
                    call eet_int_2e4c(Dbas,ish,jsh,Abas,ksh,lsh,
     $                mem2, dbl_mb(k_scr), max2e, dbl_mb(k_gJ))
c
!                     call eet_int_2e4c_ex(Dbas,ksh,jsh,Abas,ish,lsh,
!      $                mem2, dbl_mb(k_scK), max2e, dbl_mb(k_gK))
c
                    call eet_omega2_mult(
     &                ilo,ihi,jlo,jhi,
     &                klo,khi,llo,lhi,
     &                ish,jsh,ksh,lsh,
     &                G,ng,
     &                dbl_mb(k_gJ),
     &                dbl_mb(k_gK),
     &                omega2_val,
     &                omega_ex)
c
                    omega2 = omega2 + omega2_val * 4.D0
                    Vex    = Vex    + omega_ex   * 2.D0
c
                  endif  !norm screening !sij*skl.gt.tol2e
                  next = nxtask(nproc,1)
                endif   !shellnum
                shellnum = shellnum + 1
              enddo
            enddo
          endif  !sij x smax ! sij .gt. tol2e
        enddo
      enddo
c
      next = nxtask(-nproc,1)
      call ga_sync()
c
      call ga_dgop(msg_eet_2e,omega2,1,'+')
c
c     destroy arrays
c     ==============
      call schwarz_tidy()
c
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(lpb_lk))
     &  call errquit('eet_2e: pop lpb_lk', lpb_lk, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(lpb_kl))
     &  call errquit('eet_2e: pop lpb_kl', lpb_kl, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(lpa_ji))
     &  call errquit('eet_2e: pop lpa_ji', lpa_ji, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(lpa_ij))
     &  call errquit('eet_2e: pop lpa_ij', lpa_ij, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(l_scK))
     &  call errquit('eet_2e: pop l_scr', l_scr, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(l_gK))
     &  call errquit('eet_2e: pop l_gJ', l_gJ, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(l_scr))
     &  call errquit('eet_2e: pop l_scr', l_scr, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not. ma_pop_stack(l_gJ))
     &  call errquit('eet_2e: pop l_gJ', l_gJ, MA_ERR)
c ------------------------------------------------------------------------------
c
      return
      end