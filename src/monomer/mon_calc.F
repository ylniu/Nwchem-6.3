      logical function mon_calc(rtdb,Vba)
      implicit none
#include "errquit.fh"
c
c     This function calculates Vba, the electron transfer matrix element.
*
* $Id: mon_calc.F 19821 2010-12-14 07:46:49Z d3y133 $
c
*
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
#include "stdio.fh"
c
      integer rtdb, me
      integer i, j, k, l
      integer l_vals,k_vals
      integer la_mo, ka_mo, lb_mo, kb_mo
      integer g_tmp ,g_s ,g_handle(2) ,g_d ,g_temp
      integer g_atmp,g_as,g_ahandle(2),g_ad,g_atemp, ga_zai
      integer g_btmp,g_bs,g_bhandle(2),g_bd,g_btemp, gb_zai
      integer nocc(2), nocca(2), noccb(2), jroot(2)
      integer nva(2), nvb(2), iset
      integer g_u,g_vt,g_phandle(2)
      integer g_av,g_bu,g_t,g_h1
      integer g_avecsa,g_avecsb,g_bvecsa,g_bvecsb
      integer mem2, max2e, lena, lenb
      integer method_2e, anmo(2), bnmo(2), nbf_file, nset
      integer nseta, nsetb

      logical mon_movecs_read
      external mon_movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_header
      external movecs_read_header

      character*255 a_file, b_file, basis_name,
     E     a_filedir,b_filedir
      character*2 pname(2)

      double precision Vba, kappa, Hda_dipole, epsilon, Vex
      double precision Dvec(3), Avec(3), Rvec(3)
      double precision Dnorm, Anorm, Rnorm
      double precision prod(2)
      double precision Sba,omega1
      double precision detu(2),detv(2),detprod,omega2
      double precision Haa,Hbb,Hba,val, dum
      double precision Dcoc(3), Acoc(3)
      double precision Ddip(3), Adip(3), Rdip(3)
      double precision cosa , cosb, cosc
      double precision alpha, beta, theta
      data Haa,Hbb,Hba/3*0.0d0/
c     -----------------------------------------------
      integer Dgeom, Ageom, Dbas, Abas, Dnbf, Anbf, bases(2)
      integer l_occD, k_occD, l_evaD, k_evaD 
      integer l_occA, k_occA, l_evaA, k_evaA
      integer lena_cis, lenb_cis
      integer ga_vars, ga_ci, ga_cimo
      integer gb_vars, gb_ci, gb_cimo
      integer a_nci, b_nci, a_nroots, b_nroots
      character*255 a_cis_file, b_cis_file, a_cis_dir, b_cis_dir
      logical  mon_cis_read_header
      external mon_cis_read_header
      logical  mon_cis_read
      external mon_cis_read
      logical  gen_den
      external gen_den
c ------------------------------------------------------------------------------
c     Initializaion of variables: method_2e, tol2e, epsilon.
c ------------------------------------------------------------------------------
      if (.not.rtdb_get(rtdb,'mon:method_2e',MT_INT,1,method_2e)) then
          method_2e = 1
      else
          method_2e = 2
      endif
c ------------------------------------------------------------------------------
      if (.not. rtdb_get(rtdb, 'mon:tol2e', MT_DBL, 1, tol2e)) then
c     tolerance for the electron-electron interaction energy.
        tol2e = 1.d-12
      endif
c ------------------------------------------------------------------------------
      if (.not. rtdb_get(rtdb, 'mon:epsilon', MT_DBL, 1, epsilon)) then
c     dielectric constant
        epsilon = 1.d0
      endif
c ------------------------------------------------------------------------------
      me=ga_nodeid()
      jroot=1
      if (me.eq.0) then
      Call Util_Print_Centered(LuOut,'Electron Transfer Calculation'
     &     ,40,.true.)
      write(LuOut,*)
      endif
c
c -------------
c for printing:
c -------------
c
c     ! niuyingli
c       Output nothing to the mon output file
      call util_print_push
      call util_print_rtdb_load(rtdb,"mon")
c
c -------------------------------
c load geometry and symmetry info
c -------------------------------
c
      if (.not. geom_create(Dgeom, 'Dgeom'))
     $     call errquit('et_calc: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, Dgeom, 'Dgeom'))
     $     call errquit('et_calc: no Dgeom', 0, RTDB_ERR)
c ------------------------------------------------------------------------------
c     calculate charge of center of Donor for the calculation of distance between
c     donor and acceptor.
      call mon_coc(rtdb, Dgeom, 'Dgeom', Dcoc)
c -------------------------------
      if (.not. geom_create(Ageom, 'Ageom'))
     $     call errquit('et_calc: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, Ageom, 'Ageom'))
     $     call errquit('et_calc: no Ageom', 0, RTDB_ERR)
c ------------------------------------------------------------------------------
c     calculate charge of center of Acceptor for the calculation of distance between
c     donor and acceptor.
      call mon_coc(rtdb, Ageom, 'Ageom', Acoc)
c
c ----------------------------------------
c load the basis set and get info about it
c ----------------------------------------
c
      if (.not. bas_create(Dbas, 'Dbas'))
     $     call errquit('et_calc: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, Dgeom, Dbas, 'Dbas'))
     $        call errquit('et_calc: no Dbas set', 0, RTDB_ERR)
c ----------------------------------------
      if (.not. bas_create(Abas, 'Abas'))
     $     call errquit('et_calc: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, Ageom, Abas, 'Abas'))
     $        call errquit('et_calc: no Abas set', 0, RTDB_ERR)

c
c  -------
c  get nbf
c  -------
c
c     ! niuyingli
c       numbf, nbf: number of basis function
      if (.not.bas_numbf(Dbas,Dnbf))
     $     call errquit('Dbas_numbf failed',0, BASIS_ERR)
c  -------
      if (.not.bas_numbf(Abas,Anbf))
     $     call errquit('Abas_numbf failed',0, BASIS_ERR)

c
c ------------------
c initialize the API 
c ------------------
c
c     ??? niuyingli. To be analysed
      bases(1) = Dbas
      bases(2) = Abas
      call int_init(rtdb, 2, bases)
      call int_acc_std()
      call int_normalize(rtdb,Dbas)
      call int_normalize(rtdb,Abas)
c
c  -------------------------------
c  calculate the overlap integrals
c  -------------------------------
c
      if(.not.ga_create(mt_dbl,Dnbf,Dnbf,'g_as',0,0,g_as))
     $        call errquit('task_mon ga_create failed', 0, GA_ERR)
      call ga_zero(g_as)
c     ! niuyingli: calculate overlap integral of basises.
c     ! from     : src/ddscf/int_1e_ga.F
C> This routine computes 1-electron integrals of the specified kind and
C> adds them to a global array. The routine always computes a
C> rectangular matrix. It is even capable of using two different
C> basis sets. Different kinds of integrals can be computed, they are:
C>
C> * "kinetic": the kinetic energy integrals
C>
C> * "potential": the nuclear attraction
C>
C> * "overlap": the overlap integrals
      call int_1e_ga(Dbas, Dbas, g_as, 'overlap', .false.)  

      if (util_print("mon",print_debug)) then
       call ga_print(g_as)
      endif
c   ------------------------------
      if(.not.ga_create(mt_dbl,Anbf,Anbf,'g_bs',0,0,g_bs))
     $        call errquit('task_mon ga_create failed', 0, GA_ERR)
      call ga_zero(g_bs)
      call int_1e_ga(Abas, Abas, g_bs, 'overlap', .false.)  
c     ! niuyingli print global array g_bs. Source code: mon_input.F
c       could be set in mon.inp file
      if (util_print("mon",print_debug)) then
       call ga_print(g_bs)
      endif
c
c ---------------------------
c Read MO vectors from files
c ---------------------------
c
c     ! niuyingli : a_vecs -> a.mo
c                   b_vecs -> b.mo
c                   a_cis  -> a.ci_singlet
c                   b_cis  -> b.ci_singlet
c     from: src/monrans/mon_vectors_inp.F
      if (.not. rtdb_cget(rtdb, 'a_vecs', 1, a_file))
     $     call errquit('mon_calc: couldnt get a_vecs filename',0,
     &       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'b_vecs', 1, b_file))
     $     call errquit('mon_calc: couldnt get b_vecs filename',0,
     &       RTDB_ERR)
c
      lena = inp_strlen(a_file)
      lenb = inp_strlen(b_file)
c
      if (me.eq.0) then
      write(LuOut,*)'MO vectors for reactants: ',a_file(1:lena)
      write(LuOut,*)'MO vectors for products : ',b_file(1:lenb)
      write(LuOut,*)
      endif
c
c -----------------------------------------------------
c Allocate occupation and eigenvalue vectors. 
c 'occ' is used to calculate nocc, 'eval' is never used.
c ------------------------------------------------------
c
c    ! niuyingli
c      http://hpc.pnl.gov/globalarrays/ma/MA_push_get.html
      if (.not.ma_push_get(mt_dbl,Dnbf,'occa',l_occD,k_occD))
     $        call errquit('ma_push_get failed for occa',k_occD, MA_ERR)

      if (.not.ma_push_get(mt_dbl,Dnbf,'evaa',l_evaD,k_evaD))
     $        call errquit('ma_push_get failed for evaa',k_evaD, MA_ERR)
c ------------------------------------------------------
      if (.not.ma_push_get(mt_dbl,Anbf,'occb',l_occA,k_occA))
     $        call errquit('ma_push_get failed for occb',k_occA, MA_ERR)

      if (.not.ma_push_get(mt_dbl,Anbf,'evab',l_evaA,k_evaA))
     $        call errquit('ma_push_get failed for evab',k_evaA, MA_ERR)
c
c ------------------------------------------------
c Read header to check that nset=2
c ------------------------------------------------
c
c a_file:
c -------
c
      call mon_filename(a_file,a_filedir)
c    ! niuyingli src/ddscf/vectors.F:1
      if(.not.movecs_read_header(a_filedir(1:inp_strlen(a_filedir)),
     E     title,basis_name,
     &   scftype,nbf_file,nseta,anmo,2)) 
     &  call errquit('mon_calc: movecs_read_header failed',me,DISK_ERR)
c    ! niuyingli
c      Dnbf     is from mon.inp
c      nbf_file is from mo file
      if (nbf_file.ne.Dnbf) then
        call errquit('mon_calc: Dbas function mismatch',nbf_file,
     &       BASIS_ERR)
      endif

c b_file:
c -------
c
      call mon_filename(b_file,b_filedir)
      if(.not.movecs_read_header(b_filedir(1:inp_strlen(b_filedir)),
     E     title,basis_name,
     &   scftype,nbf_file,nsetb,bnmo,2)) 
     &  call errquit('mon_calc: movecs_read_header failed',me,DISK_ERR)

      if (nbf_file.ne.Anbf) then
        call errquit('mon_calc: Abas function mismatch',nbf_file,
     &       BASIS_ERR)
      endif
c -------------------------------------------------------
      if (nseta .ne. nsetb) then
        call errquit('mon_calc: nseta <> nsetb',nsetb,
     &       INPUT_ERR)
      end if
      nset=nseta
c
c -------------------------------------------------------
c allocate arrays for MO coefficients, for state A
c -------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,Dnbf,Dnbf,'MO_A_alpha',
     $ 0,0,g_avecsa))
     $ call errquit('ga_create failed for A alpha',g_avecsa,GA_ERR)
      if (nset.eq.2) then
        if (.not.ga_create(mt_dbl,Dnbf,Dnbf,'MO_A_beta',
     $ 0,0,g_avecsb))
     $ call errquit('ga_create failed for A beta',g_avecsb,GA_ERR)
      end if
      call ga_zero (g_avecsa)
      if (nset.eq.2) then
        call ga_zero (g_avecsb)
      end if
c 
c ------------------------------------------------------------------
c read occ, evals, mo coefficients for alpha orbitals from a_file
c ------------------------------------------------------------------
c
c    ! niuyingli
c      src/monrans/mon_movecs_read.F
      if(.not.mon_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     A     1,dbl_mb(k_occD)
     $   ,dbl_mb(k_evaD),g_avecsa,dum))
     $ call errquit('mon_movecs_read failed for alpha mos',0, DISK_ERR)
c
c -------------------------------------
c get nocc(1): # of occupied alpha mo's
c -------------------------------------
c
      nocc(1)=0
      do i = 1,Dnbf
       nocc(1) = nocc(1) + dbl_mb(k_occD+i-1)
      enddo
      nocca(1) = nocc(1)

      if (util_print("mon",print_high)) then
      if(me.eq.0) then
        write(LuOut,*) nocc(1), " occupied alpha orbitals"
        write(LuOut,*)
      endif
        call ga_print(g_avecsa)
      endif
c
c ------------------------------------------------------------------
c read occ, evals, mo coefficients for beta orbitals from a_file
c ------------------------------------------------------------------
c
      if (nset .eq. 2) then
        if(.not.mon_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     E     2,dbl_mb(k_occD)
     $     ,dbl_mb(k_evaD),g_avecsb, dum))
     $ call errquit('mon_movecs_read failed for beta mos',0, DISK_ERR)
      end if
c
c -------------------------------------
c get nocc(2): # of occupied alpha mo's
c -------------------------------------
c
      noccb(1) = nocca(1)
      if (nset .eq. 2) then
        nocc(2)=0
        do i = 1,Dnbf
         nocc(2) = nocc(2) + dbl_mb(k_occD+i-1)
        enddo
        nocca(2) = nocc(2)
        noccb(2) = nocca(2)

        if (util_print("mon",print_high)) then
        if(me.eq.0) then
          write(LuOut,*) nocc(2), " occupied  beta orbitals"
          write(LuOut,*)
        endif
          call ga_print(g_avecsb)
        endif
      end if
      nocca = nocca / 2
      noccb = noccb / 2
      do iset=1, nset
c     ! niuyingli
c                   total number of orbitals  - number of occupatied orbitals
        nva(iset) = anmo(iset) - nocca(iset)
        nvb(iset) = bnmo(iset) - noccb(iset)
      end do
c
c -------------------------------------------------------
c allocate arrays for MO coefficients, for state B
c -------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,Anbf,Anbf,'MO_B_alpha',
     $  0,0,g_bvecsa))
     $  call errquit('ga_create failed for B alpha',g_bvecsa,GA_ERR)

      if (nset .eq. 2) then
        if (.not.ga_create(mt_dbl,Anbf,Anbf,'MO_B_beta',
     $    0,0,g_bvecsb))
     $  call errquit('ga_create failed for B beta',g_bvecsb,GA_ERR)
      end if

      call ga_zero (g_bvecsa)
      if (nset .eq. 2) then
        call ga_zero (g_bvecsb)
      end if
c
c -----------------------------------------------------------
c read occ, evals, mo coefficients from b_file.  It is
c assumed that the number of occupied alpha and beta orbitals
c is the same as in the a_file.
c -----------------------------------------------------------
c
      if(.not.mon_movecs_read(b_filedir(1:inp_strlen(b_filedir)),
     E     1,dbl_mb(k_occA)
     $   ,dbl_mb(k_evaA) ,g_bvecsa, dum))
     $   call errquit('mon_movecs_read failed for alpha mos',0,DISK_ERR)
c
      if (nset .eq. 2) then
        if(.not.mon_movecs_read(b_filedir(1:inp_strlen(b_filedir)),
     E     2,dbl_mb(k_occA)
     $  ,dbl_mb(k_evaA) ,g_bvecsb, dum))
     $   call errquit('mon_movecs_read failed for beta mos',0,DISK_ERR)
      end if
c
      if (util_print("mon",print_high)) then
        call ga_print(g_bvecsa)
        if (nset .eq. 2) then
          call ga_print(g_bvecsb)
        end if
      endif
c ------------------------------------------------------------------------------
c
c ---------------------------
c Read CI vector filenames
c ---------------------------
c
      if (.not. rtdb_cget(rtdb, 'a_cis', 1, a_cis_file))
     $     call errquit('mon_calc: couldnt get a_cis filename',0,
     &       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'b_cis', 1, b_cis_file))
     $     call errquit('mon_calc: couldnt get b_cis filename',0,
     &       RTDB_ERR)
c
      lena_cis = inp_strlen(a_cis_file)
      lenb_cis = inp_strlen(b_cis_file)
c
      if (me.eq.0) then
      write(LuOut,*)'MO cis for reactants: ',a_cis_file(1:lena_cis)
      write(LuOut,*)'MO cis for products : ',b_cis_file(1:lenb_cis)
      write(LuOut,*)
      endif
c
c ---------------------------
c Read CI vector filename header
c ---------------------------
c
      if(.not.ga_create(mt_int,10,1,'anci',0,0,ga_vars))
     $ call errquit('task_et gci_create failed', 0, GA_ERR)
      if(.not.ga_create(mt_int,10,1,'bnci',0,0,gb_vars))
     $ call errquit('task_et gci_create failed', 0, GA_ERR)
      call ga_zero(ga_vars)
      call ga_zero(gb_vars)
c ---------------------------
      call mon_filename(a_cis_file,a_cis_dir)
c    ! niuyingli
c     src/monrans/mon_cis_read.F
      if(.not.mon_cis_read_header(a_cis_dir(1:inp_strlen(a_cis_dir)),
     & ga_vars))
     & call errquit('mon_calc: a_cis_read_header failed',me,DISK_ERR)
c ---------------------------
      call mon_filename(b_cis_file,b_cis_dir)
      if(.not.mon_cis_read_header(b_cis_dir(1:inp_strlen(b_cis_dir)),
     & gb_vars))
     & call errquit('mon_calc: b_cis_read_header failed',me,DISK_ERR)
c
c ---------------------------
c Read CI vectors
c ---------------------------
c
      call ga_get(ga_vars, 1, 1, 1, 1, a_nroots, 1)
      call ga_get(ga_vars, 4, 4, 1, 1, a_nci   , 1)
c
      call ga_get(gb_vars, 1, 1, 1, 1, b_nroots, 1)
      call ga_get(gb_vars, 4, 4, 1, 1, b_nci   , 1)
c ------------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,a_nci*2,a_nroots,'arci',0,0,ga_ci))
     $ call errquit('task_et g_a_ci_create failed', 0, GA_ERR)
c
      if(.not.ga_create(mt_dbl,b_nci*2,b_nroots,'brci',0,0,gb_ci))
     $ call errquit('task_et g_b_ci_create failed', 0, GA_ERR)
c ------------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,nva(1),nocca(1),'CI_A',0,0,ga_zai))
     $ call errquit('task_et ga_zai create failed', 0, GA_ERR)
c
      if(.not.ga_create(mt_dbl,nvb(1),noccb(1),'CI_B',0,0,gb_zai))
     $ call errquit('task_et gb_zai create failed', 0, GA_ERR)
c ------------------------------------------------------------------------------
      call ga_zero(ga_ci )
      call ga_zero(gb_ci )
      call ga_zero(ga_zai)
      call ga_zero(gb_zai)
c ------------------------------------------------------------------------------
      if(.not.mon_cis_read(a_cis_dir(1:inp_strlen(a_cis_dir)),
     & jroot(1), Dnbf, nva(1), nocca(1), ga_vars, ga_ci,
     & ga_zai,g_avecsa))
     & call errquit('mon_calc: a_cis_read_header failed',me,DISK_ERR)
c ---------------------------
      if(.not.mon_cis_read(b_cis_dir(1:inp_strlen(b_cis_dir)),
     & jroot(2), Anbf,nvb(1), noccb(1), gb_vars, gb_ci,
     & gb_zai,g_bvecsa))
     & call errquit('mon_calc: b_cis_read_header failed',me,DISK_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(ga_ci))
     &  call errquit('could not destroy ga_ci array',
     &  ga_ci, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(gb_ci))
     &  call errquit('could not destroy gb_ci array',
     &  gb_ci, GA_ERR)
c ------------------------------------------------------------------------------
c     call ga_print(g_avecsa)
      call mon_print(g_avecsa)
c     call ga_print(g_bvecsa)
      call mon_print(g_bvecsa)
c --------------------------------------------------
c     call ga_print(ga_zai)
      call mon_print(ga_zai)
c     call ga_print(gb_zai)
      call mon_print(gb_zai)
c --------------------------------------------------

c
c --------------------------------------------------
c Read converged scf or dft energy from movecs file.
c --------------------------------------------------
c
      if(.not.mon_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     $   3,dbl_mb(k_occD),dbl_mb(k_evaD),1, Haa)) 
     $   call errquit('mon_movecs_read failed for Haa',0, DISK_ERR)
c
      if(.not.mon_movecs_read(b_filedir(1:inp_strlen(b_filedir)),
     $   3,dbl_mb(k_occA),dbl_mb(k_evaA),1, Hbb))
     $   call errquit('mon_movecs_read failed for Hbb',0, DISK_ERR)
c
      if (me.eq.0) then
      write(LuOut,9995) 'Electronic energy of reactants   ','H(RR)',Haa
      write(LuOut,9995) 'Electronic energy of products    ','H(PP)',Hbb
      endif
c ------------------------------------------------------------------------------
      g_ahandle(1)=g_avecsa
      g_bhandle(1)=g_bvecsa
      if (nset .eq. 2) then
        g_ahandle(2)=g_avecsb
        g_bhandle(2)=g_bvecsb
      end if
c ------------------------------------------------------------------------------
c CI: g_ci
c MO: g_ahandle(1), g_avecsa
c ------------------------------------------------------------------------------
      if(.not.ga_create(mt_dbl,anmo(1),anmo(1),'CIMO_A',0,0,ga_cimo))
     $ call errquit('task_mon ga_cimo failed', 0, GA_ERR)
      call ga_zero(ga_cimo)
c ------------------------------------------------------------------------------
      if(.not. gen_den(nset, Dnbf, nocca, nva, anmo,
     $  ga_vars, ga_zai,g_ahandle(1), ga_cimo))
     $   call errquit('mon_movecs_read failed for alpha mos',0,DISK_ERR)
c
c ------------------------------------------------------------------------------
c
      if(.not.ga_create(mt_dbl,bnmo(1),bnmo(1),'CIMO_B',0,0,gb_cimo))
     $ call errquit('task_mon gb_cimo failed', 0, GA_ERR)
      call ga_zero(gb_cimo)
c ------------------------------------------------------------------------------
      if(.not. gen_den(nset, Anbf, noccb, nvb, bnmo,
     $  gb_vars, gb_zai,g_bhandle(1), gb_cimo))
     $   call errquit('mon_movecs_read failed for alpha mos',0,DISK_ERR)
c ------------------------------------------------------------------------------
!       call ga_print(ga_cimo)
      call mon_print(ga_cimo)
!       call ga_print(gb_cimo)
      call mon_print(gb_cimo)
c ------------------------------------------------------------------------------
      ! calculate transition dipole moment
      call mon_tran_dip(Dgeom, ga_cimo, Dbas, Dcoc, Ddip)
      call mon_tran_dip(Ageom, gb_cimo, Abas, Acoc, Adip)
      ! calculate displacement between donor and acceptor
      Rdip  = Acoc - Dcoc
      Rnorm = sqrt(Dot_product(Rdip, Rdip))
      Dnorm = sqrt(Dot_product(Ddip, Ddip))
      Anorm = sqrt(Dot_product(Adip, Adip))
      Rvec  = Rdip / Rnorm
      Dvec  = Ddip / Dnorm
      Avec  = Adip / Anorm
      !
      cosa  = Dot_product(Dvec, Rvec)
      cosb  = Dot_product(Avec, Rvec)
      cosc  = Dot_product(Dvec, Avec)
      !
      alpha = acos(cosa)
      beta  = acos(cosb)
      theta = acos(cosc)
      !
      kappa = cosc - 3.D0 * cosa * cosb
      Hda_dipole = (Dnorm * Anorm) / Rnorm**3 * kappa
c ------------------------------------------------------------------------------
      if (me==0) then
        write(LuOut,*)
        write(LuOut,'(2x,"------------------------------------------")')
        write(LuOut,'(2x,"Dipole-Dipole interaction")')
        write(LuOut,'(2x,"------------------------------------------")')
        write(LuOut,*)
        write(LuOut,'(2x,"Donor    transition dipole moment: (",
     &    f12.7 , 2(" ,",f12.7)," ),", f12.7, " a.u.")') Ddip, Dnorm
        write(LuOut,'(2x,"                                   (",
     &    f12.7 , 2(" ,",f12.7)," ),", f12.7, " Debye")') 
     &    Ddip / 0.393430307D0, Dnorm  / 0.393430307D0
        write(LuOut,'(2x,"Acceptor transition dipole moment: (",
     &    f12.7 , 2(" ,",f12.7)," ),", f12.7, " a.u.")') Adip, Anorm
        write(LuOut,'(2x,"                                   (",
     &    f12.7 , 2(" ,",f12.7)," ),", f12.7, " Debye")')
     &    Adip / 0.393430307D0, Anorm / 0.393430307D0
        write(LuOut,*)
        write(LuOut,'(2x,"Rda       ", 12x ,":", f12.7," Angstrom")')
     &    Rnorm * 0.529177D0
        write(LuOut,'(2x,"epsilon   ", 12x ,":", f12.7)') epsilon
        write(LuOut,*)
        write(LuOut,'(2x,"kappa = cos(DA) - 3 * cos(DR) * cos(AR)")')
        write(LuOut,'(2x,"kappa     ", 12x ,":", f12.7)') kappa
        write(LuOut,'(2x,"kappa^2   ", 12x ,":", f12.7)') kappa**2
        write(LuOut,*)
        write(LuOut,'(2x,"Hda = |\mu_D| * |\mu_A| / R_{DA}^3")')
        write(LuOut,'(2x,"Hda       ", 12x    ,":", f12.7," eV")')
     &    (Dnorm * Anorm) / Rnorm**3 * 27.2114d0
        write(LuOut,'(2x,"Hda * kappa  ", 9x  ,":", f12.7," eV")')
     &    (Dnorm * Anorm) / Rnorm**3 * kappa * 27.2114d0
        write(LuOut,'(2x,"Hda * kappa / epsilon :", f12.7," eV")')
     &    (Dnorm * Anorm) / Rnorm**3 * kappa/epsilon*27.2114d0
        write(LuOut,*)
      end if
c ------------------------------------------------------------------------------
c     stop
c ------------------------------------------------------------------------------
      if (method_2e.EQ.1) then
      call int_terminate()
      call mon_fock(rtdb,geom,basis,
     &  ga_cimo,gb_cimo,Dnbf,tol2e,omega2)
      endif
c ------------------------------------------------------------------------------
      if (method_2e.EQ.2) then
       call int_mem_2e4c(mem2,max2e)
       mem2 = max(mem2,max2e)
       call mon_2e(rtdb,Dnbf, Anbf, Dbas, Abas, Dgeom, Ageom,
     &   max2e,mem2,tol2e,ga_cimo,gb_cimo,omega2,Vex)
      endif
c
c ----------------------------------------------------
c destroy density arrays and calculate omega2 and Hba
c ----------------------------------------------------
c
      if (.not.ga_destroy(g_as))
     &  call errquit('could not destroy g_as array',
     &  g_as, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(g_bs))
     &  call errquit('could not destroy g_bs array',
     &  g_bs, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(g_avecsa))
     &  call errquit('could not destroy g_avecsa array',
     &  g_avecsa, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(g_bvecsa))
     &  call errquit('could not destroy g_bvecsa array',
     &  g_bvecsa, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(ga_vars))
     &  call errquit('could not destroy ga_vars array',
     &  ga_vars, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(gb_vars))
     &  call errquit('could not destroy gb_vars array',
     &  gb_vars, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(ga_zai))
     &  call errquit('could not destroy ga_zai array',
     &  ga_zai, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(gb_zai))
     &  call errquit('could not destroy gb_zai array',
     &  gb_zai, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(ga_cimo))
     &  call errquit('could not destroy ga_cimo array',
     &  ga_cimo, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(gb_cimo))
     &  call errquit('could not destroy gb_cimo array',
     &  gb_cimo, GA_ERR)
c ------------------------------------------------------------------------------
      if (nset.eq.2) then
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(g_avecsb))
     &  call errquit('could not destroy g_avecsb array',
     &  g_avecsb, GA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ga_destroy(g_bvecsb))
     &  call errquit('could not destroy g_bvecsb array',
     &  g_bvecsb, GA_ERR)
c ------------------------------------------------------------------------------
      end if
c ------------------------------------------------------------------------------
      Hba = omega2
c
       Vba =  Dabs(Hba)
c
       if (.not. rtdb_put(rtdb,'mon:V(RP)', MT_DBL, 1, Vba))
     $        call errquit('mon_calc: rtdb_put of V(RP) failed',0,
     &       RTDB_ERR)
c
c -------------
c Print results
c -------------
c
      Hda_dipole = abs(Hda_dipole)
      if(me.eq.0) then
        write(LuOut,9996)
     &    Vba,
     &    Vba * 2.194746d5,
     &    Vba * 27.2114d0,
     &    Vba*627.510d0
c
        write(LuOut,9997)
     &    Hda_dipole,
     &    Hda_dipole * 2.194746d5,
     &    Hda_dipole * 27.2114d0,
     &    Hda_dipole*627.510d0
c
        write(LuOut,8996)
     &    Vex,
     &    Vex * 2.194746d5,
     &    Vex * 27.2114d0,
     &    Vex*627.510d0
c
        write(LuOut,9998)
     &    Vba/epsilon,
     &    Vba/epsilon * 2.194746d5,
     &    Vba/epsilon * 27.2114d0,
     &    Vba/epsilon * 627.510d0
c
        write(LuOut,9998)
     &    Vex/epsilon,
     &    Vex/epsilon * 2.194746d5,
     &    Vex/epsilon * 27.2114d0,
     &    Vex/epsilon * 627.510d0
c
        write(LuOut,9999)
     &    Hda_dipole/epsilon,
     &    Hda_dipole/epsilon * 2.194746d5,
     &    Hda_dipole/epsilon * 27.2114d0,
     &    Hda_dipole/epsilon * 627.510d0
      endif
c
 9995 format (a34,1x,a6,1x,F18.10)
 9996 format( 1x,'Electron Transfer Coupling  Energy |VC(RP)|',
     &      F18.10,' Hartree',/,
     &  45x,F17.10,' cm-1',/,
     &  45x,F17.10,' eV'  ,/,
     &  45x,F17.10,' kcal/mol',/)

 9997 format( 1x,'Dipole - Dipole interaction Energy |VD(RP)|',
     &      F18.10,' Hartree',/,
     &  45x,F17.10,' cm-1',/,
     &  45x,F17.10,' eV'  ,/,
     &  45x,F17.10,' kcal/mol',/)

 8996 format( 1x,'Exchange  energy  Coupling  Energy |VX(RP)|',
     &      F18.10,' Hartree',/,
     &  45x,F17.10,' cm-1',/,
     &  45x,F17.10,' eV'  ,/,
     &  45x,F17.10,' kcal/mol',/)

 9998 format( 1x,' |VC(RP)| / epsilon                        ',
     &      F18.10,' Hartree',/,
     &  45x,F17.10,' cm-1',/,
     &  45x,F17.10,' eV'  ,/,
     &  45x,F17.10,' kcal/mol',/)

 9999 format( 1x,' |VD(RP)| / epsilon                        ',
     &      F18.10,' Hartree',/,
     &  45x,F17.10,' cm-1',/,
     &  45x,F17.10,' eV'  ,/,
     &  45x,F17.10,' kcal/mol',/)
 8998 format( 1x,' |VX(RP)| / epsilon                        ',
     &      F18.10,' Hartree',/,
     &  45x,F17.10,' cm-1',/,
     &  45x,F17.10,' eV'  ,/,
     &  45x,F17.10,' kcal/mol',/)

c

c ------------------------------------------------------------------------------
      if (.not.ma_chop_stack(l_evaA))
     &  call errquit('unable to destroy MA l_evaA', l_evaA, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ma_chop_stack(l_occA))
     &  call errquit('unable to destroy MA l_occA', l_occA, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ma_chop_stack(l_evaD))
     &  call errquit('unable to destroy MA l_evaD', l_evaD, MA_ERR)
c ------------------------------------------------------------------------------
      if (.not.ma_chop_stack(l_occD))
     &  call errquit('unable to destroy MA l_occD', l_occD, MA_ERR)
c ------------------------------------------------------------------------------

!       if (.not.ma_chop_stack(l_occA))
!      &  call errquit('unable to destroy MA l_occA',
!      &  l_occA, MA_ERR)
c ------------------------------------------------------------------------------
c ------------------------------------------------------------------------------
c ------------------------------------------------------------------------------
c
      mon_calc=.true.
c
      call util_print_pop
c
         if (.not. geom_destroy(Dgeom)) call errquit
     $        ('mon_calc: geom destroy failed', 0, GEOM_ERR)
         if (.not. geom_destroy(Ageom)) call errquit
     $        ('mon_calc: geom destroy failed', 0, GEOM_ERR)
         if (.not. bas_destroy(Dbas)) call errquit
     $        ('mon_calc: basis destroy failed',0, BASIS_ERR)
         if (.not. bas_destroy(Abas)) call errquit
     $        ('mon_calc: basis destroy failed',0, BASIS_ERR)
c
      call int_terminate()
c
c
      return
      end
c
c ------------------------------------------------------------------------------
c
      subroutine mon_filename(file,dirfile)
      implicit none
#include "inp.fh"
#include "global.fh"
      character*(*) file,dirfile
c
      character*255 dirss
      integer lenn0,lenn1,i
c

      call util_directory_name(dirss, .false., ga_nodeid())
      dirfile(1:inp_strlen(dirss))=dirss(1:inp_strlen(dirss))
      dirfile(inp_strlen(dirss)+1:inp_strlen(dirss)+1)='/'
      lenn0=inp_strlen(dirss)+2
      lenn1=inp_strlen(dirss)+inp_strlen(file)+1
      dirfile(lenn0:lenn1)=file(1:inp_strlen(file))
      do i=lenn1+1,255
      dirfile(i:i)=' '
      enddo
      
      return
      end
