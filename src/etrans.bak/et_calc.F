      logical function et_calc(rtdb,Vba)
      implicit none
#include "errquit.fh"
c
c     This function calculates Vba, the electron transfer matrix element.
*
* $Id: et_calc.F 19821 2010-12-14 07:46:49Z d3y133 $
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
      integer i, j, k
      integer l_vals,k_vals
      integer g_tmp,g_s,g_ahandle(2),g_bhandle(2),g_d,g_temp
      integer nocc(2)
      integer g_u,g_vt,g_phandle(2)
      integer g_av,g_bu,g_t,g_h1
      integer g_avecsa,g_avecsb,g_bvecsa,g_bvecsb
      integer mem2, max2e, lena, lenb, lena_cis, lenb_cis
      integer method_2e, anmo(2), bnmo(2), nbf_file, nset

      logical et_movecs_read
      external et_movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_header
      external movecs_read_header

      character*255 a_file, b_file, basis_name, a_cis_file, b_cis_file,
     E     a_filedir,b_filedir, a_cis_dir, b_cis_dir
      character*2 pname(2)

      double precision Vba
      double precision prod(2)
      double precision Sba,omega1
      double precision detu(2),detv(2),detprod,omega2
      double precision Haa,Hbb,Hba,val, dum
!     Begin adding by niuyingli
		integer Dgeom, Ageom, Dbas, Abas, Dnbf, Anbf, bases(2)
      integer l_occD, k_occD, l_occA, k_occA
      integer l_evalD,k_evalD, l_evalA,k_evalA
      integer l_vecsDi, k_vecsDi, l_vecsAk, k_vecsAk
      integer l_vecsDj, k_vecsDj, l_vecsAl, k_vecsAl
		integer ga_nci, ga_rci 
		integer gb_nci, gb_rci 
		integer a_nroots, a_ipol, a_nocc, a_nmo, a_nfc, a_nfv
		integer b_nroots, b_ipol, b_nocc, b_nmo, b_nfc, b_nfv
		integer a_nci, a_ncoeff, a_i_nc
		integer b_nci, b_ncoeff, b_i_nc
		integer iroot, inci, iii, jjj, kkk
!     End   adding by niuyingli
      data Haa,Hbb,Hba/3*0.0d0/

      me=ga_nodeid()

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
      call util_print_push
      call util_print_rtdb_load(rtdb,"et")
c
c -------------------------------
c load geometry and symmetry info
c -------------------------------
c
      if (.not. geom_create(Dgeom, 'Dgeom'))
     $     call errquit('et_calc: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, Dgeom, 'Dgeom'))
     $     call errquit('et_calc: no Dgeom', 0, RTDB_ERR)
!
      if (.not. geom_create(Ageom, 'Ageom'))
     $     call errquit('et_calc: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, Ageom, 'Ageom'))
     $     call errquit('et_calc: no Ageom', 0, RTDB_ERR)
c
c ----------------------------------------
c load the basis set and get info about it
c ----------------------------------------
c
      if (.not. bas_create(Dbas, 'Dbas'))
     $     call errquit('et_calc: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, Dgeom, Dbas, 'Dbas'))
     $        call errquit('et_calc: no Dbas set', 0, RTDB_ERR)
!
      if (.not. bas_create(Abas, 'Abas'))
     $     call errquit('et_calc: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, Ageom, Abas, 'Abas'))
     $        call errquit('et_calc: no Abas set', 0, RTDB_ERR)
c
c  -------
c  get nbf
c  -------
c
      if (.not.bas_numbf(Dbas,Dnbf))
     $     call errquit('Dbas_numbf failed',0, BASIS_ERR)
!
      if (.not.bas_numbf(Abas,Anbf))
     $     call errquit('Abas_numbf failed',0, BASIS_ERR)
c
c ------------------
c initialize the API 
c ------------------
c
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
! Begin modifying by niuyingli
!      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
!     $        call errquit('task_et ga_create failed', 0, GA_ERR)
!      call ga_zero(g_s)
!      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)  
!
!      if (util_print("et",print_debug)) then
!       call ga_print(g_s)
!      endif
! Begin modifying by niuyingli
! Begin adding by niuyingli
c
c ---------------------------------------------------------------------
c calculate donor <AO|AO> used later to see if MOs are orthonormality
c----------------------------------------------------------------------

      if(.not.ga_create(mt_dbl,Dnbf,Dnbf,'s',0,0,g_s))
     $        call errquit('task_et ga_create failed', 0, GA_ERR)
      call ga_zero(g_s)
      call int_1e_ga(Dbas, Dbas, g_s, 'overlap', .false.)

      if (util_print("et",print_debug)) then
       call ga_print(g_s)
      endif
! End   adding by niuyingli
c
c ---------------------------
c Read MO vectors from files
c ---------------------------
c
      if (.not. rtdb_cget(rtdb, 'a_vecs', 1, a_file))
     $     call errquit('et_calc: couldnt get a_vecs filename',0,
     &       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'b_vecs', 1, b_file))
     $     call errquit('et_calc: couldnt get a_vecs filename',0,
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
      if (.not. rtdb_cget(rtdb, 'a_cis', 1, a_cis_file))
     $     call errquit('et_calc: couldnt get a_cis filename',0,
     &       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'b_cis', 1, b_cis_file))
     $     call errquit('et_calc: couldnt get a_cis filename',0,
     &       RTDB_ERR)
c
      lena_cis = inp_strlen(a_cis_file)
      lenb_cis = inp_strlen(b_cis_file)
      if (me.eq.0) then
      write(LuOut,*)'CI vectors for reactants: ',a_cis_file(1:lena_cis)
      write(LuOut,*)'CI vectors for products : ',b_cis_file(1:lenb_cis)
      write(LuOut,*)
      endif
		!
		!
		!=========================================================================
		!
		!
      if(.not.ga_create(mt_int,10,1,'anci',0,0,ga_nci))
     $        call errquit('task_et gci_create failed', 0, GA_ERR)
      call ga_zero(ga_nci)
		if (me==0) then
			open(41, file=a_cis_file, status="old")
			read(41,*)
			read(41,*)
			read(41,*) a_nroots
			read(41,*) a_nocc
			read(41,*) a_nmo
			read(41,*) 
			read(41,*) 
			read(41,*) a_nci
			read(41,*)
			iii=a_nci*2*a_nroots
			if (.not.ma_push_get(mt_dbl,iii,'acis', a_ncoeff,a_i_nc))
     $		call errquit('movecs_write: ma failed', 1, MA_ERR)
			iii=-1
			do iroot=1, a_nroots
				read(41,*)
				read(41,*)
				do inci=1, a_nci * 2
					iii=iii+1
					read(41,*) dbl_mb(a_i_nc+iii)
				end do
			end do
			close(41)
			call ga_put(ga_nci, 1, 1, 1, 1, a_nroots, 1)
			call ga_put(ga_nci, 2, 2, 1, 1, a_nocc  , 1)
			call ga_put(ga_nci, 3, 3, 1, 1, a_nmo   , 1)
			call ga_put(ga_nci, 4, 4, 1, 1, a_nci   , 1)
		end if
		call ga_sync()
		if (me/=0) then
			call ga_get(ga_nci, 1, 1, 1, 1, a_nroots, 1)
			call ga_get(ga_nci, 2, 2, 1, 1, a_nocc  , 1)
			call ga_get(ga_nci, 3, 3, 1, 1, a_nmo   , 1)
			call ga_get(ga_nci, 4, 4, 1, 1, a_nci   , 1)
		end if
		call ga_sync()
		!-----------------------------------------------------------------------------
		jjj=a_nci*2
		if(.not.ga_create(mt_dbl,jjj,a_nroots,'arci',0,0,ga_rci))
     $	call errquit('task_et gci_create failed', 0, GA_ERR)
		call ga_zero(ga_rci)
		if (me==0) then
			do iroot=1, a_nroots
				kkk = a_i_nc + (iroot-1)*jjj
				call ga_put(ga_rci,1,jjj,iroot,iroot,dbl_mb(kkk),1)
			end do
		end if
		!-------------------------------------------------------------------------
		call ga_sync()
		if (me/=0) then
			iii=a_nci*2*a_nroots
			if (.not.ma_push_get(mt_dbl,iii,'acis', a_ncoeff,a_i_nc))
     $		call errquit('movecs_write: ma failed', 1, MA_ERR)
			jjj=a_nci*2
			do iroot=1, a_nroots
				kkk = a_i_nc + (iroot-1)*jjj
				call ga_get(ga_rci,1,jjj,iroot,iroot,dbl_mb(kkk),1)
			end do
		end if
		!-------------------------------------------------------------------------
		jjj=a_nci*2*a_nroots
		call ga_sync()
		!
		!
		!=========================================================================
		!
		!
      if(.not.ga_create(mt_int,10,1,'bnci',0,0,gb_nci))
     $        call errquit('task_et gci_create failed', 0, GA_ERR)
      call ga_zero(gb_nci)
		if (me==0) then
			open(41, file=b_cis_file, status="old")
			read(41,*)
			read(41,*)
			read(41,*) b_nroots
			read(41,*) b_nocc
			read(41,*) b_nmo
			read(41,*) 
			read(41,*) 
			read(41,*) b_nci
			read(41,*)
			iii=b_nci*2*b_nroots
			if (.not.ma_push_get(mt_dbl,iii,'bcis', b_ncoeff,b_i_nc))
     $		call errquit('movecs_write: ma failed', 1, MA_ERR)
			iii=-1
			do iroot=1, b_nroots
				read(41,*)
				read(41,*)
				do inci=1, b_nci * 2
					iii=iii+1
					read(41,*) dbl_mb(b_i_nc+iii)
				end do
			end do
			close(41)
			call ga_put(gb_nci, 1, 1, 1, 1, b_nroots, 1)
			call ga_put(gb_nci, 2, 2, 1, 1, b_nocc  , 1)
			call ga_put(gb_nci, 3, 3, 1, 1, b_nmo   , 1)
			call ga_put(gb_nci, 4, 4, 1, 1, b_nci   , 1)
		end if
		call ga_sync()
		if (me/=0) then
			call ga_get(gb_nci, 1, 1, 1, 1, b_nroots, 1)
			call ga_get(gb_nci, 2, 2, 1, 1, b_nocc  , 1)
			call ga_get(gb_nci, 3, 3, 1, 1, b_nmo   , 1)
			call ga_get(gb_nci, 4, 4, 1, 1, b_nci   , 1)
		end if
		call ga_sync()
		!-----------------------------------------------------------------------------
		jjj=a_nci*2
		if(.not.ga_create(mt_dbl,jjj,b_nroots,'rci',0,0,gb_rci))
     $	call errquit('task_et gci_create failed', 0, GA_ERR)
		call ga_zero(gb_rci)
		if (me==0) then
			do iroot=1, b_nroots
				kkk = b_i_nc + (iroot-1)*jjj
				call ga_put(gb_rci,1,jjj,iroot,iroot,dbl_mb(kkk),1)
			end do
		end if
		!-------------------------------------------------------------------------
		call ga_sync()
		if (me/=0) then
			iii=b_nci*2*b_nroots
			if (.not.ma_push_get(mt_dbl,iii,'cis', b_ncoeff,b_i_nc))
     $		call errquit('movecs_write: mb failed', 1, MA_ERR)
			jjj=b_nci*2
			do iroot=1, b_nroots
				kkk = b_i_nc + (iroot-1)*jjj
				call ga_get(gb_rci,1,jjj,iroot,iroot,dbl_mb(kkk),1)
			end do
		end if
		!-------------------------------------------------------------------------
		jjj=b_nci*2*b_nroots
		call ga_sync()
		!
		!
		!=========================================================================
		!
		!
c
c -----------------------------------------------------
c Allocate occupation and eigenvalue vectors. 
c 'occ' is used to calculate nocc, 'eval' is never used.
c ------------------------------------------------------
c
! Begin adding by niuyingli
      if (.not.ma_push_get(mt_dbl,Dnbf,'occD',l_occD,k_occD))
     $         call errquit('ma_push_get failed for occ',k_occD, MA_ERR)

      if (.not.ma_push_get(mt_dbl,Dnbf,'evalsD',l_evalD,k_evalD))
     $        call errquit('ma_push_get failed for evals',k_evalD,
     &       MA_ERR)
C
      if (.not.ma_push_get(mt_dbl,Anbf,'occA',l_occA,k_occA))
     $         call errquit('ma_push_get failed for occ',k_occA, MA_ERR)

      if (.not.ma_push_get(mt_dbl,Anbf,'evalsA',l_evalA,k_evalA))
     $        call errquit('ma_push_get failed for evals',k_evalA,
     &       MA_ERR)

! End   adding by niuyingli
c
c ------------------------------------------------
c Read header to check that nset=2
c ------------------------------------------------
c
c a_file:
c -------
c
      call et_filename(a_file,a_filedir)
      if(.not.movecs_read_header(a_filedir(1:inp_strlen(a_filedir)),
     E     title,basis_name,
     &   scftype,nbf_file,nset,anmo,2)) 
     &   call errquit('et_calc: movecs_read_header failed',me, DISK_ERR)

      if (nbf_file.ne.Dnbf) then
        call errquit('et_calc: basis function mismatch',nbf_file,
     &       BASIS_ERR)
      endif

c b_file:
c -------
c
      call et_filename(b_file,b_filedir)
      if(.not.movecs_read_header(b_filedir(1:inp_strlen(b_filedir)),
     E     title,basis_name,
     &   scftype,nbf_file,nset,bnmo,2)) 
     &   call errquit('et_calc: movecs_read_header failed',me, DISK_ERR)

      if (nbf_file.ne.Anbf) then
        call errquit('et_calc: basis function mismatch',nbf_file,
     &       BASIS_ERR)
      endif
c
c -------------------------------------------------------
c allocate arrays for MO coefficients, for state A
c -------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,Dnbf,Dnbf,'donor',0,0,g_avecsa))
     $         call errquit('ga_create failed for A alpha',g_avecsa,
     &       GA_ERR)

c      if (.not.ga_create(mt_dbl,Dnbf,Dnbf,'A beta',0,0,g_avecsb))
c     $         call errquit('ga_create failed for A beta',g_avecsb,
c     &       GA_ERR)

      call ga_zero (g_avecsa)
c      call ga_zero (g_avecsb)
c 
c ------------------------------------------------------------------
c read occ, evals, mo coefficients for alpha orbitals from a_file
c ------------------------------------------------------------------
c
      if(.not.et_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     A     1,dbl_mb(k_occD)
     $   ,dbl_mb(k_evalD),g_avecsa,dum))
     $ call errquit('et_movecs_read failed for alpha mos',0, DISK_ERR)
c
c -------------------------------------
c get nocc(1): # of occupied alpha mos
c -------------------------------------
c
		write(*,*) "ok1", me
      nocc(1)=0
      do i = 1,Dnbf
       nocc(1) = nocc(1) + dbl_mb(k_occ+i-1)
      enddo

      if (util_print("et",print_high)) then
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
      if(.not.et_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     E     1,dbl_mb(k_occA)
     $     ,dbl_mb(k_eval),g_avecsb, dum))
     $ call errquit('et_movecs_read failed for beta mos',0, DISK_ERR)
		call ga_print(g_avecsb)
c
c -------------------------------------
c get nocc(2): # of occupied alpha mos
c -------------------------------------
c
      nocc(2)=0
      do i = 1,nbf
       nocc(2) = nocc(2) + dbl_mb(k_occ+i-1)
      enddo

      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(LuOut,*) nocc(2), " occupied  beta orbitals"
        write(LuOut,*)
      endif
        call ga_print(g_avecsb)
      endif
		write(*,*) "okok", nocc
c
c -------------------------------------------------------
c allocate arrays for MO coefficients, for state B
c -------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,Anbf,Anbf,'B alpha',0,0,g_bvecsa))
     $         call errquit('ga_create failed for B alpha',g_bvecsa,
     &       GA_ERR)
c      if (.not.ga_create(mt_dbl,Anbf,Anbf,'B beta',0,0,g_bvecsb))
c     $         call errquit('ga_create failed for B beta',g_bvecsb,
c     &       GA_ERR)

      call ga_zero (g_bvecsa)
c      call ga_zero (g_bvecsb)
c
c -----------------------------------------------------------
c read occ, evals, mo coefficients from b_file.  It is 
c assumed that the number of occupied alpha and beta orbitals
c is the same as in the a_file.
c -----------------------------------------------------------
c
c
c -----------------------------------------------
c calculate D=B^T*S*A for alpha and beta electrons
c -----------------------------------------------
c
      g_ahandle(1)=g_avecsa
      g_ahandle(2)=g_avecsb
      g_bhandle(1)=g_bvecsa
      g_bhandle(2)=g_bvecsb
c
c
c detprod is the product of det(U)*det(V) for alpha and beta electrons,
c and is used later.
c ----------------------------------------------------------------------
c
      detprod=1.0d0
c
c -----------------------------------------------------
c  pname is needed for generalized density matrices:
c  Pa is for alpha electrons, Pb is for beta electrons.
c -----------------------------------------------------
c
      pname(1)='pa'
      pname(2)='pb'
c
c ------------------------
c begin alpha/beta loop j
c ------------------------
c
      do j = 1,2
c
      if (util_print("et",print_high)) then
      if(me.eq.0) then
       write(LuOut,*) 
       write(LuOut,*) '---------------'
       if(j.eq.1) write(LuOut,*) 'ALPHA ELECTRONS'
       if(j.eq.2) write(LuOut,*) ' BETA ELECTRONS'
       write(LuOut,*) '---------------'
      endif
      endif
c
c  -----------------------------------------------
c  Perform matrix multiplication to get D = B^T*S*A
c  -----------------------------------------------
c
      call ga_sync()
		write(*,*) me,"begin", nocc(j)
      call ga_sync()
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'d',0,0,g_d))
     $        call errquit('task_et ga_create failed', g_d, GA_ERR)
		write(*,*) me,"middle1"
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'temp',0,0,g_tmp))
     $        call errquit('task_et ga_create failed', g_d, GA_ERR)
		write(*,*) me,"middle2"
c
      call ga_dgemm('T','N',nocc(j),nbf,nbf,1.d0,g_bhandle(j)
     &               ,g_s,0.d0,g_tmp)
		write(*,*) me,"middle3"
c
      call ga_dgemm('N','N',nocc(j),nocc(j),nbf,1.d0,g_tmp
     &               ,g_ahandle(j),0.d0,g_d)
		write(*,*) me,"middle4"
c
      if (.not.ga_destroy(g_tmp))
     &  call errquit('could not destroy tmp array',g_tmp, GA_ERR)
		write(*,*) me,"middle5"
c
      if (util_print("et",print_debug)) then
       call ga_print(g_d)
      endif
		write(*,*) me,"middle6"
      call ga_sync()
		write(*,*) me,"end"
      call ga_sync()
		stop

C
      if (.not. ma_push_get(mt_dbl,Dnbf,'mo_write',l_vecsDi,k_vecsDi))
     $        call errquit('movecs_write: ma failed', Dnbf, MA_ERR)
	call ga_get(g_ahandle,1,Dnbf,1,1,dbl_mb(k_vecsDi),1)
      if (.not.ga_destroy(g_d))
     &  call errquit('could not destroy D array',g_d, GA_ERR)
		write(63,'(<Dnbf>f11.5)') dbl_mb(k_vecsDi:k_vecsDi+Dnbf-1)
		call int_2e4c_ga(Dbas,Abas,Dnbf,Anbf,g_ahandle,g_bhandle,
     & ga_nci, gb_nci, ga_rci, gb_rci, me)
		stop

c
c ----------------------------------------------
c Perform SVD to get U, V^T, and singluar values
c ----------------------------------------------
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'u',0,0,g_u))
     $        call errquit('task_et ga_create failed', g_u, GA_ERR)
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'vt',0,0,g_vt))
     $        call errquit('task_et ga_create failed', g_vt, GA_ERR)
c
      call ga_zero(g_u)
      call ga_zero(g_vt)
      if (.not.ma_push_get(mt_dbl,nocc(j),'values',l_vals,k_vals))
     $        call errquit('ma_push_get failed for values ',k_vals,
     &       MA_ERR)
c
      call ga_svd_seq(g_d,g_u,g_vt,dbl_mb(k_vals))
c
c prod is the product of singular values.
c ---------------------------------------
c
      prod(j) = 1.0d0
      do k = 1,nocc(j)
         prod(j) = prod(j)*dbl_mb(k_vals+k-1)
      enddo
c      write(LuOut,*)'prod',j,prod(j)
c      write(LuOut,*)'singular values',(dbl_mb(k_vals+i-1),i=1,nocc(j))
c
c
      if (.not.ga_destroy(g_d))
     &  call errquit('could not destroy D array',g_d, GA_ERR)
c
      if (util_print("et",print_debug)) then
      call ga_print(g_u)
      call ga_print(g_vt)
      endif
		stop
c
c
c --------------------------
c calculate det(U) and det(V)
c --------------------------
c
       call et_det('N',g_u,nocc(j),detu(j))
       call et_det('T',g_vt,nocc(j),detv(j))
       detprod = detprod*detu(j)*detv(j)
c       write(LuOut,*) j,'detu,detv',detu(j),detv(j)
c
c ------------------------------------------------------------
c                                               ^       ^
c perform matrix mult to get transformed A & B: A = AV, B = BU
c ------------------------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'av',0,0,g_av))
     $        call errquit('task_et ga_create failed', g_av, GA_ERR)
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'bu',0,0,g_bu))
     $        call errquit('task_et ga_create failed', g_bu, GA_ERR)
c
      call ga_zero(g_av)
      call ga_dgemm('N','T',nbf,nocc(j),nocc(j),1.d0,g_ahandle(j),g_vt
     &         ,0.d0,g_av)
c
      call ga_zero(g_bu)
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_bhandle(j),g_u
     &         ,0.d0,g_bu)
c
c
      if (.not.ga_destroy(g_ahandle(j))) 
     &  call errquit('could not destroy mo vectors array',g_ahandle(j),
     &       GA_ERR)
      if (.not.ga_destroy(g_bhandle(j))) 
     &  call errquit('could not destroy mo vectors array',g_bhandle(j),
     &       GA_ERR)
c
      if (.not.ga_destroy(g_u))
     &  call errquit('could not destroy U array',g_u,
     &       GA_ERR)
      if (.not.ga_destroy(g_vt))
     &  call errquit('could not destroy VT array',g_vt,
     &       GA_ERR)
c
      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(LuOut,*) 
        write(LuOut,*) 'SVD transformed MO vectors for state A'
      endif
        call ga_print (g_av)
      if(me.eq.0) then
        write(LuOut,*) 
        write(LuOut,*) 'SVD transformed MO vectors for state B'
      endif
        call ga_print (g_bu)
      endif
c
c
c    --------------------------------------------------------------
c                                   N
c                                -------
c     Calculates   T(i) =     PRODUCT  d(j) where N = nalpha or nbeta
c                                -------
c                                j.NE.i  
c    --------------------------------------------------------------
c    
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'t',0,0,g_t))
     $        call errquit('task_et ga_create failed', g_t, GA_ERR)
      call ga_zero(g_t)
c
      if(nocc(j).gt.1) then
c
        do i = 1,nocc(j)
          val = 1.0d0
            do k = 1,nocc(j)
              if(k.ne.i) val = val * dbl_mb(k_vals+k-1)
            enddo
          call ga_put(g_t,i,i,i,i,val,1)
        enddo
      else
         val = 1.0d0
         call ga_put(g_t,1,1,1,1,val,1)
c
      endif
      call ga_sync()
c
      if (util_print("et",print_debug)) then
      call ga_print(g_t) 
      endif
c
c -----------------------------------------------------
c                                            ^ ^
c  calculate generalized density matrix, P = ATB^T
c  Pa is for alpha electrons, Pb is for beta electrons.
c -----------------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,pname(j),0,0,g_phandle(j)))
     $        call errquit('task_et ga_create failed', g_phandle(j),
     &       GA_ERR)
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'temp',0,0,g_temp))
     $        call errquit('task_et ga_create failed', g_temp, GA_ERR)
c
      call ga_zero(g_phandle(j))
c
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_av,g_t
     &         ,0.d0,g_temp)
c
      call ga_dgemm('N','T',nbf,nbf,nocc(j),1.0d0,g_temp,g_bu
     &         ,0.0d0,g_phandle(j))
c
c
      if (.not.ga_destroy(g_av)) 
     &  call errquit('could not destroy transformed MO vectors',g_av,
     &       GA_ERR)
      if (.not.ga_destroy(g_bu)) 
     &  call errquit('could not destroy transformed MO vectors',g_bu,
     &       GA_ERR)
c
      if (.not.ga_destroy(g_t)) 
     &  call errquit('could not destroy T array',g_t, GA_ERR)
      if (.not.ga_destroy(g_temp)) 
     &  call errquit('could not destroy temp array',g_temp, GA_ERR)
c
c
      if (util_print("et",print_debug)) then
      call ga_print(g_phandle(j))
      endif
c
c
c --------------------
c end alpha/beta loop:
c --------------------
      enddo    
c
      if (.not.ga_destroy(g_s)) 
     &  call errquit('could not destroy overlap array',g_s, GA_ERR)

c
c Pa is scaled by prod(b), Pb is scaled by prod(a)
c ------------------------------------------------
c
       call ga_dscal(g_phandle(1),prod(2))
       call ga_dscal(g_phandle(2),prod(1))

      if (util_print("et",print_debug)) then
      if(me.eq.0) then
       write(LuOut,*) 
       write(LuOut,*) 'scaled generalized density matrices:' 
       write(LuOut,*) 
      endif
       call ga_print(g_phandle(1))
       call ga_print(g_phandle(2))
      endif
		write(*,*)me,"ok"
		stop
c
c
c -----------------------------------------
c calculate overlap of state A with state B
c -----------------------------------------
c
       Sba = prod(1)*prod(2)*detprod
c
c
      if(me.eq.0) then
       write(LuOut,*)
       write(LuOut,9997) Sba
       write(LuOut,*)
      endif
c
 9997 format(1x,'Reactants/Products overlap S(RP) : ',1Pd9.2)
c
c --------------------------------------------
c calculate the 1e contribution to Vba, omega1
c --------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'h1',0,0,g_h1))
     $        call errquit('task_et ga_create failed', g_h1, GA_ERR)
c
c get H1 integrals
c ----------------
c
      call ga_zero(g_h1)
      call int_1e_ga(basis, basis, g_h1, 'kinetic', .false.)  
      call int_1e_ga(basis, basis, g_h1, 'potential', .false.)  
c
      if (util_print("et",print_debug)) then
          call ga_print(g_h1)
      endif
c
       omega1 = 0.0d0
       omega1 = omega1 + ga_ddot (g_h1,g_phandle(1))
       omega1 = omega1 + ga_ddot (g_h1,g_phandle(2))
c
       call ga_sync()
       omega1 = omega1*detprod
c
      if(me.eq.0) then
       write(LuOut,*) 'Reactants/Products interaction energy:    '
       write(LuOut,*) '-------------------------------------           '
       write(LuOut,9995) 'One-electron contribution        ','H1(RP)',
     &    omega1
      endif
c
      if (.not.ga_destroy(g_h1))
     &  call errquit('could not destroy H1 array',g_h1, GA_ERR)
c
c --------------------------------------------
c calculate the 2e contribution to Vba, omega2
c --------------------------------------------
c
      if (.not. rtdb_get(rtdb, 'et:tol2e', MT_DBL, 1, tol2e)) then
       tol2e    =  min(1.0d-7,dabs(Sba)*1.0d-7)
       tol2e    =  max(1.0d-12,tol2e)
      endif

c ------------------------------------------------------------------------
c Two methods are available. One uses the scf fock builder (method_2e=1),
c the other uses et_2e (method_2e=2). et_2e is more accurate, but much more
c time consuming and should only be used for very small systems, or very 
c weakly interacting systems (with small overlap), or for debugging purposes.
c ------------------------------------------------------------------------
c
      if (.not.rtdb_get(rtdb,'et:method_2e',MT_INT,1,method_2e)) then
          method_2e = 1
      endif
	! debug niuyingli
	method_2e=2
c
      if (method_2e.EQ.1) then
c for some unknown reason I have to terminate and re-init the API for
c et_fock to work with ECP's
      call int_terminate()
      call et_fock(rtdb,geom,basis,
     &  g_phandle(1),g_phandle(2),nbf,tol2e,omega2)
      endif
c
      if (method_2e.EQ.2) then
       call int_mem_2e4c(mem2,max2e)
       mem2 = max(mem2,max2e)
       call et_2e(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &           g_phandle(1),g_phandle(2),omega2)
      endif
c
c ----------------------------------------------------
c destroy density arrays and calculate omega2 and Hba
c ----------------------------------------------------
c
      if (.not.ga_destroy(g_phandle(1)))
     &  call errquit('could not destroy P array',g_phandle(1), GA_ERR)
      if (.not.ga_destroy(g_phandle(2)))
     &  call errquit('could not destroy P array',g_phandle(2), GA_ERR)
c
       omega2 = omega2 * (detprod/ ( prod(1)*prod(2) ) ) 
       omega2 = omega2 / 2.0d0
c
       Hba = omega1+omega2
c
      if(me.eq.0) then
       write(LuOut,9995) 'Two-electron contribution        ','H2(RP)',
     &    omega2
       write(LuOut,9995) 'Total interaction energy         ','H(RP)',
     &    Hba
       write(LuOut,*)
      endif
c
c
c ---------------------------------------------
c Finally, calculate the ET matrix element, Vba
c ---------------------------------------------
c
c
       Vba =  Hba - Sba*(Haa+Hbb)*0.5d0
       Vba = Vba/(1.0d0-Sba*Sba)
       Vba = dabs (Vba)
c
       if (.not. rtdb_put(rtdb,'et:V(RP)', MT_DBL, 1, Vba))
     $        call errquit('et_calc: rtdb_put of V(RP) failed',0,
     &       RTDB_ERR)
c
c
c -------------
c Print results
c -------------
c
      if(me.eq.0) then
      write(LuOut,9996) Vba,Vba*2.194746d5,
     &              Vba*27.2114d0,Vba*627.510d0
      write(LuOut,*) 
      endif
c
 9995 format (a34,1x,a6,1x,F18.10)
 9996 format( 1x,'Electron Transfer Coupling Energy |V(RP)|',F18.10,/, 
     &  49x,F11.3,' cm-1',/,
     &  49x,F11.6,' eV',/,
     &  49x,F11.3,' kcal/mol')
c


      if (.not.ma_chop_stack(l_vals)) 
     &    call errquit('unable to destroy MA',l_vals, MA_ERR)
      if (.not.ma_chop_stack(l_eval)) 
     &    call errquit('unable to destroy MA',l_eval, MA_ERR)
      if (.not.ma_chop_stack(l_occ)) 
     &    call errquit('unable to destroy MA',l_occ, MA_ERR)
c
c
      et_calc=.true.
c
      call util_print_pop
c
         if (.not. geom_destroy(geom)) call errquit
     $        ('et_calc: geom destroy failed', 0, GEOM_ERR)
         if (.not. bas_destroy(basis)) call errquit
     $        ('et_calc: basis destroy failed',0, BASIS_ERR)
c
      call int_terminate()
c
c
      return
      end
      subroutine et_filename(file,dirfile)
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
c     $Id: et_2e.F,v 1.16 2003-10-17 22:54:54 carlfahl Exp $
!-------------------------------------------------------------------------------
! Added by niuyingli
      subroutine int_2e4c_ga(Dbas,Abas,Dnbf,Anbf,g_D,g_A, 
     $ ga_nci, gb_nci, ga_rci, gb_rci, me)
      implicit none
#include "errquit.fh"
c
c     calculates the two-electron contribution to Vda
c
c     Dbas, Abas  (for int_2e4c)                                [input]
c     g_s (two electron integrals for Vda)                      [output]
c
#include "schwarz.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
c

      integer Dbas, Abas, Dnbf, Anbf, g_D, g_A

c   local variable 
      integer ish, jsh, ksh, lsh, i, j, k, l, m, n
      integer max2e, mem2, l_buf, k_buf, l_scr, k_scr
      integer l_vecsDi, k_vecsDi, l_vecsAk, k_vecsAk
      integer l_vecsDj, k_vecsDj, l_vecsAl, k_vecsAl
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer nshell_d, nshell_a
      integer Dnbf_row, Anbf_row, x, y
      character(len=80) newmovecs
      logical print_vda

c      real(8), allocatable ::  Dcoef(:,:), Acoef(:,:)
      real(8), allocatable ::  DMOco(:), AMOco(:)
      real(8), allocatable ::  g(:,:,:,:)
      real(8)  sum_vda, sum_tmp, sum_dial, sum_dia, tmp
      integer  num_D_CI, num_A_CI
      integer, allocatable :: DMOi(:), DMOa(:), AMOi(:), AMOa(:)
		integer ga_nci, ga_rci
		integer gb_nci, gb_rci
		integer me
		integer a_nroots, a_ipol, a_nocc, a_nmo, a_nfc, a_nfv
		integer b_nroots, b_ipol, b_nocc, b_nmo, b_nfc, b_nfv
		integer a_nci, a_ncoeff, a_i_nc
		integer b_nci, b_ncoeff, b_i_nc
		integer iroot, inci, iii, jjj, kkk
      if (.not. bas_numcont(Dbas, nshell_d)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', Dbas,
     &       BASIS_ERR)

      if (.not. bas_numcont(Abas, nshell_a)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', Abas,
     &       BASIS_ERR)

c
c     allocate necessary local temporary arrays on the stack
c
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c
c     k_* are the offsets corrsponding to the l_* handles
       call int_acc_set(1d-30)
c
c
      call int_mem_2e4c(max2e, mem2)
c
      if (.not. MA_push_get(MT_DBL,max2e,'int_2e4c_ga:buf',l_buf,k_buf))
     $     call errquit('int_2e4c_ga: ma failed', max2e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem2,'int_2e4c_ga:scr',l_scr,k_scr))
     $     call errquit('int_2e4c_ga: ma failed', mem2, MA_ERR)
cc
c
c ----------  read_input------------------------
        open(7,file="input",status="old")
        read(7,"(14X,A40)") newmovecs
        read(7,"(14X,L10)") print_vda
        read(7,"(14X,I6)") num_D_CI
        read(7,"(14X,I6)") num_A_CI

        allocate(DMOi(num_D_CI))
        allocate(DMOa(num_D_CI))
        allocate(DMOco(num_D_CI))
        allocate(AMOi(num_A_CI))
        allocate(AMOa(num_A_CI))
        allocate(AMOco(num_A_CI))
          
       DMOi=0
       DMOa=0
       DMOco=0
       AMOi=0
       AMOa=0
       AMOco=0
       do n=1, num_D_CI
         read(7,"(14X,2I5,F18.8)") DMOi(n), DMOa(n), DMOco(n) !CI coefficients of D
       enddo

       do n=1, num_A_CI
         read(7,"(14X,2I5,F18.8)") AMOi(n), AMOa(n), AMOco(n) !CI coefficients of A
       enddo

       close(7)
       write(LuOut,*) DMOi(1),AMOa(1),AMOco(1)

      if (.not. ma_push_get(mt_dbl,Dnbf,'mo_write',l_vecsDi,k_vecsDi))
     $        call errquit('movecs_write: ma failed', Dnbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl,Dnbf,'mo_write',l_vecsDj,k_vecsDj))
     $        call errquit('movecs_write: ma failed', Dnbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl,Anbf,'mo_write',l_vecsAk,k_vecsAk))
     $        call errquit('movecs_write: ma failed', Anbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl,Anbf,'mo_write',l_vecsAl,k_vecsAl))
     $        call errquit('movecs_write: ma failed', Anbf, MA_ERR)
c---------------------------------------------------------------
      sum_vda=0.0D0
c
c     Loop thru shells with static parallel work decomposition
c

      if (print_vda) then
      open(77,file='vda.dat',status="unknown")
        write(77,*) " i  j  k  l   g(i,j,k,l) Dcoef(i,DMOi(n)) 
     & Dcoef(j,DMOa(n)) Acoef(k,AMOi(m)) Acoef(l,AMOa(m)) "
      endif

          do n=1, num_D_CI
          call ga_get(g_D,1,Dnbf,DMOi(n),DMOi(n),dbl_mb(k_vecsDi),1)
          call ga_get(g_D,1,Dnbf,DMOa(n),DMOa(n),dbl_mb(k_vecsDj),1)
            do m=1, num_A_CI
          call ga_get(g_A,1,Anbf,AMOi(m),AMOi(m),dbl_mb(k_vecsAk),1)
          call ga_get(g_A,1,Anbf,AMOa(m),AMOa(m),dbl_mb(k_vecsAl),1)

      do lsh = 1, nshell_a
        if (.not. bas_cn2bfr(Abas, lsh, llo, lhi))
     $       call errquit('int_2e4c_ga: bas_cn2bfr ?', Abas,
     &       BASIS_ERR)
        do ksh = 1, nshell_a
          if (.not. bas_cn2bfr(Abas, ksh, klo, khi))
     $         call errquit('int_2e4c_ga: bas_cn2bfr ?', Abas,
     &       BASIS_ERR)
          do jsh = 1, nshell_d
            if (.not. bas_cn2bfr(Dbas, jsh, jlo, jhi))
     $           call errquit('int_2e4c_ga: bas_cn2bfr ?', Dbas,
     &       BASIS_ERR)
            do ish = 1, nshell_d  
              if (.not. bas_cn2bfr(Dbas, ish, ilo, ihi))
     $              call errquit('int_2e4c_ga: bas_cn2bfr ?', Dbas,
     &       BASIS_ERR) 
             call int_2e4c (Dbas, ish, jsh, Abas, ksh, lsh,
     $              mem2, dbl_mb(k_scr), max2e, dbl_mb(k_buf))

       allocate(g(ilo:ihi,jlo:jhi,klo:khi,llo:lhi))
        g=0.0d0        

       call copy_2e4c_int (ilo, ihi, jlo, jhi, klo, khi,
     &           llo, lhi, dbl_mb(k_buf), g)
  
          sum_tmp=0.0D0
              do i = ilo, ihi
                do j = jlo, jhi
                  do k = klo, khi
                    do l = llo, lhi
      sum_tmp=sum_tmp+dbl_mb(k_vecsDi+i-1)*dbl_mb(k_vecsDj+j-1)
     & *dbl_mb(k_vecsAk+k-1)*dbl_mb(k_vecsAl+l-1)*g(i,j,k,l)
     & *DMOco(n)*AMOco(m)
				write(81,'(4i6, es20.12)') i,j,k,l, g(i,j,k,l)
                      if (print_vda) then
                      endif   
                    enddo
                  enddo
                enddo
             enddo

 
            sum_vda=sum_vda+sum_tmp
 
        deallocate(g)

                end do
            end do
         end do
       end do
       enddo
       enddo
		stop


!             call print_2e4c_int (Dnbf, Anbf, dbl_mb(k_buf))

c      if (.not. MA_pop_stack(l_vecsD)) call errquit
c     $     ('int_2e4c_ga: pop failed', 0, GA_ERR)
c      if (.not. MA_pop_stack(l_vecsA)) call errquit
c     $     ('int_2e4c_ga: pop failed', 0, GA_ERR)

              sum_vda=2*sum_vda

      write(LuOut,*) "coulomb coupling between D&A:"
      write(LuOut,*) sum_vda,"a.u.", sum_vda*27.2114,"eV"
      if (print_vda) then
      write(77,*) "coulomb coupling between D&A:"
      write(77,*) sum_vda,"a.u.", sum_vda*27.2114,"eV"
      endif
c      write(LuOut,*) "when i=j=k=l,Vda:"
c      write(LuOut,*) sum_dial,"a.u.", sum_dial*27.2114,"eV"
c      close(77)
c
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_2e4c_ga: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_2e4c_ga: pop failed', 0, GA_ERR)
c
      deallocate(DMOi)
      deallocate(DMOa)
      deallocate(DMOco)
      deallocate(AMOi)
      deallocate(AMOa)
      deallocate(AMOco)
c      deallocate(Dcoef)
c      deallocate(Acoef)

   
      end

!      subroutine print_2e4c_int (Dnbf, Anbf, g)
!      implicit none
!      integer Dnbf, Anbf

       subroutine copy_2e4c_int (ilo, ihi, jlo, jhi, klo, khi,
     &           llo, lhi, g_s, g)
       implicit none
c
       integer ilo, ihi, jlo, jhi, klo,
     &        khi, llo, lhi, Dnbf, Anbf
c     &        num_D_CI, num_A_CI
c     &        DMOi(num_D_CI), DMOa(num_D_CI),
c     &        AMOi(num_A_CI), AMOa(num_A_CI)
c      double precision  DMOco(num_D_CI), AMOco(num_A_CI),
c     &                  Dcoef(Dnbf,Dnbf), Acoef(Anbf,Anbf),
      double precision  g(ilo:ihi,jlo:jhi,klo:khi,llo:lhi)
c      double precision  sum_vdan, sum_tmp
      double precision g_s(ilo:ihi,jlo:jhi,klo:khi,llo:lhi)  
      integer i, j, k, l
c     
          do i = ilo, ihi
            do j = jlo, jhi
              do k = klo, khi
                do l = llo, lhi
             g(i,j,k,l)=g_s(i,j,k,l)
                 enddo
               enddo
            enddo
          enddo
cc 
      end
