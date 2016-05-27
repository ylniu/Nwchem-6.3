		logical function ebt_calc(rtdb,Vba)
		implicit none
		!-------------------------------------------------------------
		!
		! $Id: ebt_calc.F 19821 2016-04-16 08:46:49Z d3y133 $
		!
		! Program Author
		!   Ying-Li Niu
		!   niuyl@nanoctr.cn
		!   National Center for Nanoscience and Technology,
		!   Beijing, Peoples republic of China.
		!
		! This function calculates energy transfer matrix element.
		!
		! Paper
		!   J. Chern. Phys. 101, 10521, (1994).
		!
		! Authors
		!   Richard D. Harcourt
		!   Gregory D. Scholes
		!   Kenneth P. Ghiggino
		!
		! Title
		!   Rate expressions for excitation transfer. 
		!   II. Electronic considerations of direct and 
		!   through-configuration exciton resonance interactions
		!
		!-------------------------------------------------------------
#include "errquit.fh"
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
		!-------------------------------------------------------------
		integer  rtdb, me
		integer  i, j, k, l, fid
		integer  g_tmp,g_s,g_d
		integer  nocc(2)
		integer  g_h1
		integer  mem2, max2e, maxl
		integer  method_2e, output_mo, nmo2(2), nbf_file, nset
		integer  rotate_homo
		integer  rotate_lumo
		!-------------------------------------------------------------
		logical  ebt_movecs_read
		external ebt_movecs_read
		logical  movecs_write
		external movecs_write
		logical  movecs_read_header
		external movecs_read_header
		!-------------------------------------------------------------
		integer              :: g_mo
		integer              :: g_omo4
		integer              :: g_mo4
		integer              :: g_s4
		integer              :: g_h4
		integer              :: g_den(4,4)
		integer              :: ga_den
		integer              :: gb_den
		integer              :: natom1
		integer              :: natom2
		integer, allocatable :: mol1(:)
		integer, allocatable :: mol2(:)
		real(8)              :: au2a  = 0.529177D0
		real(8)              :: au2cm = 219474.6307D0
		real(8)              :: coc(3,2)
		real(8)              :: dip(3,2)
		real(8)              :: theta(3)
		real(8)              :: rtmp
		real(8)              :: dr
		real(8)              :: coup12
		real(8)              :: coup14
		real(8)              :: const
		real(8)              :: H, V, kappa
		real(8)              :: time(11)
		real(8)              :: dtime(10)
		real(8)              :: hab(4,4)
		real(8)              :: sab(4,4)
		real(8)              :: Jab(4,4,4,4)
		real(8), allocatable :: den_a(:,:)
		real(8), allocatable :: den_b(:,:)
		real(8), allocatable :: omo(:,:)
		real(8), allocatable :: den(:,:,:,:)
		real(8), allocatable :: mo(:,:)
		character(200)       :: gname
		character(200)       :: a_file
		character(200)       :: basis_name
		character(200)       :: a_filedir
		character(200)       :: f_movecs_out
		character(200)       :: outfile
		!-------------------------------------------------------------
		double precision Vba
		double precision Haa, dum
		data Haa/0.0d0/
		!-------------------------------------------------------------
		fid=88
		!-------------------------------------------------------------
		me=ga_nodeid()
		if (me.eq.0) then
			Call Util_Print_Centered(LuOut,
     &	'Electron Bridge Transfer Calculation',40,.true.)
			write(LuOut,*)
		end if
		!-------------------------------------------------------------
		rotate_homo=1
		if (.not.rtdb_get(rtdb,'ebt:rotate_homo',MT_INT,1,rotate_homo)) then
		endif
		!-------------------------------------------------------------
		rotate_lumo=1
		if (.not.rtdb_get(rtdb,'ebt:rotate_lumo',MT_INT,1,rotate_lumo)) then
		endif
		!-------------------------------------------------------------
		output_mo=0
		if (.not.rtdb_get(rtdb,'ebt:output_mo',MT_INT,1,output_mo)) then
		endif
		!-------------------------------------------------------------
		method_2e=1
		if (.not.rtdb_get(rtdb,'ebt:method_2e',MT_INT,1,method_2e)) then
		endif
		!-------------------------------------------------------------
		outfile="ebridge.ebt"
		if (.not.rtdb_cget(rtdb,'outfile',1,outfile)) then
		endif
		!-------------------------------------------------------------
		! for printing:
		!--------------
		call util_print_push
		call util_print_rtdb_load(rtdb,"ebt")
		!-------------------------------------------------------------
		! load geometry and symmetry info
		!--------------------------------
		if (.not. geom_create(geom, 'geometry'))
     $	call errquit('ebt_calc: geom_create?', 0, GEOM_ERR)
		if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $	call errquit('ebt_calc: no geometry ', 0, RTDB_ERR)
		!-------------------------------------------------------------
		! load the basis set and get info about it
		!-----------------------------------------
		if (.not. bas_create(basis, 'ao basis'))
     $	call errquit('ebt_calc: bas_create?', 0, BASIS_ERR)
		if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $	call errquit('ebt_calc: no ao basis set', 0, RTDB_ERR)
		!-------------------------------------------------------------
		! get nbf
		!--------
		if (.not.bas_numbf(basis,nbf))
     $	call errquit('bas_numbf failed',0, BASIS_ERR)
		!-------------------------------------------------------------
		! initialize the API 
		!-------------------
		call int_init(rtdb, 1, basis)
		call int_acc_std()
		call int_normalize(rtdb,basis)
		!-------------------------------------------------------------
		! calculate the overlap integrals and H1 integrals
		!-------------------------------------------------
		if (.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
     $	call errquit('task_et ga_create failed', 0, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_create(mt_dbl,nbf,nbf,'h1',0,0,g_h1))
     $	call errquit('task_et ga_create failed', g_h1, GA_ERR)
		!-------------------------------------------------------------
		call ga_zero(g_s )
		call ga_zero(g_h1)
		!-------------------------------------------------------------
		call int_1e_ga(basis, basis, g_s,  'overlap'  , .false.)  
		call int_1e_ga(basis, basis, g_h1, 'kinetic'  , .false.)
		call int_1e_ga(basis, basis, g_h1, 'potential', .false.)  
		!-------------------------------------------------------------
		! Read MO vectors from files
		!---------------------------
		if (.not. rtdb_cget(rtdb, 'a_vecs', 1, a_file))
     &	call errquit('ebt_calc: couldnt get a_vecs filename',0,RTDB_ERR)
		if (.not. rtdb_cget(rtdb, 'b_vecs', 1, f_movecs_out))
     &	call errquit('ebt_calc: couldnt get b_vecs filename',0,RTDB_ERR)
		
		if (me.eq.0) then
			write(LuOut,*)'MO vectors from file: ',trim(a_file)
		end if
		!-------------------------------------------------------------
		! Allocate occupation and eigenvalue vectors. 
		! 'occ' is used to calculate nocc, 'eval' is never used.
		!-------------------------------------------------------
		if (.not.ma_push_get(mt_dbl,nbf,'occ',l_occ,k_occ))
     &	call errquit('ma_push_get failed for occ',k_occ,MA_ERR)
		!-------------------------------------------------------------
		if (.not.ma_push_get(mt_dbl,nbf,'evals',l_eval,k_eval))
     &	call errquit('ma_push_get failed for evals',k_eval,MA_ERR)
		!-------------------------------------------------------------
		! read from file:
		!----------------
		call ebt_filename(a_file,a_filedir)
		if(.not.movecs_read_header(a_filedir(1:inp_strlen(a_filedir)),
     &	title,basis_name,scftype,nbf_file,nset,nmo2,2)) 
     &	call errquit('ebt_calc: movecs_read_header failed',me, DISK_ERR)
		if (nbf_file.ne.nbf) then
			call errquit('ebt_calc: basis function mismatch',nbf_file,
     &	BASIS_ERR)
		endif
		!-------------------------------------------------------------
		! allocate arrays for MO coefficients, for state A
		!-------------------------------------------------
		if (.not.ga_create(mt_dbl,nbf,nbf,'ga_den',0,0,ga_den))
     $	call errquit('ga_create failed for ga_den',ga_den, GA_ERR)
		if (.not.ga_create(mt_dbl,nbf,nbf,'gb_den',0,0,gb_den))
     $	call errquit('ga_create failed for gb_den',gb_den, GA_ERR)
		if (.not.ga_create(mt_dbl,nbf,nbf,'A alpha',0,0,g_mo))
     $	call errquit('ga_create failed for A alpha',g_mo, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_create(mt_dbl,nbf,4,'omo4',0,0,g_omo4))
     $	call errquit('ga_create failed for A alpha',g_omo4, GA_ERR)
		!-------------------------------------------------------------
			if(.not.ga_create(mt_dbl,4,nbf,'temp',0,0,g_tmp))
     $		call errquit('task_et ga_create failed', g_tmp, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_create(mt_dbl,nbf,4,'mo4',0,0,g_mo4))
     $	call errquit('ga_create failed for A alpha',g_mo4, GA_ERR)
		!-------------------------------------------------------------
		if(.not.ga_create(mt_dbl,4,4,'s4',0,0,g_s4))
     $	call errquit('task_et ga_create failed', g_s4, GA_ERR)
		!-------------------------------------------------------------
		if(.not.ga_create(mt_dbl,4,4,'h4',0,0,g_h4))
     $	call errquit('task_et ga_create failed', g_h4, GA_ERR)
		!-------------------------------------------------------------
		call ga_zero (g_s4   )
		call ga_zero (g_h4   )
		call ga_zero (g_mo   )
		call ga_zero (g_mo4  )
		call ga_zero (g_omo4 )
		!-------------------------------------------------------------
		! read occ, evals, mo coefficients for alpha orbitals from a_file
		!-------------------------------------------------------------
		if(.not.ebt_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     &	1,dbl_mb(k_occ),dbl_mb(k_eval),g_mo,dum))
     &	call errquit('ebt_movecs_read failed for alpha mos',0, DISK_ERR)
		!-------------------------------------------------------------
		! Read converged scf or dft energy from movecs file.
		!-------------------------------------------------------------
		if(.not.ebt_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     &	2,dbl_mb(k_occ),dbl_mb(k_eval),1, Haa)) 
     &	call errquit('ebt_movecs_read failed for Haa',0, DISK_ERR)
		!-------------------------------------------------------------
		if (me.eq.0) then
			write(LuOut,9995) 'Electronic energy ','H(RR)',Haa
		endif
		!-------------------------------------------------------------
		! get nocc(1): # of occupied alpha mo's
		!--------------------------------------
		nocc(1)=0
   	do i = 1,nbf
    		nocc(1) = nocc(1) + dbl_mb(k_occ+i-1)
   	enddo
		nocc=nocc/2
    	if(me.eq.0) then
    		write(LuOut,*) nmo2(1), " occupied orbitals"
    		write(LuOut,*) nocc(1), " occupied alpha orbitals"
    		write(LuOut,*)
    	endif
		!-------------------------------------------------------------
		allocate(omo(nbf,          4))
		allocate(mo (nbf,          4))
		allocate(den(nbf, nbf, 4,  4))
		!-------------------------------------------------------------
		do i = 1, 4
			k=nocc(1)-2+i
			call ga_get(g_mo  , 1, nbf, k, k, omo(1,i), 1)
			call ga_put(g_omo4, 1, nbf, i, i, omo(1,i), 1)
		end do
		!-------------------------------------------------------------
		if (rotate_homo==1) then
			const   = 1.D0 / sqrt(2.D0)
			!--------------------------
			! a1
			mo(:,1) = (omo(:,1) + omo(:,2)) * const
			!--------------------------
			! b1
			mo(:,3) = (omo(:,1) - omo(:,2)) * const
			!--------------------------
		else
			!--------------------------
			! a1
			mo(:,1) = omo(:,1)
			!--------------------------
			! b1
			mo(:,3) = omo(:,2)
			!--------------------------
		end if
		!-------------------------------------------------------------
		if (rotate_lumo==1) then
			const   = 1.D0 / sqrt(2.D0)
			!--------------------------
			! a2
			mo(:,2) = (omo(:,3) + omo(:,4)) * const
			!--------------------------
			! b2
			mo(:,4) = (omo(:,3) - omo(:,4)) * const
			!--------------------------
		else
			!--------------------------
			! a2
			mo(:,2) = omo(:,3)
			!--------------------------
			! b2
			mo(:,4) = omo(:,4)
			!--------------------------
		end if
		!-------------------------------------------------------------
		allocate(den_a(nbf,nbf))
		allocate(den_b(nbf,nbf))
		!-------------------------------------------------------------
		do i=1, nbf
			do j=1, nbf
				den_a(i,j) = mo(i,1) * mo(j,2)
				den_b(i,j) = mo(i,3) * mo(j,4)
			end do
		end do
		do i = 1, nbf
			call ga_put(ga_den, 1, nbf, i, i, den_a(1,i), 1)
			call ga_put(gb_den, 1, nbf, i, i, den_b(1,i), 1)
		end do
		!-------------------------------------------------------------
		coup12=0.D0
		coup14=0.D0
		!-------------------------------------------------------------
		coc=0.D0
		!-------------------------------------------------------------
		if (.not.rtdb_get(rtdb,'ebt:natom1',MT_INT,1,natom1)) then
		endif
		!
		if (.not.rtdb_get(rtdb,'ebt:natom2',MT_INT,1,natom2)) then
		endif
		!-------------------------------------------------------------
		allocate(mol1(natom1))
		allocate(mol2(natom2))
		!-------------------------------------------------------------
		if (.not.rtdb_get(rtdb,'ebt:mol1',MT_INT,natom1,mol1)) then
		endif
		!
		if (.not.rtdb_get(rtdb,'ebt:mol2',MT_INT,natom2,mol2)) then
		endif
		!-------------------------------------------------------------
		call ebt_coc(rtdb, geom, 'geom', coc(:,1), natom1, mol1)
		call ebt_coc(rtdb, geom, 'geom', coc(:,2), natom2, mol2)
		!-------------------------------------------------------------
		deallocate(mol1)
		deallocate(mol2)
		!-------------------------------------------------------------
		call ebt_tran_dip (geom, ga_den   , basis, coc(1,1), dip(1,1))
		call ebt_tran_dip (geom, gb_den   , basis, coc(1,2), dip(1,2))
		!call ebt_tran_dip2(geom, den_a,nbf, basis, coc(1,1), dip(1,1))
		!call ebt_tran_dip2(geom, den_b,nbf, basis, coc(1,2), dip(1,2))
		!-------------------------------------------------------------
		call ebt_draw_xyz(rtdb, geom, 'geom', coc, dip)
		!-------------------------------------------------------------
		call ebt_dip_coupling(coc, dip, H, V, kappa, dr, theta)
		if (me==0) then
			write(LuOut,'(2x,"Dipole-dipole interaction"         )')
			write(LuOut,'(2x,"kappa = cos c - 3 cos a * cos b"   )')
			write(LuOut,*)
			write(LuOut,'(2x,"Center of charge    (a) = ",
     &		3f18.7, x, "Angstrom")') coc(:,1) * au2a
			write(LuOut,'(2x,"Center of charge    (b) = ",
     &		3f18.7, x, "Angstrom")') coc(:,2) * au2a
			!----------------------------------------------------------
			write(LuOut,'(2x,"R vector                = ",
     &		3f18.7, x, "Angstrom")') (coc(:,2)-coc(:,1))*au2a
			!----------------------------------------------------------
			write(LuOut,'(2x,"Trans dipole moment (a) = ",
     &		3f18.7, x, "a.u.")') dip(:,1)
			write(LuOut,'(2x,"Trans dipole moment (b) = ",
     &		3f18.7, x, "a.u.")') dip(:,2)
			!----------------------------------------------------------
			write(LuOut,'(2x,"alpha, beta, theta      = ",
     &		3f18.7, x, "Degree")')
     &		theta
			!----------------------------------------------------------
			write(LuOut,'(2x,"kappa                   = ",
     &		 f18.7               )') kappa
			write(LuOut,'(2x,"Distance                = ",
     &		 f18.7, x, "Angstrom")') dr*au2a
			write(LuOut,'(2x,"V ( |mu1||mu2|/R^3 )    = ",
     &		 f18.7, x, "cm^-1"   )') V * au2cm
			write(LuOut,'(2x,"H ( V.kappa)            = ",
     &		 f18.7, x, "cm^-1"   )') H * au2cm
			!----------------------------------------------------------
		end if
		!-------------------------------------------------------------
		deallocate(den_a)
		deallocate(den_b)
		!-------------------------------------------------------------
		do i=1, nbf
			coup12 = coup12 + abs(mo(i,1)) * abs(mo(i,2))
			coup14 = coup14 + abs(mo(i,1)) * abs(mo(i,4))
		end do
		if (coup12 < coup14) then
			do i=1, nbf
				rtmp    = mo(i,2)
				mo(i,2) = mo(i,4)
				mo(i,4) = rtmp
			end do
		end if
		!-------------------------------------------------------------
		do i=1, 4
			do j=1, 4
				!-------------------------------------------------------
				call ga_sync()
				write(gname,'("g_den_",i0,"_",i0)') i,j
		if(.not.ga_create(mt_dbl,nbf,nbf,trim(gname),0,0,g_den(j,i)))
     &	call errquit('task_et ga_create failed', g_den(i,j), GA_ERR)
				!-------------------------------------------------------
				do k=1, nbf
					do l=1, nbf
						den(l,k,j,i) = mo(l,j) * mo(k,i)
					end do ! l
				end do ! k
				!-------------------------------------------------------
				do k=1, nbf
					call ga_put(g_den(j,i),1,nbf,k,k,den(1,k,j,i),1)
				end do ! k
				!-------------------------------------------------------
			end do ! j
		end do ! i
		!-------------------------------------------------------------
		do i = 1, 4
			k=nocc(1)-2+i
			call ga_put(g_mo4, 1, nbf, i, i, mo(1,i), 1)
			call ga_put(g_mo , 1, nbf, k, k, mo(1,i), 1)
		end do
		!-------------------------------------------------------------
		! Write new mo
		!
		if (.not. movecs_write(rtdb, basis, f_movecs_out,
     $	'scf', title, nbf, nset, nmo2, dbl_mb(k_occ),
     $ nbf, dbl_mb(k_eval), nbf, g_mo)) call errquit
     $	('scf_movecs_write: movecs_write failed', 0, DISK_ERR)
		if(me==0) write(LuOut,'(2x, "Only output new mo file")')
		if(me==0) write(LuOut,'(2x, "Generate coupling matrix")')
		!-------------------------------------------------------------
		! Perform matrix multiplication to get Sab = Cma'*Smn*Cnb
		! Perform matrix multiplication to get Hab = Cma'*Hmn*Cnb
		!-------------------------------------------------------------
		call ga_sync()
		!-------------------------------------------------------------
		call ga_dgemm('T','N',4,nbf,nbf,1.d0,g_mo4,  g_s,0.d0,g_tmp)
		call ga_dgemm('N','N',4,  4,nbf,1.d0,g_tmp,g_mo4,0.d0,g_s4)
		!-------------------------------------------------------------
		call ga_sync()
		call ga_dgemm('T','N',4,nbf,nbf,1.d0,g_mo4,  g_h1,0.d0,g_tmp)
		call ga_dgemm('N','N',4,  4,nbf,1.d0,g_tmp,g_mo4,0.d0,g_h4)
		!-------------------------------------------------------------
		call ga_sync()
		!-------------------------------------------------------------
		do i = 1, 4
			call ga_get(g_s4,1, 4, i, i, sab(1,i), 1)
			call ga_get(g_h4,1, 4, i, i, hab(1,i), 1)
		end do
		!-------------------------------------------------------------
		if (util_print("ebt",print_debug)) then
			call ga_print(g_s    )
  			call ga_print(g_h1   )
  			call ga_print(g_mo   )
			call ga_print(g_omo4 )
			call ga_print(g_mo4  )
			call ga_print(g_s4   )
			call ga_print(g_h4   )
		endif
		!-------------------------------------------------------------
		time(1)= util_cpusec()
		if (output_mo/=1) then
			tol2e = 1.0d-12
			call int_mem_2e4c(mem2,max2e)
			mem2 = max(mem2,max2e)
			!----------------------------------------------------------
			if (me==0) write(LuOut,'(2x, "Two-Electron Integration")')
			if (method_2e==1) then
				call ebt_2e(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &			g_den, Jab)
			else
				call ebt_2e2(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &			g_den, den, Jab)
			end if
			!----------------------------------------------------------
			call ga_sync()
			!----------------------------------------------------------
			call full_Jab(Jab)
			!-------------------------------------------------------------
			call ga_sync()
			!-------------------------------------------------------------
			call ebt_trans(sab,hab,Jab)
			!-------------------------------------------------------------
			if(me==0) then
				call ebt_output(Jab)
				open(fid, file=outfile)
					write(fid, '("Sab")')
					do i=1, 4
						write(fid, '(4es20.10)') (sab(i,j), j=1, 4)
					end do
					write(fid, '("Hab")')
					do i=1, 4
						write(fid, '(4es20.10)') (hab(i,j), j=1, 4)
					end do
					write(fid, '("Jab")')
					do i=1, 4
					do j=1, 4
					do k=1, 4
					do l=1, 4
						write(fid, '(4i5, es20.10)') i,j,k,l,Jab(i,j,k,l)
					end do
					end do
					end do
					end do
				close(fid)
			end if
		end if
		!----------------------------------------------------------------
		time (2) = util_cpusec()
		dtime(1) = time(2) - time(1)
		if(me==0) then
			write(LuOut, '("jtime", f15.7)') dtime(1)
		end if
		!----------------------------------------------------------------
		deallocate(omo)
		deallocate(mo )
		deallocate(den)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_tmp))
     &	call errquit('could not destroy tmp array',g_tmp, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(ga_den))
     &	call errquit('could not destroy ga_den array',ga_den, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(gb_den))
     &	call errquit('could not destroy gb_den array',gb_den, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_mo))
     &	call errquit('could not destroy g_mo array',g_mo, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_omo4))
     &	call errquit('could not destroy g_omo4 array',g_omo4, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_mo4))
     &	call errquit('could not destroy g_mo4 array',g_mo4, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_s)) 
     &	call errquit('could not destroy overlap array',g_s, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_s4))
     &	call errquit('could not destroy g_s4 array',g_s4, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ga_destroy(g_h4))
     &	call errquit('could not destroy g_h4 array',g_h4, GA_ERR)
		!-------------------------------------------------------------
		if (.not.ma_chop_stack(l_eval)) 
     &	call errquit('unable to destroy MA',l_eval, MA_ERR)
		if (.not.ma_chop_stack(l_occ)) 
     &	call errquit('unable to destroy MA',l_occ , MA_ERR)
		!-------------------------------------------------------------
		!
		ebt_calc=.true.
		call util_print_pop
		!
		!-------------------------------------------------------------
		if (.not. geom_destroy(geom))
     &	call errquit('ebt_calc: geom destroy failed', 0, GEOM_ERR)
		if (.not. bas_destroy(basis))
     &	call errquit('ebt_calc: basis destroy failed',0, BASIS_ERR)
		!-------------------------------------------------------------
		call int_terminate()
 9995 format (a34,1x,a6,1x,F18.10)
		!-------------------------------------------------------------
		return
		end
