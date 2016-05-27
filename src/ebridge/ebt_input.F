c $Id: ebt_input.F 19708 2010-10-29 18:04:21Z d3y133 $
      subroutine ebt_input(rtdb)
      implicit none
#include "errquit.fh"
c
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c     
      integer rtdb, method_2e, output_mo, rotate_homo, rotate_lumo
      character*255 test, outfile
      double precision tol2e
c
 100  if (.not. inp_read()) call errquit('ebt_input: inp_read failed',0,
     &       INPUT_ERR)
c
       if (.not. inp_a(test))
     $     call errquit('scf_input: failed to read keyword', 0,
     &       INPUT_ERR)
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "ebt")
c
      goto 100
c
c     MO vector filenames are read
c     ----------------------------
      else if (inp_compare(.false.,'vectors', test)) then
         call ebt_vecs_input(rtdb,'ebt')
      goto 100

      else if (inp_compare(.false.,'molecules', test)) then
         call ebt_mols_input(rtdb,'ebt')
      goto 100
c
c     screening tolerances are read
c     -----------------------------
      else if (inp_compare(.false.,'tol2e', test)) then
         if (.not. inp_f(tol2e)) call errquit
     $        ('ebt_input: tol2e <thresh>',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'ebt:tol2e', MT_DBL, 1, tol2e))
     $        call errquit('ebt_input: rtdb_put of tol2e failed',0,
     &       RTDB_ERR)
c 
      goto 100
c
c     output
c     -----------------------------
      else if (inp_compare(.false.,'output', test)) then
			if (inp_a(test)) then
				if (.not. rtdb_cput(rtdb, 'outfile', 1, outfile)) call errquit
     $			('ebt_outfile: failed to store outfile', 0,
     &			RTDB_ERR)
			end if
      goto 100
c
c     rotate
c     ---------------------------
      else if (inp_compare(.false.,'rotate', test)) then
			rotate_homo=1
			rotate_lumo=1
			if (inp_a(test)) then
				if (trim(test)=="homo") then
					rotate_homo=1
					rotate_lumo=0
				else if (trim(test)=="lumo") then
					rotate_homo=1
					rotate_lumo=0
				else if (trim(test)=="all" ) then
					rotate_homo=1
					rotate_lumo=1
				else if (trim(test)=="none") then
					rotate_homo=0
					rotate_lumo=0
				else
					write(LuOut,'("bet_input:
     $				 rotate should be homo, lumo, all, none!")')
				end if
			end if
			if (.not. rtdb_put(rtdb, 'ebt:rotate_homo', MT_INT, 1, rotate_homo))
     $		call errquit('ebt_input: failed to store rotate_homo', 0,
     &		RTDB_ERR)
				
			if (.not. rtdb_put(rtdb, 'ebt:rotate_lumo', MT_INT, 1, rotate_lumo))
     $		call errquit('ebt_input: failed to store rotate_lumo', 0,
     &		RTDB_ERR)
				
c
      goto 100
c
c     output_mo_only
c     ---------------------------
      else if (inp_compare(.false.,'output_mo_only', test)) then
			output_mo=1
         if (.not. rtdb_put(rtdb,'ebt:output_mo', 
     $     MT_INT, 1, output_mo))
     $     call errquit('ebt_input: rtdb_put of
     $     output_mo failed',0,RTDB_ERR)
c
      goto 100
c
c     choice of 2e method is read
c     ---------------------------
      else if (inp_compare(.false.,'method_2e', test)) then
        method_2e = 1
         if (.not. inp_i(method_2e)) call errquit
     $        ('ebt_input: method_2e ',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'ebt:method_2e', MT_INT, 1, method_2e))
     $        call errquit('ebt_input: rtdb_put of method_2e failed',0,
     &       RTDB_ERR)
c
      goto 100
c
      else if (inp_compare(.false.,'method_b', test)) then
        method_2e = 2
c
         if (.not. rtdb_put(rtdb,'ebt:method_2e', MT_INT, 1, method_2e))
     $        call errquit('ebt_input: rtdb_put of method_2e failed',0,
     &       RTDB_ERR)
c
      goto 100
c
c     end of input block
c     -------------------
      else if (.not. inp_compare(.false.,'end', test)) then
       continue
      endif
c     
c
      return
      end
