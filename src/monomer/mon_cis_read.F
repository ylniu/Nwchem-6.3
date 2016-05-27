      logical function mon_cis_read_header(filename, g_vars)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
c
      character*(*) filename
      integer g_vars, unitno
      integer nroots, nocc, nmo, nci
      parameter (unitno = 67)
c
		mon_cis_read_header = .true.
      call ga_sync()
      if (ga_nodeid() .eq. 0) then
			open(unitno, file=filename, status="old", err=1000)
			read(unitno,*, err=2000) ! ltransden F
			read(unitno,*, err=2000) ! ipol = 1
			read(unitno,*, err=2000) nroots
			read(unitno,*, err=2000) nocc
			read(unitno,*, err=2000) nmo
			read(unitno,*, err=2000) ! nfc =  0
			read(unitno,*, err=2000) ! nfv =  0
			read(unitno,*, err=2000) nci
			read(unitno,*, err=2000) ! enter (space line)
			close(unitno)
			call ga_put(g_vars, 1, 1, 1, 1, nroots, 1)
			call ga_put(g_vars, 2, 2, 1, 1, nocc  , 1)
			call ga_put(g_vars, 3, 3, 1, 1, nmo   , 1)
			call ga_put(g_vars, 4, 4, 1, 1, nci   , 1)
		end if
		call ga_sync()
		return
		
 1000 write(6,*) ' mon_cis_read_header: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      mon_cis_read_header = .false.
      close(unitno)
      return
 2000 write(6,*) ' mon_cis_read_header: failed to read ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      mon_cis_read_header = .false.
      close(unitno)
      return
      end
!
!
!===============================================================================
!
!
      logical function mon_cis_read(filename, jroot,
     $  nmo, nv, nocc, g_vars, g_ci, g_zai, g_mo)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
c
      character*(*) filename
      integer jroot, g_vars, g_ci, g_zai, unitno
      integer nci, iroot, nroots, inci, nocc, nv
      integer i, j, k, me
      integer ncoeff, idx_nc
      integer ii, aa, nmo, g_mo
      real(8) :: mos(13,13)
      character*255 tmp, line
      logical :: read_gaussian
c ------------------------------------------------------------------
		read_gaussian = .false.
		mon_cis_read = .true.
		me = ga_nodeid()
      if (me .eq. 0) then
			open(unitno, file=filename, status="old", err=1000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000)
			read(unitno,*, err=2000) nci
			read(unitno,*, err=2000)
			call ga_get(g_vars, 1, 1, 1, 1, nroots, 1)
			call ga_get(g_vars, 4, 4, 1, 1, nci   , 1)
			i=nci * 2 * nroots
			if (.not.ma_push_get(mt_dbl,i,'cis', ncoeff,idx_nc))
     $		call errquit('cis_write: cis failed', 1, MA_ERR)
			i=-1
			do iroot=1, nroots
				read(unitno,*, err=2000)
				read(unitno,*, err=2000)
				do inci=1, nci * 2
					i=i+1
					read(unitno,*, err=2000) dbl_mb(idx_nc+i)
C					??? need to be checked
					dbl_mb(idx_nc+i) = dbl_mb(idx_nc+i) / sqrt(2.D0)
				end do
			end do
			close(unitno)
			! ---------------------------------------------------------------------
			do iroot=1, nroots
				j = idx_nc + (iroot-1)*nci*2
				call ga_put(g_ci,1,nci*2,iroot,iroot,dbl_mb(j),1)
				if (iroot==jroot) then
					do i=1, nocc
						k=j + (i-1) * nv
C					??? need to be checked
						call ga_put(g_zai,1,nv,i,i,dbl_mb(k),1)
					end do
				end if
			end do
			! ---------------------------------------------------------------------
			if (read_gaussian) then
			dbl_mb(idx_nc:idx_nc+nocc*nv-1) = 0.D0
			open(97,file="coupling-10.0.out")
				do i=1,7
					read(97,*)
				end do
				do i=1, 32
					read(97,'(a)') line
					read(line,*) ii, aa
					aa = aa - 5
					j=idx_nc+(ii-1) * nv + aa -1
					read(line,*) (tmp, k=1, 13), dbl_mb(j)
				end do
				do i=1, nocc
					j=idx_nc+(i-1) * nv
					call ga_put(g_zai,1,nv,i,i,dbl_mb(j),1)
				end do
				!-------------------------------------------------------------------
				do while (line(2:4)/="MOs")
					read(97,'(a)') line
				end do
				do i=1, 13
					read(97,*) tmp, (mos(i,j), j=1, 13)
				end do
				do i=1, 13
					call ga_put(g_mo,1,13,i,i,mos(1,i),1)
				end do
			close(97)
			end if
c ------------------------------------------------------------------------------
      if (.not.ma_chop_stack(ncoeff))
     &  call errquit('unable to destroy MA ncoeff',
     &  ncoeff, MA_ERR)
c ------------------------------------------------------------------------------
			! ---------------------------------------------------------------------
		end if
		call ga_sync()
		return
		
 1000 write(6,*) ' mon_cis_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      mon_cis_read = .false.
      close(unitno)
      return
 2000 write(6,*) ' mon_cis_read: failed to read ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      mon_cis_read = .false.
      close(unitno)
      return
      end
