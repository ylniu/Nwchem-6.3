      subroutine ebt_draw_xyz(rtdb, geom, name, coc, dip)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      !-------------------------------------------------------------------------
      integer rtdb              !< [Input] the RTDB handle
      integer geom
      integer fid
      double precision coc(3,2)
      double precision dip(3,2)
      real(8) :: au2a=0.529177D0
      character*(*) name        !< [Input] the name of the geometry
      character*256 tmp  , file_prefix, fxyz
      !-------------------------------------------------------------------------
      integer i, j, k
      logical s, geom_rtdb_load
      double precision tcha
      !-------------------------------------------------------------------------
      fid=38
      names(geom) = name
      lenn(geom) = inp_strlen(name)
      trans(geom) = 'junk'
      call util_file_prefix(' ',file_prefix)
      if (.not. rtdb_cget(rtdb, name, 1, trans(geom)))
     $     trans(geom) = name
      !-------------------------------------------------------------------------
      lent(geom) = inp_strlen(trans(geom))
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      s = .true.
c     
      tmp(k:) = ' ' 
      tmp(k:) = ':ncenter'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, ncenter(geom))
      geom_rtdb_load = s

      !-------------------------------------------------------------------------
      tmp(k:) = ' '
      tmp(k:) = ':coords'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent*3, 
     $     coords(1,1,geom))
      !-------------------------------------------------------------------------
      tmp(k:) = ' ' 
      tmp(k:) = ':tags'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent, tags(1,geom))
      !-------------------------------------------------------------------------
		fxyz=trim(adjustl(file_prefix))//"dip.xyz"
		open(fid, file=fxyz)
			write(fid,'(i4)') ncenter(geom) + 2
			write(fid,'("title")')
			do i=1,ncenter(geom)
				write(fid,'(a4,6f16.10)')
     &			trim(tags(i,geom)), (coords(j,i,geom)*au2a,j=1,3),(0.D0,j=1,3)
			end do
			do i=1, 2
				write(fid,'(a4,6f16.10)')
     &			"X", (coc(j,i)*au2a,j=1,3), (dip(j,i),j=1,3)
			end do
		close(fid)
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      return
      end subroutine
