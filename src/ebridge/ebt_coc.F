      subroutine ebt_coc(rtdb, geom, name, coc, natom, mol)
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
      integer natom
      integer mol(natom)
      double precision coc(3)
      character*(*) name        !< [Input] the name of the geometry
      character*256 tmp
      !-------------------------------------------------------------------------
      integer i, j, k
      logical s, geom_rtdb_load
      double precision tcha
      !-------------------------------------------------------------------------
      names(geom) = name
      lenn(geom) = inp_strlen(name)
      trans(geom) = 'junk'
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
      tmp(k:) = ' ' 
      tmp(k:) = ':vel'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent*3, 
     $     velocities(1,1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':charges'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent, charge(1,geom))
      !-------------------------------------------------------------------------
      tcha = 0.D0
      coc  = 0.D0
      do k=1, natom
         i=mol(k)
        tcha = tcha + charge(i,geom)
        do j=1, 3
          coc(j) = coc(j) + charge(i,geom) * coords(j,i,geom)
        end do
      end do
      coc = coc / tcha
      !-------------------------------------------------------------------------
      return
      end subroutine
