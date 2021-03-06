      subroutine argos_diana_getref
c
c $Id: argos_diana_getref.F 19708 2010-10-29 18:04:21Z d3y133 $
c
c     read the reference coordinates
c     ------------------------------
c
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      filref=card(8:80)
c
c     get the size from the reference file
c
      call argos_diana_sizref(nsai,msai,nwmi,mwmi,nwai,mwai)
c
c     allocate memory and initialize
c
      call argos_diana_alloc(.true.)
c
c     read the reference coordinates
c
      call argos_diana_rdref(dbl_mb(i_xref))
c
      call argos_diana_all(int_mb(i_ssel),1,dbl_mb(i_swt),nsel)
      lsel=.false.
      call argos_diana_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      call argos_diana_wall(int_mb(i_wsel),1,int_mb(i_ndxw))
c
      if(me.eq.0) then
      write(*,'(a,a)') ' Reference coordinates read from ',
     + filref(1:index(filref,' ')-1)
      write(*,'(/,a,i5,/)') ' Number of atoms is ',nsai
      endif
      if(.not.ma_verify_allocator_stuff()) print*,'Oopsi 2'
c
      return
      end
