      subroutine eet_print(g)
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
#include "stdio.fh"
c ------------------------------------------------------------------------------
      character(255) g_name
      integer g, g_type
      integer ni, nv, mp
      integer me
      integer i, j, k, l
      integer la, ka
      integer p, q
      character*20 b(6)
c ------------------------------------------------------------------------------
      me =  ga_nodeid()
      call ga_sync()
      if (me==0) then
        mp = 6
        call ga_inquire(g, g_type, ni, nv)
        call ga_inquire_name(g, g_name)
        ! ----------------------------------------------------------------------
        if (.not.ma_push_get(mt_dbl,ni*nv,'a',la,ka))
     $    call errquit('ma_push_get failed for ka',ka, MA_ERR)
        ! ----------------------------------------------------------------------
        write(LuOut,*)
        write(LuOut,*) "----------------------------------------"
        write(LuOut,*) "PRINT: ", trim(g_name)
        write(LuOut,*) "----------------------------------------"
        p = nv/mp
        q = mod(nv,mp)
        b = "-------------------"

        do i=1, nv
          j=ka + (i-1) * ni
          call ga_get(g, 1, ni, i, i, dbl_mb(j), ni)
        end do
        do i=1, p
          k=(i-1) * mp
          write(LuOut,'(10x, 6i20    )') (k + j, j=1, mp)
          write(LuOut,'(10x, 6(x,a19))') (b(j) , j=1, mp)
          do l=1, ni
            write(LuOut,'(i8, 2x, 6es20.12)') l,
     $        (dbl_mb(ka+( (i-1)*mp+(j-1) )*ni + l-1), j=1, mp)
          end do
          write(LuOut,*)
        end do
        if (q.ne.0) then
          i=p+1
          k=(i-1) * mp
          write(LuOut,'(10x, 6i20    )') (k + j, j=1, q)
          write(LuOut,'(10x, 6(x,a19))') (b(j) , j=1, q)
          do l=1, ni
            write(LuOut,'(i8, 2x, 6es20.12)') l,
     $        (dbl_mb(ka+( (i-1)*mp+(j-1) )*ni + l-1), j=1, q)
          end do
          write(LuOut,*)
        end if
c ------------------------------------------------------------------------------
      if (.not.ma_pop_stack(la))
     &  call errquit('unable to destroy MA la', la, MA_ERR)
c ------------------------------------------------------------------------------
      end if
c ------------------------------------------------------------------------------
      call ga_sync()
c ------------------------------------------------------------------------------
      return
      end subroutine
