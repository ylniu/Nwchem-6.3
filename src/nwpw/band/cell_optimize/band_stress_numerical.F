*
* $Id: band_stress_numerical.F 19707 2010-10-29 17:59:36Z d3y133 $
*

      logical function band_stress_numerical(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "mafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "rtdb.fh"
#include "geom.fh"

*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,value2

      real*8 stress(3,3)
      real*8 unita(3,3)

*     **** external functions ****
      logical  band_minimizer,band_stress_fd
      external band_minimizer,band_stress_fd

      taskid = ga_nodeid()


      if (taskid.eq.MASTER) then
        write(*,*) '  ******************************************** '
        write(*,*) '  *                                          * '
        write(*,*) '  *    BAND Stress Numerical - version 1.0   * '
        write(*,*) '  *                                          * '
        write(*,*) '  *      version last updated 4/24/02        * '
        write(*,*) '  *                                          * '
        write(*,*) '  ******************************************** '
        call nwpw_message(1)
      end if


*     **** minimize psi ***
      value = band_minimizer(rtdb,1)
      if(.not. value) 
     > call errquit(
     > 'pspw_stress_numerical:band_gradient did not converge: error',
     >  0,CALC_ERR)

*     **** calculate stress tensor ***
      value = .true.
      value = value.and.rtdb_parallel(.true.)
      value = value.and.band_stress_fd(rtdb)
      if(.not.value) call errquit(
     >   'band_stress_numerical failed: error',0,CALC_ERR)
 
      call ga_sync()

      band_stress_numerical = value
      return
      end
 

