      subroutine argos_cafe_print_dera(lfnout,number)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfnout,number
c
      call argos_cafe_print_der(lfnout,number,byte_mb(i_snam),
     + dbl_mb(i_sti))
c
      return
      end
c $Id: argos_cafe_print_dera.F 21176 2011-10-10 06:35:49Z d3y133 $
