      subroutine tce_hash_ind(k_offset_a,key,indexc)
cccx      subroutine tce_hash_ind(hash,key,indexc)
cccx      subroutine tce_hash_ind(hash,key,offset,indexc)
c
c $Id: tce_hash_ind.F 19706 2010-10-29 17:52:31Z d3y133 $
c
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "tce.fh"
      integer key
      integer length
      integer i
      integer less,middle,more
      integer indexc
      integer k_offset_a
c
      length=int_mb(k_offset_a)
      less = 1
      more = length
      indexc = 1
      return
      end
