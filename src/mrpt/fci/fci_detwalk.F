*
* $Id: fci_detwalk.F 19708 2010-10-29 18:04:21Z d3y133 $
*
      subroutine fci_det_to_walk(idet, owalk, mwalk)
      implicit none
#include "fciinfo.fh"
c
c     Given the index of a determinant, idet, return the
c     corresponding lexical weights of orbital and primitive
c     spin functions.  Without having an index vector the
c     length of the CI expansion the only way to do this
c     is with a binary search, so this operation costs 
c     about log2(ndets).
c
      integer idet    ! [input] Index into determinant CI vector
      integer owalk   ! [output] Corresponding orbital walk
      integer mwalk   ! [output] Corresponding primitive spin function walk
c
      integer left, right, middle
c
      left = 1
      right = nwalko
c
 10   middle = (left+right)/2
      if (indexo(middle) .ge. idet) then
         right = middle
      else 
         left = middle
      endif
      if ((right-left) .gt. 1) goto 10
c
      if (indexo(right) .lt. idet) then
         owalk = right
      else
         owalk = left
      endif
c
c     Adjust for walks of incorrect symmetry
c
 20   if (indexo(owalk).eq.indexo(owalk+1)) then
         owalk = owalk + 1
         goto 20
      endif
c
      mwalk = idet - indexo(owalk)
c
      if (mwalk.le.0 .or. mwalk.gt.(indexo(owalk+1)-indexo(owalk)))
     $     call fci_err('fci_det_to_walk: bad mwalk ', mwalk)
c
      end

