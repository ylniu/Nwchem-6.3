      subroutine ebt_vecs_input(rtdb,module)
C     $Id: ebt_vectors_inp.F 19708 2010-10-29 18:04:21Z d3y133 $
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer rtdb
      character*(*) module
c     
c     Copied from vectors_input. This program is used to read
c     the movecs filename for the MO vectors for the outfilename.
c     
c     Input vectors must be specified as:
c     vectors outfilename filename products filename
c     
      integer mlen, lena, lenb
      character*255 test, vectors_in, vectors_out
      character*255 a_vecs, b_vecs
      integer num_dirs, ind
      parameter (num_dirs = 2)
      character*9 dirs(num_dirs)
      data dirs /'outfilename'/
c     
      vectors_in = ' '
      vectors_out = ' '
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 10000
c     
c  read 'outfilename' keyword signifying that A state vectors are to be read next:
c
      if (inp_compare(.false.,'output', test)) then
c
      if (.not. inp_a(test)) goto 10000
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
c     
c     Does not match a keyword ... assume it is the name of the
c     file with the input vectors
c     
c         if (vectors_in .ne. ' ') goto 10000
         lena = inp_strlen(test)
         a_vecs = test(1:lena)
         if (.not. rtdb_cput(rtdb, 'a_vecs', 1, a_vecs)) call errquit
     $        ('ebt_vectors_input: failed to store vectors_input', 0,
     &       RTDB_ERR)
      endif
c
      else 
         write(6,10001)
         write(6,10003)
         a_vecs = 'movecsa.mo'
         if (.not. rtdb_cput(rtdb, 'a_vecs', 1, a_vecs)) call errquit
     $        ('ebt_vectors_input: failed to store vectors_input', 0,
     &       RTDB_ERR)
         b_vecs = 'movecsb.mo'
         if (.not. rtdb_cput(rtdb, 'b_vecs', 1, b_vecs)) call errquit
     $        ('ebt_vectors_input: failed to store vectors_input', 0,
     &       RTDB_ERR)
      endif
      return
c
c     
10001 format(/' vectors [outfilename] filename [products] filename  '/)
10003 format(/' default outfilename and products MO vector filenames: '
     $     /'         movecsa.mo and movecsb.mo    '/)
c
c
10000 call errquit('ebt_vectors_input: invalid format', 0, INPUT_ERR)
 999  continue
      call util_flush(6)
c     
      return
      end

