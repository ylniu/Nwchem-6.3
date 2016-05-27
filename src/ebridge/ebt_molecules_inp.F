      subroutine ebt_mols_input(rtdb,module)
C     $Id: ebt_vectors_inp.F 19708 2010-10-29 18:04:21Z d3y133 $
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      character*(*) module
c     
c     Copied from vectors_input. This program is used to read
c     the movecs filename for the MO vectors for the reactants.
c     
c     Input vectors must be specified as:
c     vectors reactants filename products filename
c     
      integer mlen, lena, lenb
      character*255 test, vectors_in, vectors_out
      character*255 a_vecs, b_vecs
      integer num_dirs, natom1, natom2, nd, ind
      integer, allocatable :: mol1(:), mol2(:)
      parameter (num_dirs = 2)
      character*9 dirs(num_dirs)
      data dirs /'reactants', 'products'/
c     
      vectors_in = ' '
      vectors_out = ' '
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 10000
c-------------------------------------------------------------------------------
c  read 'reactants' keyword signifying that A state vectors are to be read next:
c
      if (inp_compare(.false.,'reactants', test)) then
c-------------------------------------------------------------------------------
      if (.not. inp_a(test)) goto 10000
      call get_natom (natom1,nd,trim(test)      )
      allocate(mol1(natom1))
      call get_mol   (natom1,nd,trim(test), mol1)
      if (.not. rtdb_put(rtdb, 'ebt:natom1', MT_INT, 1, natom1))
     $  call errquit('ebt_input: failed to store natom1',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'ebt:mol1', MT_INT, natom1, mol1))
     $  call errquit('ebt_input: failed to store nmol1',0,RTDB_ERR)
      deallocate(mol1)
c-------------------------------------------------------------------------------
c  read movecs filename for the products
c
      else if (inp_compare(.false.,'products', test)) then
c-------------------------------------------------------------------------------
      if (.not. inp_a(test)) goto 10000
      call get_natom (natom2,nd,trim(test)      )
      allocate(mol2(natom2))
      call get_mol   (natom2,nd,trim(test), mol2)
      if (.not. rtdb_put(rtdb, 'ebt:natom2', MT_INT, 1, natom2))
     $  call errquit('ebt_input: failed to store natom2',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'ebt:mol2', MT_INT, natom2, mol2))
     $  call errquit('ebt_input: failed to store nmol2',0,RTDB_ERR)
      deallocate(mol2)
c-------------------------------------------------------------------------------
      else 
         write(6,10001)
         goto 10000
      endif
      return
c
c     
10001 format(/' molecules [reactants] filename [products] filename  '/)
c
c
10000 call errquit('ebt_vectors_input: invalid format', 0, INPUT_ERR)
 999  continue
      call util_flush(6)
c     
      return
      end

