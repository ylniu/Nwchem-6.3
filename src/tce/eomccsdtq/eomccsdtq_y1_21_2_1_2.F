      SUBROUTINE eomccsdtq_y1_21_2_1_2(d_a,k_a_offset,d_c,k_c_offset)
C     $Id: eomccsdtq_y1_21_2_1_2.F 24183 2013-05-03 23:47:08Z edo $
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i3 ( h5 h6 h9 h14 p1 p3 p7 p8 )_y + = 4/5 * y ( h5 h6 h9 h14 p1 p3 p7 p8 )_y
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER h5b
      INTEGER h6b
      INTEGER h9b
      INTEGER h14b
      INTEGER p1b
      INTEGER p3b
      INTEGER p7b
      INTEGER p8b
      INTEGER dimc
      INTEGER h5b_1
      INTEGER h6b_1
      INTEGER h9b_1
      INTEGER h14b_1
      INTEGER p1b_1
      INTEGER p3b_1
      INTEGER p7b_1
      INTEGER p8b_1
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_c
      INTEGER k_c
      EXTERNAL NXTASK
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs,1)
      DO h5b = 1,noab
      DO h6b = h5b,noab
      DO h9b = h6b,noab
      DO h14b = h9b,noab
      DO p1b = noab+1,noab+nvab
      DO p3b = noab+1,noab+nvab
      DO p7b = p3b,noab+nvab
      DO p8b = p7b,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-1
     &)+int_mb(k_spin+h9b-1)+int_mb(k_spin+h14b-1)+int_mb(k_spin+p1b-1)+
     &int_mb(k_spin+p3b-1)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1).ne.
     &16)) THEN
      IF (int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-1)+int_mb(k_spin+h9b-1)
     &+int_mb(k_spin+h14b-1) .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p3b
     &-1)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)) THEN
      IF (ieor(int_mb(k_sym+h5b-1),ieor(int_mb(k_sym+h6b-1),ieor(int_mb(
     &k_sym+h9b-1),ieor(int_mb(k_sym+h14b-1),ieor(int_mb(k_sym+p1b-1),ie
     &or(int_mb(k_sym+p3b-1),ieor(int_mb(k_sym+p7b-1),int_mb(k_sym+p8b-1
     &)))))))) .eq. irrep_y) THEN
      dimc = int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1) * int_mb(k_ra
     &nge+h9b-1) * int_mb(k_range+h14b-1) * int_mb(k_range+p1b-1) * int_
     &mb(k_range+p3b-1) * int_mb(k_range+p7b-1) * int_mb(k_range+p8b-1)
      CALL TCE_RESTRICTED_8(h5b,h6b,h9b,h14b,p1b,p3b,p7b,p8b,h5b_1,h6b_1
     &,h9b_1,h14b_1,p1b_1,p3b_1,p7b_1,p8b_1)
      dim_common = 1
      dima_sort = int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1) * int_mb
     &(k_range+h9b-1) * int_mb(k_range+h14b-1) * int_mb(k_range+p1b-1) *
     & int_mb(k_range+p3b-1) * int_mb(k_range+p7b-1) * int_mb(k_range+p8
     &b-1)
      dima = dim_common * dima_sort
      IF (dima .gt. 0) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('eomccsdtq_y1_21_2_1_2',0,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &eomccsdtq_y1_21_2_1_2',1,MA_ERR)
      IF ((p8b .lt. p1b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(p1b_1
     & - noab - 1 + nvab * (p8b_1 - noab - 1 + nvab * (p7b_1 - noab - 1 
     &+ nvab * (p3b_1 - noab - 1 + nvab * (h14b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h9b-1),int_mb(k_range+h14b-1
     &),int_mb(k_range+p3b-1),int_mb(k_range+p7b-1),int_mb(k_range+p8b-1
     &),int_mb(k_range+p1b-1),7,6,5,8,4,3,2,1,-1.0d0)
      END IF
      IF ((p7b .lt. p1b) .and. (p1b .le. p8b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(p8b_1
     & - noab - 1 + nvab * (p1b_1 - noab - 1 + nvab * (p7b_1 - noab - 1 
     &+ nvab * (p3b_1 - noab - 1 + nvab * (h14b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h9b-1),int_mb(k_range+h14b-1
     &),int_mb(k_range+p3b-1),int_mb(k_range+p7b-1),int_mb(k_range+p1b-1
     &),int_mb(k_range+p8b-1),8,6,5,7,4,3,2,1,1.0d0)
      END IF
      IF ((p3b .lt. p1b) .and. (p1b .le. p7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(p8b_1
     & - noab - 1 + nvab * (p7b_1 - noab - 1 + nvab * (p1b_1 - noab - 1 
     &+ nvab * (p3b_1 - noab - 1 + nvab * (h14b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h9b-1),int_mb(k_range+h14b-1
     &),int_mb(k_range+p3b-1),int_mb(k_range+p1b-1),int_mb(k_range+p7b-1
     &),int_mb(k_range+p8b-1),8,7,5,6,4,3,2,1,-1.0d0)
      END IF
      IF ((p1b .le. p3b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(p8b_1
     & - noab - 1 + nvab * (p7b_1 - noab - 1 + nvab * (p3b_1 - noab - 1 
     &+ nvab * (p1b_1 - noab - 1 + nvab * (h14b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h9b-1),int_mb(k_range+h14b-1
     &),int_mb(k_range+p1b-1),int_mb(k_range+p3b-1),int_mb(k_range+p7b-1
     &),int_mb(k_range+p8b-1),8,7,6,5,4,3,2,1,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('eomccsdtq_y1_21_2_1_2',2
     &,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &eomccsdtq_y1_21_2_1_2',3,MA_ERR)
      CALL TCE_SORT_8(dbl_mb(k_a_sort),dbl_mb(k_c),int_mb(k_range+p8b-1)
     &,int_mb(k_range+p7b-1),int_mb(k_range+p3b-1),int_mb(k_range+p1b-1)
     &,int_mb(k_range+h14b-1),int_mb(k_range+h9b-1),int_mb(k_range+h6b-1
     &),int_mb(k_range+h5b-1),8,7,6,5,4,3,2,1,4.0d0/5.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p8b -
     & noab - 1 + nvab * (p7b - noab - 1 + nvab * (p3b - noab - 1 + nvab
     & * (p1b - noab - 1 + nvab * (h14b - 1 + noab * (h9b - 1 + noab * (
     &h6b - 1 + noab * (h5b - 1)))))))))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('eomccsdtq_y1_21_2_1_2',4
     &,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('eomccsdtq_y1_21_2_1
     &_2',5,MA_ERR)
      END IF
      END IF
      END IF
      END IF
      next = NXTASK(nprocs,1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs,1)
      call GA_SYNC()
      RETURN
      END
