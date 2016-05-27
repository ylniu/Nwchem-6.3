		subroutine full_Jab(Jab)
			implicit none
			real(8) :: Jab(4,4,4,4)
			integer :: i, j, k, l, maxl
			real(8) :: a
			!-------------------------------------------------------------
			do i=4,1,-1
			do j=i,1,-1
			do k=i,1,-1
				maxl=k
				if (k==i) maxl=j
			do l=maxl,1,-1
				a=Jab(i,j,k,l)
				Jab(i,j,l,k) = a
				Jab(j,i,k,l) = a
				Jab(j,i,l,k) = a
				Jab(k,l,i,j) = a
				Jab(l,k,i,j) = a
				Jab(k,l,j,i) = a
				Jab(l,k,j,i) = a
			end do
			end do
			end do
			end do
			!-------------------------------------------------------------
			return
		end subroutine full_Jab
