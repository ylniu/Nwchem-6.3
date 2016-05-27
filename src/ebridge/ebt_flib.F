	subroutine get_natom(natom_of_list,nd,word)
	implicit none
	integer                    :: natom_of_list
	character(*)               :: word
	character(len(word))       :: word1
	integer                    :: natom, ia, ib, nd
	integer      , allocatable :: n(:)
	character(80)              :: ca, cb
	character(80), allocatable :: b(:)
	!----------------------------------------------------------------------------
	integer              :: la, lw, i, ind, j, k
	!----------------------------------------------------------------------------
	word1=trim(adjustl(word))
	!----------------------------------------------------------------------------
	i=1
	do while (.true.)
		ind = index(word1,",")
		lw  = len(trim(word1))
		if(ind>0) then
			i=i  + 1
			j=lw - ind
			word1(1:j)  = word1(ind+1:lw)
			word1(j+1:) = ""
		else
			exit
		end if
	end do
	!----------------------------------------------------------------------------
	nd = i
	!----------------------------------------------------------------------------
	allocate(b(i))
	allocate(n(i))
	!----------------------------------------------------------------------------
	word1=trim(adjustl(word))
	!----------------------------------------------------------------------------
	do k=1, i-1
		ind  = index(word1,",")
		lw   = len(trim(word1))
		b(k) = trim(word1(1:ind-1))
		j    = lw - ind
		word1(1:j)  = word1(ind+1:lw)
		word1(j+1:) = ""
	end do
	!----------------------------------------------------------------------------
	b(i) = trim(adjustl(word1))
	!----------------------------------------------------------------------------
	natom=0
	do k=1, i
		ind  = index(b(k),"-")
		if (ind==0) then
			n(k)=1
		else
			ca=trim(b(k)(1:ind-1))
			cb=trim(b(k)(ind+1:len(b(k))))
			read(ca,*) ia
			read(cb,*) ib
			n(k) = ib - ia + 1
		end if
		natom = natom + n(k)
	end do
	!----------------------------------------------------------------------------
	natom_of_list = natom
	!----------------------------------------------------------------------------
	deallocate(b)
	deallocate(n)
	!----------------------------------------------------------------------------
	return
	end subroutine
	!----------------------------------------------------------------------------
	subroutine get_mol(natom,nd,word, mol)
	implicit none
	integer                    :: natom, nd
	integer                    :: mol(natom)
	character(*)               :: word
	character(len(word))       :: word1
	integer                    :: k, ind, lw, j, l, i, ib, ia
	!----------------------------------------------------------------------------
	character(80)              :: ca, cb
	character(80) ,allocatable :: b(:)
	!----------------------------------------------------------------------------
	word1=trim(adjustl(word))
	!----------------------------------------------------------------------------
	allocate(b(nd))
	!----------------------------------------------------------------------------
	do k=1, nd-1
		ind  = index(word1,",")
		lw   = len(trim(word1))
		b(k) = trim(word1(1:ind-1))
		j    = lw - ind
		word1(1:j)  = word1(ind+1:lw)
		word1(j+1:) = ""
	end do
	!----------------------------------------------------------------------------
	b(nd) = trim(adjustl(word1))
	!----------------------------------------------------------------------------
	l=0
	do k=1, nd
		ind  = index(b(k),"-")
		if (ind==0) then
			l = l + 1
			read(b(k),*) mol(l)
		else
			ca=trim(b(k)(1:ind-1))
			cb=trim(b(k)(ind+1:len(b(k))))
			read(ca,*) ia
			read(cb,*) ib
			do i=1,ib-ia+1
				l=l+1
				mol(l) = ia-1 + i
			end do
		end if
	end do
	!----------------------------------------------------------------------------
	deallocate(b)
	!----------------------------------------------------------------------------
	return
	end subroutine
