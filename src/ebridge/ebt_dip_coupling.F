	subroutine ebt_dip_coupling(coc, dip, H, V, kappa, r_norm, theta)
	implicit none
	!----------------------------------------------------------------------------
	real(8) :: coc(3,2)
	real(8) :: dip(3,2)
	!----------------------------------------------------------------------------
	real(8) :: r_norm
	real(8) :: theta(3)
	real(8) :: dip_a_norm
	real(8) :: dip_b_norm
	real(8) :: r_unit(3)
	real(8) :: dip_a_unit(3)
	real(8) :: dip_b_unit(3)
	!----------------------------------------------------------------------------
	real(8) :: a, b, c
	real(8) :: cosa
	real(8) :: cosb
	real(8) :: cosc
	real(8) :: V
	real(8) :: H
	real(8) :: kappa
	real(8) :: PI
	real(8) :: eps
	!----------------------------------------------------------------------------
	PI          = acos(-1.D0)
	eps         = 1.D-8
	!----------------------------------------------------------------------------
	r_unit      = coc(:,2) - coc(:,1)
	!----------------------------------------------------------------------------
	dip_a_unit  = dip(:,1)
	dip_b_unit  = dip(:,2)
	!----------------------------------------------------------------------------
	r_norm      = sqrt(dot_product(    r_unit,     r_unit))
	dip_a_norm  = sqrt(dot_product(dip_a_unit, dip_a_unit))
	dip_b_norm  = sqrt(dot_product(dip_b_unit, dip_b_unit))
	!----------------------------------------------------------------------------
	r_unit      = r_unit     / r_norm
	dip_a_unit  = dip_a_unit / dip_a_norm
	dip_b_unit  = dip_b_unit / dip_b_norm
	!----------------------------------------------------------------------------
	cosa        = dot_product(dip_a_unit,     r_unit)
	cosb        = dot_product(dip_b_unit,     r_unit)
	cosc        = dot_product(dip_a_unit, dip_b_unit)
	kappa       = cosc - 3.D0 * cosa * cosb
	!----------------------------------------------------------------------------
	if (abs(cosa-1.D0) < eps) then
		theta(1) = 0.D0
	else
		theta(1) = acos(cosa) * 180.D0 / PI
	end if
	!----------------------------------------------------------------------------
	if (abs(cosb-1.D0) < eps) then
		theta(2) = 0.D0
	else
		theta(2) = acos(cosb) * 180.D0 / PI
	end if
	!----------------------------------------------------------------------------
	if (abs(cosc-1.D0) < eps) then
		theta(3) = 0.D0
	else
		theta(3) = acos(cosc) * 180.D0 / PI
	end if
	!----------------------------------------------------------------------------
	V           = dip_a_norm * dip_a_norm / r_norm**3
	H           = V * kappa
	!----------------------------------------------------------------------------
	return
	end subroutine
