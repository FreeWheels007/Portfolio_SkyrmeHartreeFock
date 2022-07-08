subroutine coulomb_potential(nPoints, zVec, rho, d2rho, chargeFactor, boxLength, dz, vPotz)
        implicit none
        integer, parameter :: dp = kind(0.d0), xyPoints = 100
        real(dp), parameter :: numConst = (3.d0/3.14159d0)**(1.d0/3.d0)
        integer, intent(in) :: nPoints
	real(dp), dimension(nPoints), intent(in) :: zVec, rho, d2rho
        real(dp), intent(in) :: chargeFactor, boxLength, dz
        real(dp), dimension(nPoints), intent(out) :: vPotz
	real(dp), dimension(nPoints) :: zVec2, zWeights
        real(dp), dimension(:), allocatable :: integrandzp
        real(dp), dimension(:, :), allocatable :: integrand
	real(dp), dimension(xyPoints) :: xyVec, xyWeights
        real(dp) :: term1Factor, term2Factor
        integer :: k, kp
        
        !Set x y axis using gaulegint
        call gauleg(0.d0, boxLength, xyVec, xyWeights, xyPoints)
        !set z axis for gaulegint, zVec2 is just necassary dummy array for gaulegint
        call gauleg(0.d0, boxLength, zVec2, zWeights, nPoints)

        allocate(integrandzp(nPoints), integrand(nPoints, nPoints))

	!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,kp)
        do k = 1, nPoints
                do kp = 1, xyPoints
                        integrand(kp, k) = d2rho(kp)*integralOverXYXpYp(xyPoints, zVec(k), zVec(kp), xyVec, xyWeights)
                enddo

                ! integrate over z', leaving vector over z
                integrandzp(k) = gaulegint(nPoints, integrand(:, k), zVec2, zWeights)
        enddo
	!$OMP END PARALLEL DO
        
        ! combine potential with exhange tern
        term1Factor = chargeFactor/(4.d0*boxLength**2)
        term2Factor = 3.d0*chargeFactor*numConst/4.d0
        vPotz(:) = term1Factor*integrandzp(:) - term2Factor*(rho(:))**(1.d0/3.d0)

        deallocate(integrandzp, integrand)

        ! end of subroutine body, now declaring helper functions not to be used by python script
contains

real(dp) function integralOverXYXpYp(xyPoints, z, zp, xyVec, xyWeights)
        implicit none
        integer, intent(in) :: xyPoints
	real(dp), dimension(xyPoints), intent(in) :: xyVec, xyWeights
	real(dp), intent(in) :: z, zp
        real(dp), dimension(:), allocatable :: integrandxpypx
        real(dp), dimension(:, :), allocatable :: integrandxpyp
        real(dp), dimension(:, :, :), allocatable :: integrandxp
        real(dp), dimension(:, :, :, :), allocatable :: integrand
        real(dp) :: integralTemp
        integer :: i, j, ip, jp

        allocate(integrandxpypx(xyPoints), integrandxpyp(xyPoints, xyPoints), &
                integrandxp(xyPoints, xyPoints, xyPoints), integrand(xyPoints, xyPoints, xyPoints, xyPoints))

        do j = 1, xyPoints
                do i = 1, xyPoints
                        do jp = 1, xyPoints
                                do ip = 1, xyPoints
                                        integrand(ip, jp, i, j) = sqrt((z - zp)**2 + (xyVec(j) - xyVec(jp))**2 + &
                                                                  (xyVec(i) -xyVec(ip))**2)
                                end do

                                integrandxp(jp, i, j) = gaulegintxy(xyPoints, integrand(:, jp, i, j), xyWeights)
                        end do

                        integrandxpyp(i, j) = gaulegintxy(xyPoints, integrandxp(:, i, j), xyWeights)
                end do

                integrandxpypx(j) = gaulegintxy(xyPoints, integrandxpyp(:, j), xyWeights)
        end do

        ! final integral resulting in value only dependant on z and z'
        integralTemp = gaulegintxy(xyPoints, integrandxpypx(:), xyWeights)

        deallocate(integrandxpypx, integrandxpyp, integrandxp, integrand)

        integralOverXYXpYp = integralTemp
end function

real(dp) function gaulegint(accuracy, A, x_array, w_array)
        implicit none
        integer, intent(in) :: accuracy
	real(dp), dimension(accuracy), intent(in) :: A, x_array, w_array
        integer :: interlow, interhigh, a_index
	real(dp) :: interweight, stepsize
        
        stepsize = dz
        
        gaulegint = 0.0d0
        do a_index = 1, accuracy
                !linear interpolation between points is used for better accuracy
                interlow = floor(x_array(a_index)/stepsize)+1
                interhigh = interlow + 1
                interweight = (x_array(a_index)/stepsize) - floor(x_array(a_index)/stepsize)
                
                if (interlow == accuracy) then
                        gaulegint = gaulegint +  w_array(a_index)*(A(accuracy)*(1.0d0 - interweight) + A(1)*interweight)
                else
                        gaulegint = gaulegint +  w_array(a_index)*(A(interlow)*(1.0d0 - interweight) + A(interhigh)*interweight)
                end if
        end do
end function

real(dp) function gaulegintxy(accuracy, A, weights)
        implicit none
        integer, intent(in) :: accuracy
        real(dp), dimension(accuracy), intent(in) :: A, weights
        integer :: xy_index

        gaulegintxy = 0.d0

        do xy_index = 1, accuracy
                gaulegintxy = gaulegintxy + weights(xy_index)*A(xy_index)
        end do
end function

SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER, intent(in) :: n
      real(dp), intent(in) :: x1,x2
      real(dp), intent(out) :: x(n),w(n)
      real(dp), PARAMETER :: EPS=1.d-14
      INTEGER :: i,j,m
      real(dp) :: p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
END SUBROUTINE gauleg

end subroutine coulomb_potential
