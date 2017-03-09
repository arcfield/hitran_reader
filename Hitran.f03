program Hitran_Calc
 implicit none
 
 double precision :: nu_ini = 100.0, nu_end = 2000.0, nu_step = 1.0e-3
 double precision, allocatable :: nu(:), inten(:), gam_air(:), n_air(:), del_air(:), gam(:), delnu(:)
 double precision, allocatable :: f(:)
 double precision :: A, E_pprime, gam_self
 real, parameter :: prs = 950.0, pi = 4*atan(1.0)
 integer :: i, j, delnu_siz, lines, mole_ID, iso_ID


 ! Opening HITRAN data file and storing to an array
 open (unit=1, file='h2o_HITRAN.par')
 lines= 0
 do
  read(1,*,IOSTAT=delnu_siz) nu
  if(delnu_siz /= 0) exit
  lines = lines + 1
 end do
 delnu_siz = 0
 close(unit=1)

 delnu_siz = ((nu_end - nu_ini) / nu_step) + 1

 ! Allocating size of the variables
 allocate(nu(lines))
 allocate(inten(lines))
 allocate(gam_air(lines))
 allocate(n_air(lines))
 allocate(del_air(lines))
 allocate(gam(lines))
 allocate(delnu(delnu_siz))
 allocate(f(delnu_siz))


 ! Storing data to the variables
 open(unit=1, file='h2o_HITRAN.par')
 do i = 1,lines
  read(1,'(I2, I1, F12.6, E10.3, E10.3, F5.4, F5.3, F10.4, F4.2, F8.6)') &
  mole_ID, iso_ID, nu(i), inten(i), A, gam_air(i), gam_self, E_pprime, n_air(i), del_air(i)
 end do


! do i = 1,lines
!  print '(F14.6, E12.3, F7.4, F6.2, F10.6)', nu(i), inten(i), gam_air(i), n_air(i), del_air(i)
! end do

 ! 'arange' routine to discretize the wave numbers
 call arange(nu_ini, nu_end, nu_step, delnu)
 

 ! Calculation of 'gamma'
 do i = 1,lines
  gam(i) = ((296.0/300.0)**n_air(i))*(gam_air(i)*(prs/1013.0))
 end do



 ! 'f': Absorption coefficient calculation
 do i = 1,lines
  f = f + inten(i) * ((1.0/pi)*(gam(i)/(gam(i)**2.0 + (delnu - (nu(i) + del_air(i)*prs/1013.0))**2.0)))
 end do


! do i = 1,delnu_siz 
!  print *, f(i)
! end do

 deallocate(nu)
 deallocate(inten)
 deallocate(gam_air)
 deallocate(n_air)
 deallocate(del_air)
 deallocate(gam)
 deallocate(delnu)
 deallocate(f)

 close(unit=1)
end program Hitran_Calc


subroutine arange(ini, last, step, output)
 integer :: i, siz
 double precision, intent(in) :: ini, last, step
 double precision, dimension(:), intent(out) :: output(siz)

 siz = ((last - ini) / step) + 1 

 do i = 1,siz
  output(i) = ini + (step * (i - 1))
 end do

end subroutine arange

