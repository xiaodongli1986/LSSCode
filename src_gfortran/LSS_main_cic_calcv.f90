program main_cic_calcv

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, printstr1
        integer :: i,j,k, block1,block2,nc
        real*4, allocatable:: rho(:,:,:)
	
        printstr1 = "calcluate velocity from density, based on linear perturbation theory"
	printstr = "LSS_cic_calcv   -inputfile ...   -outputfile ...   -nc 0"
	if(iargc().le.1) then
		print *, printstr
		stop
	endif

        nc = -1

	outputfile = ""
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
		elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
			read(tmpstr2,"(A)") outputfile
		elseif(trim(adjustl(tmpstr1)).eq."-nc") then
			read(tmpstr2,"(A)") nc
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

	if(trim(adjustl(outputfile)).eq."") then
		
	endif

        !open(unit=1000,file=trim(adjustl(inputfile)),action='read',form='binary',access='sequential') ! gfortran fmt
        open(unit=1000,file=trim(adjustl(inputfile)),action='read',form='unformatted',access='stream')
        read(1000) block1
        print *, 'block1 = ', block1
        if(nc == -1) then
                nc = int((block1/4.)**(1./3.) + 0.5)
                print *, ' reset nc as ', nc
        endif
        allocate(rho(nc,nc,nc))
        read(1000) rho
        read(1000) block2
        print *, 'block2 = ', block2

   



end program

subroutine get_velocity(rho, vel)
	use, intrinsic :: iso_c_binding
  	!include 'fftw3.f03'
	implicit none
        real :: ahf
	real :: z, om, h0

	! Sizes of 3D transform
  	integer(C_INTPTR_T), parameter :: N1 = 256
  	integer(C_INTPTR_T), parameter :: N2 = 256
  	integer(C_INTPTR_T), parameter :: N3 = 256

	real(C_DOUBLE), INTENT(IN), pointer :: rho(:,:,:)
	real(C_DOUBLE), INTENT(OUT), pointer :: vel(:,:,:)
	!complex(C_DOUBLE_COMPLEX), pointer :: rho_kspace(:,:,:), vel_kspace(:,:,:)
	real(C_DOUBLE), pointer :: rho_kspace(:,:,:), vel_kspace(:,:,:)
	type(C_PTR) :: plan_rho, plan_vel, data_rho, data_vel

	
	! allocating array
  	!data_rho = fftw_alloc_complex(int((N3/2+1) * N2 * N1, C_SIZE_T))
  	!data_vel = fftw_alloc_complex(int((N3/2+1) * N2 * N1, C_SIZE_T))
	!call c_f_pointer(data_rho, rho, [2*(N3/2+1), N2, N1])
	!call c_f_pointer(data_rho, rho_kspace, [N3/2+1, N2, N1])
	!call c_f_pointer(data_vel, vel_kspace, [N3/2+1, N2, N1])
	!call c_f_pointer(data_vel, rho, [2*(N3/2+1), N2, N1])


	!plan_rho = fftw_plan_r2r_3d(N1,N2,N3, rho, rho_kspace, FFFTW_R2HC, FFTW_R2HC, FFTW_R2HC, FTW_ESTIMATE)
	!call fftw_execute_r2r(plan_rho, rho, rho_kspace)


	!!if()
	!!k2 = twopi/N1 + twopi/N2 + twopi/N3
	!!vel_kspace = rho_kspace*ahf*1.0/sqrt(k2)


	!plan_vel = fftw_plan_r2r_3d(N1,N2,N3, vel_kspace, vel, FFTW_HC2R, FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE)
	!call fftw_execute_r2r(plan_vel, vel_kspace, vel)


	!call fftw_destroy_plan(plan_rho)
	!call fftw_free(data_rho)
	!call fftw_destroy_plan(plan_vel)
	!call fftw_free(data_vel)

end subroutine get_velocity


real function ahf(z,om,h0)
	real :: z, om, f, h
	real :: h0!, ahf
	h =h0*sqrt(om*(1.0+z)**3 + (1.0-om))
	f = om**0.55
	ahf = h*f*1.0/(1.0 + z)

end function ahf
