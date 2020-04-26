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

        open(unit=1000,file=trim(adjustl(inputfile)),action='read',form='binary',access='sequential')
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


subroutine get_velocity_field(n , inputfile, outputfile)        

