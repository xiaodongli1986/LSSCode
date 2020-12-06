program main_converttoRhct

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr
	integer :: i
        real(dl) :: weight, r2, redshift, x,y,z,r, DA2
	
	printstr = "Now it is empty!"
	if(iargc().le.1) then
		print *, printstr
	endif

	outputfile = ""
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
		elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
			read(tmpstr2,"(A)") outputfile
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

	if(trim(adjustl(outputfile)).eq."") then
		outputfile=trim(adjustl(inputfile))//'.Rhct'
	endif

	print *, 'Convert the sample from om=0.26 WMAP5 comoving coordinates to the coordinates in Rh = ct cosmology. Assuming the format is x,y,z,w '
        print *, gb_h
        call cosmo_funs_init(.true.)
        gb_omegam = 0.26; gb_w = -1; gb_h = 0.73_dl;
        call de_calc_comovr()

        open(unit=1000,file=inputfile)
        open(unit=1001,file=outputfile)
        do while(.true.)
                read(1000,*,end=101) x,y,z,weight
                r = sqrt(x*x+y*y+z*z)
                redshift = get_z(r, 0.0_dl, 5.0_dl) !de_zfromintpl(r)
                print *, r, redshift, Hz(redshift)
                print *, '((100.0*gb_h) * (1.0+z)), Hz(redshift) = ', ((100.0*gb_h) * (1.0+redshift)), Hz(redshift)
                DA2 = 1.0 / ((100.0*gb_h) * (1.0+redshift)) * log(1.0+redshift) * const_c*gb_h
                r2 = DA2 * (1.0+redshift)
                print *,'redshift, r, DA2, r2, r2/r = ', redshift, r, DA2, r2, r2/r
                stop
                write(1001,'(4e14.5)') x*r2/r, y*r2/r, z*r2/r, weight
                cycle
101             exit
        enddo
        close(1000)
        close(1001)

end program
