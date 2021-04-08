program main_add_gass_error_periodic_box

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, sigma, mu
	integer :: i
	
	printstr = "Add "
	if(iargc().le.1) then
		print *, printstr
		stop
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
		
	endif

	print *, 'This is an empty program add_gass_error_periodic_box!!'

end program
