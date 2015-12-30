program main_2pcf__numline_sumweight

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr
	integer :: i, icol_weight, num_line
	real(dl) :: tmpX(1000), sumweight
	
	printstr = "EXE -inputfile inputfile -outputfile outputfile -icol_weight icol_weight"
	if(iargc().le.1) then
		print *, printstr
		stop
	endif

	icol_weight = 4
	outputfile = ""
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
		elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
			read(tmpstr2,"(A)") outputfile
		elseif(trim(adjustl(tmpstr1)).eq."-icol_weight") then
			read(tmpstr2,*) icol_weight
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

	if(trim(adjustl(outputfile)).eq."") then
		outputfile = trim(adjustl(inputfile))//'.numline_sumweight'
	endif

	open(unit=100,file=inputfile)
	num_line = 0
	sumweight = 0.0_dl
	do while(.true.)
		read(100,*,end=100) tmpX(1:icol_weight)
		num_line = num_line + 1
		sumweight = sumweight + tmpX(icol_weight)
		cycle
100 		exit
	enddo
	close(100)
	
	open(unit=101,file=outputfile)
	write(101,'(i10,e40.30)') num_line, sumweight
	close(101)

	write(*,'(A,i10,e40.30)') trim(adjustl(inputfile)), num_line, sumweight

end program
