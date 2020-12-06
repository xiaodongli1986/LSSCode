
! What treatment we have adopted now:


program LSS_main_TakePartofFile

use LSS_tools

	implicit none
	character(len=char_len) :: printstr, inputfile, outputfile, tmpstr1, tmpstr2
	real(dl) :: ratio, x
	integer :: numarg, i, j


	printstr = '### EXE -inputfile inputfile -outputfile outputfile -ratio ratio'

	ratio = 0.01
	outputfile = 'NONE'
	inputfile = 'NONE'

	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-inputfile') then
			read(tmpstr2,'(A)') inputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-outputfile') then
			read(tmpstr2,'(A)') outputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-ratio') then
			read(tmpstr2, *)  ratio
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	if(ratio .ge. 1.0) then
		print *, 'ratio larger than 1, no need to run the programe!!!', ratio
		stop
	endif

	if(trim(adjustl(outputfile)) .eq. 'NONE') then
		write(tmpstr2, '(f8.2)') 100*ratio
		outputfile = trim(adjustl(inputfile))//'.'//trim(adjustl(tmpstr2))//'percent'
	endif

	i = 0; j=0
	open(unit=1,file=inputfile)
	open(unit=2,file=outputfile)
	do while(.true.)
		read(1,'(A)',end=100) tmpstr1
		i = i+1
		call random_number(x)
		if(x < ratio) then
			write(2,'(A)') trim(adjustl(tmpstr1))
			j = j+1
		endif
		cycle
100		exit
	enddo
	close(1)
	close(2)

	write(*,'(A,i9,A,A,A,i9,A,A)') 'Finishing processing ', i, ' lines from ', trim(adjustl(inputfile)), &
		'; ', j, ' lines written to ', trim(adjustl(outputfile))
end program LSS_main_TakePartofFile
