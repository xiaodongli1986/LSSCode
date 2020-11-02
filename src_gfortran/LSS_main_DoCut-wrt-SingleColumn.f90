
program main

use LSS_cosmo_funs

implicit none

	integer :: i, skiprow, icol, maxcol
	real(dl) :: tmp(10000), bound, x
	integer(8) :: numarg, nlines, nlines_remain
	character(len=char_len) :: inputfile, outputfile, printstr, tmpstr1, tmpstr2, suffix
	logical :: cutabove
	
	
	printstr = 'Usage: EXE -inputfile intpufile -icol icol -cutabove cutabove -bound bound -skiprow skiprow '//&
		'-outputfile outputfile'//&
		'### Cut the file; key value is a single column; cut above or cut below; given bound'

	! Default values
	cutabove = .true.
	icol = 1
	skiprow = 0; maxcol = 3
	inputfile = 'NONE'
	outputfile = 'NONE'
	suffix = ''
	
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
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,'(A)') suffix
		elseif(trim(adjustl(tmpstr1)).eq.'-cutabove') then
			read(tmpstr2,*) cutabove
		elseif(trim(adjustl(tmpstr1)).eq.'-icol') then
			read(tmpstr2,*) icol
		elseif(trim(adjustl(tmpstr1)).eq.'-bound') then
			read(tmpstr2,*) bound
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	maxcol = icol
	if(maxcol > 10000) then
		print *, 'maxcol larger than 10000; increase size of tmp!: ', maxcol
		stop
	endif

	if(trim(adjustl(outputfile)).eq.'NONE') then
		write(tmpstr1,*) icol
		outputfile = trim(adjustl(inputfile))//'.cutted-icol'//trim(adjustl(tmpstr1))//trim(adjustl(suffix))
	endif

	open(unit=1,file=inputfile,action='read')
	open(unit=2,file=outputfile)
	do i = 1, skiprow
		read(1,'(A)') tmpstr1
		write(2,'(A)') trim(adjustl(tmpstr1))
	enddo
	nlines = 0
	nlines_remain = 0
	do while(.true.)
		read(1,'(A)',end=101) tmpstr1
		read(tmpstr1,*) tmp(1:maxcol)
		nlines = nlines + 1
		x = tmp(icol) 
		if( (cutabove .and. x > bound) .or. ((.not.cutabove).and.x<bound)) then
			write(2,'(A)') trim(tmpstr1)
			nlines_remain = nlines_remain +1
		endif
		cycle
101		exit
	enddo
	close(1)
	close(2)
	write(*,'(A,i10,A,i8,A,i4,e14.7,A,L3,A,i10,A,A)') ' Finishing processing ', nlines, ' lines (skipping first ', &
		skiprow,' rows). Cut column, bound =  ', &
		icol, bound, ', cutabove = ', cutabove, '. ',nlines_remain, ' lines remain and written to ' , trim(adjustl(outputfile)), '.'

end program main
