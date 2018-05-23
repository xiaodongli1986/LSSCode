
program main

use LSS_cosmo_funs

implicit none

	integer :: i, skiprow, maxcol
	real(dl) :: om,w, ra,dec,x,y,z,r,redshift, tmp(10000)
	integer(8) :: numarg, nlines
	character(len=char_len) :: inputfile, outputfile, printstr, tmpstr1, tmpstr2
	logical :: addw1
	
	printstr = 'Please set maxcol >= 4 to have more columns. Usage: EXE -inputfile intpufile -om -om -w -w -skiprow skiprow'//&
		'-outputfile outputfile -maxcol maxcol -addw1 T/F '//&
		'### fmt of input: x, y, z, ***'//&
		'### fmt of output: ra, dec, z, *** ### Converting the ra, dec, z to x, y, z; based on some cosmology'

	! Default values
	om = 0.26; w = -1.0; skiprow = 0; maxcol = 3
	inputfile = 'NONE'
	outputfile = 'NONE'
	
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
		elseif(trim(adjustl(tmpstr1)).eq.'-om') then
			read(tmpstr2,*) om
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		elseif(trim(adjustl(tmpstr1)).eq.'-maxcol') then
			read(tmpstr2,*) maxcol
		elseif(trim(adjustl(tmpstr1)).eq.'-addw1') then
			read(tmpstr2,*) addw1
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	if(maxcol > 10000) then
		print *, 'maxcol larger than 10000; increase size of tmp!: ', maxcol
		stop
	endif

	! intput cosmology
	call cosmo_funs_init(.true.)
	gb_omegam = om; gb_w = w; gb_h = 0.73_dl;
	call de_calc_comovr()
	
	if(trim(adjustl(outputfile)).eq.'NONE') then
		outputfile = trim(adjustl(inputfile))//'.'//trim(adjustl(get_omwstr(om,w)))
	endif

	open(unit=1,file=inputfile,action='read')
	open(unit=2,file=outputfile)
	do i = 1, skiprow
		read(1,'(A)') tmpstr1
		write(2,'(A)') trim(adjustl(tmpstr1))
	enddo
	nlines = 0
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines + 1
		x=tmp(1); y=tmp(2); z=tmp(3)
		call xyz_to_radecr(x,y,z,ra,dec,r)
                redshift = de_zfromintpl(r)
		if(.not.addw1) then
			write(2,'(<maxcol>(e15.7,1x))') ra, dec, redshift, tmp(4:maxcol)
		else
			write(2,'(<maxcol>(e15.7,1x)," 1")') ra, dec, redshift, tmp(4:maxcol)
		endif
		cycle
101		exit
	enddo
	close(1)
	close(2)
	write(*,'(A,i10,A,i8,A,2f10.4,A,i5,A,A,A)') ' Finishing processing ', nlines, ' lines (skipping first ', &
		skiprow,' rows). Converted to cosmology ', &
		om, w, ', ', maxcol, ' columns. Result written to ' , trim(adjustl(outputfile)), '.'

end program main
