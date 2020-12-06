
! July 23: This is an old code. We can:
!	# test it;
!	# update it (distance rather than redshift split; no consideration of ra/dec split and just redshift split; ...)

!####################################
! scan chisqs
!####################################

program main

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: inputfilename, suffix
	integer :: xcol, ycol, zcol, masscol
	real(dl) :: deltar, nbar, rmin, rmax, frac_surface

	character(len=char_len) :: outputfilename, tmpstr1, tmpstr2, tmpstr3
	character(len=3000) :: printstr
	integer :: i, maxcol, numarg, linenumber, irbin,numrbin
	real(dl) :: tmp(1000), x,y,z,r, mass, rmin_real, rmax_real, massmin, massmax, rminthisbin,rmaxthisbin
	real(dl), allocatable :: masses(:), rs(:)


	print *, '	This code is not finished yet!!!'
	stop

	printstr = 'Usage:  EXE -input inputfilename -xcol xcol -ycol ycol -zcol zcol -masscol masscol -deltar deltar -frac_surface frac_surface -deltar deltar -suffix suffix -nbar nbar -rmin rmin -rmax rmax ### if all-sky sample then frac_surfact=1 if half-sky then 0.5 by default 1 ###'

	! Default values
	xcol=1; ycol=2; zcol=3; masscol=4;
	frac_surface = 1.0
	deltar = 30.0
	nbar = 0.0
	suffix = '.constant-nbar-down-sampling'

	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif

	tmpstr1 = '   haha it is a test   '

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-xcol') then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq.'-ycol') then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq.'-zcol') then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq.'-masscol') then
			read(tmpstr2,*) masscol
		elseif(trim(adjustl(tmpstr1)).eq.'-deltar') then
			read(tmpstr2,*) deltar
		elseif(trim(adjustl(tmpstr1)).eq.'-rmin') then
			read(tmpstr2,*) rmin
		elseif(trim(adjustl(tmpstr1)).eq.'-rmax') then
			read(tmpstr2,*) rmax
		elseif(trim(adjustl(tmpstr1)).eq.'-frac_surface') then
			read(tmpstr2,*) frac_surface
		elseif(trim(adjustl(tmpstr1)).eq.'-nbar') then
			read(tmpstr2,*) nbar
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,'(A)') suffix
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	if(nbar.eq.0.0) then
		print *, 'ERROR!! Must input nbar!: nbar = ', nbar
		stop
	endif

      	maxcol = max(xcol,ycol,zcol,masscol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif

	write(tmpstr2, '(f10.7)') nbar
	outputfilename = trim(adjustl(inputfilename))//trim(adjustl(suffix))//'.nbar'//trim(adjustl(tmpstr2))
	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,A)') 	'   outputfilename:  ', trim(adjustl(outputfilename))
	write(*,'(A,4i3)')  	'   cols of x,y,z,mass: ', xcol,ycol,zcol,masscol
	write(*,'(A,e14.7)')  	'   nbar: ', nbar
	write(*,'(A,e14.7)')  	'   deltar: ', deltar
	write(*,'(A,e14.7)')  	'   frac_surface: ', frac_surface
	write(*,'(A,2e15.7)')  	'   rmin, rmax (down sampling confined in this region): ', rmin, rmax
	write(*,'(A,A)')	'   suffix = ', trim(adjustl(suffix))
	write(*,'(A,i3)')	'   maxcol = ', maxcol
	print *, '#####################################'
	print *

	call count_line_number(inputfilename, linenumber)
	print *, 'linenumber = ', linenumber
	allocate(rs(linenumber),masses(linenumber))

	open(unit=1,file=inputfilename)
	i = 1
	rmin_real = 1.0e20; rmax_real = -rmin_real;
	massmin = rmin_real; massmax = -massmin;
	do while(.true.)
		read(1,*,end=100) tmp(1:maxcol)
		r = sqrt(tmp(xcol)**2.0 + tmp(ycol)**2.0 + tmp(zcol)**2.0)
		mass = tmp(masscol)
		rs(i) = r
		masses(i) = mass
		rmin_real = min(rmin_real, r)
		rmax_real = max(rmax_real, r)
		massmin = min(massmin, mass)
		massmax = max(massmax, mass)
		i = i+1
		cycle
100		exit
	enddo

!	This code is not finished yet!!!
!	num
!	do ibin = 1, 

	! range of r
	print *, 'After scanning file we find rmin_real, rmax_real = ', rmin_real, rmax_real
	if(rmin_real>rmin .or. rmax_real < rmax) then	
		print *, 'WARNING!!! Find inputed rmin, rmax does not covered by the real range of data rmin_real, rmax_real :'
		print *, '  rmin, rmax (input) = ', rmin, rmax
		print *, '  rmin, rmax (real)  = ', rmin_real, rmax_real
	endif

	! range of mass
	print *, 'After scanning file we find massmin, massmax = ', real(massmin), real(massmax)

	! line number check
	if(i-1.ne.linenumber) then
		print *, 'ERROR! #-line mismatch: i-1, linenumber = ', i-1, linenumber
		stop
	endif


end program main
