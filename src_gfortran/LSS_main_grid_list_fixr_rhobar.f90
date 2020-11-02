program main_grid_list_fixr_rhobar

use LSS_grid_tools

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, printstr, outputfmtstr, inputfile, xyzfile, outputfile, fixrstr
	type(gridtype) :: grid
	real(dl), allocatable :: xyzw(:,:)
	integer :: numarg, numNB, n0NB, i, npos, nbd
	logical :: touchboundary
	real(dl) :: x,y,z,fixr,rhomean,rho,drhox,drhoy,drhoz,nbarmean,nbar,dnbarx,dnbary,dnbarz, &
		numNBmean
	

	printstr = 'EXE -inputfile inputfile -xyzfile xyzfile -fixr fixr -outputfile outputfile '//&
		'# intputfile is the data in fmt of x,y,z,mass; xyzfile is list of positions '//&
		'where you want to compute rho, in fmt of x,y,z; '//&
		'fixr is the fixed distance within which you search for the neighbors and calculate rho'

	if(iargc().le.1) then
		print *, printstr
		stop
	endif
	
	outputfile = ''
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-inputfile') then
			read(tmpstr2,'(A)') inputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-outputfile') then
			read(tmpstr2,'(A)') outputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzfile') then
			read(tmpstr2,'(A)') xyzfile
		elseif(trim(adjustl(tmpstr1)).eq.'-fixr') then
			read(tmpstr2,*) fixr
			fixrstr = tmpstr2
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo
	

	if(trim(adjustl(outputfile)).eq.'') then
		outputfile = trim(adjustl(xyzfile))//'.rhodrholist.fixr'//trim(adjustl(fixrstr))
	endif

	write(*,'(A)') ' Computing mass/number density: '
	write(*,'(A,A)') '   Input datafile (fmt: x,y,z,weight): ', trim(adjustl(inputfile))
	write(*,'(A,A)') '   Input position file (fmt: x,y,z)  : ', trim(adjustl(xyzfile))
	write(*,'(A,A)') '   Fixed radius of smoothing sphere  : ', trim(adjustl(fixrstr))

	
	! initialize the grid from file with fmt x,y,z,weight
	call grid_init_fromfile(grid, inputfile, .true.)

	! process the file with positions, compute various densities...
	open(unit=1,file=xyzfile,action='read')
	open(unit=2,file=outputfile)

	!gfortran fmt
	outputfmtstr = '# result using fixed distance '//trim(adjustl(fixrstr))//&
		'. fmt:   x,y,z,  numNB (number of neighbour found),  '//&
		'touchboundary (if True, the sphere near this point touches'//&
		' the boundary, results may be not reliable),   rhomean (mean density '//&
		'within the volume),   rho,drhox,drhoy,drhoz (density and gradients '//&
		'estimated using spline kernel),   nbarmean,nbar,dnbarx,dnbary,dnbarz (same as rho but number density)' 
	write(2,'(A)') trim(adjustl(outputfmtstr))

	numNBmean = 0
	npos = 0
	nbd = 0
	n0NB = 0
	do while(.true.)
		read(1,*,end=100) x,y,z
		npos = npos+1
		call grid_fixr_rho_drho(grid, x,y,z, fixr, &
			numNB,touchboundary, &
			rhomean,rho,drhox,drhoy,drhoz,&
			nbarmean,nbar,dnbarx,dnbary,dnbarz)
		numNBmean = numNBmean+numNB
		if(numNB.eq.0) n0NB = n0NB+1
		if(touchboundary) nbd = nbd+1
		write(2,'(3e15.7,i10,L2,2e15.7,3x,3e15.7,5x,2e15.7,3x,3e15.7)') x,y,z,numNB,touchboundary, &
			rhomean,rho,drhox,drhoy,drhoz,&
			nbarmean,nbar,dnbarx,dnbary,dnbarz
		cycle
100		exit
	enddo
	close(1)
	close(2)
	numNBmean = numNBmean/dble(npos)
	write(*,'(A,f10.3,A,i11)') ' Mean # of in-sphere point = ', numNBmean, ';    empty sphere: ', n0NB
	write(*,'(A,i11,A,f6.3,A)'), ' Touching boundary event: ', nbd, '(', nbd/dble(npos)*100, '%)'
	write(*,'(A,i11,A,i11,A)') ' Finishing processing ', npos, ' positions from ', grid%np, ' data points.'
	write(*,'(A,A)') ' Result written to: ', trim(adjustl(outputfile))

end program
