
program main

use LSS_tools

implicit none

	real(dl) :: rmin, rmax, massmin, massmax, x,y,z,mass,r, tmp(1000), rmin2, rmax2, massmin2, massmax2, mass1,mass2,masscut, &
		xmin,ymin,zmin,xmax,ymax,zmax, vol, frac_surface
	real(dl), allocatable :: redges(:), massedges(:), massmins(:), massmaxs(:)
	integer :: numarg, numrbin,nummassbin, xcol,ycol,zcol,masscol,maxcol, i,j,nowrbin,nowmassbin, skiprow
	integer(8) :: nlines, nlines_degrade, numoutrange, numbigmass, numdegrade, num,num1,num2
	integer(8), allocatable :: binnednum(:,:)
	logical :: logmass, dodegrade
	character(len=2000) :: inputfilename, outputfilename, outputfilemassedges, outputfilename1, outputfilename2, printstr, tmpstr1,tmpstr2, degraded_filename
	
	print *
	printstr = "Usage: EXE -input intpufilename -rmin rmin "//&
		'-rmax rmax -massmin massmin -massmax massmax '//&
		'-logmass logmass -numrbin numrbin -nummassbin nummassbin '//&
		'-skiprow skiprow'//&
		'-xcol xcol -ycol ycol -zcol zcol -masscol masscol '//&
		'-dodegrade dodegrade -degradedfile degraded_filename -numdegrade numdegrade -frac_surface frac_surface'//&
		'### dodegrade will choose a suitable masscut and '//&
		'degrade the data to a subsample with number numdegrade'

	! Default values
	rmin = 0.0; rmax = 100000.0d0; 
	massmin = 1.0d10; massmax = 1.0d16; logmass = .true.
	xcol=1; ycol=2; zcol=3; masscol=4
	numrbin = 1; nummassbin = 100
	frac_surface=1.0
        skiprow=0
        degraded_filename = ''
	
	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif
	
	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-rmin') then
			read(tmpstr2,*) rmin
		elseif(trim(adjustl(tmpstr1)).eq.'-rmax') then
			read(tmpstr2,*) rmax
		elseif(trim(adjustl(tmpstr1)).eq.'-numrbin') then
			read(tmpstr2,*) numrbin
		elseif(trim(adjustl(tmpstr1)).eq.'-massmin') then
			read(tmpstr2,*) massmin
		elseif(trim(adjustl(tmpstr1)).eq.'-massmax') then
			read(tmpstr2,*) massmax
		elseif(trim(adjustl(tmpstr1)).eq.'-logmass') then
			read(tmpstr2,*) logmass
		elseif(trim(adjustl(tmpstr1)).eq.'-nummassbin') then
			read(tmpstr2,*) nummassbin
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		elseif(trim(adjustl(tmpstr1)).eq.'-xcol') then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq.'-ycol') then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq.'-zcol') then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq.'-masscol') then
			read(tmpstr2,*) masscol
		elseif(trim(adjustl(tmpstr1)).eq.'-dodegrade') then
			read(tmpstr2,*) dodegrade
		elseif(trim(adjustl(tmpstr1)).eq.'-numdegrade') then
			read(tmpstr2,*) numdegrade
		elseif(trim(adjustl(tmpstr1)).eq.'-degradedfile') then
			read(tmpstr2,'(A)') degraded_filename
		elseif(trim(adjustl(tmpstr1)).eq.'-frac_surface') then
			read(tmpstr2,*) frac_surface
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	maxcol = max(xcol,ycol,zcol,masscol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif
	outputfilename = trim(adjustl(inputfilename))//'.rmassinfo'
	outputfilemassedges = trim(adjustl(inputfilename))//'.massedges'
	
	print *, '#####################################'
	print *, ' Settings'
	print *, '  r: range, numrbin = ', real(rmin), real(rmax), numrbin
	print *, '  mass: range, numbin, logmass = ', real(massmin), real(massmax), nummassbin, logmass
	print *, '  cols of x,y,z,mass: ', xcol,ycol,zcol,masscol
	write(*,'(A,A)') '   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,A)') '   outputfilename (r-mass stat): ', trim(adjustl(outputfilename))
	write(*,'(A,A)') '   outputfilename (mass edges): ', trim(adjustl(outputfilemassedges))
	if(dodegrade) then
		print *, '  dodegrade with #-goal: ', numdegrade
	endif
	print *, '  frac_surface (used to calculate nbar) = ', frac_surface
	print *, '#####################################'


	allocate(redges(numrbin+1),massedges(nummassbin+1),binnednum(numrbin,nummassbin),massmins(numrbin),massmaxs(numrbin))
	
	do i = 1, numrbin+1
		redges(i) = rmin + (rmax-rmin)/(dble(numrbin))*(i-1)
	enddo

	do i = 1, numrbin
		massmins(i) =  1.0e30
		massmaxs(i) = -1.0e30
	enddo
	
	do i = 1, nummassbin+1
		if(logmass) then
			massedges(i) = exp(log(massmin) + (log(massmax)-log(massmin))/(dble(nummassbin))*(i-1))
		else
			massedges(i) = massmin + (massmax-massmin)/(dble(nummassbin))*(i-1)
		endif
	enddo
	
	write(*,*)     '  Edges of r:    ', redges(1:numrbin+1)
	write(*,*) '  Edges of mass: ', massedges(1:nummassbin+1)
	
	binnednum = 0
	
	nlines =0 
	numoutrange =0
	numbigmass =0
	rmin2=1.0e30;rmax2=-rmin2;
	xmin = 1.0e30; ymin = 1.0e30; zmin = 1.0e30
	xmax = -xmin; ymax = -ymin; zmax = -zmin
	massmin2=rmin2;massmax2=-massmin2
	open(unit=1,file=inputfilename,action='read')
	do while(.true.)
                if(nlines+1<=skiprow) then
                        read(1,*,end=101) tmpstr1
                        nlines = nlines+1
                        cycle
                endif
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines+1
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); mass=tmp(masscol)
		r = sqrt(x*x+y*y+z*z)
		rmin2=min(r,rmin2); rmax2=max(r,rmax2)
		xmin = min(x,xmin); xmax = max(x,xmax)
		ymin = min(y,ymin); ymax = max(y,ymax)
		zmin = min(z,zmin); zmax = max(z,zmax)
		massmin2=min(mass,massmin2); massmax2=max(mass,massmax2)
		
		if(numrbin<100) then
			nowrbin = find_ibin(r,redges,numrbin+1)
		else
			nowrbin = find_ibin_2split(r,redges,numrbin+1)
		endif
		if(nummassbin<100) then
			nowmassbin = find_ibin(mass,massedges,nummassbin+1)
		else
			nowmassbin = find_ibin_2split(mass,massedges,nummassbin+1)
		endif

		if(nowrbin.ge.1.and.nowrbin.le.numrbin.and.nowmassbin.ge.1.and.nowmassbin.le.nummassbin) then
			binnednum(nowrbin,nowmassbin) = binnednum(nowrbin,nowmassbin)+1
			massmins(nowrbin) = min(massmins(nowrbin),mass)
			massmaxs(nowrbin) = max(massmaxs(nowrbin),mass)
		else
			numoutrange = numoutrange+1
		endif
		
		if(mass>massedges(nummassbin+1)) numbigmass=numbigmass+1
			
		cycle
101		exit
	enddo		
	close(1)
	
	write(*,'(A,i10,A,i10)') ' Finishing processing ', nlines, ' lines; numoutrange = ', numoutrange
	write(*,'(A,3("(",2f15.7,"); "), A,2f15.7,A,2e15.7,A,i10,A,i10)') &
		'Done. range of x,y,z: ', xmin,xmax,ymin,ymax,zmin,zmax, '; min/max of r:', real(rmin2),real(rmax2), &
		';  min/max of mass:', real(massmin2),real(massmax2), &
		';  numoutrange = ', numoutrange, '; numbigmass = ', numbigmass
	open(unit=2,file=outputfilename)
	open(unit=3,file=outputfilemassedges)
	print *, 'Result: rmin, rmax, massmin, massmax, num'
	write(2,'(A,3("(",2f15.7,"); "), A,2f15.7,A,2e15.7,A,i12,A,i12,A,i12)') &
		'#Result: rmin, rmax, massmin, massmax, num, nbar. p.s. range of x,y,z: ',&
		xmin,xmax,ymin,ymax,zmin,zmax, 'min/max of r:', real(rmin2),real(rmax2), &
		';  min/max of mass:', real(massmin2),real(massmax2), &
		';  numoutrange = ', numoutrange, '; numbigmass = ', numbigmass, &
		';  total num = ', nlines
	write(3,'(A,3("(",2f15.7,"); "), A,2f15.7,A,2e15.7,A,i12,A,i12,A,i12)') &
		'#Result: rmin, rmax, minimal & maximal mass in this bin . p.s. range of x,y,z: ',&
		xmin,xmax,ymin,ymax,zmin,zmax, 'min/max of r:', real(rmin2),real(rmax2), &
		';  min/max of mass:', real(massmin2),real(massmax2), &
		';  numoutrange = ', numoutrange, '; numbigmass = ', numbigmass, &
		';  total num = ', nlines
	do i = 1, numrbin
		write(3, '(2f15.7,4x,2e15.7)') redges(i), redges(i+1), massmins(i), massmaxs(i)
		do j = 1, nummassbin
			vol = 4.0/3.0*const_pi * (redges(i+1)**3.0-redges(i)**3.0) * frac_surface
			write(*,'(2f15.7,2e15.7, i10, e15.7)') redges(i), redges(i+1), massedges(j), massedges(j+1), binnednum(i,j), binnednum(i,j)/vol
			write(2,'(2f15.7,2e15.7, i10, e15.7)') redges(i), redges(i+1), massedges(j), massedges(j+1), binnednum(i,j), binnednum(i,j)/vol
		enddo
	enddo
	close(1)
	close(3)
	
	if(.not.dodegrade) stop

        if(trim(adjustl(degraded_filename)) .eq. '') then
        	outputfilename1 = trim(adjustl(inputfilename))//'.degraded'
        	outputfilename2 = trim(adjustl(inputfilename))//'.degraded.info'	
        else
        	outputfilename1 = trim(adjustl(degraded_filename))
        	outputfilename2 = trim(adjustl(outputfilename1))//'.info'	
        endif

	if(numdegrade.gt.nlines) then
		print *, 'numdegrade larger than nline: no need to do degrade! Programme will just copy file!'
		open(unit=3,file=outputfilename2)
		write(3,'(A,A,A,i10,A,i10,A,A)') '# This is for file ', trim(adjustl(outputfilename)), '. This file has ',nlines,&
			' lines. We want to degrade it to ',numdegrade,' halos. New sample written to: ', trim(adjustl(outputfilename1))
		write(3,'(A)') '# numdegrade larger than nline: no need to do degrade! Programme will just copy file!'
		close(3)
		call system('cp '//trim(adjustl(inputfilename))//' '//trim(adjustl(outputfilename1)))
		
		stop
	endif
	
	do j = nummassbin, 1, -1
		num = sum(binnednum(1:numrbin,j:nummassbin)) + numbigmass
		if(num>numdegrade) exit
	enddo
	
	mass2 = massedges(j); num2 = sum(binnednum(1:numrbin,j:nummassbin)) + numbigmass
	
	do i = j+1, nummassbin
		num = sum(binnednum(1:numrbin,i:nummassbin)) + numbigmass
		if(num.ne.num2) exit
	enddo
	mass1 = massedges(i); num1 = sum(binnednum(1:numrbin,i:nummassbin)) + numbigmass
	masscut = mass1 + dble(numdegrade-num1)/dble(num2-num1)*(mass2-mass1)
	
	write(*,'(1x,A,e15.7)') 'Do degrade: Applying masscut ', masscut
	
	! TBD...
	nlines_degrade =0
        nlines = 0
	open(unit=1,file=inputfilename,action='read')
	open(unit=2,file=outputfilename1)
	open(unit=3,file=outputfilename2)
	do while(.true.)
                if(nlines+1<=skiprow) then
                        read(1,*,end=102) tmpstr1
                        nlines = nlines+1
                        cycle
                endif
		read(1,'(A)',end=102) tmpstr1
		read(tmpstr1,*) tmp(1:maxcol)
                nlines = nlines+1
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); mass=tmp(masscol)
		if(mass > masscut) then
			write(2,'(A)') trim(adjustl(tmpstr1))
			nlines_degrade = nlines_degrade+1
		endif
		cycle
102		exit
	enddo
	close(2);

	write(*,'(A,A,A,i10,A,i10,A,A)') '# This is for file ', trim(adjustl(outputfilename)), '. This file has ',nlines,&
		' lines. We want to degrade it to ',numdegrade,' halos. New sample written to: ', trim(adjustl(outputfilename1))
	write(3,'(A,A,A,i10,A,i10,A,A)') '# This is for file ', trim(adjustl(outputfilename)), '. This file has ',nlines,&
		' lines. We want to degrade it to ',numdegrade,' halos. New sample written to: ', trim(adjustl(outputfilename1))

	write(*,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass1, ', there are ', num1, ' halos'
	write(3,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass1, ', there are ', num1, ' halos'

	write(*,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass2, ', there are ', num2, ' halos'
	write(3,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass2, ', there are ', num2, ' halos'

	write(*,'(A,e15.7,A,i10,A,f15.7)') '# We interploate and adopt  mass>masscut=', masscut, ', we obtain ', nlines_degrade, &
		 ' halos. Rato (num-obtaned/num-goal): ', dble(nlines_degrade)/dble(numdegrade)
	write(3,'(A,e15.7,A,i10,A,f15.7)') '# We interploate and adopt  mass>masscut=', masscut, ', we obtain ', nlines_degrade, &
		 ' halos. Rato (num-obtaned/num-goal): ', dble(nlines_degrade)/dble(numdegrade)
	close(3)
end program main
