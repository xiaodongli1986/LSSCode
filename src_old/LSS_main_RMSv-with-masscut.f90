
! We calculate the r.m.s velocity of the sample, with a given mass cut imposed

program main
 
use LSS_cosmo_funs

implicit none

	integer :: i, skiprow, xcol,ycol,zcol,vxcol,vycol,vzcol,masscol, maxcol
	real(dl) :: tmp(10000), massbound, x,y,z,vx,vy,vz,mass, rmin,rmax, deltar, r,v,vsq, vlos,vlossq, time0,time1,time2
	real(dl), allocatable :: meanvinbin(:), meanvsqinbin(:), meanvlosinbin(:), meanvlossqinbin(:)
	integer, allocatable :: numinbin(:)
	integer(8) :: numarg, nlines, nlines_remain, numrbin,irbin
	character(len=char_len) :: inputfile, outputfile, printstr, tmpstr1, tmpstr2, suffix,quanname
	logical :: cutabove
	
	
	printstr = 'Usage: EXE -inputfile intpufile -cutabove cutabove -massbound massbound-skiprow skiprow '//&
		'-xcol xcol -ycol ycol -zcol zcol -vxcol vxcol -vycol vycol -vzcol vzcol -masscol masscol '//&
		' -rmin rmin -rmax rmax -deltar deltar '//&
		' -quanname quanname -suffix suffix'//&
		'-outputfile outputfile'//&
		'### Compute RMS v (3-d and los component) above given mass cut; '

	! Default values
	cutabove = .true.
	xcol=1;ycol=2;zcol=3;vxcol=4;vycol=5;vzcol=6;masscol=7
	skiprow = 0; 
	inputfile = 'NONE'
	outputfile = 'NONE'
	quanname = 'RMSvel'
	
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
		elseif(trim(adjustl(tmpstr1)).eq.'-quanname') then
			read(tmpstr2,'(A)') quanname
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,'(A)') suffix
		elseif(trim(adjustl(tmpstr1)).eq.'-cutabove') then
			read(tmpstr2,*) cutabove
		elseif(trim(adjustl(tmpstr1)).eq.'-xcol') then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq.'-ycol') then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq.'-zcol') then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq.'-vxcol') then
			read(tmpstr2,*) vxcol
		elseif(trim(adjustl(tmpstr1)).eq.'-vycol') then
			read(tmpstr2,*) vycol
		elseif(trim(adjustl(tmpstr1)).eq.'-vzcol') then
			read(tmpstr2,*) vzcol
		elseif(trim(adjustl(tmpstr1)).eq.'-masscol') then
			read(tmpstr2,*) masscol
		elseif(trim(adjustl(tmpstr1)).eq.'-massbound') then
			read(tmpstr2,*) massbound
		elseif(trim(adjustl(tmpstr1)).eq.'-rmin') then
			read(tmpstr2,*) rmin
		elseif(trim(adjustl(tmpstr1)).eq.'-rmax') then
			read(tmpstr2,*) rmax
		elseif(trim(adjustl(tmpstr1)).eq.'-deltar') then
			read(tmpstr2,*) deltar
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	maxcol = max(xcol,ycol,zcol,vxcol,vycol,vzcol,masscol)
	if(maxcol > 10000) then
		print *, 'maxcol larger than 10000; increase size of tmp!: ', maxcol
		stop
	endif

	if(trim(adjustl(outputfile)).eq.'NONE') then
		write(tmpstr1,'(e17.4)') massbound
		outputfile = trim(adjustl(inputfile))//'.cutted-massgt'//trim(adjustl(tmpstr1))//'e11.'// &
			trim(adjustl(quanname))//'.info'//trim(adjustl(suffix))
	endif

	numrbin = ceiling((rmax - rmin) / deltar)
	write(*,'(A)') 'Inputfile = ', trim(adjustl(inputfile))
	write(*,'(A)') 'Outputfile = ', trim(adjustl(outputfile))
	write(*,'(i4,A,2f13.3,A,f10.3)'), numrbin, ' bins within r = ', rmin, rmax, '; deltar = ', deltar
	write(*,'(A,7i4)'), 'Cols of x,y,z,vx,vy,vz,mass: ', xcol,ycol,zcol,vxcol,vycol,vzcol,masscol
	allocate(numinbin(numrbin), meanvsqinbin(numrbin),  meanvinbin(numrbin), meanvlossqinbin(numrbin), meanvlosinbin(numrbin))

	! go through the file; calculate mean v, mean vsq...
	open(unit=1,file=inputfile,action='read')
	do i = 1, skiprow
		read(1,'(A)') tmpstr1
		write(2,'(A)') trim(adjustl(tmpstr1))
	enddo
	nlines = 0
	nlines_remain = 0
	do irbin = 1, numrbin
		numinbin(irbin) = 0
		meanvinbin(irbin) = 0.0
		meanvsqinbin(irbin) = 0.0
		meanvlosinbin(irbin) = 0.0
		meanvlossqinbin(irbin) = 0.0
	enddo
	call cpu_time(time0); time1=time0
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines + 1
		! read in mass, x,y,z,vx,vy,vz; calculate r
		mass=tmp(masscol)
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol)
		vx=tmp(vxcol); vy=tmp(vycol); vz=tmp(vzcol)
		r = sqrt(x*x+y*y+z*z);
		! check satisfying condition
		if((r>rmin.and.r<rmax) .and. ((cutabove.and. mass > massbound) .or. ((.not.cutabove).and.x<massbound)) ) then
			nlines_remain = nlines_remain +1			
			vsq = vx*vx+vy*vy+vz*vz
			v = sqrt(vsq);
			vlos = (vx*x + vy*y + vz*z) / r
			vlossq = vlos**2.0
			irbin = ceiling((r-rmin)/deltar)
			if(irbin .gt. numrbin .or. irbin .lt. 1) then
				print *, 'WARNING!!! Rare irbin: ', irbin
			endif
			irbin = max(irbin,1)
			irbin = min(irbin,numrbin)
			! count; sum...
			numinbin(irbin) = numinbin(irbin) + 1
			meanvinbin(irbin) = meanvinbin(irbin) + v
			meanvsqinbin(irbin) = meanvsqinbin(irbin) + vsq
			meanvlosinbin(irbin) = meanvlosinbin(irbin) + vlos
			meanvlossqinbin(irbin) = meanvlossqinbin(irbin) + vlossq
		endif
		call cpu_time(time2)
		if(time2-time1.ge.10.0) then
			write(*,'(f13.3,A,i10,A,i10)'), time2-time1, 'sec ellapes. ', &
				nlines, ' lines processed; in-condition #-line = ', nlines_remain
			time1=time2
		endif
		cycle
101		exit
	enddo
	close(1)

	do irbin = 1, numrbin
		meanvinbin(irbin) = meanvinbin(irbin) / dble(numinbin(irbin))
		meanvsqinbin(irbin) = meanvsqinbin(irbin) / dble(numinbin(irbin))
		meanvlosinbin(irbin) = meanvlosinbin(irbin) / dble(numinbin(irbin))
		meanvlossqinbin(irbin) = meanvlossqinbin(irbin) / dble(numinbin(irbin))
	enddo

	open(unit=2,file=outputfile)
	write(*,'(A,i10,A,i8,A,i4,e14.7,A,L3,A,i10,A,A,A)') ' Finishing processing ', nlines, ' lines (skipping first ', &
		skiprow,' rows). Cut column, massbound =  ', masscol, massbound, &
		', cutabove = ', cutabove, '. ',nlines_remain, ' lines participate statistics. Result written to ' , &
		trim(adjustl(outputfile)), &
		'; fmt: rmin, rmax, numinbin, meanvinbin, meanvsqinbin, Varvinbin, meanvlos, meanvsqlos, Varvlos'
	write(2,'(A,i10,A,i8,A,i4,e14.7,A,L3,A,i10,A,A,A)') ' Finishing processing ', nlines, ' lines (skipping first ', &
		skiprow,' rows). Cut column, massbound =  ', masscol, massbound, &
		', cutabove = ', cutabove, '. ',nlines_remain, ' lines participate statistics', '###', &
		' fmt: rmin, rmax, numinbin, meanvinbin, meanvsqinbin, Varvinbin, meanvlos, meanvsqlos, Varvlos'
	do irbin = 1, numrbin
		write(*,'(2f13.3,i10,6e15.7)') rmin+deltar*(irbin-1), min(rmin+deltar*irbin, rmax), &
			numinbin(irbin), meanvinbin(irbin), meanvsqinbin(irbin), &
			sqrt(meanvsqinbin(irbin) - meanvinbin(irbin)**2.0), &
			meanvlosinbin(irbin), meanvlossqinbin(irbin), &
			sqrt(meanvlossqinbin(irbin) - meanvlosinbin(irbin)**2.0)
		write(2,'(2f13.3,i10,5e15.7)') rmin+deltar*(irbin-1), min(rmin+deltar*irbin, rmax), &
			numinbin(irbin), meanvinbin(irbin), meanvsqinbin(irbin), &
			sqrt(meanvsqinbin(irbin) - meanvinbin(irbin)**2.0), &
			meanvlosinbin(irbin), meanvlossqinbin(irbin), &
			sqrt(meanvlossqinbin(irbin) - meanvlosinbin(irbin)**2.0)
	enddo
	close(2)
	deallocate(meanvsqinbin, numinbin, meanvinbin, meanvlosinbin, meanvlossqinbin)
end program main
