! Aug 29: We consider rescale of volume due to change of the cosmology!!! (dovoleff)


program main

use LSS_cosmo_funs

implicit none

	real(dl) :: rmin,rmax,dr, redshift,omegam,w,dist,boxdistance, x,y,z,r,vx,vy,vz,v,vr,costheta, xmin,xmax,ymin,ymax,zmin,zmax,vdotr,&
                radiusmin, radiusmax, tmp(1000), &
		zobs, shiftdist,shiftrat, shiftedcord, volscale, voldft, volnow, rescalefac
	integer :: xcol,ycol,zcol,vxcol,vycol,vzcol,maxcol, i,j, numshiftback, numshiftbackorig, skiprow
	integer(8) :: numarg, nlines
	character(len=char_len) :: printstr, inputfilename, outputname, &
		outputfilenamer,outputfilenameinfo, tmpstr1,tmpstr2, suffix
	
	printstr = 'Usage: EXE -input intpufilename -rmin rmin -rmax rmax'//&
		' -xcol xcol -ycol ycol -zcol zcol '//&
		'-vxcol vxcol -vycol vycol -vzcol vzcol '//&
		'-omegam omegam -w w '//&
		'-skiprow skiprow -suffix suffix'

	! Default values
	rmin = -1.0e30; rmax = 1.0e30;
	xcol=1; ycol=2; zcol=3; vxcol=4; vycol=5; vzcol=6; skiprow=0
	omegam = 0.26; w = -1.0; redshift = 0.0
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
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-rmin') then
			read(tmpstr2,*) rmin
		elseif(trim(adjustl(tmpstr1)).eq.'-rmax') then
			read(tmpstr2,*) rmax
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
		elseif(trim(adjustl(tmpstr1)).eq.'-omegam' .or. trim(adjustl(tmpstr1)).eq.'-om') then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,*) suffix
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo


	maxcol = max(xcol,ycol,zcol,vxcol,vycol,vzcol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif

	dr = rmax-rmin
	
	gb_h = 0.73_dl;

	call cosmo_funs_init(.true.)
	gb_omegam = omegam; gb_w = w
	call de_calc_comovr()
	
	if(trim(adjustl(suffix)).ne.'') then
		outputname = trim(adjustl(inputfilename))//'.shiftr'
	else
		outputname = inputfilename
	endif
	
	outputfilenamer = trim(adjustl(outputname))//'.shiftr'
	outputfilenameinfo = trim(adjustl(outputname))//'.shiftr.info'

	open(unit=100,file=outputfilenameinfo)

	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,2f16.7)')   '   range of shell radius = ', rmin, rmax
	write(*,'(A,4f16.7)') 	'   omegam, w = ', real(omegam), real(w)
	write(*,'(A,6i3)')  	'   cols of x,y,z,vx,vy,vz: ', xcol,ycol,zcol,vxcol,vycol,vzcol
	write(*,'(A,i5)')	'   skip rows ', skiprow

	write(100,'(A)')   	'  Settings:'
	write(100,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(100,'(A,2f16.7)')   '   range of shell radius = ', rmin, rmax
	write(100,'(A,4f16.7)') 	'   omegam, w = ', real(omegam), real(w)
	write(100,'(A,6i3)')  	'   cols of x,y,z,vx,vy,vz: ', xcol,ycol,zcol,vxcol,vycol,vzcol
	write(100,'(A,i5)')	'   skip rows ', skiprow
	print *, '#####################################'

	open(unit=40,file=outputfilenamer)
	open(unit=5,file=trim(adjustl(outputfilenamer))//'.rs_zs_zobss_shifts')


	nlines =0; numshiftback=0; numshiftbackorig=0
	xmin = 1.0e30; ymin = 1.0e30; zmin = 1.0e30; radiusmin = 1.0e30 
	xmax = -xmin; ymax = -ymin; zmax = -zmin; radiusmax = -radiusmin
	open(unit=1,file=inputfilename,action='read')
	do i = 1, skiprow
		read(1,*) tmpstr1
	enddo
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines+1
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); r= sqrt(x*x+y*y+z*z);
		vx=tmp(vxcol); vy=tmp(vycol); vz=tmp(vzcol); 
		xmin = min(x,xmin); xmax = max(x,xmax)
		ymin = min(y,ymin); ymax = max(y,ymax)
		zmin = min(z,zmin); zmax = max(z,zmax)
		radiusmin = min(r,radiusmin); radiusmax = max(r,radiusmax)
		if(.true.) then
			!print *, 'velocity at z direction is really larger than velocity at a single direction... please check it!'; stop
                        v = sqrt(vx**2 + vy**2 + vz**2)
                        !print *, 'v = ', v
                        vdotr = vx*x + vy*y + vz*z
                        !print *, 'vdotr = ', vdotr
                        !print *, 'r = ', r
			costheta = vdotr / v / r
                        !print *, 'costheta = ', costheta
			vr = v * costheta
                        redshift = de_zfromintpl(r)
			zobs = redshift + vr * (1.0_dl+redshift) / const_c 
			shiftrat = de_get_comovr(zobs) / r 
			!print *, 'r_obs, dist, shiftr = ', de_get_comovr(zobs), dist, shiftdist

			write(40,'(3e15.7)') x*shiftrat, y*shiftrat, z*shiftrat
		endif
		write(5,'(4e15.7)') r, redshift, zobs, r*(shiftrat-1.)
		!write(*,'(9e15.7)') r, redshift, vx, vy, vz, costheta, sqrt(vx**2+vy**2+vz**2),  vr, zobs, r*(shiftrat-1.)
		cycle
101		exit
	enddo
	
	close(40)

	write(*,'(A,i10,A,f12.3,f12.3,A)') ' Finishing processing ', nlines, ' lines. Radius within range of (', radiusmin, radiusmax, ')'
	write(100,'(A,i10,A,f12.3,f12.3,A)') ' Finishing processing ', nlines, &
		' lines. Radius within range of (', radiusmin, radiusmax, ')'

	close(100)

end program main
