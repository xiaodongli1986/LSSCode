!### Python code for a quick check of the all-sky random program
!!LSS_make-allsky-random 'ran.allsky' 0 200000 0 1
!X, Y, Z = XYZfromdata(np.loadtxt('ran.allsky0'))
!R = list_xyz_to_r(X, Y, Z)
!X, Y, Z = plt.hist(R, bins=100)
!X = [X[row]/((Y[row+1]+Y[row])*0.5)**2.0 for row in range(len(X))]
!plt.plot(Y[0:len(Y)-1], X)

program main 

use LSS_cosmo_funs

	implicit none
	! HR3 file, data
        integer :: i, iran, nran, ranseed, puts(33) ! gfortran fmt puts
        real(dl) :: rmin, rmax, r,x,y,z
        character(len=char_len) :: randatafile, tmpstr

	print *, '## Generating randoms uniformly distributed within all sky'
        i = iargc()
        if(i.ne.5) then
        	print *, 'Useage:'
        	print *, '  EXE randatafile ranseed nran rmin rmax'
        	stop
        endif
        
        call getarg(1,randatafile)
        call getarg(2,tmpstr)
        read(tmpstr,*) ranseed
        randatafile = trim(adjustl(randatafile))//trim(adjustl(tmpstr))
        call getarg(3,tmpstr)
        read(tmpstr,*) nran
        call getarg(4,tmpstr)
        read(tmpstr,*) rmin
        call getarg(5,tmpstr)
        read(tmpstr,*) rmax
        
        write(*,'(i9,A,f10.2,A,f10.2,A,A)') nran, ' randoms in all sky, r ', rmin, &
		' to', rmax, ', file: ', trim(adjustl(randatafile))
        
	puts = 0
	puts(1) = 0
	puts(2) = ranseed
        call random_seed(put=puts)

        iran = 0
        
        open(unit=1,file=randatafile)
!        open(unit=2,file=trim(adjustl(randatafile))//'.1percent')
        do while(iran.lt.nran)
	        call random_number(x)
	        call random_number(y)
	        call random_number(z)
	        x = x*2*rmax - rmax
	        y = y*2*rmax - rmax
	        z = z*2*rmax - rmax
	        r = sqrt(x*x+y*y+z*z)
	        if(r<rmax .and. r>rmin) then
	        	write(1,'(3(e15.7,1x))') x,y,z
	        	iran = iran+1
!	                if(mod(iran,100).eq.1) then
!	        	        write(2,'(3(e15.7,1x))') x,y,z
 !       		endif
	        endif
	enddo
	close(1); 
        !close(2)
end program main
