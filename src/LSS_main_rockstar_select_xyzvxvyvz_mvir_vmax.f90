program main_rockstar_select_xyzvxvyvz_mvir_vmax

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, headfile
	integer :: iline
        real :: tmpreal(15)
	
        printstr = 'Usage: LSS_rockstar_select_vxvyvz_mvir_vmax inputfile'
        write(*,'(A)') "(LSS_rockstar_select_xyzvxvyvz_mvir_vmax) select x,y,z, vx,vy,vz, mvir, vmax out of the rockstar halo sample.  print the head" 
        if(iargc().le.0) then

		print *, printstr
		stop
	endif

        call getarg(1,inputfile)
        outputfile = trim(adjustl(inputfile))//'.xyzvxvyvz_mvir_vmax'
        headfile = trim(adjustl(inputfile))//'.head'


        open(unit=1001,file=trim(adjustl(inputfile)),action='read')
        open(unit=1002,file=trim(adjustl(outputfile)),action='write')

        iline = 1
        do while(.true.)
                read(1001,'(A)',end=100) tmpstr1
                if(tmpstr1(1:1)=='#') cycle
                !print *, tmpstr1
                read(tmpstr1,*) tmpreal
                write(1002,'(8e15.7)') tmpreal(9:14), tmpreal(3),tmpreal(6)
                iline = iline+1
                cycle
100             exit
        enddo
        close(1001); close(1002)
        write(*,'(i12,A,A)') iline, ' lines written to ', trim(adjustl(outputfile))
        write(*,'(A)') 'head written to ...'

end program
