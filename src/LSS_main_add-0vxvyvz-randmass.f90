program main

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, printstr, inputfile, outputfile
        real(dl) :: log10massmin, log10massmax, nowlog10mass
        integer :: i, nlines

        printstr = 'Usgae: EXE inputfile log10massmin log10massmax'
        if(iarg().ne. 3) then
                print *, trim(adjustl(printstr))
                stop
        endif
        call getarg(1, inputfile)
        call getarg(2, tmpstr1)
        read(tmpstr1, *) log10massmin
        call getarg(3, tmpstr2)
        read(tmpstr2, *) log10massmax
        outputfile = trim(adjustl(inputfile))//'.0vxvyvz.log10massfrom'//trim(adjustl(tmpstr1))//'to'//trim(adjustl(tmpstr2))
        print *, '####################################################'
        print *, 'Add zero velocity, random mass to ', trim(adjustl(inputfile)), ', output to ', trim(adjustl(outputfile))
        if(.true.) then
                print *, 'Calling random seed...'
                call random_seed()
        endif

        ! processing files...
        nlines = 0
        print *, 'Processing files...'
        open(unit=100,file=inputfile,action='read')
        open(unit=200,file=outputfile,action='write')
        do while(.true.)
                read(100,'(A)',end=100) tmpstr1
                nlines = nlines+1
                call random_number(nowlog10mass)
                nowlog10mass = nowlog10mass * (log10massmax - log10massmin) + log10massmin
                write(200,'(A,3i2,e15.7)') trim(adjustl(tmpstr1)), 0,0,0, 10**nowlog10mass
                cycle
100             exit                
        enddo
                
        print *, 'Finishing processing ', nlines, 'lines.'



end program
