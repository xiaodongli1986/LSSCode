program main_mpi_random_select

use mpi

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, filelist, suffixstr = '', printstr, file1, file2
        character(len=char_len), allocatable :: files(:)
        integer :: i, nfile, ifile, iline, iwriteline
        real(dl) :: randrat = 0.001, x

        ! mpi variables
	integer :: ierr, nproc, myid

       	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)
        call random_seed()
        write(*,*) 'myid, nproc = ', myid, nproc
	
        printstr = "Randomly select a part of the input file. Usage: LSS_mpi_random_select -inputfile inputfilename -randrat your_rat -suffixstr your_suffix"
	if(iargc().le.1) then
		print *, printstr
		stop
	endif

	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
                elseif(trim(adjustl(tmpstr1)).eq."-randrat ") then
			read(tmpstr2,"(A)") randrat
		elseif(trim(adjustl(tmpstr1)).eq."-suffixstr") then
			read(tmpstr2,"(A)") suffixstr
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

        if(trim(adjustl(suffixstr)) .eq. '') then
                write(suffixstr,'(f14.7)') randrat
                suffixstr = '.randrat'//trim(adjustl(suffixstr))
        endif

        write(*,'(A)') 'inputfile = '//trim(adjustl(inputfile))
        write(*,'(A,f14.7)') 'radrat    = ', randrat
        write(*,'(A,A)')'suffixstr = ', trim(adjustl(suffixstr))

        filelist = 'LSS_mpi_random_select.filelist'
        write(*,'(A,A)') 'Files to be processed stored in: ', trim(adjustl(filelist))

        call system('ls '//trim(adjustl(inputfile))//' > '//trim(adjustl(filelist)))

        open(file = filelist, action='read', unit = 10000)
        nfile = 0
        do while(.true.)
                read(10000,*,end=100) tmpstr1
                nfile = nfile + 1
        enddo
100     close(10000)

        allocate(files(nfile))
        open(file = filelist, action='read', unit = 10000)
        do ifile = 1, nfile
                read(10000,*,end=100) files(ifile)
        enddo
        close(10000)


        do ifile = myid+1, nfile, nproc
                file1 = files(ifile)
                file2 = trim(adjustl(file1))//trim(adjustl(suffixstr))
                write(*, '(A,i4,i4,A,A)')'myid, nproc = ', myid, nproc, '; processing ', trim(adjustl(file1))
                open(file = file1, action='read', unit = 12903480)
                open(file = file2, action='write', unit = 12903481)
                iline = 0; iwriteline = 0
                do while(.true.)
                        read(12903480,'(A)',end=101) tmpstr1
                        iline = iline + 1
                        call random_number(x)
                        if (x < randrat) then
                                write(12903481, '(A)') trim(adjustl(tmpstr1))
                                iwriteline = iwriteline + 1
                        endif
                enddo
101             close(10000)
                write(*, '(A,i4,i4,";",i12,A,A,";",i12,A,A)')'myid, nproc = ', myid, nproc, iline, ' lines read from ', trim(adjustl(file1)), iwriteline, ' lines write to ', trim(adjustl(file2))
        enddo



	call mpi_barrier(mpi_comm_world,ierr)
	call mpi_finalize(ierr)




end program
