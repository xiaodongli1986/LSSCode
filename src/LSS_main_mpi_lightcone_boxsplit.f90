program main_mpi_lightcone_boxsplit

use mpi

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, inputfilelist, outputfile, printstr, suffixstr, infofile, outputname
        character(len=char_len), allocatable :: inputfiles(:), outputfiles(:,:)
	integer :: i, nbox, nfile
        real(dl) :: overlap_distance, xyzmin, xyzmax
        real(dl), allocatable :: xyz_ranges(:,:,:)

	! mpi variables
	integer :: ierr, nproc, myid

	printstr = '### Cut a sample into small samples with overlapping. Usage: mpirun -np ??? EXE -inputfile ? -outputfile ? -nbox ? -overlap_distance ? -xyzmin ? -xyzmax ?'

	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)
	
	printstr = "Now it is empty!"
	if(iargc().le.1) then
		print *, printstr
		stop
	endif

	suffix = ""
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfilelist") then
			read(tmpstr2,"(A)") inputfilelist
		elseif(trim(adjustl(tmpstr1)).eq."-suffix") then
			read(tmpstr2,"(A)") suffix
		elseif(trim(adjustl(tmpstr1)).eq."-outputname") then
			read(tmpstr2,"(A)") outputname
		elseif(trim(adjustl(tmpstr1)).eq."-nbox") then
			read(tmpstr2,*) nbox
		elseif(trim(adjustl(tmpstr1)).eq."-overlap_distance") then
			read(tmpstr2,*) overlap_distance
		elseif(trim(adjustl(tmpstr1)).eq."-xyzmin") then
			read(tmpstr2,*) xyzmin
		elseif(trim(adjustl(tmpstr1)).eq."-xyzmax") then
			read(tmpstr2,*) xyzmax
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

        write(tmpstr2,*) nbox
        write(tmpstr3,'(5.1f)') overlap_distance
        write(tmpstr4,'(5.1f)') xyzmin
        write(tmpstr5,'(5.1f)') xyzmax




	if(trim(adjustl(suffix)).eq."") then
                suffix = '.nbox'//trim(adjustl(tmpstr2))//'_overlap'//trim(adjustl(tmpstr3))//'_xyz'//trim(adjustl(xyzmin))//'to'//trim(adjustl(xyzmax))
	endif
        infofile = trim(adjustl(intputfile))//trim(adjustl(suffix))//'.info'

        call count_line_number(inputfile, nfile)
        allocate(inputfiles(nfile),outputfiles(nbox**3,nproc))

        open(unit=100,file=inputfilelist)
        do ifile = 1, nfile
                read(100,'(A)') inputfiles(ifile)
        enddo

        do ibox = 1, nbox*nbox*nbox
                do iproc = 0, nproc-1
                        write(tmpstr2, *) ibox
                        write(tmpstr3, *) iproc
                        outputfiles(ifile,iproc) = trim(adjustl(outputname))//trim(adjustl(suffix))//'.ibox'//trim(adjustl(tmpstr2))//'.iproc'//trim(adjustl(tmpstr3))
                enddo
        enddo


	print *, 'This is an empty program mpi_lightcone_boxsplit!!'

end program

! What treatment we have adopted now:

