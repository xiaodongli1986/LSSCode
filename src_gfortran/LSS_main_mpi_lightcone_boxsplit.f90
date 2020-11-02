program main_mpi_lightcone_boxsplit

use mpi

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3, tmpstr4, &
		tmpstr5, inputfilelist, outputfile, printstr, suffix, infofile, outputname
        character(len=char_len), allocatable :: inputfiles(:), outputfiles(:)
        integer :: i, nbox, nfile, ifile, ibox, iproc, ixs(2), iys(2), izs(2), nwrites(3), i1,i2,i3, ix,iy,iz, iwrite
        real(dl) :: overlap_distance, xyzmin, xyzmax, x,y,z
        real(dl), allocatable :: xyz_ranges(:,:,:)
        logical :: add1

	! mpi variables
	integer :: ierr, nproc, myid

	printstr = '### Cut a sample into small samples with overlapping. Usage: mpirun '//&
		'-np ??? EXE -inputfilelist ? -outputname ?  -nbox ? -overlap_distance ? -xyzmin ? -xyzmax ? -add1 T/F'

	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)
        iproc = myid + 1
        print *, 'myid, iproc = ', myid, iproc
	
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
		elseif(trim(adjustl(tmpstr1)).eq."-add1") then
			read(tmpstr2,*) add1
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

        write(tmpstr2,*) nbox
        write(tmpstr3,'(f10.1)') overlap_distance
        write(tmpstr4,'(f10.1)') xyzmin
        write(tmpstr5,'(f10.1)') xyzmax




	if(trim(adjustl(suffix)).eq."") then
                suffix = '.nbox'//trim(adjustl(tmpstr2))//'_overlap'//&
			trim(adjustl(tmpstr3))//'_xyz'//trim(adjustl(tmpstr4))//'to'//trim(adjustl(tmpstr5))
	endif
        infofile = trim(adjustl(outputname))//trim(adjustl(suffix))//'.info'

        call count_line_number(inputfilelist, nfile)
        print *, 'In total ', nfile, 'files'
        allocate(inputfiles(nfile),outputfiles(nbox*nbox*nbox))

        open(unit=100,file=inputfilelist)
        do ifile = 1, nfile
                read(100,'(A)') inputfiles(ifile)
        enddo
        close(100)
        print *, 'Finishing read in ', nfile, 'filenames'

        write(tmpstr3, *) iproc
        do ibox = 1, nbox*nbox*nbox
               write(tmpstr2, *) ibox
               outputfiles(ibox) = trim(adjustl(outputname))//trim(adjustl(suffix))//'.ibox'//&
	           trim(adjustl(tmpstr2))//'.iproc'//trim(adjustl(tmpstr3))
               print *, trim(adjustl(outputfiles(ibox)))
               open(unit=ibox+200000, file = trim(adjustl(outputfiles(ibox))), action='write')
        enddo
 

        do ifile = iproc, nfile, nproc
                print *, 'opening ', trim(adjustl(inputfiles(ifile))), ' for read...'
                open(unit=nbox*nbox*nbox+100000, file=trim(adjustl(inputfiles(ifile))), action='read')
                print *, 'file opened...'
                do while (.true.)
                        read(nbox*nbox*nbox+100000, '(A)', end=100) tmpstr1
                        if(add1) tmpstr1 = trim(adjustl(tmpstr1))//' '//'1'
                        read(tmpstr1, *) x,y,z
                        ixs(1) = int((x-overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1 
                        ixs(2) = int((x+overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        iys(1) = int((y-overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        iys(2) = int((y+overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        izs(1) = int((z-overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        izs(2) = int((z+overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        ixs(1) = max(ixs(1),1); ixs(2) = min(ixs(2),nbox)
                        iys(1) = max(iys(1),1); iys(2) = min(iys(2),nbox)
                        izs(1) = max(izs(1),1); izs(2) = min(izs(2),nbox)
                        nwrites = 1
                        if(ixs(2) .ne. ixs(1)) nwrites(1) = 2
                        if(iys(2) .ne. iys(1)) nwrites(2) = 2
                        if(izs(2) .ne. izs(1)) nwrites(3) = 2
                        do i1 = 1, nwrites(1)
                        do i2 = 1, nwrites(2)
                        do i3 = 1, nwrites(3)
                                ix = ixs(i1)
                                iy = iys(i2)
                                iz = izs(i3)
                                iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                write(iwrite, '(A)') trim(adjustl(tmpstr1))
                        enddo
                        enddo
                        enddo
                        cycle
100                     exit
                enddo
                close(nbox*nbox*nbox+100000)
        enddo
        do ibox = 1, nbox*nbox*nbox
               close(ibox + 200000)
        enddo

        
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_finalize(ierr)



end program

! What treatment we have adopted now:

