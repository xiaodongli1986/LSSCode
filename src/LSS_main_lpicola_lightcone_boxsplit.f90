module lightcone_boxsplit_tools
use LSS_cosmo_funs
implicit none
contains

        subroutine determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                real(dl) :: x, y, z, overlap_distance, xyzmin, xyzmax
                integer :: ixs(2), iys(2), izs(2), nwrites(3), nbox, flag
                if(x<xyzmin .or. x>xyzmax .or. y<xyzmin .or. y>xyzmax .or. z<xyzmin .or. z>xyzmax) then
                        flag = 0; return
                else
                        flag = 1
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
                endif
        end subroutine determine_iwrites
end module lightcone_boxsplit_tools



program main_mpi_lightcone_boxsplit

!use mpi
use lightcone_boxsplit_tools


implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, inputfilelist, outputfile, printstr, suffix, infofile, outputname, headfile=''
        character(len=char_len), allocatable :: inputfiles(:), outputfiles(:)
        integer :: i, nbox, nfile, ifile, ibox, ixs(2), iys(2), izs(2), nwrites(3), i1,i2,i3, ix,iy,iz, iwrite, flag, nown
        real(dl) :: overlap_distance, xyzmin, xyzmax, x,y,z, boxsize,parmass,omegam, redshift,hubble,w
        real(dl), allocatable :: xyz_ranges(:,:,:)
        logical :: add1, binary_IO = .false., block_write = .true.


        !gadget head
        type head
                integer :: npart(6)
                double precision :: mass(6)
                double precision :: time
                double precision :: redshift
                integer :: flag_sfr, flag_feedback, npartTotal(6), flag_cooling, num_files
                double precision :: BoxSize, Omega0, OmegaLambda, HubbleParam;
                integer :: fill(24)
        end type head
        type(head) :: headinfo

        integer*4 :: ntotal, EOF, block1, nowunit
        real, allocatable::  xyzvs(:,:), tmp_xyzvs(:,:)
        integer, allocatable :: iwrite_counts(:,:,:)
        !integer, allocatable :: nparticles(:)

        type data_arrays
                integer :: n
                integer, allocatable :: ids(:)
                real, allocatable :: xyzs(:,:)
                real, allocatable :: vxvyvzs(:,:)
        end type data_arrays
        type(data_arrays), allocatable :: all_data(:,:,:)

                

	! mpi variables
!	integer :: ierr, nproc, myid

	printstr = ' ## Usage: EXE   -inputfilelist ?   -outputname ?   -nbox ?   -overlap_distance ?   -xyzmin ?   -xyzmax ?  -add1 T/F   -binary_IO F    -headfile headfile     ### Cut the lpicola lightcone sample into small boxes with overlapping region.  Convenient for halo-finding.    '

!	call mpi_init(ierr)
!	call mpi_comm_size(mpi_comm_world,nproc,ierr)
!	call mpi_comm_rank(mpi_comm_world,myid,ierr)
!        iproc = myid + 1
!        print *, 'myid, iproc = ', myid, iproc
	
	if(iargc().le.1) then
		write(*,'(A)') trim(adjustl(printstr))
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
		elseif(trim(adjustl(tmpstr1)).eq."-binary_IO") then
			read(tmpstr2,*) binary_IO
                elseif(trim(adjustl(tmpstr1)).eq.'-headfile' .or. trim(adjustl(tmpstr1)).eq.'-lpicola_headfile')  then
                        read(tmpstr2,'(A)') headfile
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo


        if(trim(adjustl(headfile)).ne.'') then
                print *
                write(*,'(A)') ' Read in parameters from lpicola headfile  (This will ignore hand-given values of omegam, w, boxsize, parmass, hubble!)'
                write(*,'(A,A)') '        ',trim(adjustl(headfile))
                print *
                open(unit=10002,file=headfile,action='read')
                read(10002,*) tmpstr1
                read(10002,*) ntotal, boxsize, parmass, redshift, omegam, hubble, w;
                if(w.ge.-0.001) then
                        print *, '   set w as -1 (found w = ',w,')'
                        w = -1.
                endif
                close(10002)
        endif

        write(tmpstr2,*) nbox
        write(tmpstr3,'(f10.1)') overlap_distance
        write(tmpstr4,'(f10.1)') xyzmin
        write(tmpstr5,'(f10.1)') xyzmax



	if(trim(adjustl(suffix)).eq."") then
                suffix = '.nbox'//trim(adjustl(tmpstr2))//'_overlap'//trim(adjustl(tmpstr3))//'_xyz'//trim(adjustl(tmpstr4))//'to'//trim(adjustl(tmpstr5))
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
        print *, 'Finishing read in ', nfile, 'filenames. We will read-in them.'

!        write(tmpstr3, *) iproc
        
        do ibox = 1, nbox*nbox*nbox
               write(tmpstr2, *) ibox
               outputfiles(ibox) = trim(adjustl(outputname))//trim(adjustl(suffix))//'.ibox'//trim(adjustl(tmpstr2)) !//'.iproc'//trim(adjustl(tmpstr3))
        enddo
        allocate(iwrite_counts(nbox,nbox,nbox))
        iwrite_counts = 0

        ! ascii format for input and output files!
        if(.not.binary_IO) then

           print *, 'IO format = ASCII'
           do ibox = 1, nbox*nbox*nbox
               print *, trim(adjustl(outputfiles(ibox)))
               open(unit=ibox+200000, file = trim(adjustl(outputfiles(ibox))), action='write')
           enddo
 

           print *, 'Start read-in files.'
           do ifile = 1, nfile
                write(*,'(A,A,A)') '     opening ', trim(adjustl(inputfiles(ifile))), ' for read...'
                open(unit=nbox*nbox*nbox+100000, file=trim(adjustl(inputfiles(ifile))), action='read')
                do while (.true.)
                        read(nbox*nbox*nbox+100000, '(A)', end=100) tmpstr1
                        if(add1) tmpstr1 = trim(adjustl(tmpstr1))//' '//'1'
                        read(tmpstr1, *) x,y,z

                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                        if(flag.eq.0) cycle
                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                write(iwrite, '(A)') trim(adjustl(tmpstr1))
                        enddo; enddo; enddo
                        cycle


                        !if(x<xyzmin .or. x>xyzmax .or. y<xyzmin .or. y>xyzmax .or. z<xyzmin .or. z>xyzmax) cycle
                        !ixs(1) = int((x-overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1 
                        !ixs(2) = int((x+overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        !iys(1) = int((y-overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        !iys(2) = int((y+overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        !izs(1) = int((z-overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        !izs(2) = int((z+overlap_distance/2.-xyzmin) / ((xyzmax-xyzmin)/float(nbox))) +1
                        !ixs(1) = max(ixs(1),1); ixs(2) = min(ixs(2),nbox)
                        !iys(1) = max(iys(1),1); iys(2) = min(iys(2),nbox)
                        !izs(1) = max(izs(1),1); izs(2) = min(izs(2),nbox)
                        !nwrites = 1
                        !if(ixs(2) .ne. ixs(1)) nwrites(1) = 2
                        !if(iys(2) .ne. iys(1)) nwrites(2) = 2
                        !if(izs(2) .ne. izs(1)) nwrites(3) = 2
                        !do i1 = 1, nwrites(1)
                        !do i2 = 1, nwrites(2)
                        !do i3 = 1, nwrites(3)
                        !        ix = ixs(i1)
                        !        iy = iys(i2)
                        !        iz = izs(i3)
                        !        iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                        !        write(iwrite, '(A)') trim(adjustl(tmpstr1))
                        !enddo
                        !enddo
                        !enddo
                        !cycle
100                     exit
                enddo
                close(nbox*nbox*nbox+100000)
           enddo
           do ibox = 1, nbox*nbox*nbox
               close(ibox + 200000)
           enddo
        else
           print *, 'IO_format = binary; will output following files...'
                                !open(unit=nowunit+100, file='test.binary', action='write', access='stream', form='binary')
           print *, 'Start read-in files.'

           !#############################################################
           ! count #-lines
           ntotal=0
           do ifile = 1, nfile
                write(*,'(A,A,A,"(",i4," of",i4,")")') '     opening ', trim(adjustl(inputfiles(ifile))), ' for countline...', ifile, nfile
                nowunit = nbox*nbox*nbox+100000
                open(unit=nowunit, file=trim(adjustl(inputfiles(ifile))), action='read',form='unformatted')
                do
                        read(nowunit,IOSTAT=EOF) block1
                        !print *, 'block1 = ', block1
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program'
                        else if (EOF .lt. 0) then
                               ! If EOF < 0 then we have read all the chunks
                                exit
                        else
                                !print *, 'before allocating'
                                !call system('sleep 1')
                                allocate(tmp_xyzvs(6,block1)); read(nowunit) tmp_xyzvs;
                                !print *, 'after allocating'
                                !call system('sleep 1')
                                do i = 1, block1
                                        x=tmp_xyzvs(1,i); y=tmp_xyzvs(2,i); z=tmp_xyzvs(3,i); 
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                                        if(flag.eq.0) cycle
                                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                                iwrite_counts(ix,iy,iz) = iwrite_counts(ix,iy,iz) + 1
                                                !iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                !write(iwrite, '(A)') trim(adjustl(tmpstr1))
                                        enddo; enddo; enddo
                                enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                enddo
                close(nowunit)
           enddo

           print *; print *, '#########################'
           write(*,'(A,i15,A,i8,A)') ' Finishing count-line. In total we have ', ntotal, ' particles in ', nfile, ' files.'

           ! write everything once-for-all
           if(block_write)then
             allocate(all_data(nbox,nbox,nbox))
             do i1 = 1,nbox; do i2=1,nbox; do i3=1,nbox
                        all_data(i1,i2,i3).n = 0
                        allocate(all_data(i1,i2,i3).ids(iwrite_counts(i1,i2,i3)))
                        allocate(all_data(i1,i2,i3).xyzs(3,iwrite_counts(i1,i2,i3)))
                        allocate(all_data(i1,i2,i3).vxvyvzs(3,iwrite_counts(i1,i2,i3)))
             enddo; enddo; enddo
             ! x,y,z,vx,vy,vz info written to all_data storage
             ntotal=0
             do ifile = 1, nfile
                write(*,'(A,A,A,"(",i4," of",i4,")")') '     opening ', trim(adjustl(inputfiles(ifile))), ' for block read-in...', ifile, nfile
                nowunit = nbox*nbox*nbox+100000
                open(unit=nowunit, file=trim(adjustl(inputfiles(ifile))), action='read',form='unformatted')
                do
                        read(nowunit,IOSTAT=EOF) block1
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program'
                        else if (EOF .lt. 0) then
                               ! If EOF < 0 then we have read all the chunks
                                exit
                        else
                                allocate(tmp_xyzvs(6,block1)); read(nowunit) tmp_xyzvs;
                                do i = 1, block1
                                        x=tmp_xyzvs(1,i); y=tmp_xyzvs(2,i); z=tmp_xyzvs(3,i); 
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                                        if(flag.eq.0) cycle
                                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                                all_data(ix,iy,iz).n = all_data(ix,iy,iz).n + 1; nown = all_data(ix,iy,iz).n 
                                                all_data(ix,iy,iz).xyzs(1:3,nown) = tmp_xyzvs(1:3,i)
                                                all_data(ix,iy,iz).vxvyvzs(1:3,nown) = tmp_xyzvs(4:6,i)
                                                all_data(ix,iy,iz).ids(nown) = ntotal+i
                                                !iwrite_counts(ix,iy,iz) = iwrite_counts(ix,iy,iz) + 1
                                                !iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                !write(iwrite, '(A)') trim(adjustl(tmpstr1))
                                        enddo; enddo; enddo
                                enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                enddo
                close(nowunit)
             enddo

             do i1 = 1,nbox; do i2=1,nbox; do i3=1,nbox
                        write(*,'(A,i3,i3,i3, A,i10,i10,A,i10,i10)') '  Finishing all_data storage: i1,i2,i3 = ', i1,i2,i3, '; #-par, iwrite(i1,i2,i3) = ', all_data(i1,i2,i3).n, iwrite_counts(i1,i2,i3), '; check difference (should be zero): ', all_data(i1,i2,i3).n - iwrite_counts(i1,i2,i3) 
             enddo; enddo; enddo
           endif

           do ibox = 1, nbox*nbox*nbox
               write(*,'(5x, A)') trim(adjustl(outputfiles(ibox)))
               open(unit=ibox+200000, file = trim(adjustl(outputfiles(ibox))), action='write', access='stream', form='binary')
           enddo


           !#############################################################
           ! write head to files
           headinfo.npart = 0; headinfo.npart(2) = ntotal
           headinfo.mass = 0.; headinfo.mass(2) = parmass/1.0e10
           headinfo.time = 1.; headinfo.redshift = 0.;
           headinfo.flag_sfr = 0; headinfo.flag_feedback=0;
           headinfo.npartTotal = 0; headinfo.npartTotal(2) = ntotal
           headinfo.flag_cooling = 0; headinfo.num_files = 1
           headinfo.BoxSize = boxsize;
           headinfo.Omega0 = omegam;
           headinfo.OmegaLambda = 1. - omegam ; 
           print *, 'WARNING! set omegal as 1-omegam: ', headinfo.OmegaLambda
           print *, 'WARNING! set omegal as 1-omegam: ', headinfo.OmegaLambda
           print *, 'WARNING! set omegal as 1-omegam: ', headinfo.OmegaLambda
           headinfo.HubbleParam = hubble
           headinfo.fill = 0
           
           write(*,'(A)') ' * writing heads... '
           do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                print *, 'i1,i2,i3, counts = ', i1, i2,i3, iwrite_counts(i1,i2,i3)
                headinfo.npart(2) = iwrite_counts(i1,i2,i3)
                iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                write(iwrite) 256
                write(iwrite) headinfo
                write(iwrite) 256
                write(iwrite) iwrite_counts(i1,i2,i3)*3*4
           enddo; enddo; enddo
           if(block_write) then
                write(*,'(A)') ' * writing x,y,z, vx,vy,vz (block_write mode)... '
                do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                        iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                        write(*,'(A,A)') ' ** writing to file', trim(adjustl(outputfiles(iwrite-200000)))
                        write(iwrite) all_data(i1,i2,i3).xyzs
                        write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                        write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                        write(iwrite) all_data(i1,i2,i3).vxvyvzs
                        write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                        write(iwrite) iwrite_counts(i1,i2,i3)*4
                        write(iwrite) all_data(i1,i2,i3).ids
                        write(iwrite) iwrite_counts(i1,i2,i3)*4
                        close(iwrite)
                enddo; enddo; enddo
                stop
           endif

           !#############################################################
           ! write x,y,z to files
           write(*,'(A)') ' * writing x,y,z (one-by-one mode)... '
           ntotal=0
           do ifile = 1, nfile
                !write(*,'(A,A,A)') '     opening ', trim(adjustl(inputfiles(ifile))), ' ...'
                nowunit = nbox*nbox*nbox+100000
                open(unit=nowunit, file=trim(adjustl(inputfiles(ifile))), action='read',form='unformatted')
                !print *, 'file opened...'
                do
                        read(nowunit,IOSTAT=EOF) block1
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program'
                        else if (EOF .lt. 0) then
                                exit
                        else
                                allocate(tmp_xyzvs(6,block1)); read(nowunit) tmp_xyzvs;
                                do i = 1, block1
                                        x=tmp_xyzvs(1,i); y=tmp_xyzvs(2,i); z=tmp_xyzvs(3,i); 
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                                        if(flag.eq.0) cycle
                                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                                iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                write(iwrite) tmp_xyzvs(1:3,i); !write(iwrite) tmp_xyzvs(2,i); write(iwrite) tmp_xyzvs(3,i)
                                        enddo; enddo; enddo
                                enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                enddo
                close(nowunit)
           enddo
           do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                write(iwrite) iwrite_counts(i1,i2,i3)*3*4
           enddo; enddo; enddo

           !#############################################################
           ! write vx,vy,vz to files
           write(*,'(A)') ' * writing vx,vy,vz (one-by-one mode)... '
           ntotal=0
           do ifile = 1, nfile
                !write(*,'(A,A,A)') '     opening ', trim(adjustl(inputfiles(ifile))), ' ...'
                nowunit = nbox*nbox*nbox+100000
                open(unit=nowunit, file=trim(adjustl(inputfiles(ifile))), action='read',form='unformatted')
                !print *, 'file opened...'
                do
                        read(nowunit,IOSTAT=EOF) block1
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program'
                        else if (EOF .lt. 0) then
                                exit
                        else
                                allocate(tmp_xyzvs(6,block1)); read(nowunit) tmp_xyzvs;
                                do i = 1, block1
                                        x=tmp_xyzvs(1,i); y=tmp_xyzvs(2,i); z=tmp_xyzvs(3,i); 
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                                        if(flag.eq.0) cycle
                                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                                iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                write(iwrite) tmp_xyzvs(4:6,i); !write(iwrite) tmp_xyzvs(5,i); write(iwrite) tmp_xyzvs(6,i)
                                        enddo; enddo; enddo
                                enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                enddo
                close(nowunit)
           enddo
           do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                write(iwrite) iwrite_counts(i1,i2,i3)*4
           enddo; enddo; enddo

           !#############################################################
           ! write ids to files
           write(*,'(A)') ' * writing ids (one-by-one mode)... '
           ntotal=0
           do ifile = 1, nfile
                !write(*,'(A,A,A)') '     opening ', trim(adjustl(inputfiles(ifile))), ' ...'
                nowunit = nbox*nbox*nbox+100000
                open(unit=nowunit, file=trim(adjustl(inputfiles(ifile))), action='read',form='unformatted')
                !print *, 'file opened...'
                do
                        read(nowunit,IOSTAT=EOF) block1
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program'
                        else if (EOF .lt. 0) then
                                exit
                        else
                                allocate(tmp_xyzvs(6,block1)); read(nowunit) tmp_xyzvs;
                                do i = 1, block1
                                        x=tmp_xyzvs(1,i); y=tmp_xyzvs(2,i); z=tmp_xyzvs(3,i); 
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                                        if(flag.eq.0) cycle
                                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                                iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                write(iwrite) ntotal + i
                                        enddo; enddo; enddo
                                enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                enddo
                close(nowunit)
           enddo
           do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                write(iwrite) iwrite_counts(i1,i2,i3)*4
           enddo; enddo; enddo
        endif

        
        do ibox = 1, nbox*nbox*nbox
               close(ibox+200000)
        enddo
!        call mpi_barrier(mpi_comm_world,ierr)
!        call mpi_finalize(ierr))



end program

! What treatment we have adopted now:

