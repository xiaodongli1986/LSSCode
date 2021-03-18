module lightcone_boxsplit_tools_mpi
use LSS_cosmo_funs
implicit none
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


        subroutine read_in_gadget_npar(filename, npar)
                character(len=char_len) :: filename
                integer *8 :: npar
                type(head) :: head_info
                open(file=trim(adjustl(filename)), unit=1000,action='read', form='unformatted')
                read(1000) head_info
                npar = head_info%npart(2)
        end subroutine read_in_gadget_npar

       	subroutine read_in_colanpar(filename, npar)
                character(len=char_len) :: filename, filenamenpar
                integer *8 :: npar
                filenamenpar = trim(adjustl(filename))//".npar"
		!gfortran fmt
                !open(file=trim(adjustl(filenamenpar)),unit=1000,action='read',form='binary',access='stream'); read(1000) npar; close(1000)
                open(file=trim(adjustl(filenamenpar)),unit=1000,action='read',form='unformatted',access='stream') 
                read(1000) npar; close(1000)
                print*, 'inputfile.npar is: ', npar, trim(adjustl(filenamenpar))
                !print*, 'npar is: ', npar
        end subroutine read_in_colanpar

	subroutine read_in_coladata_slow(filename, npar, rmid, ids, xyzvs)
                character(len=char_len) :: filename
                integer *8 :: npar,  ids(npar), i 
                real :: xyzvs(6,npar), rmid


                !open(file=trim(adjustl(filename)),unit=100000001, action='read', form='binary', access='stream') ! gfortran fmt
                open(file=trim(adjustl(filename)),unit=100000001, action='read', form='unformatted', access='stream')
                read(100000001) rmid
                do i =1, npar
                        read(100000001) ids(i); read(100000001) xyzvs(1:6,i)
                enddo
                close(100000001) 
                !write(*, '(A,f12.3)')       '   read-in rmid = ', rmid
                !write(*, '(A,i15,i15)')     '   ids1 begin/end with ' , ids(1), ids(npar)
                !write(*, '(A,6(f12.3))')    '   xyzvs1 begin with ', xyzvs(:,1)
                !write(*, '(A,6(f12.3))')    '   xyzvs1  end  with ', xyzvs(:,npar)
	end subroutine read_in_coladata_slow

	subroutine read_in_coladata(filename, npar, rmid, ids, xyzvs)
                character(len=char_len) :: filename
                integer *8 :: npar,  ids(npar), i
                real :: xyzvs(6,npar), id_xyzvs(8,npar), rmid

                !print *, ' TEST BEGIN of read_in_coladata'

                !open(file=trim(adjustl(filename)),unit=1001, action='read', form='binary', access='stream') ! gfortran fmt
                open(file=trim(adjustl(filename)),unit=1001, action='read', form='unformatted', access='stream')
                read(1001) rmid
                !print *, 'TEST read in rmin = ', rmid
                read(1001) id_xyzvs ! cheat: read-in integer*8 as two real numbers ! since we do not use ids, it is OK
                xyzvs = id_xyzvs(3:8,:)
                write(*, '(A,f12.3)')       '   read-in rmid = ', rmid
                !write(*, '(A,6(f12.4))')    '   xyzvs begin with ', xyzvs(:,220:240)
                write(*, '(A,6(f12.3))')    '   xyzvs  end  with ', xyzvs(:,npar)
                close(1001)
                !print *, ' TEST END  of read_in_coladata'
        end subroutine read_in_coladata

end module lightcone_boxsplit_tools_mpi



program main_mpi_lightcone_boxsplit

use mpi
use lightcone_boxsplit_tools_mpi


implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, &
		inputfilelist, outputfile, printstr, suffix, infofile, outputname, headfile='', inputtypestr, par_random_output_file
        character(len=char_len), allocatable :: inputfiles(:), outputfiles(:)
        integer, parameter:: type_lpicola=0, type_cola = 1, type_gadget = 2, par_random_output_unit = 10000000
        integer :: i, nbox, nfile, ifile, ibox, ixs(2), iys(2), izs(2), nwrites(3), &
		i1,i2,i3, ix,iy,iz, iwrite, flag, nown, inputtype = type_lpicola, n_par_random_output=0, icycle
        real(dl) :: overlap_distance, xyzmin, xyzmax, x,y,z, boxsize,parmass,omegam, &
		redshift,hubble,w,  par_random_output_rat=0.01, tmpx
        real(dl), allocatable :: xyz_ranges(:,:,:)
        logical :: add1, binary_IO = .false., block_write = .true.
        double precision :: tmpdoubles(64)


        type(head) :: headinfo

        integer*4 :: EOF, block1, nowunit
        integer*8 :: ntotal
        real, allocatable::  xyzvs(:,:), tmp_xyzvs(:,:)
        integer, allocatable :: iwrite_counts(:,:,:)
        !integer, allocatable :: nparticles(:)

        real :: cola_rmid
        integer *8 :: cola_npar
        integer *8, allocatable :: cola_ids(:)
        real, allocatable :: cola_xyzvs(:,:)
        integer, parameter :: data_max_len = 3000000

        type data_arrays
                integer :: n
                integer, allocatable :: ids(:)
                real, allocatable :: xyzs(:,:)
                real, allocatable :: vxvyvzs(:,:)
                integer :: ntotal
        end type data_arrays
        type(data_arrays), allocatable :: all_data(:,:,:)

                

!	 mpi variables
	integer :: ierr, nproc, myid, iproc

	printstr = ' ## Usage: EXE   -inputfilelist ?   -outputname ?   -nbox ?   -overlap_distance ?   '//&
		'-xyzmin ?   -xyzmax ?  -add1 T/F   -binary_IO F    -headfile headfile     -inputtype '//&
		'cola/lpicola (cola or lpicola?)    -par_random_output_rat 0.01  ### Cut the lpicola lightcone '//&
		'sample into small boxes with overlapping region.  Convenient for halo-finding.    '

	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)
        iproc = myid + 1
        print *, 'myid, iproc, nproc = ', myid, iproc, nproc
	
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
            print *, 'now read in nbox... tmpstr2 = ', trim(adjustl(tmpstr2))
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
		elseif(trim(adjustl(tmpstr1)).eq."-inputtype") then
			read(tmpstr2,*) inputtypestr
		elseif(trim(adjustl(tmpstr1)).eq."-par_random_output_rat") then
			read(tmpstr2,*) par_random_output_rat
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo


        if(trim(adjustl(headfile)).ne.'') then
                print *
                write(*,'(A)') ' Read in parameters from lpicola headfile  (This will ignore hand-given values of '//&
			'omegam, w, boxsize, parmass, hubble!)'
                write(*,'(A,A)') '        ',trim(adjustl(headfile))
                print *
                open(unit=10002,file=headfile,action='read')
                read(10002,*) tmpstr1
                read(10002,*) ntotal, boxsize, parmass, redshift, omegam, hubble, w;
                !write(*, '(A,f12.3)')       '   read-in omegam = ', omegam
                !write(*, '(A,f12.3)')       '   read-in hubble = ', hubble
                !write(*, '(A,f12.3)')       '   read-in ntotal = ', ntotal
                if(w.ge.-0.001) then
                        print *, '   set w as -1 (found w = ',w,')'
                        w = -1.
                endif
                close(10002)
        endif
        print *, nbox, xyzmin, xyzmax
        write(tmpstr2,*) nbox
        write(tmpstr3,'(f10.1)') overlap_distance
        write(tmpstr4,'(f10.1)') xyzmin
        write(tmpstr5,'(f10.1)') xyzmax


        if(trim(adjustl(inputtypestr)).eq.'lpicola') then
                inputtype = type_lpicola
                write(*,'(A,A)') ' set inputtype as lpicola: inputtypestr= ', trim(adjustl(inputtypestr))
        elseif(trim(adjustl(inputtypestr)).eq.'cola') then
                inputtype = type_cola
                write(*,'(A,A)') ' set inputtype as cola: inputtypestr= ', trim(adjustl(inputtypestr))
        elseif(trim(adjustl(inputtypestr)).eq.'gadget') then
                inputtype = type_gadget
                write(*,'(A,A)') ' set inputtype as gadget: inputtypestr= ', trim(adjustl(inputtypestr))
        else
                inputtype = type_lpicola
                write(*,'(A,A)') ' set inputtype as lpicola: inputtypestr= ', trim(adjustl(inputtypestr))
        endif



	if(trim(adjustl(suffix)).eq."") then
                suffix = '.nbox'//trim(adjustl(tmpstr2))//'_overlap'//trim(adjustl(tmpstr3))&
			//'_xyz'//trim(adjustl(tmpstr4))//'to'//trim(adjustl(tmpstr5))
	endif
        infofile = trim(adjustl(outputname))//trim(adjustl(suffix))//'.info'

        ! par_random_output_file
        write(tmpstr1, '(f10.3)') par_random_output_rat
        !par_random_output_file = trim(adjustl(outputname))//'_'//trim(adjustl(tmpstr1))//'particles'
        write(tmpstr2, *) iproc
        par_random_output_file = trim(adjustl(outputname))//'_'//trim(adjustl(tmpstr1))//'particles.iproc'//trim(adjustl(tmpstr2)) ! mpi
        write(*,'(A)') ' output '//trim(adjustl(tmpstr1))//' of all particles to '//trim(adjustl(par_random_output_file))

        call count_line_number(inputfilelist, nfile)
        print *, 'In total ', nfile, 'files'
        allocate(inputfiles(nfile),outputfiles(nbox*nbox*nbox))

        open(unit=100,file=inputfilelist)
        do ifile = 1, nfile
                read(100,'(A)') inputfiles(ifile)
        enddo
        close(100)
        print *, 'Finishing read in ', nfile, 'filenames. We will read-in them.'

        write(tmpstr3, *) iproc ! mpi
        do ibox = 1, nbox*nbox*nbox
               write(tmpstr2, *) ibox
               outputfiles(ibox) = trim(adjustl(outputname))//trim(adjustl(suffix))//'.ibox'//&
                   trim(adjustl(tmpstr2)) //'.iproc'//trim(adjustl(tmpstr3)) ! mpi
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
           open(unit=par_random_output_unit, file=trim(adjustl(par_random_output_file)), action='write') !par_random_output
 

           print *, 'Start read-in files.'
           !do ifile = 1, nfile 
           do ifile = iproc, nfile, nproc  ! mpi
                !write(*,'(A,A,A)') '     opening ', trim(adjustl(inputfiles(ifile))), ' for read...'
                write(*,'(A,i5,A,A,A)') '(iproc ',iproc, ')     opening ', trim(adjustl(inputfiles(ifile))), ' for read...'
                open(unit=nbox*nbox*nbox+100000, file=trim(adjustl(inputfiles(ifile))), action='read')
                do while (.true.)
                        read(nbox*nbox*nbox+100000, '(A)', end=100) tmpstr1
                        if(add1) tmpstr1 = trim(adjustl(tmpstr1))//' '//'1'
                        read(tmpstr1, *) x,y,z

           !open(unit=par_random_output_unit, file=trim(adjustl(par_random_output_file)), action='write') !par_random_output
                        call random_number(tmpx) !par_random_output
                        if(tmpx<par_random_output_rat .and. x>0 .and. y>0 .and. z>0) then
                                write(par_random_output_unit,'(A)')  trim(adjustl(tmpstr1)) !par_random_output
                                n_par_random_output = n_par_random_output + 1
                        endif
                        
           !            close(par_random_output_unit)

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
           close(par_random_output_unit)
        else
           print *, 'IO_format = binary; will output following files...'
                                !open(unit=nowunit+100, file='test.binary', action='write', access='stream', form='binary')
           print *, 'Start read-in files.'

           !#############################################################
           ! count #-lines
           ntotal=0
           !open(unit=par_random_output_unit,file=trim(adjustl(par_random_output_file)), action='write', access='stream',form='binary') !par_random_output ! gfortran fmt
           open(unit=par_random_output_unit,file=trim(adjustl(par_random_output_file)), &
	   	action='write', access='stream',form='unformatted') !par_random_output
           !do ifile = 1, nfile !mpi
           do ifile = iproc, nfile, nproc
              write(*,'(A,i5,A,A,A,"(",i4," of",i4,")")') '(iproc',iproc,')     opening ', &
	           trim(adjustl(inputfiles(ifile))), ' for countline...', ifile, nfile
              if(inputtype .eq. type_lpicola) then
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

                                        call random_number(tmpx) !par_random_output
                                        if(tmpx<par_random_output_rat .and. x>0 .and. y>0 .and. z>0) then
                                                write(par_random_output_unit)  tmp_xyzvs(:,i)
                                                n_par_random_output = n_par_random_output + 1
                                        endif

                                        call determine_iwrites(x, y, z, overlap_distance, &
						xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
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
             elseif(inputtype .eq. type_cola) then
                call read_in_colanpar(inputfiles(ifile), cola_npar)
                !print*, 'input_file is', inputfiles(ifile)
                if(cola_npar.eq.0) then
                        print *, '             cycle because of npar =0: ', int(cola_npar); cycle
                endif
                !print *, 'TEST A'
                allocate(cola_ids(cola_npar), cola_xyzvs(6,cola_npar))
                call read_in_coladata(inputfiles(ifile), cola_npar, cola_rmid, cola_ids, cola_xyzvs)
                !stop !! maqlin debug
                print*, 'cola_npar is', cola_npar
                do i = 1, cola_npar
                        x=cola_xyzvs(1,i); y=cola_xyzvs(2,i); z=cola_xyzvs(3,i)

                        call random_number(tmpx) !par_random_output
                        if(tmpx<par_random_output_rat .and. x>0 .and. y>0 .and. z>0) then
                                write(par_random_output_unit)  cola_xyzvs(:,i)
                                n_par_random_output = n_par_random_output + 1
                        endif

                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                        if(flag.eq.0) cycle
                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                iwrite_counts(ix,iy,iz) = iwrite_counts(ix,iy,iz) + 1
                                !iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                !write(iwrite, '(A)') trim(adjustl(tmpstr1))
                        enddo; enddo; enddo
                enddo
                deallocate(cola_ids, cola_xyzvs) ! maqlin debug
                !print *, 'TEST B'
                ntotal  = ntotal + cola_npar
             endif
           enddo
           close(par_random_output_unit) ! par_random_output

           print *; print *, '#########################'
           write(*,'(A,i5,A,i15,A,i8,A)') '(', iproc, ') Finishing count-line. In total we have ', ntotal, &
                ' particles in ', nfile, ' files.'
           write(*,'(A,i5,A,i15,A,A)') '(', iproc, ' Finishing particle output (stage 1). In total we have ', &
	   	n_par_random_output, ' particles write to ', trim(adjustl(par_random_output_file))

           do ibox = 1, nbox*nbox*nbox
               write(*,'(A,i5,": ",5x, A)') '* iproc = ',iproc, trim(adjustl(outputfiles(ibox)))
               !open(unit=ibox+200000, file = trim(adjustl(outputfiles(ibox))), action='write', access='stream', form='binary') ! gfortran fmt
               open(unit=ibox+200000, file = trim(adjustl(outputfiles(ibox))), action='write', access='stream', form='unformatted')
           enddo
           !#############################################################
           ! write head to files
           headinfo%npart = 0; headinfo%npart(2) = ntotal
           headinfo%mass = 0.; headinfo%mass(2) = parmass/1.0e10
           headinfo%time = 1.; headinfo%redshift = 0.;
           headinfo%flag_sfr = 0; headinfo%flag_feedback=0;
           headinfo%npartTotal = 0; headinfo%npartTotal(2) = ntotal
           headinfo%flag_cooling = 0; headinfo%num_files = 1
           headinfo%BoxSize = boxsize;
           headinfo%Omega0 = omegam;
           headinfo%OmegaLambda = 1. - omegam ; 
           print *, 'WARNING! set omegal as 1-omegam: ', headinfo%OmegaLambda
           print *, 'WARNING! set omegal as 1-omegam: ', headinfo%OmegaLambda
           print *, 'WARNING! set omegal as 1-omegam: ', headinfo%OmegaLambda
           headinfo%HubbleParam = hubble
           headinfo%fill = 0
           
           write(*,'(A,i5,A)') ' * iproc = ', iproc, ':  writing heads... '
           do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                print *, '(iproc',iproc,') i1,i2,i3, counts = ', i1, i2,i3, iwrite_counts(i1,i2,i3)
                headinfo%npart(2) = iwrite_counts(i1,i2,i3)
                iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                write(iwrite) 256
                write(iwrite) headinfo
                write(iwrite) 256
                write(iwrite) iwrite_counts(i1,i2,i3)*3*4
           enddo; enddo; enddo

           ! write everything once-for-all
           if(block_write)then
             allocate(all_data(nbox,nbox,nbox))
             do i1 = 1,nbox; do i2=1,nbox; do i3=1,nbox
                        all_data(i1,i2,i3)%n = 0
                        all_data(i1,i2,i3)%ntotal = 0
                        !allocate(all_data(i1,i2,i3).ids(iwrite_counts(i1,i2,i3)))
                        !allocate(all_data(i1,i2,i3).xyzs(3,iwrite_counts(i1,i2,i3)))
                        !allocate(all_data(i1,i2,i3).vxvyvzs(3,iwrite_counts(i1,i2,i3)))
                        allocate(all_data(i1,i2,i3)%ids(data_max_len))
                        allocate(all_data(i1,i2,i3)%xyzs(3,data_max_len))
                        allocate(all_data(i1,i2,i3)%vxvyvzs(3,data_max_len))
             enddo; enddo; enddo
             ! x,y,z,vx,vy,vz info written to all_data storage
             ! icycle == 1: write xyz; icycle==2: write v; icycle==3: write id
             do icycle = 1, 3 !doicycle
              write(*,*) '(iproc ',iproc,")  block write start!!!: icycle = ", icycle
              write(*,*) '(iproc ',iproc,")  block write start!!!: icycle = ", icycle
              write(*,*) '(iproc ',iproc,")  block write start!!!: icycle = ", icycle
              ntotal=0
              !do ifile = 1, nfile !doifile
              do ifile = iproc, nfile, nproc !doifile ! mpi
                write(*,'(A,i5,A,A,A,"(",i4," of",i4,")")') '(iproc',iproc,')     opening ', &
			trim(adjustl(inputfiles(ifile))), ' for block read-in...', ifile, nfile
                if(inputtype.eq.type_lpicola) then
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
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, &
						nwrites, ixs, iys, izs, nbox, flag)
                                        if(flag.eq.0) cycle
                                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                                ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                                all_data(ix,iy,iz)%n = all_data(ix,iy,iz)%n + 1; nown = all_data(ix,iy,iz)%n 
                                                all_data(ix,iy,iz)%ntotal = all_data(ix,iy,iz)%ntotal + 1
                                                all_data(ix,iy,iz)%xyzs(1:3,nown) = tmp_xyzvs(1:3,i)
                                                all_data(ix,iy,iz)%vxvyvzs(1:3,nown) = tmp_xyzvs(4:6,i)
                                                all_data(ix,iy,iz)%ids(nown) = ntotal+i
                                                if(all_data(ix,iy,iz)%n .eq. data_max_len) then
                                                        iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                        if(icycle.eq.1) then
                                                                write(iwrite) all_data(ix,iy,iz)%xyzs
                                                        elseif(icycle.eq.2) then
                                                                write(iwrite) all_data(ix,iy,iz)%vxvyvzs
                                                        elseif(icycle.eq.3) then
                                                                write(iwrite) all_data(ix,iy,iz)%ids
                                                        endif
                                                        all_data(ix,iy,iz)%n = 0
                                                endif
                                        enddo; enddo; enddo
                                enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                  enddo
                  close(nowunit)
               elseif(inputtype.eq.type_cola) then
                 call read_in_colanpar(inputfiles(ifile), cola_npar)
                 if(cola_npar.eq.0) then
                        print *, 'cycle because of npar =0: ', cola_npar; cycle
                 endif
                 !print*, '#####################################'
                 allocate(cola_ids(cola_npar), cola_xyzvs(6,cola_npar))
                 !print*, 'cola_npar is :', cola_npar ! maqlin debug
                 !print *, 'TEST C'
                 call read_in_coladata(inputfiles(ifile), cola_npar, cola_rmid, cola_ids, cola_xyzvs)
                 do i = 1, cola_npar
                        x=cola_xyzvs(1,i); y=cola_xyzvs(2,i); z=cola_xyzvs(3,i)
                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, nwrites, ixs, iys, izs, nbox, flag)
                        if(flag.eq.0) cycle
                        do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
                                        ix = ixs(i1); iy = iys(i2); iz = izs(i3)
                                        all_data(ix,iy,iz)%n = all_data(ix,iy,iz)%n + 1; nown = all_data(ix,iy,iz)%n 
                                        all_data(ix,iy,iz)%ntotal = all_data(ix,iy,iz)%ntotal + 1
                                        all_data(ix,iy,iz)%xyzs(1:3,nown) = cola_xyzvs(1:3,i)
                                        all_data(ix,iy,iz)%vxvyvzs(1:3,nown) = cola_xyzvs(4:6,i)
                                        all_data(ix,iy,iz)%ids(nown) = ntotal+i
                                        if(all_data(ix,iy,iz)%n .eq. data_max_len) then
                                                        iwrite = (ix-1)*nbox*nbox + (iy-1)*nbox + iz + 200000
                                                        if(icycle.eq.1) then
                                                                write(iwrite) all_data(ix,iy,iz)%xyzs
                                                        elseif(icycle.eq.2) then
                                                                write(iwrite) all_data(ix,iy,iz)%vxvyvzs
                                                        elseif(icycle.eq.3) then
                                                                write(iwrite) all_data(ix,iy,iz)%ids
                                                        endif
                                                        all_data(ix,iy,iz)%n = 0
                                        endif
                        enddo; enddo; enddo
                        ntotal  = ntotal + cola_npar
                 enddo
                 deallocate(cola_ids, cola_xyzvs)
                 ! need to store cola particles to this structure!!! all_data...
                 ! need to : generate headfile for cola;; and run for test... 
                 ! need to read this code carefully...
                 ! need to add inputs to cola_halo; 
                 !print *, 'TEST D'
               endif

              enddo !enddoifile

              ! write residual particles in all_data
              do i1=1,nbox; do i2=1,nbox; do i3=1,nbox
                        ix=i1; iy=i2; iz=i3;
                        iwrite = (i1-1)*nbox*nbox + (i2-1)*nbox + i3 + 200000
                        if(icycle.eq.1) then
                          if(all_data(ix,iy,iz)%n > 0) write(iwrite) all_data(ix,iy,iz)%xyzs(:,1:all_data(ix,iy,iz)%n)
                          write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                          write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                        elseif(icycle.eq.2) then
                          if(all_data(ix,iy,iz)%n > 0) write(iwrite) all_data(ix,iy,iz)%vxvyvzs(:,1:all_data(ix,iy,iz)%n)
                          write(iwrite) iwrite_counts(i1,i2,i3)*3*4
                          write(iwrite) iwrite_counts(i1,i2,i3)*4
                        else
                          if(all_data(ix,iy,iz)%n > 0) write(iwrite) all_data(ix,iy,iz)%ids(1:all_data(ix,iy,iz)%n)
                          write(iwrite) iwrite_counts(i1,i2,i3)*4
                          close(iwrite)
                        endif
                        all_data(ix,iy,iz)%n = 0 ! set as zero, be ready for next time output
              enddo; enddo; enddo

             enddo !endodicycle
             do i1 = 1,nbox; do i2=1,nbox; do i3=1,nbox
                        write(*,'(A,i5,A,i3,i3,i3, A,i10,i10,A,i10,i10)') '(iproc',iproc,&
                                ')  Finishing all_data storage: i1,i2,i3 = ', &
				i1,i2,i3, '; #-par, iwrite(i1,i2,i3) = ', all_data(i1,i2,i3)%ntotal/3, &
				iwrite_counts(i1,i2,i3), '; check difference (should be zero): ', &
				all_data(i1,i2,i3)%ntotal - iwrite_counts(i1,i2,i3)*3
             enddo; enddo; enddo
           endif



           if(block_write) then
                goto 200 !stop
           endif

           !#############################################################
           ! write x,y,z to files
           write(*,'(A,i5,A)') ' * (iproc ',iproc,') writing x,y,z (one-by-one mode)... '
           ntotal=0
           !do ifile = 1, nfile
           do ifile = iproc, nfile, nproc ! mpi
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
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, &
						nwrites, ixs, iys, izs, nbox, flag)
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
           print*, '##########################################'
           write(*,'(A,i5,A)') ' * (iproc',iproc,') writing vx,vy,vz (one-by-one mode)... '
           ntotal=0
           !do ifile = 1, nfile
           do ifile = iproc, nfile, nproc ! mpi
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
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, &
						nwrites, ixs, iys, izs, nbox, flag)
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
           write(*,'(A,i5,A)') ' * (iproc',iproc,')writing ids (one-by-one mode)... '
           ntotal=0
           !do ifile = 1, nfile
           do ifile = iproc, nfile, nproc ! mpi
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
                                        call determine_iwrites(x, y, z, overlap_distance, xyzmin, xyzmax, &
						nwrites, ixs, iys, izs, nbox, flag)
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

        ! output particles
200     allocate(tmp_xyzvs(7,n_par_random_output))
        !open(unit=par_random_output_unit,file=trim(adjustl(par_random_output_file)),action='read', access='stream',form='binary') ! gfortran fmt
        open(unit=par_random_output_unit,file=trim(adjustl(par_random_output_file)),action='read', &
		access='stream',form='unformatted')
        read(par_random_output_unit) tmp_xyzvs(1:6,:)
        tmp_xyzvs(7,:) = 1.0
        close(par_random_output_unit)

        open(file=trim(adjustl(par_random_output_file)),unit=par_random_output_unit,form='unformatted',action='write')
        tmpdoubles = 0.
        tmpdoubles(1) = n_par_random_output
        tmpdoubles(2) = headinfo%boxsize
        tmpdoubles(3) = parmass 
        tmpdoubles(4) = headinfo%redshift
        tmpdoubles(5) = headinfo%Omega0
        tmpdoubles(6) = headinfo%HubbleParam


        write(par_random_output_unit) tmpdoubles

        write(*,*)
        write(*,'(A,i5,A,i15,A)') '(iproc ',iproc,')', n_par_random_output, ' particles selected for output ...'
        write(*,'(A,i5,A)') ' (iproc',iproc,')  staring write particles to '//trim(adjustl(par_random_output_file))//&
		'  format = block+[x,y,z,vx,vy,vz,1.0]*noutput+block ...'
        write(par_random_output_unit) tmp_xyzvs(1:7,:)
        close(par_random_output_unit)




        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_finalize(ierr)



end program

! What treatment we have adopted now:

