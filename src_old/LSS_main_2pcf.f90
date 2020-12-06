

program LSS_main

!use mpi
use LSS_2pcf
implicit none

	character(len=char_len) :: inputfilename, outputfilename, tmpstr1, tmpstr2
        character(len=1000) :: printstr, decomp
	real(dl) :: omegam, w, rmin,rmax
	logical :: printinfo, has_outputfilename
        integer  :: numtbin,numrbin, numarg, i, dec
        real(dl), allocatable ::  counts(:, :)
	! Initialize MPI
	!call mpi_init(ierr)
	!call mpi_comm_size(mpi_comm_world,nproc,ierr)
	!call mpi_comm_rank(mpi_comm_world,myid,ierr)

        printstr = ' Usage: ./LSS_2pcf '//&
                '-input intpufilename -output outputfilname -rmin rmin '//&
                '-rmax rmax -numrbin numrbin -numtbin numtbin -printinfo printinfo -decomp SIGPI/SMU'

	! datatype: -1 means x,y,z,w (see LSS_settings_init.f90)
	gb_i_datatype = gb_dt_xyzw
	! rantype: no ran
	gb_i_rantype = gb_noran
	! read in settings
	omegam=om_dft; w=w_dft
	print *, 'Default cosmology: omegam, w = ', om_dft, w_dft
	numarg = iargc()
	if(numarg .le. 0) then
		print *, 'small numarg: ', numarg
		write(*,'(A)') printstr
		stop
	endif
	
        printinfo = .true.
        has_outputfilename = .false.
        decomp = 'SMU'
	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-decomp') then
			read(tmpstr2,'(A)') decomp
		elseif(trim(adjustl(tmpstr1)).eq.'-output') then
	        		read(tmpstr2,'(A)') outputfilename
			        has_outputfilename = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-printinfo') then
			read(tmpstr2,*) printinfo
		elseif(trim(adjustl(tmpstr1)).eq.'-rmin') then
			read(tmpstr2,*) rmin
		elseif(trim(adjustl(tmpstr1)).eq.'-rmax') then
			read(tmpstr2,*) rmax
		elseif(trim(adjustl(tmpstr1)).eq.'-numrbin') then
			read(tmpstr2,*) numrbin 
		elseif(trim(adjustl(tmpstr1)).eq.'-numtbin') then
			read(tmpstr2,*) numtbin 
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') printstr
			stop
		endif
	enddo
	
	if(.not.has_outputfilename) then
		outputfilename = trim(adjustl(inputfilename))
	        write(tmpstr1,'(f30.3)') rmin
		write(tmpstr2,'(f30.3)') rmax
	        outputfilename = trim(adjustl(outputfilename))//'.rmin'//trim(adjustl(tmpstr1))//'.rmax'//trim(adjustl(tmpstr2))
		write(tmpstr1,*) numrbin
		write(tmpstr2,*) numtbin
	        outputfilename = trim(adjustl(outputfilename))//'.rbin'//trim(adjustl(tmpstr1))//'.tbin'//trim(adjustl(tmpstr2))
	endif

	print *, '#####################################'
	print *, ' Settings'
	write(*,'(A,A)') '   inputfilename: ', trim(adjustl(inputfilename))
	write(*,'(A,A)') '   outputfilename: ', trim(adjustl(outputfilename))
	print *, '  printinfo = ', printinfo
	print *, '#####################################'

        if (trim(adjustl(decomp)) .eq. 'SMU' ) then
		dec = 0
                allocate(counts(numrbin,numtbin))
		call Tpcf(inputfilename, rmin,rmax,numrbin, numtbin, counts, dec, printinfo)
	elseif (trim(adjustl(decomp)) .eq. 'SIGPI' ) then
		dec = 1
                allocate(counts(numrbin,numrbin))
		call Tpcf(inputfilename, rmin,rmax,numrbin, numrbin, counts, dec, printinfo)
	else
		print *, 'ERROR! decomp must be SIGPI or SMU! : decomp = ', trim(adjustl(decomp))
		stop
	endif

	!call mpi_barrier(mpi_comm_world,ierr)
	!call mpi_finalize(ierr)
end program LSS_main
