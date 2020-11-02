

program LSS_main

!use mpi
use LSS_BSK

	character(len=char_len) :: inputfilename, outputfilename, tmpstr1, tmpstr2, &
		ngldir = '~/SparseFilaments/code/ngl-beta/bin/getNeighborGraph'
	real(dl) :: omegam, w, beta = 5.0, minimalrcut = -1.0e30, maximalrcut = 1.0e30, dropstep=0.1
	logical :: printinfo = .true., do_nglcrosscheck = .false., has_outputfilename = .false.
	integer :: numNNB = 100, numarg, nproc, ierr, myid, i, numdrop=1

	! Initialize MPI
	!call mpi_init(ierr)
	!call mpi_comm_size(mpi_comm_world,nproc,ierr)
	!call mpi_comm_rank(mpi_comm_world,myid,ierr)

	! testing format I/O

	!write(*, '(A)') format_string(3, '(e14.7,1x)')
	!write(*, format_string(3, '(e14.7,1x)')) 1.32, 3.45, 5.67
	!write(*, '(A)') format_string(3, '(e14.7,1x)', '  " 1" ')
	!write(*, format_string(3, '(e14.7,1x)', '" 1"')) 1.32, 3.45, 5.67
	!stop


	! datatype: -1 means x,y,z (see LSS_settings_init.f90)
	gb_i_datatype = gb_dt_xyz

	! rantype: no ran
	gb_i_rantype = gb_noran

	! read in settings
	omegam=100.0; w=100.0;
	numarg = iargc()
	if(numarg .le. 0) then
		print *, 'small numarg: ', numarg
		write(*,'(A)') ' Usage: ./BSK_calc -omdft omdft -wdft wdft -om omegam -w w '//&
			'-input intpufilename -output outputfilname -beta beta '//&
			'-printinfo printinfo -numNNB numNNB -minimalrcut minimalrcut '//&
			'-maximalrcut maximalrcut -nglcrosscheck do_nglcrosscheck '//&
			'-ngldir ngldir -numdrop numdrop -dropstep dropstep'
                write(*,'(A)') 'omdft, wdft 样本的真实宇宙学，你的 xyz 对应的宇宙学'
                !write(*,'(A)') 'om w 可以把数据转到其他宇宙学背景下（相当于加AP），如果 no AP 就跟 omdft wdft 一样'
                write(*,'(A)') 'om w  | convert to cosmology with these parameters'
                write(*,'(A)') 'inputfilename 输入数据文件'
                write(*,'(A)') 'outputfilename 应当可以不给，按照这个名字输出一个文件夹，里面有各种结果'
                write(*,'(A)') 'beta 取3比较中庸'
                write(*,'(A)') 'numNNB 越大搜索范围越大越精确但是会慢很多，取到结果不变就可以'
                write(*,'(A)') 'minimalrct 除非你要截取一个 shell，设成0'
                write(*,'(A)') 'maximalrcut 除非你要截取一个 shell，设成很大比如 1e30'
                write(*,'(A)') 'nglcrosscheck False'
                write(*,'(A)') 'numdrop dropstep 都不用管'
                write(*,'(A)') ' Note for BSKInfo: gb_BSKs(i)%r, de_zfromintpl(gb_BSKs(i)%r), gb_BSKs(i)%length, gb_BSKs(i)%mu'
                write(*,'(A)')
		stop
	endif
	
	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-omdft') then
			read(tmpstr2,*) om_dft
		elseif(trim(adjustl(tmpstr1)).eq.'-wdft') then
			read(tmpstr2,*) w_dft
                elseif(trim(adjustl(tmpstr1)).eq.'-om') then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-output') then
			read(tmpstr2,'(A)') outputfilename
			has_outputfilename = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-beta') then
			read(tmpstr2,*) beta
		elseif(trim(adjustl(tmpstr1)).eq.'-printinfo') then
			read(tmpstr2,*) printinfo
		elseif(trim(adjustl(tmpstr1)).eq.'-numNNB') then
			read(tmpstr2,*) numNNB
		elseif(trim(adjustl(tmpstr1)).eq.'-minimalrcut') then
			read(tmpstr2,*) minimalrcut
		elseif(trim(adjustl(tmpstr1)).eq.'-maximalrcut') then
			read(tmpstr2,*) maximalrcut
		elseif(trim(adjustl(tmpstr1)).eq.'-nglcrosscheck') then
			read(tmpstr2,*) do_nglcrosscheck 
		elseif(trim(adjustl(tmpstr1)).eq.'-ngldir') then
			read(tmpstr2,*) ngldir
		elseif(trim(adjustl(tmpstr1)).eq.'-numdrop') then
			read(tmpstr2,*) numdrop
		elseif(trim(adjustl(tmpstr1)).eq.'-dropstep') then
			read(tmpstr2,*) dropstep
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') ' Usage: ./BSK_calc -om omegam -w w '//&
			 '-input intpufilename -output outputfilname -beta beta '//&
			 '-printinfo printinfo -numNNB numNNB -minimalrcut minimalrcut '//&
			 '-maximalrcut maximalrcut -nglcrosscheck do_nglcrosscheck '//&
			 '-ngldir ngldir -numdrop numdrop -dropstep dropstep'
			stop
		endif
	enddo
        if(omegam.eq.100.0) omegam = om_dft
        if(w.eq.100.0) w = om_dft
	print *, 'Default cosmology: omegam, w = ', om_dft, w_dft
	print *, 'Using   cosmology: omegam, w = ', omegam, w
	
	if(.not.has_outputfilename) then
		outputfilename = trim(adjustl(inputfilename))//'_output'
		if(minimalrcut .ne. -1.0e30) then
			write(tmpstr1,'(f30.3)') minimalrcut
			outputfilename = trim(adjustl(outputfilename))//'_rmin'//trim(adjustl(tmpstr1))
		endif
		if(maximalrcut .ne. 1.0e30) then
			write(tmpstr1,'(f30.3)') maximalrcut
			outputfilename = trim(adjustl(outputfilename))//'_rmax'//trim(adjustl(tmpstr1))
		endif
		write(tmpstr1, *) numNNB; outputfilename = trim(adjustl(outputfilename))//'_numNNB'//trim(adjustl(tmpstr1));
		write(tmpstr1, '(f9.4)') beta; outputfilename = trim(adjustl(outputfilename))//'_beta'//trim(adjustl(tmpstr1));
	endif

	print *, '#####################################'
	print *, ' Settings'
	write(*,'(A,2f8.3)') '   omegam, w = ', omegam, w
	write(*,'(A,A)') '   inputfilename: ', trim(adjustl(inputfilename))
	write(*,'(A,A)') '   outputfilename: ', trim(adjustl(outputfilename))
	print *, '  beta = ', beta
	print *, '  printinfo = ', printinfo
	print *, '  #-NNB searching for BSK: ', numNNB
	print *, '  minimalrcut: ', minimalrcut
	print *, '  maximalrcut: ', maximalrcut
	print *, '#####################################'

	call BSK(inputfilename, omegam, w, beta, outputfilename, printinfo, numNNB, minimalrcut, maximalrcut, &
		do_nglcrosscheck, ngldir, outputBSKindex=.true., outputBSKinfo=.true., outputxyz=.true., do_init=.true.)

	call BSK_stat(omegam, w, beta, outputfilename, numdrop, dropstep, printinfo)

	!call mpi_barrier(mpi_comm_world,ierr)
	!call mpi_finalize(ierr)
end program LSS_main
