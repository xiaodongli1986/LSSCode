module gadget_to_fmt3_tools
implicit none        
contains
        integer function check_blocks(block1,block2)
                integer :: block1, block2
                if(block1.eq.block2) then
                      print *, '  (check_blocks) consistent blocks:', block1, block2
                      check_blocks=1
                else
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      check_blocks=0
                endif
        end function check_blocks
end module gadget_to_fmt3_tools



program main_gadget_to_fmt3

use LSS_cosmo_funs
use gadget_to_fmt3_tools

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, output_suffix, printstr, inputfiles(10000), headfile, outputfile
        integer :: i,j,k,  ifile, nfile,  block1, block2, blockg,ntotal, nnow, numarg, chbk, nparttotal, n1,n2, ng,noutput
        real :: x,y,z,xyzmin, xyzmax, xyz_rescale,  xcut_min, xcut_max, ycut_min, &
                ycut_max, zcut_min, zcut_max, Omega0, OmegaLambda, HubbleParam, redshift, Boxsize
        real, allocatable::  xyzvs(:,:), xyzvsg(:,:)
        integer, allocatable :: is(:), ids(:)
        logical :: do_xcut1, do_ycut1, do_zcut1, do_xcut2, do_ycut2, do_zcut2, write_an_ascii_copy =.false., &
                write_a_gadget_copy =.false., just_output_head = .false., print_all_head = .false.
        double precision :: tmpdoubles(64), nonzeromass=0., randrat = 10., randx

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

	
	printstr = "Usage: LSS_gadget_to_fmt3 -input inputfile -output output_suffix -redshift redshift "//&
                "-Omega0 Omega0 -OmegaLambda OmegaLambda -HubbleParam HubbleParam  -Boxsize Boxsize  -xyz_rescale xyz_rescale"
        numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') printstr
		stop
	endif

        xyz_rescale = 1;  output_suffix=''
        do_xcut1 = .false.; do_ycut1=.false.; do_zcut1 = .false.
        do_xcut2 = .false.; do_ycut2=.false.; do_zcut2 = .false.
        Omega0 = -1.; redshift = -1.; OmegaLambda = -1.; HubbleParam = -1.; Boxsize = -1.;

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-output') then
			read(tmpstr2,*) output_suffix
		elseif(trim(adjustl(tmpstr1)).eq.'-xyz_rescale') then
			read(tmpstr2,*) xyz_rescale
		elseif(trim(adjustl(tmpstr1)).eq.'-Omega0') then
			read(tmpstr2,*) Omega0
		elseif(trim(adjustl(tmpstr1)).eq.'-redshift') then
			read(tmpstr2,*) redshift
		elseif(trim(adjustl(tmpstr1)).eq.'-OmegaLambda') then
			read(tmpstr2,*) OmegaLambda
		elseif(trim(adjustl(tmpstr1)).eq.'-HubbleParam') then
			read(tmpstr2,*) HubbleParam
		elseif(trim(adjustl(tmpstr1)).eq.'-Boxsize') then
			read(tmpstr2,*) Boxsize 
	!	elseif(trim(adjustl(tmpstr1)).eq.'-headfile') then
        !                read(tmpstr2,'(A)') headfile
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo


        if(trim(adjustl(output_suffix)).eq.'') then
                print *, 'ERROR! Find non-existed output_suffix : ', trim(adjustl(output_suffix))
                stop
        endif
                      

        write(*,'(A)') '#############################################################'//&
                '#####################################################'
        write(*,'(A)') '(LSS_gadget_change_head) Changing the head of the gadget file'
        write(*,'(A,A)') ' inputfile  = ', trim(adjustl(inputfile))
        write(*,'(A,A)') ' output_suffix = ', trim(adjustl(output_suffix))
        write(*,'(A,f10.3)') ' xyz_rescale = ', xyz_rescale
        write(*,'(A,f10.3)') ' HubbleParam = ', randrat
        write(*,'(A,f10.3)') ' Omega0 = ', Omega0
        write(*,'(A,f10.3)') ' OmegaLambda = ', OmegaLambda
        write(*,'(A,f10.3)') ' redshift = ', redshift
        write(*,'(A,f10.3)') ' Boxsize = ', Boxsize
        !write(*,'(A,L5)') ' write_an_ascii_copy = ', write_an_ascii_copy
        !write(*,'(A,L5)') ' write_a_gadget_copy = ', write_a_gadget_copy
        !write(*,'(A,L5)') ' just_output_head    = ', just_output_head
        !if (just_output_head) then
        !        write(*,'(A)') '  WARNING!! (LSS_gadget_to_fmt3) Will just output the head! fmt = 6*[double] (noutput, boxsize, mass, redshfit, Omega0, HubbleParam)'
        !endif

        if(Omega0.eq.-1. .and. OmegaLambda.eq.-1. .and. redshift.eq. -1. .and.HubbleParam.eq.-1.) then
                print *, 'ERROR! Found Omega0, OmegaLambda, redshift, Hubble Param all being -1!'
                print *, 'Will do nothing. Finish.'
                stop
        endif


        if(do_xcut1) print *, ' outputcut x >', xcut_min
        if(do_xcut2) print *, ' outputcut x <', xcut_max
        if(do_ycut1) print *, ' outputcut y >', ycut_min
        if(do_ycut2) print *, ' outputcut y <', ycut_max
        if(do_zcut1) print *, ' outputcut z >', zcut_min
        if(do_zcut2) print *, ' outputcut z <', zcut_max


        call system("ls "//trim(adjustl(inputfile))//' > gadget_to_fmt3_filelists.tmp ');

        open(file="gadget_to_fmt3_filelists.tmp",action="read",unit=100)
        nfile = 0
        do while(.true.)
                nfile = nfile+1; read(100,'(A)',end=100) inputfiles(nfile)
!                print *, nfile, trim(adjustl(inputfiles(nfile)))
                cycle
100             exit
        enddo
        close(100); nfile = nfile-1
        print *, ' found ', nfile, 'files to read-in...'
        do i = 1, nfile
                print *, i, trim(adjustl(inputfiles(i)))
        enddo

        ntotal=0
        n1 = 1
        ng = 1
        do ifile = 1, nfile
                print *
                write(*,'(A)') ' ##########################'
                write(*,'(A,A,A)'), ' open ', trim(adjustl(inputfiles(ifile))), '...'
                !open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='binary',access='sequential') ! gfortran fmt
                open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='unformatted',access='stream')
                read(10001) block1, headinfo, block2; headinfo%boxsize = headinfo%boxsize*xyz_rescale
                chbk = check_blocks(block1,block2)

                if(Omega0 .ne. -1.) then
                        if(ifile.eq.1) print *, '  ************** change Omega0 from ', headinfo%Omega0, ' to ', Omega0
                        headinfo%Omega0 = Omega0
                endif
                if(OmegaLambda.ne. -1.) then
                        if(ifile.eq.1) print *, '  ************** change OmegaLambda from ', &
                               headinfo%OmegaLambda, ' to ', OmegaLambda
                        headinfo%OmegaLambda = OmegaLambda
                endif
                if(redshift .ne. -1.) then
                        if(ifile.eq.1) print *, '  ************** change redshift from ', headinfo%redshift , ' to ', redshift
                        headinfo%redshift = redshift
                endif
                if(HubbleParam .ne. -1.) then
                        if(ifile.eq.1) print *, '  ************** change HubbleParam from ', &
                                headinfo%HubbleParam , ' to ', HubbleParam 
                        headinfo%HubbleParam = HubbleParam
                endif
                if(Boxsize.ne. -1.) then
                        if(ifile.eq.1) print *, '  ************** change HubbleParam from ', headinfo%Boxsize, ' to ', Boxsize 
                        headinfo%HubbleParam = Boxsize
                endif
                if(ifile.eq.1 ) then
                        write(*,'(A,6i12)')       ' npart(6) = ', headinfo%npart 
                        write(*,'(A,6e15.7)')     ' mass(6)  = ', headinfo%mass
                        write(*,'(A,e15.7)')      ' time     = ', headinfo%time
                        write(*,'(A,e15.7)')      ' redshift = ', headinfo%redshift
                        write(*,'(A,2i10)')       ' flag_sfr, flag_feedback  = ', headinfo%flag_sfr, headinfo%flag_feedback
                        write(*,'(A,6i12)')       ' npartTotal(6) = ', headinfo%npartTotal
                        write(*,'(A,2i12)')       ' flag_cooling, num_files = ', headinfo%flag_cooling, headinfo%num_files
                        write(*,'(A,4e15.7)')     ' BoxSize, Omega0, OmegaLambda, HubbleParam =',  &
                                headinfo%BoxSize, headinfo%Omega0, headinfo%OmegaLambda, headinfo%HubbleParam
                endif

                if(ifile.eq.1) then
                        write(*,'(A,3f14.5)')'   om /ol / h         = ', headinfo%Omega0, headinfo%OmegaLambda, headinfo%HubbleParam
                        write(*,'(A,3f14.5)')'   redshift / boxsize = ', headinfo%redshift, headinfo%boxsize
                        if(xyzmax.eq.0.) then
                                xyzmax = headinfo%boxsize
                                print *, '  xyz range (reset)  = ', xyzmin, xyzmax
                        endif
                endif
                nnow = 0;
                do i = 1, 6
                        if(headinfo%npart(i).ne.0) then
                                write(*,'("   In total",i12,A,i3,A,e14.7,A,i12,A)') headinfo%npart(i), &
                                        ' particles with type ', i, ', mass',&
                                        headinfo%mass(i),' (',headinfo%nparttotal(i),') in total'
                                nonzeromass = headinfo%mass(i)
                                nnow = nnow+headinfo%npart(i)
                        endif
                enddo
                ntotal = ntotal+nnow
                !if(ifile .eq. 1) then
                !       print *, '  output head file: ', trim(adjustl(output_suffix)) //'.head'
                !       open(file=trim(adjustl(output_suffix))//'.head',unit=999)
                !       write(*,'("  ",6e25.17)') ntotal+0., headinfo%boxsize, nonzeromass, headinfo%redshift, headinfo%Omega0, headinfo%HubbleParam
                !       write(999,'(6e25.17)') ntotal+0., headinfo%boxsize, nonzeromass, headinfo%redshift, headinfo%Omega0, headinfo%HubbleParam
                !       close(999)
                !       if(just_output_head) then
                !              print *, ' (LSS_gadget_to_fmt3) Stop as just_output_head set as True'; stop
                !       endif       
                !endif

                n2 = n1+nnow-1

                nparttotal = sum(headinfo%nparttotal)
                if(ifile.eq.1) allocate(xyzvs(7,nnow), ids(nnow))

                read(10001) block1
                read(10001) xyzvs(1:3,:); xyzvs(1:3,:)=xyzvs(1:3,:)*xyz_rescale
                read(10001) block2
                chbk = check_blocks(block1,block2)
                write(*,'(A,3f10.3,5x,3f10.3)') '     begin and end of xyzs: ', xyzvs(1:3, 1), xyzvs(1:3,nnow)

                read(10001) block1
                read(10001) xyzvs(4:6,:)
                read(10001) block2
                chbk = check_blocks(block1,block2)
                write(*,'(A,3f10.3,5x,3f10.3)') '     begin and end of   vs: ', xyzvs(4:6, 1), xyzvs(4:6,nnow)

                read(10001) block1
                read(10001) ids(:)
                read(10001) block2
                chbk = check_blocks(block1,block2)
                write(*,'(A,i15,5x,i15)') '     begin and end of   ids: ', ids(1), ids(nnow)

                close(10001)

                outputfile = trim(adjustl(inputfiles(ifile)))//trim(adjustl(output_suffix)) !//'.gadget_copy'
                write(*,'(A,A,A)'), ' write to ', trim(adjustl(outputfile)), '...'
                open(file=trim(adjustl(outputfile)),unit=2001,form='unformatted',action='write')
                write(2001)  headinfo
                write(2001)  xyzvs(1:3,:)
                write(2001)  xyzvs(4:6,:)
                write(2001)  ids
                deallocate(xyzvs, ids)
                close(2001)
        enddo

        write(*,*) 'Done.'


end program
