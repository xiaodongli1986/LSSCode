module gadget_printhead_tools
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
end module gadget_printhead_tools



program main_gadget_printhead

use LSS_cosmo_funs
use gadget_printhead_tools

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputname, printstr
        integer :: i,j,k,  ifile, nfile,  block1, block2, ntotal, nnow, numarg, chbk, nparttotal, n1,n2, noutput
        real :: x,y,z,xyzmin, xyzmax, xyz_rescale,  xcut_min, xcut_max, ycut_min, ycut_max, zcut_min, zcut_max
        real, allocatable::  xyzvs(:,:)
        integer, allocatable :: is(:)
        logical :: do_xcut1, do_ycut1, do_zcut1, do_xcut2, do_ycut2, do_zcut2, &
                write_an_ascii_copy =.false., just_output_head = .false., print_all_head = .false.


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

	
	printstr = "Usage: LSS_gadget_to_fmt3 inputfile "
        numarg = iargc()
	if(numarg.le.0) then
		write(*,'(A)') printstr
		stop
	endif


        call getarg(1,tmpstr1)
        read(tmpstr1,'(A)') inputfile

        xyz_rescale = 1.



        write(*,'(A)') '##############################################'//&
                '####################################################################'
        write(*,'(A)') '(LSS_gadget_printhead) Print the head of a gadget file'
        write(*,'(A,A)') ' inputfile  = ', trim(adjustl(inputfile))


        write(*,'(A,A,A)'), ' open ', trim(adjustl(inputfile)), '...'
        !open(file=trim(adjustl(inputfile)),unit=10001,action='read',form='binary',access='sequential') ! gfortran fmt
        open(file=trim(adjustl(inputfile)),unit=10001,action='read',form='unformatted',access='stream')
        read(10001) block1, headinfo, block2; headinfo%boxsize = headinfo%boxsize*xyz_rescale
        chbk = check_blocks(block1,block2)

        write(*,'(A,6i12)')       ' npart(6) = ', headinfo%npart 
        write(*,'(A,6e15.7)')     ' mass(6)  = ', headinfo%mass
        write(*,'(A,e15.7)')      ' time     = ', headinfo%time
        write(*,'(A,e15.7)')      ' redshift = ', headinfo%redshift
        write(*,'(A,2i10)')       ' flag_sfr, flag_feedback  = ', headinfo%flag_sfr, headinfo%flag_feedback
        write(*,'(A,6i12)')       ' npartTotal(6) = ', headinfo%npartTotal
        write(*,'(A,2i12)')       ' flag_cooling, num_files = ', headinfo%flag_cooling, headinfo%num_files
        write(*,'(A,4e15.7)')     ' BoxSize, Omega0, OmegaLambda, HubbleParam =',  headinfo%BoxSize, &
                headinfo%Omega0, headinfo%OmegaLambda, headinfo%HubbleParam
        write(*,'(A,24i12)')     ' fill(24) = ',  headinfo%fill

        !write(*,'(A,3f14.5)')'   om /ol / h         = ', headinfo%Omega0, headinfo%OmegaLambda, headinfo%HubbleParam
        !write(*,'(A,3f14.5)')'   redshift / boxsize = ', headinfo%redshift, headinfo%boxsize
        nnow = 0;
        do i = 1, 6
                if(headinfo%npart(i).ne.0) then
                        write(*,'("   In total",i12,A,i3,A,e14.7,A,i12,A)') &
                                headinfo%npart(i), ' particles with type ', i, ', mass',&
                                headinfo%mass(i),' (',headinfo%nparttotal(i),') in total'
                        nonzeromass = headinfo%mass(i)
                        nnow = nnow+headinfo%npart(i)
                endif
        enddo
        ntotal = ntotal+nnow
        print *, '  output head file: ', trim(adjustl(inputfile)) //'.printhead'
        open(file=trim(adjustl(inputfile))//'.printhead',unit=999)
        write(999,'(A)') ' # ntotal+0., headinfo%boxsize, nonzeromass, headinfo%redshift, headinfo%Omega0, headinfo%HubbleParam'
        write(999,'(6e25.17)') ntotal+0., headinfo%boxsize, nonzeromass, headinfo%redshift, headinfo%Omega0, headinfo%HubbleParam
        close(999)
end program
