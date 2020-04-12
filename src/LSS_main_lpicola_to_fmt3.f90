module lpicola_to_fmt3_tools
implicit none        
contains
        integer function check_blocks(block1,block2)
                integer :: block1, block2
                if(block1.eq.block2) then
                      !print *, '  (check_blocks) consistent blocks:', block1, block2
                      check_blocks=1
                else
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      print *, '  (check_blocks) WARNING! Inconsistent blocks: ', block1, block2
                      check_blocks=0
                endif
        end function check_blocks
end module lpicola_to_fmt3_tools



program main_lpicola_to_fmt3

use LSS_cosmo_funs
use lpicola_to_fmt3_tools

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputname, printstr, inputfiles(10000), headfile
        integer :: i,j,k,  ifile, nfile,  ntotal, nnow, numarg, chbk, nparttotal, n1,n2, noutput
        integer*4 :: block1, block2, EOF
        real :: x,y,z,xyzmin, xyzmax, xyz_rescale,  xcut_min, xcut_max, ycut_min, ycut_max, zcut_min, zcut_max, w, redshift
        real, allocatable::  xyzvs(:,:), tmp_xyzvs(:,:)
        integer, allocatable :: is(:)
        logical :: do_xcut1, do_ycut1, do_zcut1, do_xcut2, do_ycut2, do_zcut2, write_an_ascii_copy =.false., just_output_head = .false.
        double precision :: tmpdoubles(64), nonzeromass=0., randrat = 10., randx, boxsize, omegam, hubble, parmass


	
	printstr = "Usage: LSS_lpicola_to_fmt3 -input inputfile -output outputname -xyz_rescale xyz_rescale -omegam omegam -w w -hubble hubble -boxsize boxsize -parmass parmass   -randrat the_rate_you_want_to_output -xcut_min xcut_min -xcut_max xcut_max -ycut_min ycut_min -ycut_max ycut_max -zcut_min zcut_min -zcut_max zcut_max -write_an_ascii_copy T/F -just_output_head F -headfile T #### Example: LSS_gadget_to_fmt3 -input snp04000e.\? -xyz_rescale 0.001 -xyzcut_min 0 -xyzcut_max 100 -write_an_ascii_copy F"
        numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') printstr
		stop
	endif

        xyz_rescale = 1;  outputname=''
        do_xcut1 = .false.; do_ycut1=.false.; do_zcut1 = .false.
        do_xcut2 = .false.; do_ycut2=.false.; do_zcut2 = .false.

        omegam = -1.; hubble = -1.; boxsize=-1.; parmass=1.;

        headfile = ''

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-output') then
			read(tmpstr2,*) outputname
		elseif(trim(adjustl(tmpstr1)).eq.'-xyz_rescale') then
			read(tmpstr2,*) xyz_rescale
		elseif(trim(adjustl(tmpstr1)).eq.'-xcut_min') then
                        read(tmpstr2,*) xcut_min; do_xcut1 = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-xcut_max') then
                        read(tmpstr2,*) xcut_max; do_xcut2 = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-ycut_min') then
                        read(tmpstr2,*) ycut_min; do_ycut1 = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-ycut_max') then
                        read(tmpstr2,*) ycut_max; do_ycut2 = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-zcut_min') then
                        read(tmpstr2,*) zcut_min; do_zcut1 = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-zcut_max') then
                        read(tmpstr2,*) zcut_max; do_zcut2 = .true.
		elseif(trim(adjustl(tmpstr1)).eq.'-randrat') then
                        read(tmpstr2,*) randrat
		elseif(trim(adjustl(tmpstr1)).eq.'-write_an_ascii_copy') then
                        read(tmpstr2,*) write_an_ascii_copy
		elseif(trim(adjustl(tmpstr1)).eq.'-just_output_head') then
                        read(tmpstr2,*) just_output_head
		elseif(trim(adjustl(tmpstr1)).eq.'-omegam') then
                        read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
                        read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-boxsize') then
                        read(tmpstr2,*) boxsize
		elseif(trim(adjustl(tmpstr1)).eq.'-hubble') then
                        read(tmpstr2,*) hubble
		elseif(trim(adjustl(tmpstr1)).eq.'-parmass') then
                        read(tmpstr2,*) parmass
		elseif(trim(adjustl(tmpstr1)).eq.'-headfile' .or. trim(adjustl(tmpstr1)).eq.'-lpicola_headfile')  then
                        read(tmpstr2,'(A)') headfile
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

        write(*,'(A)') '##################################################################################################################'
        write(*,'(A)') '(LSS_gadget_to_fmt3) Converting gadget files to block+head+block+blcok+[x,y,z,vx,vy,vz,1.0]*npar+blcok format !'
        if(trim(adjustl(headfile)).ne.'') then
                print *
                write(*,'(A)') ' Read in parameters from lpicola headfile:  (This will ignore hand-given values of omegam, w, boxsize, parmass, hubble!)'
                write(*,'(A,A)') '        ',trim(adjustl(headfile))
                print *
                open(unit=10002,file=headfile,action='read')
                read(10002,*) tmpstr1
                read(10002,*) ntotal, boxsize, parmass, redshift, omegam, hubble, w; 
                if(w.ge.-0.1) then
                        print *, '   set w as -1 (found w = ',w,')'
                        w = -1.
                endif
                close(10002)
        endif

        write(*,*) ' boxsize   = ', boxsize
        write(*,*) ' parmass   = ', parmass
        write(*,*) ' omegam    = ', omegam 
        write(*,*) ' w         = ', w
        write(*,*) ' hubble    = ', hubble 
        

        if(omegam<0. .or. boxsize<0. .or. hubble<0.) then
                print *, '(lpicola_to_fmt3) ERROR! need value of omegam, boxsize, hubble!:', omegam, boxsize, hubble
                write(*,'(A)') trim(adjustl(printstr))
                stop
        endif



        if(trim(adjustl(outputname)).eq.'') then
                print *, 'ERROR! Find non-existed outputname : ', trim(adjustl(outputname))
                stop
        endif
                      

        write(*,'(A,A)') ' inputfile  = ', trim(adjustl(inputfile))
        write(*,'(A,A)') ' outputname = ', trim(adjustl(outputname))
        write(*,'(A,f10.3)') ' xyz_rescale= ', xyz_rescale
        write(*,'(A,f10.3)') ' randrat    = ', randrat
        write(*,'(A,L5)') ' write_an_ascii_copy = ', write_an_ascii_copy
        write(*,'(A,L5)') ' just_output_head    = ', just_output_head
        if (just_output_head) then
                write(*,'(A)') '  WARNING!! (LSS_gadget_to_fmt3) Will just output the head! fmt = 6*[double] (noutput, boxsize, mass, redshfit, Omega0, HubbleParam)'
        endif


        if(do_xcut1) print *, ' outputcut x >', xcut_min
        if(do_xcut2) print *, ' outputcut x <', xcut_max
        if(do_ycut1) print *, ' outputcut y >', ycut_min
        if(do_ycut2) print *, ' outputcut y <', ycut_max
        if(do_zcut1) print *, ' outputcut z >', zcut_min
        if(do_zcut2) print *, ' outputcut z <', zcut_max


        call system("ls "//trim(adjustl(inputfile))//' > lpicola_to_fmt3_filelists.tmp ');

        open(file="lpicola_to_fmt3_filelists.tmp",action="read",unit=100)
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

        ! count number of particles...
        ntotal=0
        do ifile = 1, nfile
                write(*,'(A,A,A)'), ' open ', trim(adjustl(inputfiles(ifile))), ' for count ...'
                open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='unformatted')
                !open(unit=10002,file=trim(adjustl(inputfiles(ifile)))//'_reform')
                !read(10001) block1, headinfo, block2; headinfo.boxsize = headinfo.boxsize*xyz_rescale
                do
                        read(10001,IOSTAT=EOF) block1
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program' 
                        else if (EOF .lt. 0) then
                               ! If EOF < 0 then we have read all the chunks
                                exit
                        else
                                allocate(tmp_xyzvs(6,block1)); read(10001) tmp_xyzvs; 
                                !do i = 1, block1
                                !        write(10002,'(6e15.7)') tmp_xyzvs(:,i)
                                !enddo
                                deallocate(tmp_xyzvs)
                                ntotal  = ntotal + block1
                        endif
                enddo
                close(10001)
                !close(10002)
        enddo
        print * 
        print *, '#########################'
        write(*,'(A,i15,A,i8,A)') ' Finishing count-line. In total we have ', ntotal, ' particles in ', nfile, ' files.'

        open(file=trim(adjustl(outputname))//'.head',unit=999)
        write(*,'("  ",7e25.17)') ntotal+0., boxsize, parmass, -1., omegam, hubble, w
        write(999,'(7e25.17)') ntotal+0., boxsize, parmass, -1., omegam, hubble, w
        close(999)
        if(just_output_head) then
                print *, ' (LSS_lpicola_to_fmt3) Stop as just_output_head set as True'; stop
        endif

        ! read in
        allocate(xyzvs(7,ntotal))
        n1 = 1
        do ifile = 1, nfile
                write(*,'(A,A,A)'), ' open ', trim(adjustl(inputfiles(ifile))), ' for read-in ...'
                open(file=trim(adjustl(inputfiles(ifile))),unit=10001,action='read',form='unformatted')
                do
                        read(10001,IOSTAT=EOF) block1
                        
                        if(EOF.gt.0) then
                                print*, 'Read error in file: ', trim(adjustl(inputfiles(ifile)))
                                print*, 'Exiting program' 
                        else if (EOF .lt. 0) then
                               ! If EOF < 0 then we have read all the chunks
                                exit
                        else
                                read(10001) xyzvs(1:6,n1:n1+block1-1)
                                n1 = n1 + block1
                        endif
                enddo
                close(10001)
        enddo

        print * 
        print *, '#########################'
        write(*,'(A,i15,A,i8,A)') ' Finishing read-in.  In total we read-in', ntotal, ' particles in ', nfile, ' files.'


        print *, '#########################'
        write(*,*) 'output...'

        !if(do_xcut1.or.do_xcut2.or.do_ycut1.or.do_ycut2.or.do_zcut1.or.do_zcut2) then
        call random_seed()
        if(.true.) then
                write(*,*) '  selecting particles in xyz range for writing...'
                allocate(is(ntotal))
                noutput = 0
                do i =1, ntotal
                        call random_number(randx)
                        if(randx>randrat) cycle
                        x=xyzvs(1,i); y = xyzvs(2,i); z = xyzvs(3,i)
                        if(do_xcut1 .and. x<xcut_min) cycle
                        if(do_xcut2 .and. x>xcut_max) cycle
                        if(do_ycut1 .and. y<ycut_min) cycle
                        if(do_ycut2 .and. y>ycut_max) cycle
                        if(do_zcut1 .and. z<zcut_min) cycle
                        if(do_zcut2 .and. z>zcut_max) cycle
                        xyzvs(7,i) = 1
                        noutput=noutput+1; is(noutput) =i 
                enddo
        endif

        write(*,*) '  start writing head. format =  block+64*double_precision+block...'
        write(*,*) '    64*double = noutput/boxsize/part-mass/redshift/omegam/h/...nuisance'
        open(file=trim(adjustl(outputname)),unit=2001,form='unformatted',action='write')
        tmpdoubles = 0.
        tmpdoubles(1) = noutput
        tmpdoubles(2) = boxsize
        tmpdoubles(3) = parmass
        tmpdoubles(4) = -1
        tmpdoubles(5) = omegam
        tmpdoubles(6) = hubble
        tmpdoubles(7) = w

        write(2001) tmpdoubles

        write(*,*)
        write(*,'(i15,A,i15,A,f12.6)') noutput, ' particles selected from ', ntotal, '       rat = ', real(noutput)/real(ntotal)
        write(*,'(A)') '   staring write particles. format = block+[x,y,z,vx,vy,vz,1.0]*noutput+block ...'
        write(2001) xyzvs(1:7,is(1:noutput))
        close(2001)

        if(write_an_ascii_copy) then
                print *
                write(*,'(A)') '   output an ascii file as copy... file = '//trim(adjustl(outputname))//'.ascii_copy'
                open(file=trim(adjustl(outputname))//'.ascii_copy',unit=2001,action='write')
                do i = 1, noutput
                        write(2001,'(6e15.7," 1")') xyzvs(1:6,is(i))
                enddo
                close(2001)
        endif


        write(*,*) 'Done.'


end program
