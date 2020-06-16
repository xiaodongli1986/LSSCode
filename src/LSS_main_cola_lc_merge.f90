module cola_lc_merge_tools
use LSS_cosmo_funs
implicit none
contains

	subroutine read_in_npar(basename, istep, num_threads, npar)
                character(len=char_len) :: basename, filename, filename_tmp, filenamenpar, tmpstr
                integer :: istep, num_threads,  ifile
                integer *8 :: npar, npar_tmp

                write(tmpstr, *) istep
                filename = trim(adjustl(basename))//trim(adjustl(tmpstr))

                if(num_threads.eq.0) then
                        filenamenpar = trim(adjustl(filename))//".npar"
                        open(file=trim(adjustl(filenamenpar)),unit=1000,action='read',form='binary',access='stream'); read(1000) npar; close(1000)
                        print *, '     get npar=',npar,' from ', trim(adjustl(filenamenpar)); 
                else
                        npar = 0
                        do ifile = 0, num_threads-1
                          write(tmpstr, *) ifile
                          filename_tmp = trim(adjustl(filename))//"_thread"//trim(adjustl(tmpstr))
                          filenamenpar = trim(adjustl(filename_tmp))//".npar"
                          open(file=trim(adjustl(filenamenpar)),unit=1000,action='read',form='binary',access='stream'); read(1000) npar_tmp; close(1000)
                          npar = npar + npar_tmp
                        enddo
                endif
	end subroutine read_in_npar

	subroutine read_in_data(basename, istep, num_threads, npar, rmid, ids, xyzvs)
	use LSS_cosmo_funs
	implicit none
                character(len=char_len) :: basename, filename, filename_tmp, filenamenpar, tmpstr
                integer :: istep, num_threads, ifile, i, istart
                integer *8 :: npar, npar_tmp
                integer*8 :: ids(npar)
                real :: xyzvs(6,npar), rmid

                write(tmpstr, *) istep
                filename = trim(adjustl(basename))//trim(adjustl(tmpstr))

                if(num_threads.eq.0) then
                        open(file=trim(adjustl(filename)),unit=1001, action='read', form='binary', access='stream')
                        read(1001) rmid
                        do i =1, npar
                                read(1001) ids(i); read(1001) xyzvs(1:6,i)
                        enddo
                        close(1001)
                else
                        istart = 1
                        do ifile = 0, num_threads-1
                                write(tmpstr, *) ifile
                                filename_tmp = trim(adjustl(filename))//"_thread"//trim(adjustl(tmpstr))
                                filenamenpar = trim(adjustl(filename_tmp))//".npar"
                                open(file=trim(adjustl(filenamenpar)),unit=1000, action='read', form='binary', access='stream')
                                read(1000) npar_tmp; close(1000)
                                write(*,'(A,i15,A,10x,A)') '     read-in ',npar_tmp,'       particles from ', trim(adjustl(filenamenpar)); 
                                open(file=trim(adjustl(filename_tmp)),unit=1001, action='read', form='binary', access='stream')
                                read(1001) rmid
                                do i =istart, istart+npar_tmp-1
                                        read(1001) ids(i); read(1001) xyzvs(1:6,i)
                                enddo
                                close(1001)
                                istart = istart + npar_tmp
                        enddo
                endif
                !write(*, '(A,f12.3)')       '   read-in rmid = ', rmid1
                !write(*, '(A,i15,i15)')     '   ids1 begin/end with ' , ids1(1), ids1(npar1)
                !write(*, '(A,6(f12.3))')    '   xyzvs1 begin with ', xyzvs1(:,1)
                !write(*, '(A,6(f12.3))')    '   xyzvs1  end  with ', xyzvs1(:,npar1)
	end subroutine read_in_data

	subroutine print_stat(xyzvs, npar)
		integer *8, intent(in) :: npar
		real, intent(in) :: xyzvs(6, npar)

		real :: xmin, xmax, xmean, ymin, ymax, ymean, zmin, zmax, zmean, rmin, rmax, rmean
		real :: x,y,z,r
		integer *8 :: i

		xmin=1e30;ymin=1e30;zmin=1e30;rmin = 1e30;
		xmax=-1e30;ymax=-1e30;zmax=-1e30;rmax = -1e30;
		xmean=0;ymean=0;zmean=0;rmean = 0;

		do i = 1, npar
		        x = xyzvs(1,i);  y = xyzvs(2,i);  z = xyzvs(3,i)
		        r = sqrt(x*x + y*y + z*z)
		        xmin = min(x,xmin); xmax = max(x,xmax); xmean = xmean + x
		        ymin = min(y,ymin); ymax = max(y,ymax); ymean = ymean + y
		        zmin = min(z,zmin); zmax = max(z,zmax); zmean = zmean + z
		        rmin = min(r,rmin); rmax = max(r,rmax); rmean = rmean + r
		enddo
		xmean = xmean / float(npar)
		ymean = ymean / float(npar)
		zmean = zmean / float(npar)
		rmean = rmean / float(npar)

		write(*, '("Finishing read in  ", i11, "particles")') npar
		write(*, '("min, max, mean of x  = ", 3f10.3)') xmin, xmax, xmean
		write(*, '("min, max, mean of y  = ", 3f10.3)') ymin, ymax, ymean
		write(*, '("min, max, mean of z  = ", 3f10.3)') zmin, zmax, zmean
		write(*, '("min, max, mean of r  = ", 3f10.3)') rmin, rmax, rmean
		        
	end subroutine print_stat

	  


	subroutine remove_some_ids(ids1, ids2, xyzvs1, xyzvs2, rmid1, rmid2, ids1_rlt, ids2_rlt, npar1, npar2, boxsize, test_flag)
		! <<dummy>>
		integer*8, intent(in) :: npar1, npar2
		integer*8, intent(in) :: ids1(npar1), ids2(npar2)
		real, intent(in) :: xyzvs1(6,npar1), xyzvs2(6,npar2), rmid1, rmid2, boxsize
		integer*2, intent(out) :: ids1_rlt(npar1), ids2_rlt(npar2)
		! <<local>>
		integer*8 :: i1, i1a,i1b, i2a,i2b, i2, id1, id2, n1_remove, n2_remove, ipar_print
		real ::  r1, r2, dr, dr1, dr2, dr1_mean, dr2_mean
		logical :: exit_loop = .false., test_flag, test_output = .false.

		ids1_rlt = 1
		ids2_rlt = 1
		n1_remove = 0; n2_remove = 0
		dr1_mean = 0; dr2_mean = 0
                write(*,'(A,\)') '    (remove_some_ids) start: 0% '

		if(test_output) &
		        open(unit=1234567,file='testdr.txt')
		i1 = 1; id1 = ids1(i1)
		i2 = 1; id2 = ids2(i2)
                ipar_print = npar1/20
		do while(.true.)
		        if(test_flag) print *, 'mark remove A'

		        ! keep increase id1 if id1<id2
		        do while(id1<id2)
		                        i1=i1+1; 
	!                                        print *, id1, id2 ! hahahahaaha
		                        if(i1>npar1) then
		                                exit_loop = .true.; exit
		                        endif
		                        id1=ids1(i1)
                                        if(mod(i1,ipar_print).eq.0) then
                                                write(*,'(A,f4.1,A,\)') '==>', i1 / (npar1+0.) * 100, '%'
                                        endif
		        enddo
		        if(exit_loop) exit;
	!                if(id1.eq.id2) goto 100;

		        if(test_flag) print *, 'mark remove B'
		        ! keep increase id2 if id2<id1
		        do while(id1>id2)
		                        i2=i2+1; 
	 !                                       print *, id1, id2 ! hahahahaaha
		                        if(i2>npar2) then
		                                exit_loop=.true.; exit
		                        endif
		                        id2=ids2(i2)
		        enddo
		        if(exit_loop) exit;
		        if(test_flag) print *, 'mark remove C'


		        if (id1.eq.id2)  then
		           ! find out region of id1, id2

		           i1a = i1; i2a = i2
		           i1b = i1a; i2b = i2a
		           if(test_flag) then
		                   print *, 'i1, i2 = ', i1, i2
		                   print *, 'mark remove D'
		           endif
		           do while(i1b+1<=npar1 .and. ids1(i1b+1) .eq. ids1(i1a)) 
		                i1b = i1b + 1
		           enddo
	!                   if(test_flag)  print *, 'mark remove E'
		           do while(i2b+1<=npar2 .and. ids2(i2b+1) .eq. ids2(i2a)) 
		                i2b = i2b + 1
	!                        if(mod(i2b,1000).eq.0) then
	!                                print *, 'i2b, npar2, ids2(i2a), ids2(i2b) = ', int(i2b), int(npar2), int(ids2(i2a)), int(ids2(i2b))
	!                        endif
		           enddo
		           if(test_flag)  then
		                   print *, 'mark remove F'
		                   print *, 'npar1, npar2 = ', int(npar1), int(npar2)
		                   print *, 'i1a, i1b, ids1(i1a), ids1(i1b) = ', int(i1a), int(i1b), int(ids1(i1a)), int(ids1(i1b))
		                   print *, 'i2a, i2b, ids2(i2a), ids2(i2b) = ', int(i2a), int(i2b), int(ids2(i2a)), int(ids2(i2b))
		           endif
		           
		           do i1 = i1a, i1b
		           do i2 = i2a, i2b
	!                                        print *, 'check: ', id1, id2 ! hahahahaaha
		                r1 = (xyzvs1(1,i1) **2 + xyzvs1(2,i1)**2 + xyzvs1(3,i1)**2) ** 0.5
		                dr1 = abs(r1 - rmid1)
		                r2 = (xyzvs2(1,i2) **2 + xyzvs2(2,i2)**2 + xyzvs2(3,i2)**2) ** 0.5
		                dr2 = abs(r2 - rmid2)
		                !print *, id1, id2, dr1, dr2
		                dr = ((xyzvs1(1,i1)-xyzvs2(1,i2))**2 + (xyzvs1(2,i1)-xyzvs2(2,i2))**2 + (xyzvs1(3,i1)-xyzvs2(3,i2))**2)**0.5
		                if (dr < boxsize/2.) then
		                  !write(*,'(A,3f10.3)') 'x,y,z of the first file:  ', xyzvs1(1:3,i1)
		                  !write(*,'(A,3f10.3, A, f10.3)') 'x,y,z of the second file: ', xyzvs2(1:3,i2), '      dr = ', dr
		                  if(dr1>dr2) then
		                        ids1_rlt(i1) = 0
		                        n1_remove = n1_remove + 1
		                        if(test_output) &
		                                write(1234567,'(A, 3f10.3,A, 6f10.3, A, 6f10.3 )') &
		                                'dr, dr1, dr2 = ', dr, dr1, dr2, ' remove ', xyzvs1(1:6,i1), ' keep ', xyzvs2(1:6,i2)
		                  else
		                        ids2_rlt(i2) = 0
		                        n2_remove = n2_remove + 1
		                        if(test_output) &
		                                write(1234567,'(A, 3f10.3,A, 6f10.3, A, 6f10.3 )') &
		                                'dr, dr1, dr2 = ', dr, dr1, dr2, ' remove ', xyzvs2(1:6,i2), ' keep ', xyzvs1(1:6,i1)
		                  endif
		                else
	!                          if(test_output) &
	!                                write(1234567,'(A, 3f10.3,A, 3f10.3, A, 3f10.3 )') &
	!                                        'dr, dr1, dr2 = ', dr, dr1, dr2, ' keep ', xyzvs1(1:3,i1), ' keep ', xyzvs2(1:3,i2)
		                endif
		           enddo
		           enddo
		           if(test_flag)  print *, 'mark remove G'
		           i1 = i1b+1; 
	!                                        print *, id1, id2 ! hahahahaaha
		           if(i1>npar1) exit
		           if(test_flag)  print *, 'mark remove H'
		           id1 = ids1(i1)
		           if(test_flag)  print *, 'mark remove I'
		        endif
		enddo
		
                print *
		write(*, '(A, 2i15,A,2f15.7)') ' (remove_some_ids) n1_remove, n2_remove = ', n1_remove, n2_remove, &
		        '; rats = ', n1_remove / real(npar1), n2_remove / real(npar2)
		if(test_output) &
		        stop; close(1234567)

	end subroutine remove_some_ids

	subroutine read_cola_lc2(filename, npar, ids, xyzvs, rmid)!, xyzvs, npar, rmid)
		use LSS_cosmo_funs 
		implicit none
		! << dummy >>
		character (len=char_len), intent(in) :: filename
		integer*8, intent(in) :: npar
		integer*8, intent(out):: ids(:)
		real, intent(out) :: xyzvs(:,:), rmid
		! << local >>
		integer :: i

		! << read-in main data>>
		open(file=trim(adjustl(filename)),unit=1001, action='read', form='binary', access='stream')
		read(1001) rmid
		ids = 1; xyzvs = 0
		!do i =1, 1
		!        read(1001) ids(i:i)
		!        read(1001) xyzvs(1:6,i)
		!enddo
		close(1001)
	end subroutine read_cola_lc2        

	subroutine read_cola_lc(filename, ids, xyzvs, npar, rmid)
		use LSS_cosmo_funs 
		implicit none
		! << dummy >>
		character (len=char_len), intent(in) :: filename
		integer*8, allocatable, intent(out) :: ids(:)
		real, allocatable, intent(out) :: xyzvs(:,:)
		real, intent(out) :: rmid
		integer*8, intent(out) :: npar
		! << local >>
		character(len=char_len) :: filename2
		integer :: i

		return
		! << read-in npar >>
		filename2 = trim(adjustl(filename))//'.npar'
		print *, trim(adjustl(filename2))

		open(file=trim(adjustl(filename2)),unit=1000, action='read', form='binary', access='stream')
		read(1000) npar
		close(1000)
		write(*, '(A,A,A,i15)')    '  (read in npar from file ', trim(adjustl(filename2)), ')   npar = ', npar

		! << read-in main data>>
		allocate(ids(npar), xyzvs(6,npar))
		open(file=trim(adjustl(filename)),unit=1001, action='read', form='binary', access='stream')
		read(1001) rmid
		do i =1, npar
		        read(1001) ids(i)
		        read(1001) xyzvs(1:6,i)
		enddo
		close(1001)
		print *, 'Finishing read in ', trim(adjustl(filename))
	 
		
	end subroutine read_cola_lc
end module cola_lc_merge_tools




program main_cola_lc_merge

use LSS_cosmo_funs
use cola_lc_merge_tools

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, basename, outputfile, printstr, filename, filenamenpar,filename1,filename2, rltfile1,rltfile2, outputname, tmpstr, gadgetfile, cmdstr1, cmdstr2
	integer :: i, istep, istep1, istep2 , status1, status2, access, num_threads = 0, npar_tmp
        integer*8 :: npar, npar1,npar2, ipar, ipar_final
        real :: boxsize = 100.0, rmid, rmid1,rmid2, dr,dr1,dr2
        integer*8, allocatable :: ids(:), ids1(:), ids2(:), ids_final(:)
        integer*4, allocatable :: ids_for_write(:)
        integer*4 :: now_id_for_write = 1  ! will write fake id (to foolish rockstar); start with this number
        real, allocatable :: xyzvs(:,:), xyzvs1(:,:),xyzvs2(:,:)
        integer*2, allocatable :: ids1_rlt(:), ids2_rlt(:), ids_rlt1(:), ids_rlt2(:), ids_final_rlt(:)
        logical :: has_rlt1, has_rlt2, write_an_ascii_copy=.false., test_flag =.false.
        integer*8 :: n_output
        real :: parmass, omegam,  nparticle, omegal, hubble, tmpx

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

        printstr = "Usage: LSS_cola_lc_merge   -basename BASE_FILE_NAME    -istep1 BEGIN_OF_ISTEP    -istep2 END_OF_ISTEP -omegam omegam -h h   -nparticle nparticle    -boxsize BOXSIZE   -num_threads num_threads    -write_an_ascii_copy ... # will output a series of gadget-format files"
	if(iargc().le.1) then
		write(*,'(A)') printstr
		stop
	endif



	outputfile = ""
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-basename") then
			read(tmpstr2,"(A)") basename
		elseif(trim(adjustl(tmpstr1)).eq."-istep1") then
			read(tmpstr2,*) istep1
		elseif(trim(adjustl(tmpstr1)).eq."-istep2") then
			read(tmpstr2,*) istep2
		elseif(trim(adjustl(tmpstr1)).eq."-boxsize") then
			read(tmpstr2,*) boxsize
		elseif(trim(adjustl(tmpstr1)).eq."-omegam") then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq."-h") then
			read(tmpstr2,*) hubble
		elseif(trim(adjustl(tmpstr1)).eq."-nparticle") then
			read(tmpstr2,*) nparticle
		elseif(trim(adjustl(tmpstr1)).eq."-num_threads") then
			read(tmpstr2,*) num_threads
		elseif(trim(adjustl(tmpstr1)).eq."-write_an_ascii_copy") then
			read(tmpstr2,*) write_an_ascii_copy
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

        omegal = 1-omegam
        parmass = particle_mass(dble(omegam), dble(hubble), dble(boxsize), dble(nparticle))

        write(*, '(A,A)') ' (cola_lc_merge) will merge many cola lightcone-slices into one lightcone.'
        write(*, '(A, i4, A, i4, A, i4)') ' (cola_lc_merge) will scan files index from ', istep1, ' to', istep2, '; num_threads = ', num_threads
        write(*, '(A, f20.3)') ' (cola_lc_merge) set boxsize as ', boxsize
        write(*, '(A, 4f10.5 )') ' (cola_lc_merge) omegam, omegal, h, nparticle^(1/3) = ', omegam, omegal, hubble, nparticle
        write(*, '(A, e15.7)') ' (cola_lc_merge) particle mass  = ', parmass

        write(*, '(A)') '====================================================='
        write(*, '(A)') ' (cola_lc_merge) find overlapped particles of the files'


        ! ================================================================================
        ! <<   1. loop check the overlap   >>
        ! << find out the overlapped ids, store them to files 
        ! ================================================================================

        write(*, '(A)') '====================================================='

        do istep = istep1, istep2-1

                print *, '*  processing istep = ', istep


                ! << some filenames >>
                write(tmpstr, *) istep
                filename1 = trim(adjustl(basename))//trim(adjustl(tmpstr))
                write(tmpstr, *) istep+1
                filename2 = trim(adjustl(basename))//trim(adjustl(tmpstr))


                ! << read-in file1>>
                call read_in_npar(basename, istep, num_threads, npar1)
                if(npar1.eq.0) then
                        print *, 'cycle because of npar1 =0: ', npar1; cycle
                endif
                allocate(ids1(npar1), xyzvs1(6,npar1))
                call read_in_data(basename, istep, num_threads, npar1, rmid1, ids1, xyzvs1)


                ! << read-in file2>>
                call read_in_npar(basename, istep+1, num_threads, npar2 )
                if(npar2.eq.0) then
                        print *, 'cycle because of npar2 =0: ', npar2; cycle
                endif
                allocate(ids2(npar2), xyzvs2(6,npar2))
                call read_in_data(basename, istep+1, num_threads, npar2, rmid2, ids2, xyzvs2)


                ! << remove identitical ids >>
                allocate(ids1_rlt(npar1), ids2_rlt(npar2))
                call remove_some_ids(ids1, ids2, xyzvs1, xyzvs2, rmid1, rmid2, ids1_rlt, ids2_rlt, npar1, npar2, real(boxsize), test_flag)


                deallocate(ids1,xyzvs1, ids2,xyzvs2)


                ! << output the ids_rlts >>
                rltfile1 = trim(adjustl(filename1)) // '.exclude_rlt_1'
                rltfile2 = trim(adjustl(filename2)) // '.exclude_rlt_2'
                open(file=trim(adjustl(rltfile1)),unit=10001, action='write', form = 'binary', access = 'stream')
                open(file=trim(adjustl(rltfile2)),unit=10002, action='write', form = 'binary', access = 'stream')

                write(10001) ids1_rlt
                write(10002) ids2_rlt
                close(10001); close(10002)

                deallocate(ids1_rlt, ids2_rlt)

        enddo



        ! ================================================================================
        ! <<   2. loop remove the overlap   >>
        ! << generating 'pure' samples without identitial ids >>
        ! ================================================================================

        write(*, '(A)') '====================================================='


        open(file=trim(adjustl(basename))//'.rmids', unit=876, action='write')
        open(unit=192038,file=trim(adjustl(basename))//'_rand0.01.txt')

        now_id_for_write = 1

        do istep = istep1, istep2-1

                write(*, '(A)') '====================================================='
                print *, '** processing istep = ', istep

                write(tmpstr, *) istep
                filename = trim(adjustl(basename))//trim(adjustl(tmpstr))

                ! << read-in file1>>
                call read_in_npar(basename, istep, num_threads, npar)
                if(npar.eq.0) then
                        print *, 'cycle because of npar =0: ', npar; cycle
                endif
                allocate(ids(npar), xyzvs(6,npar))
                call read_in_data(basename, istep, num_threads, npar, rmid, ids, xyzvs)


                ! << 1% particles output for test >>
                do i = 1, npar
                        call random_number(tmpx)
                        if(tmpx>0.99) then
                                write(192038, '(6f15.3)') xyzvs(1:6,i)
                        endif
                enddo

                write(*, '(A,f12.3)')       '   read-in rmid = ', rmid
                write(*, '(A,i15,i15)')     '   ids begin/end with ', ids(1), ids(npar)
                write(*, '(A,6(f12.3))')    '   xyzvs begin with ', xyzvs(:,1)
                write(*, '(A,6(f12.3))')    '   xyzvs  end  with ', xyzvs(:,npar)
                close(1001)


                ! << read-in ids1_rlt, ids2_rlt >> 

                allocate(ids_rlt1(npar), ids_rlt2(npar))

                rltfile1 = trim(adjustl(filename)) // '.exclude_rlt_1'
                status1 = access(trim(adjustl(rltfile1)), 'r')
                if(status1 .eq. 0) print *, '  file   exit  : ', trim(adjustl(rltfile1))
                if(status1 .eq. 1) print *, '  file not exit: ', trim(adjustl(rltfile1))
                if(status1 .eq. 0) then
                        open(file=trim(adjustl(rltfile1)),unit=2000, action='read', form='binary', access='stream')
                        read(2000) ids_rlt1
!                        print *, '   finishing read in ', trim(adjustl(rltfile1))
                        close(2000)
                        cmdstr1 = 'rm '//trim(adjustl(rltfile1))
                        call system(cmdstr1)
                endif

                rltfile2 = trim(adjustl(filename)) // '.exclude_rlt_2'
                status2 = access(trim(adjustl(rltfile2)), 'r')
                if(status2 .eq. 0) print *, '  file   exit  : ', trim(adjustl(rltfile2))
                if(status2 .eq. 1) print *, '  file not exit: ', trim(adjustl(rltfile2))
                if(status2 .eq. 0) then
                        open(file=trim(adjustl(rltfile2)),unit=2000, action='read', form='binary', access='stream')
                        read(2000) ids_rlt2
!                        print *, '   finishing read in ', trim(adjustl(rltfile2))
                        close(2000)
                        cmdstr2 = 'rm '//trim(adjustl(rltfile2))
                        call system(cmdstr2)
                endif


                ! << merge rltfile1 and rltfile2 >>
                allocate(ids_final_rlt(npar))
                if(status1 .eq. 1 .and. status2 .eq. 2) then
                        ids_final_rlt = 1
                else
                        if(status1 .eq.1 ) then
                                ids_final_rlt = ids_rlt2
                        elseif(status2 .eq. 1) then
                                ids_final_rlt = ids_rlt1
                        else
                                ids_final_rlt = ids_rlt1 * ids_rlt2
                        endif
                endif
                deallocate(ids_rlt1, ids_rlt2)


                n_output = 0
                do ipar = 1, npar
                        if(ids_final_rlt(ipar) .eq. 1) then
                                n_output = n_output + 1
                        endif
                enddo

                if(n_output.eq.0) then
                        deallocate( ids, xyzvs, ids_final_rlt )
                        cycle
                endif



                allocate(ids_final(n_output))
                ipar_final = 1
                do ipar = 1, npar
                        if(ids_final_rlt(ipar) .eq. 1) then
                                ids_final(ipar_final) = ipar
                                ipar_final = ipar_final + 1
                        endif
                enddo

                print *, ' before remove:'
                call print_stat(xyzvs, npar)
                print *, ' after remove:'
                call print_stat(xyzvs(:,ids_final(1:n_output)), n_output)


                ! compute particle mass for cola simulation

                headinfo.npart = 0; headinfo.npart(2) = n_output
                headinfo.mass = 0.; headinfo.mass(2) = parmass/1.e10
                headinfo.time = 1.; headinfo.redshift = 0; !headinfo.redshift = 0.1;
                headinfo.flag_sfr = 0; headinfo.flag_feedback=0;
                headinfo.npartTotal = 0; headinfo.npartTotal(2) = n_output
                headinfo.flag_cooling = 0; headinfo.num_files = 1
                headinfo.BoxSize = int(boxsize+0.5);
                headinfo.Omega0 = omegam;
                headinfo.OmegaLambda = omegal;
                headinfo.HubbleParam = hubble
                headinfo.fill = 0

                gadgetfile = trim(adjustl(filename))//'.filtered_gadget_file'
                open(file=trim(adjustl(gadgetfile)),unit=2001,form='unformatted',action='write')
                write(2001) headinfo
                write(2001) xyzvs(1:3,ids_final(1:n_output))
                write(2001) xyzvs(4:6,ids_final(1:n_output))

                ! fake ids
                allocate(ids_for_write(n_output))
                do i = 1, n_output
                        ids_for_write(i) = now_id_for_write + i
                enddo
                now_id_for_write = now_id_for_write + n_output + 1

                write(2001) ids_for_write
                close(2001)
                deallocate(ids_for_write)


                if(write_an_ascii_copy) then
                        print *
                        write(*,'(A)') '   output an ascii file as copy... file = '//trim(adjustl(gadgetfile))//'.ascii_copy'
                        open(file=trim(adjustl(gadgetfile))//'.ascii_copy',unit=2001,action='write')
                        do i = 1, n_output
                                write(2001,'(6f10.3, i10)') xyzvs(1:6,ids_final(i)), ids(ids_final(i))
!                        do i = 1, npar
!                                write(2001,'(6e15.7," 1")') xyzvs(1:6,i)
                        enddo
                        close(2001)
                endif


                deallocate( ids, xyzvs, ids_final_rlt, ids_final)
        enddo
        close(876)
        close(192038)

	if(trim(adjustl(outputfile)).eq."") then
		
	endif


end program
