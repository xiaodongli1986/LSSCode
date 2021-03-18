
module illustris_boxsplit_tools

use LSS_cosmo_funs       
implicit none
contains
    
  !---------------------------------------------------------------
  ! npar objects, each object having nfeature features (3 positions, 3 velocities, 1 mass, ...
  !---------------------------------------------------------------
    subroutine readin_shape(filename, npar, nfeature, boxsize)
        character(len=char_len), intent(in) :: filename
        integer, intent(out):: npar, nfeature
        double precision, intent(out) :: boxsize
        double precision :: npar_db, nfeature_db
        open(unit=258035365,file=trim(adjustl(filename)), form='unformatted', access='stream', action='read')
        read(258035365) npar_db, nfeature_db, boxsize
        close(258035365)
        npar = int(npar_db+0.1); nfeature = int(nfeature_db+0.1)
    end subroutine readin_shape
    subroutine readin_binary(filename, A, npar, nfeature)
        character(len=char_len), intent(in) :: filename
        integer, intent(in) ::  npar, nfeature
        real :: A(nfeature, npar)
        double precision :: npar_db, nfeature_db, boxsize
        open(unit=258035365,file=trim(adjustl(filename)), form='unformatted', access='stream', action='read')
        read(258035365) npar_db, nfeature_db, boxsize
        read(258035365) A
        close(258035365)
    end subroutine readin_binary

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
            ixs(1) = max(ixs(1),1); ixs(1) = min(ixs(1),nbox); ixs(2) = max(ixs(2),1); ixs(2) = min(ixs(2),nbox)
            iys(1) = max(iys(1),1); iys(1) = min(iys(1),nbox); iys(2) = max(iys(2),1); iys(2) = min(iys(2),nbox)
            izs(1) = max(izs(1),1); izs(1) = min(izs(1),nbox); izs(2) = max(izs(2),1); izs(2) = min(izs(2),nbox)
            !iys(1) = max(iys(1),1); iys(2) = min(iys(2),nbox)
            !izs(1) = max(izs(1),1); izs(2) = min(izs(2),nbox)
            nwrites = 1
            if(ixs(2) .ne. ixs(1)) nwrites(1) = 2
            if(iys(2) .ne. iys(1)) nwrites(2) = 2
            if(izs(2) .ne. izs(1)) nwrites(3) = 2
        endif
    end subroutine determine_iwrites

end module illustris_boxsplit_tools



program main_illustris_boxsplit

!use LSS_cosmo_funs
use illustris_boxsplit_tools

implicit none

    character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, &
        inputfile, outputdir, printstr, filelistfile, filename, output_infofile, output_logfile
    character(len=char_len), allocatable :: inputfiles(:), outputfiles(:,:,:)
    integer :: i, nfile, tmpi, ifile, nsplit=1, npar, nfeature, size=1, nowunit, ix,iy,iz, ixs(2),iys(2),&
      izs(2), nwrites(3), nbox, flag, ipar, i1,i2,i3, iwrite, nown, ntotal_read, ntotal_write, &
      output_infounit, output_logunit
    integer, allocatable :: seed(:), outputfiles_units(:,:,:)
    real(dl) :: tmpx, tmpy, overlaprat, overlap_distance, x,y,z, xyzmin
    real, allocatable :: input_data(:,:)
    double precision :: npar_db, nfeature_db, boxsize

    type data_arrays
       integer :: n
       real, allocatable :: features(:,:)
       integer :: ntotal
   end type data_arrays
   type(data_arrays), allocatable :: all_data(:,:,:)
   integer, parameter :: data_array_max_len = 100000
   logical :: testprint = .false.
    
    if(iargc().le.1) then
        write(*,'(A)') "#Example: "
        write(*,'(A)') 'snapnum=135'
        write(*,'(A)') 'overlaprat=0.05'
        write(*,'(A)') 'nsplit=5'
        write(*,'(A)') 'nowdir=Illustris-3/Snapshot/'
        write(*,*) 
        write(*,'(A)') "for parttype in 0 1 4 5 "
        write(*,'(A)') "do"
        write(*,'(A)') "    LSS_illustris_boxsplit    -inputfile ${nowdir}snap_${snapnum}.PartType${parttype}.\*    "// &
            '-outputdir ${nowdir}snap_$snapnum.PartType$parttype    -overlaprat $overlaprat   -nsplit $nsplit '
        write(*,'(A)') "done"
        write(*,'(A)') '# for illustris-1, using  '//&
          'LSS_illustris_boxsplit    -inputfile'// &
          '${nowdir}snap_${snapnum}.PartType${parttype}/snap_${snapnum}.PartType${parttype}.\*    -outputdir'// &
          '${nowdir}snap_$snapnum.PartType$parttype/snap_${snapnum}.PartType$parttype    -overlaprat $overlaprat  '// &
          ' -nsplit $nsplit '
        stop
    endif

    outputdir = ""
    do i = 1, iargc()
        if(mod(i,2).eq.0) cycle
        call getarg(i,tmpstr1)
        call getarg(i+1,tmpstr2)
        if(trim(adjustl(tmpstr1)).eq."-inputfile") then
            read(tmpstr2,"(A)") inputfile
        elseif(trim(adjustl(tmpstr1)).eq."-outputdir") then
            read(tmpstr2,"(A)") outputdir
        elseif(trim(adjustl(tmpstr1)).eq."-nsplit") then
            read(tmpstr2,*) nsplit
        elseif(trim(adjustl(tmpstr1)).eq."-overlaprat") then
            read(tmpstr2,*) overlaprat
        else
            write(*,*) "Unkown argument: ", trim(adjustl(tmpstr1))
            write(*,"(A)") trim(adjustl(printstr))
            stop
        endif
    enddo

    if(trim(adjustl(outputdir)).eq."") then
        write(*,*) 'ERROR! (LSS_illustris_boxsplit): must give a outputdir!'; stop
    endif
    call system("mkdir -p "//trim(adjustl(outputdir)))

    ntotal_read = 0
    ntotal_write = 0
    !----------------------------------
    ! 1. names of inputfiles 
    call random_seed(size=size); allocate(seed(size)); call random_seed(put=seed)
    call random_number(tmpx); call random_number(tmpy); tmpi = int(tmpx*100000000 + tmpy * 100000)
    write(filelistfile, *) tmpi; filelistfile= 'illustris_boxsplit_filelist_'//trim(adjustl(filelistfile))//'.tmp'
    call system("ls "//trim(adjustl(inputfile))//' > '//trim(adjustl(filelistfile)));

    call count_line_number(filelistfile, nfile)
    allocate(inputfiles(nfile),outputfiles(nsplit,nsplit,nsplit),outputfiles_units(nsplit, nsplit, nsplit),&
            all_data(nsplit,nsplit,nsplit))

    open(file=trim(adjustl(filelistfile)),action="read",unit=100); nfile = 0
    do while(.true.)
      nfile = nfile+1; read(100,'(A)',end=100) inputfiles(nfile); cycle
  100 exit
    enddo
    close(100); nfile = nfile-1


    !----------------------------------
    ! 2. names of outputfiles & allocate all_data
    write(tmpstr1, *) nsplit; write(tmpstr5, '(f10.3)') overlaprat;  nowunit = 100000;

    output_infofile = trim(adjustl(outputdir))//'/nsplit'//trim(adjustl(tmpstr1))// &
        '_overlaprat'//trim(adjustl(tmpstr5))//'.info'
    output_logfile = trim(adjustl(outputdir))//'/nsplit'//trim(adjustl(tmpstr1))// &
        '_overlaprat'//trim(adjustl(tmpstr5))//'.log'
    output_infounit = 10; open(unit=output_infounit, file=trim(adjustl(output_infofile)))
    output_logunit = 11; open(unit=output_logunit, file=trim(adjustl(output_logfile)))

    write(*,*) '(LSS_illustris_boxsplit) will generate ', nsplit*nsplit*nsplit, 'files'
    write(output_logunit,*) '(LSS_illustris_boxsplit) will generate ', nsplit*nsplit*nsplit, 'files'

    do ix=1,nsplit; do iy=1,nsplit; do iz=1,nsplit
      !----------------------------------
      ! 2.1 give it a name
      write(tmpstr2,*) ix; write(tmpstr3,*) iy; write(tmpstr4,*) iz; 
      outputfiles(ix,iy,iz) = trim(adjustl(outputdir))//'/nsplit'//trim(adjustl(tmpstr1))//'_overlaprat'//&
        trim(adjustl(tmpstr5))//'_subbox'//trim(adjustl(tmpstr2))//'_'//trim(adjustl(tmpstr3))//'_'//trim(adjustl(tmpstr4))
      write(*,'(5x, A)') trim(adjustl(outputfiles(ix,iy,iz)))
      write(output_logunit,'(5x, A)') trim(adjustl(outputfiles(ix,iy,iz)))
      !----------------------------------
      ! 2.2 give it a unit
      outputfiles_units(ix,iy,iz) = nowunit; nowunit=nowunit+1
      !----------------------------------
      ! 2.3 open this file
      open(unit=outputfiles_units(ix,iy,iz), file = trim(adjustl(outputfiles(ix,iy,iz))), &
         action='write', access='stream', form='unformatted')
    enddo; enddo; enddo

    

!    type data_arrays
!       integer :: n
!       integer, allocatable :: ids(:)
!       real, allocatable :: features(:,:)
!       integer :: ntotal
!   end type data_arrays
!   type(data_arrays), allocatable :: all_data(:,:,:)
!   integer, parameter :: data_array_max_len = 1000000
    !----------------------------------
    ! 3. write data to outputfiles
    write(*,*) '(LSS_illustris_boxsplit) In total ', nfile, 'files for read-in'
    write(output_logunit,*) '(LSS_illustris_boxsplit) In total ', nfile, 'files for read-in'
    write(output_infounit,'(A)') '# outputfilename  npar  nfeature  boxsize  ix iy iz'
    do ifile = 1, nfile
      filename = inputfiles(ifile); 
      write(*,'(A,i4,A,A)') '      ', ifile, '-th file:     ', trim(adjustl(filename))
      write(output_logunit,'(A,i4,A,A)') '      ', ifile, '-th file:     ', trim(adjustl(filename))
      if(testprint) &
              write(*,*) '      (test) will read in data with npar, nfeature =', npar, nfeature
      call readin_shape(filename, npar, nfeature, boxsize); allocate(input_data(nfeature, npar))
      !write(outout_infofile,'(A,i12,i4,f15.3)') trim(adjustl(filename)), npar, nfeature, boxsize
      if(ifile.eq.1) then
          do ix=1,nsplit; do iy=1,nsplit; do iz=1,nsplit
              !----------------------------------
              ! allocate all_data
              allocate(all_data(ix,iy,iz)%features(nfeature, data_array_max_len))
              all_data(ix,iy,iz)%n = 0; all_data(ix,iy,iz)%ntotal=0
          enddo; enddo; enddo
      endif

      overlap_distance = boxsize / dble(nsplit) * overlaprat
      if(testprint) &
            write(*,*) '      (test) overlap_distance = ', overlap_distance
      if(ifile.eq.1) then
              write(*,'(12x,A,f10.3,i4,f10.3)') 'boxsize, nsplit, overlap_distance = ', &
                 boxsize,nsplit,overlap_distance
              write(output_logunit,'(12x,A,f10.3,i4,f10.3)') 'boxsize, nsplit, overlap_distance = ', &
                 boxsize,nsplit,overlap_distance
      endif
      if(testprint) &
              write(*,'(A,A)') '      (test) begin read in data: ', trim(adjustl(filename))
      call readin_binary(filename, input_data, npar, nfeature); ntotal_read=ntotal_read+npar
      if(testprint) &
              write(*,*) '      (test) finish read in data: ', npar, nfeature
      do ipar = 1,npar     
        x=input_data(1,ipar); y=input_data(2,ipar); z=input_data(3,ipar);
        call determine_iwrites(x, y, z, overlap_distance, 0._dl, boxsize, nwrites, ixs, iys, izs, nsplit, flag)
        !if (flag.eq.1) then
        if(testprint .and. mod(ipar,100000).eq.0) print *, 'processing ', ipar, '-th particle...'
          do i1 = 1, nwrites(1); do i2 = 1, nwrites(2); do i3 = 1, nwrites(3)
             ix = ixs(i1); iy = iys(i2); iz = izs(i3)
             if(ix>nsplit.or.iy>nsplit.or.iz>nsplit) then
                     print *, 'ERROR (LSS_illustris_boxsplit)! nsplit overflows: ipar, ix,iy,iz, x,y,z = ', &
                         ipar, ix,iy,iz, x,y,z
             endif
             all_data(ix,iy,iz)%n = all_data(ix,iy,iz)%n + 1; nown = all_data(ix,iy,iz)%n
             all_data(ix,iy,iz)%ntotal = all_data(ix,iy,iz)%ntotal + 1
             all_data(ix,iy,iz)%features(:,nown) = input_data(:,ipar)
             if(all_data(ix,iy,iz)%n .eq. data_array_max_len) then
                 iwrite = outputfiles_units(ix,iy,iz)
                 write(iwrite) all_data(ix,iy,iz)%features; ntotal_write=ntotal_write+data_array_max_len
                 all_data(ix,iy,iz)%n = 0
             endif
          enddo; enddo; enddo
        !endif
      enddo
      deallocate(input_data)
    enddo
    do ix=1,nsplit; do iy=1,nsplit; do iz=1,nsplit
        iwrite = outputfiles_units(ix,iy,iz)
        if(all_data(ix,iy,iz)%n > 0) write(iwrite) all_data(ix,iy,iz)%features(:,1:all_data(ix,iy,iz)%n)
        ntotal_write = ntotal_write+all_data(ix,iy,iz)%n
        all_data(ix,iy,iz)%n = 0 ! set as zero, be ready for next time output
        write(output_infounit,'(A,i12,i4,f15.3,i4,i4,i4)') trim(adjustl(outputfiles(ix,iy,iz))), &
          all_data(ix,iy,iz)%ntotal, nfeature, boxsize, ix,iy,iz
        close(iwrite)
    enddo; enddo; enddo

    write(*,'(A,i12,i12)') ' (LSS_illustris_boxsplit) total objects read/write = ', ntotal_read, ntotal_write
    write(output_logunit,'(A,i12,i12)') ' (LSS_illustris_boxsplit) total objects read/write = ', ntotal_read, ntotal_write
    write(*,'(A,A)') ' (LSS_illustris_boxsplit) information of outputfiles stored in ', trim(adjustl(output_infofile))
    write(output_logunit,'(A,A)') ' (LSS_illustris_boxsplit) information of outputfiles stored in ', trim(adjustl(output_infofile))
    write(*,'(A,A)') ' (LSS_illustris_boxsplit) detailed information stored in ', trim(adjustl(output_logfile))
    write(output_logunit,'(A,A)') ' (LSS_illustris_boxsplit) detailed information stored in ', trim(adjustl(output_logfile))
    close(output_infounit); close(output_logunit)
    !----------------------------------
    call system("rm "//trim(adjustl(filelistfile)))

end program
