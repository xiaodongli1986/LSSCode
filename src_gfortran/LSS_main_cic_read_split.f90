

program main_cic

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputname, printstr, inputfiles(10000), readin_fmt_str
        integer :: i,j,k, i1,i2,j1,j2,k1,k2, ifile, nfile, nc, &
		tmpint(64), block1, block2, ntotal, nnow, numarg, chbk, nsplit, n, readin_fmt
        real :: xyzmin, xyzmax, tmpfloats(64), xyz_rescale, xyz(3), vxvyvz(3)  !, xshift=0., yshift=0., zshift=0.
        real, allocatable:: rhogrid(:,:,:), vxgrid(:,:,:), vygrid(:,:,:), vzgrid(:,:,:) , xyzs(:,:), vs(:,:), rgrid(:,:,:)
        integer, parameter :: fmt_gadget = 1, fmt_ascii = 2
        logical :: do_rgrid = .false.

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

	!headinfo%BoxSize=1.0

	
	printstr = "Usage: LSS_cic_read_split -input inputfile   -nc n-cells       "//&
		"-nsplit nsplit    "//&
		"#### Example: LSS_cic_read_split -input subgrids.rhogrid -nc 900 -nsplit 20 "
        numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') printstr
		stop
	endif


	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfile
		elseif(trim(adjustl(tmpstr1)).eq.'-nc') then
			read(tmpstr2,*) nc
		elseif(trim(adjustl(tmpstr1)).eq.'-output') then
			read(tmpstr2,'(A)') outputname
		elseif(trim(adjustl(tmpstr1)).eq.'-nsplit') then
			read(tmpstr2,*) nsplit
        endif
!		elseif(trim(adjustl(tmpstr1)).eq.'-xshift') then
!			read(tmpstr2,*) xshift
!		elseif(trim(adjustl(tmpstr1)).eq.'-yshift') then
!			read(tmpstr2,*) yshift
!		elseif(trim(adjustl(tmpstr1)).eq.'-zshift') then
!			read(tmpstr2,*) zshift
	enddo

    if(nsplit.eq.1) then
            print *, 'nothig to do! (found nsplit=1)'
            stop
    endif

    write(*,'(A,A)')  '  inputfile  = ', trim(adjustl(inputfile))
    print *, ' nc (celss) = ', nc
    write(*,'(A,i3,A,i8,A)'), '  nsplit^3   = ', nsplit, '^3 = ', nsplit**3, ' sets of outupt files.'


    write(*,*) 'allocating fields...'
    allocate(rhogrid(nc,nc,nc))

    open(file=trim(adjustl(inputfile)),unit=2000,form='unformatted',action='read')
    read(2000) rhogrid

    write(tmpstr1, *) nsplit
    outputname = trim(adjustl(inputfile))//'_nsplit'//trim(adjustl(tmpstr1))
    tmpstr1 = 'mkdir -p '//trim(adjustl(outputname))
    write(*,'(A,A)') 'Executing cmd: ', trim(adjustl(tmpstr1))
    call system(trim(adjustl(tmpstr1)))


    print *, '#########################'
    write(*,*) 'output...'
    n = nc / nsplit
    print *, ' nc, nsplit, n = ', nc, nsplit, n
    do i=1,nsplit
    do j=1,nsplit
    do k=1,nsplit
         ifile = k+(j-1)*nsplit+(i-1)*nsplit*nsplit
         i1 = (i-1)*n+1; i2=i1+n-1
         j1 = (j-1)*n+1; j2=j1+n-1
         k1 = (k-1)*n+1; k2=k1+n-1
         write(tmpstr1,*) ifile
         write(*,'(3x,A,A,A,2i4,3x,2i4,3x,2i4)') 'write to ', trim(adjustl(outputname)), '. i/j/k range = ', i1,i2, j1,j2, k1,k2
         open(file=trim(adjustl(outputname))//'/ifile'//trim(adjustl(tmpstr1)),&
	         unit=2001,form='unformatted',action='write')
         write(2001) rhogrid(k1:k2,j1:j2,i1:i2)
         close(2001); 
    enddo
    enddo
    enddo

            


end program
