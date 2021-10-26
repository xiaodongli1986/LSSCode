program main_box_split_shift

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr, tmpstr1, tmpstr2, tmpstr3, inputfile, outputfile, printstr
	integer :: i, nsplit, i1, i2, i3, handle, xcol,ycol,zcol, ix,iy,iz, maxcol, iline
	real(dl) :: boxsize, subboxsize, x,y,z

	character(len=char_len), allocatable :: outputfile_array(:,:,:)
	integer, allocatable :: outputfile_handles(:,:,:)
	real(dl), allocatable :: values(:)
	logical :: add1 = .false.


	printstr = " (LSS_box_split_shift) Split a box into n^3 sub-boxes, and shift these sub-boxes back to the origin point (0,0,0). "//&
		" x,y,z should be the 1,2,3-th column"//&
	   "  Usage: EXE -inputfile intputfile -boxsize boxsize -nsplit nsplit "//&
	   "[-outputfile outputfile -maxcol maxcol -xcol xcol -ycol ycol -zcol zcol]"
	if(iargc().le.1) then
		print *, printstr
		stop
	endif

	outputfile = ""
	maxcol = -1
	xcol = 1
	ycol = 2
	zcol = 3
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
		elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
			read(tmpstr2,"(A)") outputfile
		elseif(trim(adjustl(tmpstr1)).eq."-boxsize") then
			read(tmpstr2,*) boxsize
		elseif(trim(adjustl(tmpstr1)).eq."-nsplit") then
			read(tmpstr2,*) nsplit
		elseif(trim(adjustl(tmpstr1)).eq."-maxcol") then
			read(tmpstr2,*) maxcol
		elseif(trim(adjustl(tmpstr1)).eq."-xcol") then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq."-ycol") then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq."-zcol") then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq."-add1") then
			read(tmpstr2,*) add1
		else

			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

	write(tmpstr1, *) nsplit
	if(trim(adjustl(outputfile)).eq."") then
		outputfile = trim(adjustl(inputfile))//'.nsplit'//trim(adjustl(tmpstr1))
	endif

	print *, 'Will read data from  ', trim(adjustl(inputfile))
	print *, 'Will output data to  ', trim(adjustl(outputfile))
	print *, 'Boxsize = ', boxsize
	print *, 'nsplit  = ', nsplit

	maxcol = max(maxcol, xcol)
	maxcol = max(maxcol, ycol)
	maxcol = max(maxcol, zcol)

	print *, 'xcol, ycol, zcol, maxcol  = ', xcol, ycol, zcol, maxcol

	call system('mkdir -p '//trim(adjustl(outputfile)))
	call system('sleep 2')

	allocate(outputfile_array(nsplit,nsplit,nsplit),outputfile_handles(nsplit,nsplit,nsplit))
	handle = 77583
	do i1 = 1, nsplit
	do i2 = 1, nsplit
	do i3 = 1, nsplit
		write(tmpstr1, *) i1
		write(tmpstr2, *) i2
		write(tmpstr3, *) i3
        outputfile_array(i1,i2,i3) = trim(adjustl(outputfile))//'/ifile_'//trim(adjustl(tmpstr1))//'_'&
			//trim(adjustl(tmpstr2))//'_'//trim(adjustl(tmpstr3))
		outputfile_handles(i1,i2,i3) = handle
		open(unit=handle,file=trim(adjustl(outputfile_array(i1,i2,i3))),action='write')
		handle = handle + 1
	enddo
	enddo
	enddo

	handle = handle + 1
	open(unit=handle, file=trim(adjustl(inputfile)),action='read')

	subboxsize = boxsize / real(nsplit)

	allocate(values(maxcol+1))

	iline = 0
	do while(.true.)
		read(handle, *,end=100)   values(1:maxcol)
		x = values(xcol); y = values(ycol); z = values(zcol)
		ix = int(x / subboxsize) 
		iy = int(y / subboxsize) 
		iz = int(z / subboxsize) 
		ix = min(ix, nsplit-1); iy = min(iy, nsplit-1); iz = min(iz, nsplit-1); 
		ix = max(ix, 0); iy = max(iy, 0); iz = max(iz, 0);
		
		values(xcol) = x - ix*subboxsize
		values(ycol) = y - iy*subboxsize
		values(zcol) = z - iz*subboxsize

		if(.not. add1) then
			write(outputfile_handles(ix+1,iy+1,iz+1), format_string(maxcol, '(e14.7,1x)')) values(1:maxcol)
		else
			values(maxcol+1) = 1
			write(outputfile_handles(ix+1,iy+1,iz+1), format_string(maxcol+1, '(e14.7,1x)')) values(1:maxcol+1)
		endif
		iline = iline + 1
		cycle 
100		exit
	enddo

	print *, '(LSS_box_split_shift) Finishing read and write ', iline, 'lines.'

	deallocate(values,outputfile_handles,outputfile_array)


end program
