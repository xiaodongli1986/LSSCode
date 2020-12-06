program main_boss2pcf_fmtxyzwfile

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, datatype
	logical :: ignorecp, ignoremask
	integer :: i, maxcol, dt, nlines, nlinesw
	real(dl) :: tmp(100), ra,dec,r, x,y,z, weight, wfkp,wcp, AW,VW
	
	printstr = "EXE -datatype datatype -inputfile inputfile [-outputfile outputfile] [-ignorecp ignorecp] [-ignoremask ignoremask] ### if ignorecp, will not take CP into consideration; only applicable for mock or mockran; by default ignorecp=False"
	if(iargc().le.1) then
		print *, printstr
		stop
	endif

	ignorecp = .false.
        ignoremask = .false.
	outputfile = ""

	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
		elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
			read(tmpstr2,"(A)") outputfile
		elseif(trim(adjustl(tmpstr1)).eq."-datatype") then
			read(tmpstr2,"(A)") datatype
		elseif(trim(adjustl(tmpstr1)).eq."-ignorecp") then
			read(tmpstr2,*) ignorecp
		elseif(trim(adjustl(tmpstr1)).eq."-ignoremask") then
			read(tmpstr2,*) ignoremask
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

	if(trim(adjustl(outputfile)).eq."") then
		if(ignorecp.or.ignoremask) then
			tmpstr1 = '.xyzw'
			if(ignorecp) tmpstr1 = '.cpignored'//trim(adjustl(tmpstr1))
			if(ignoremask) tmpstr1 = '.maskignored'//trim(adjustl(tmpstr1))
			outputfile = trim(adjustl(inputfile))//trim(adjustl(tmpstr1))
		else
			outputfile = trim(adjustl(inputfile))//'.xyzw'
		endif
	endif

	print *, 'Creating formatted x,y,z,w files'
	print *, '  input: ', trim(adjustl(inputfile))
	print *, '  output: ', trim(adjustl(outputfile))
	print *, '  datatype: ', trim(adjustl(datatype))
	print *, '  ignorecp: ', ignorecp
	print *, '  ignoremask: ', ignoremask

        ! intput cosmology
        gb_omegam = 0.26d0; gb_w = -1.0; gb_h = 0.73_dl;
	print *
	print *, '   Setting cosmology as om, w = ', gb_omegam, gb_w
        call cosmo_funs_init(.true.)
        call de_calc_comovr()

	if(trim(adjustl(datatype)).eq.'data') then
		dt = 1 ! 1 means data; 2,3,4 means dataran, mock, mockran
		maxcol = 12
	elseif(trim(adjustl(datatype)).eq.'dataran') then
		dt = 2 ! 
		maxcol = 4
	elseif(trim(adjustl(datatype)).eq.'mock') then
		dt = 3 !
		maxcol = 14
	elseif(trim(adjustl(datatype)).eq.'mockran') then
		dt = 4 !
		maxcol = 14
	elseif(trim(adjustl(datatype)).eq.'compact') then
		dt = 5 !
		maxcol = 3
	else
		write(*,'(A,A)') ' ERROR (LSS_boss2pcf_fmtxyzwfile): wrong datatype: must be data, dataran, mock or mockran!', trim(adjustl(datatype))
		stop
	endif

	if(ignorecp .and. (dt.eq.1 .or. dt.eq.2 .or. dt.eq.5)) then
		write(*,'(A,A)') ' ERROR (LSS_boss2pcf_fmtxyzwfile): ignorecp only possible when datatype is mock, mockran!: ', trim(adjustl(datatype))
		stop
	endif
	if(ignoremask .and. (dt.eq.1 .or. dt.eq.2 .or. dt.eq.5)) then
		write(*,'(A,A)') ' ERROR (LSS_boss2pcf_fmtxyzwfile): ignoremask only possible when datatype is mock, mockran!: ', trim(adjustl(datatype))
		stop
	endif

	
	! Some information copied from boss2pcf.py
	open(unit=1,file=inputfile,action='read')
	open(unit=2,file=outputfile)
        nlines = 0; nlinesw = 0
        do while(.true.)
                read(1,*,end=101) tmp(1:maxcol)
                nlines = nlines + 1

		! read in x,y,z
		if(dt.eq.2) then
                	ra=tmp(1); dec=tmp(2); r=de_get_comovr(tmp(3))
                	call radecr_to_xyz(ra,dec,r, x,y,z)
		else
			x=tmp(1); y=tmp(2); z=tmp(3)
		endif

		! read in weight
		if(dt.eq.1) then
			weight = tmp(12) * tmp(7) * (tmp(10) + tmp(11) - 1.0_dl) ! data
		elseif(dt.eq.2) then
			weight = tmp(4) ! data ran
		elseif(dt.eq.3 .or. dt.eq.4) then ! mock, mockran
			AW = tmp(6+1); VW = tmp(7+1); wcp = tmp(8+1); wfkp = tmp(13+1)
			if(ignoremask) then
				AW = 1; VW = 1
			endif
			if(AW<0.5 .or. VW<0.5) cycle
			if(ignorecp) then
				weight = wfkp
			else
				if(wcp<0.5) cycle
				weight = wfkp * wcp
			endif
		elseif(dt.eq.5) then
                        weight = 1.0
		endif

		write(2,'(4e15.7)') x,y,z,weight
		nlinesw = nlinesw + 1

		cycle
101		exit
	enddo
	close(1); close(2);
	write(*,'(A,i8,A,i8,A,A)') 'Finishing processing ', nlines, ' lines;', nlinesw, ' lines written to ', trim(adjustl(outputfile))
end program
