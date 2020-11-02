program main_scpfrom

use LSS_cosmo_funs

implicit none

	character(len=char_len) :: cmd, clustername, dir1, filename, dir2

        print *, '###################################################'
	print *, 'This code will scp files from KIAS cluster to local directory'
	print *, 'Warning: local files maybe overwritten; be careful!'
	print *
	! gfortran fmt
	!if (iarg() .ne. 4) then
	if (iargc() .ne. 4) then
		print *, 'Usage: EXE clustername clusterdir filename localdir'
		stop
	endif

	call getarg(1, clustername)
	call getarg(2, dir1)
	call getarg(3, filename)
	call getarg(4, dir2)

	cmd = 'scp -r xiaodongli@'//trim(adjustl(clustername))//'.kias.re.kr:'//trim(adjustl(dir1))//'/'//&
		trim(adjustl(filename))//' '//trim(adjustl(dir2))
	print *, 'We will execute the command ', trim(adjustl(cmd)), '; Press any key to continue...'
	read(*,*)
	call system(cmd)

end program
