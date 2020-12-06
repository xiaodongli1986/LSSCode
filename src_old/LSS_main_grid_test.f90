program main_grid_test

use LSS_grid_tools

implicit none

	integer, parameter :: np = 200000, nx = 100
	integer :: i,j,l, ix,iy,iz, ip, numNB,iplist(maxnumip),iNB
	real(dl) :: x,y,z,xs(np), ys(np), zs(np), ws(np), rlist(maxnumip), r0(3),r1(3),fixr,&
		rhomean,rho,drhox,drhoy,drhoz,nbarmean,nbar,dnbarx,dnbary,dnbarz
	logical :: touchboundary
	type(gridtype) :: grid

	! random list of x,y,z,weight
	call random_seed()
	call random_number(xs)
	call random_number(ys)
	call random_number(zs)
	call random_number(ws)
	do i = 1, np
!		ws(i) = 0.5 !exp(10.0 + 5.0*ws(i))
	enddo

	if(.true.) then
		open(unit=1,file='grid_rho_drho_testdata.txt')
		open(unit=2,file='grid_rho_drho_testxyz.txt')
		do i = 1, np
			write(1,'(4e15.7)') xs(i), ys(i), zs(i), ws(i)
		enddo
		do j = 1, 3000
			call random_number(x)
			call random_number(y)
			call random_number(z)
			write(2,'(3e15.7)') x,y,z
		enddo
		close(1); close(2)
	endif

	! build cell
	call grid_init(grid,xs,ys,zs,ws,np,nx,.true.)

	! print elements in some cell
	print *, 'Print elements in some cells...'
	do ix=10,10
	do iy=20,22
	do iz=30,30
		if(grid%cellmat(ix,iy,iz)%np>0) then
			print *, 'Cell-index: ', ix,iy,iz 
			do ip = 1, grid%cellmat(ix,iy,iz)%np
				i = grid%cellmat(ix,iy,iz)%iplist(ip)
				print *, real(xs(i)),real(ys(i)),real(zs(i))
			enddo
		endif
	enddo
	enddo
	enddo

	! fixed distance search
	print *, 'Do a fixed distance search: print list of ip and r...'

	r0 = 0.5_dl
	fixr = 0.05_dl
	print *, 'Searching near position ', real(r0)
	write(*,'(A,f10.3)') ' Searching with maximal distance  ', real(fixr)
	call grid_fixr_search(grid, r0(1), r0(2), r0(3), fixr, numNB,touchboundary,iplist, rlist)

	print *, 'Check touchbounday: ', touchboundary
	print *, 'Number of found neighbors ', numNB
	do iNB = 1, numNB
		ip = iplist(iNB)
		r1(1)=grid%xs(ip)
		r1(2)=grid%ys(ip)
		r1(3)=grid%zs(ip)
		write(*,'(A,i11,3e15.7,3x,f10.5)') 'ip, position, distance =',ip, real(r1), real(distance(r0,r1,3))
	enddo

	! compute rho, drho
	call grid_fixr_rho_drho(grid, r0(1), r0(2), r0(3), fixr, &
		numNB,touchboundary, &
		rhomean,rho,drhox,drhoy,drhoz,&
		nbarmean,nbar,dnbarx,dnbary,dnbarz)
	print *, 'numNB = ', numNB
	print *, 'rhomean = ', 	rhomean
	print *, 'rho = ', rho
	print *, 'drhos = ', real(drhox),real(drhoy),real(drhoz)

	print *, 'nbarmean = ', 	nbarmean
	print *, 'nbar = ', nbar
	print *, 'dnbars = ', real(dnbarx),real(dnbary),real(dnbarz)
	! Check there is no diff between the arrays input and arrays in grid
!	write(*,'(<np>(f4.1))'), grid%xs - xs!
!	write(*,'(<np>(f4.1))'), grid%ys - ys
	!write(*,'(<np>(f4.1))'), grid%zs - zs
!	write(*,'(<np>(f4.1))'), grid%ws - ws

end program
