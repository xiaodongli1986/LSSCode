module periodic_funs

use LSS_constants_types
implicit none
contains

subroutine periodic_shift_ps(x,y,z, xsize, ysize, zsize, ps)
        real(dl) :: x,y,z, xsize, ysize, zsize, ps(3,26)
        integer :: i 
        i = 1
        ps(1,i) = x+xsize; ps(2,i) = y; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y+ysize; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y; ps(3,i) = z+zsize; i=i+1

        ps(1,i) = x-xsize; ps(2,i) = y; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y-ysize; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y; ps(3,i) = z-zsize; i=i+1

        ps(1,i) = x+xsize; ps(2,i) = y+ysize; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y+ysize; ps(3,i) = z+zsize; i=i+1
        ps(1,i) = x+xsize; ps(2,i) = y; ps(3,i) = z+zsize; i=i+1

        ps(1,i) = x+xsize; ps(2,i) = y-ysize; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y+ysize; ps(3,i) = z-zsize; i=i+1
        ps(1,i) = x+xsize; ps(2,i) = y; ps(3,i) = z-zsize; i=i+1

        ps(1,i) = x-xsize; ps(2,i) = y+ysize; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y-ysize; ps(3,i) = z+zsize; i=i+1
        ps(1,i) = x-xsize; ps(2,i) = y; ps(3,i) = z+zsize; i=i+1

        ps(1,i) = x-xsize; ps(2,i) = y-ysize; ps(3,i) = z; i=i+1
        ps(1,i) = x; ps(2,i) = y-ysize; ps(3,i) = z-zsize; i=i+1
        ps(1,i) = x-xsize; ps(2,i) = y; ps(3,i) = z-zsize; i=i+1

        ps(1,i) = x+xsize; ps(2,i) = y+ysize; ps(3,i) = z+zsize; i=i+1
        ps(1,i) = x-xsize; ps(2,i) = y+ysize; ps(3,i) = z+zsize; i=i+1
        ps(1,i) = x+xsize; ps(2,i) = y-ysize; ps(3,i) = z+zsize; i=i+1
        ps(1,i) = x+xsize; ps(2,i) = y+ysize; ps(3,i) = z-zsize; i=i+1
        ps(1,i) = x-xsize; ps(2,i) = y-ysize; ps(3,i) = z+zsize; i=i+1
        ps(1,i) = x-xsize; ps(2,i) = y+ysize; ps(3,i) = z-zsize; i=i+1
        ps(1,i) = x+xsize; ps(2,i) = y-ysize; ps(3,i) = z-zsize; i=i+1
        ps(1,i) = x-xsize; ps(2,i) = y-ysize; ps(3,i) = z-zsize; i=i+1
end subroutine periodic_shift_ps

subroutine ps_within_margin(ps1, n1, ps2, n2, xsize, ysize, zsize,  margin)
        integer :: i, n1, n2
        real(dl) :: ps1(3,n1), ps2(3,n1), xsize, ysize, zsize, margin, x,y,z
        n2 = 1
        do i = 1, n1
                x=ps1(1,i); y=ps1(2,i); z=ps1(3,i)
                if(-margin<x .and. x<xsize+margin .and. -margin<y .and. y<ysize+margin .and. -margin<z .and. z<zsize+margin) then
                        ps2(1,n2) = x; ps2(2,n2) = y; ps2(3,n2) = z; n2=n2+1
                endif
        enddo
        n2 = n2-1
end subroutine ps_within_margin
end module periodic_funs

program main_rho_knn_periodic_cuboid

use LSS_chisq
use periodic_funs

implicit none

	character(len=char_len) :: inputfile, outputfile, printstr, optionstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5
	integer :: i,j,k,di, xcol, ycol, zcol, numNB, maxcol, npar1, npar2, n1,n2
        integer, allocatable :: selected_list(:)
        real(dl) :: margin, xsize, ysize, zsize, tmpA(1000), x,y,z, ps1(3,26), ps2(3,26), rho,drhox,drhoy,drhoz,max_dist
        real(dl), allocatable :: par2s(:,:)
        logical :: printinfo =.True., erflag
        logical, allocatable :: compute_list(:)
	
	printstr = " Computing density around each particle. Using k nearest neighbours. Periodical cuboid. \\n"//&
                "      Usage::: EXE -inputfile inputfile -xcol 1 -ycol 2 -zcol 3 -margin 30 -numNB 30 -outputfile yourfilename -xsize 1024 -ysize 1024 -zsize 1024 \\n"//&
                " numNB: number of NNBs; outputfile name can be automatically generated (if not given)."
	if(iargc().le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif

	outputfile = ""
        xcol = 1; ycol = 2; zcol = 3; numNB = 30;
	do i = 1, iargc()
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq."-inputfile") then
			read(tmpstr2,"(A)") inputfile
		elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
			read(tmpstr2,"(A)") outputfile
		elseif(trim(adjustl(tmpstr1)).eq."-xcol") then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq."-ycol") then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq."-zcol") then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq."-numNB") then
			read(tmpstr2,*) numNB
		elseif(trim(adjustl(tmpstr1)).eq."-margin") then
			read(tmpstr2,*) margin
		elseif(trim(adjustl(tmpstr1)).eq."-xsize") then
			read(tmpstr2,*) xsize
		elseif(trim(adjustl(tmpstr1)).eq."-ysize") then
			read(tmpstr2,*) ysize
		elseif(trim(adjustl(tmpstr1)).eq."-zsize") then
			read(tmpstr2,*) zsize
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

        write(tmpstr1, *) numNB
        write(tmpstr2, '(f15.2)') margin
        write(tmpstr3, '(f15.2)') xsize
        write(tmpstr4, '(f15.2)') ysize
        write(tmpstr5, '(f15.2)') zsize
        optionstr = ".numNB"//trim(adjustl(tmpstr1))//".margin"//trim(adjustl(tmpstr2))//".xsize"//trim(adjustl(tmpstr3))//".ysize"//trim(adjustl(tmpstr4))//".zsize"//trim(adjustl(tmpstr5))
        write(*,'(A,A)') ' (rho_knn_periodic_cuboid) options: ', trim(adjustl(optionstr))
	if(trim(adjustl(outputfile)).eq."") then
                outputfile = trim(adjustl(inputfile))//trim(adjustl(optionstr))
	endif
        write(*,'(A,A)') ' (rho_knn_periodic_cuboid) inputfile: ', trim(adjustl(inputfile))
        write(*,'(A,A)') ' (rho_knn_periodic_cuboid) outputfile: ', trim(adjustl(outputfile))

        maxcol = xcol; maxcol = max(maxcol, ycol, zcol);

        ! read in the file 
        npar1 = 0; npar2 = 0
        open(file=inputfile,unit=90812743)
        open(file=outputfile,unit=981)
        do while(.true.)
                read(90812743, *, end=100) tmpA(1:maxcol)
                x = tmpA(xcol); y = tmpA(ycol); z= tmpA(zcol); npar1 = npar1+1;
                call periodic_shift_ps(x,y,z, xsize, ysize, zsize, ps1)
                call ps_within_margin(ps1, 26, ps2, n2, xsize, ysize, zsize, margin)
                npar2 = npar2+n2+1
                write(981,'(3e15.7)') x,y,z
                do j = 1, n2
                        write(981,'(3e15.7)') ps2(:,j)
                enddo
                cycle
100             exit                
        enddo
        close(90812743); close(981)
        write(*,'(A,i10,A,i10,A)') ' (rho_knn_periodic_cuboid) processing ', npar1, ' particles; ', npar2, ' within box+margin.'

        i=1; 
        allocate(compute_list(npar2)); compute_list = .false.
        open(file=inputfile,unit=90812743)
        open(file=outputfile,unit=981)
        do while(.true.)
                read(90812743, *, end=101) tmpA(1:maxcol)
                x = tmpA(xcol); y = tmpA(ycol); z= tmpA(zcol); npar1 = npar1+1;
                call periodic_shift_ps(x,y,z, xsize, ysize, zsize, ps1); compute_list(i)=.true.; i=i+1;
                call ps_within_margin(ps1, 26, ps2, n2, xsize, ysize, zsize, margin); i=i+n2; 
                write(981,'(3e15.7)') x,y,z
                do j = 1, n2
                        write(981,'(3e15.7)') ps2(:,j)
                enddo
                cycle
101             exit
        enddo
        close(90812743); close(981)





        ! expand the data to have a margin
        !allocate(par2s(3,npar2))
        !open(file=inputfile,unit=90812743)
        !npar2 = 1
        !do while(.true.)
        !        read(90812743, *, end=101) tmpA(1:maxcol)
        !        x = tmpA(xcol); y = tmpA(ycol); z= tmpA(zcol); 
        !        call periodic_shift_ps(x,y,z, boxsize, ps1)
        !        par2s(1,npar2) = x; par2s(2,npar2) = y; par2s(3,npar2) = z; npar2 = npar2+1
        !        call ps_within_margin(ps1, 26, ps2, n2, boxsize, margin)
        !        do j = 1, n2
        !                 par2s(:,npar2) = ps2(:,j); npar2=npar2+1
        !        enddo
        !        cycle
!101             exit                
!        enddo
!        close(90812743); npar2= npar2-1
!        write(*,'(A,i10,A,i10,A)') ' (rho_knn_periodic_cuboid) ', npar2, ' particles stored in par2s.'



        ! build the grid for NNB search...
        write(*,'(A,i10,A,i10,A)') ' (rho_knn_periodic_cuboid) start to building cells...'
        gb_datafile = outputfile
        gb_i_datatype = gb_dt_xyz
        gb_i_rantype = gb_noran
        gb_usenumdensity = .true.
        gb_dodensitynorm = .false.
        call cosmo_funs_init(printinfo)
        call readin_dataran(printinfo)
        call init_cosmo(om_dft,w_dft,h_dft,printinfo)
        call init_mult_lists(printinfo)
        call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)

        ! compute rho near each particle
        open(file=outputfile,unit=981)
        n1 = 0; n2=0; npar1 = 0;
        allocate(selected_list(numNB))
        di = gb_numdata / 20
        do i = 1, gb_numdata
                if(mod(i,di).eq.di-1) then
                        write(*,'(i2,A,$)') int(real(i)/real(gb_numdata)*100), '%==>'
                endif
                x = gb_xyz_list(1,i); y = gb_xyz_list(2,i); z = gb_xyz_list(3,i)
                !if(0.<x .and. x<xsize .and. 0.<y .and. y<ysize .and. 0.<z .and. z<zsize) then
                if(compute_list(i)) then
                        !call NNBExactSearch_data(x,y,z, numNB,selected_list,erflag)
                        call nb_list0(x,y,z,numNB,rho,drhox,drhoy,drhoz,max_dist,erflag)
                        !print *, rho, drhox, drhoy, drhoz, max_dist
                        if(max_dist .gt. margin) n1=n1+1
                        if(x-max_dist < -margin .or. x+max_dist > xsize+margin .or.y-max_dist < -margin .or. y+max_dist > &
                               ysize+margin .or. z-max_dist < -margin .or. z+max_dist > zsize+margin) n2=n2+1
                        npar1 = npar1 + 1
                        write(981,'(8e15.7)') x,y,z,rho,drhox,drhoy,drhoz,max_dist
                endif
        enddo
        print *
        write(*,'(A,i10,A,i10,A,f5.2,A,i10,A,f5.2,A)') ' (rho_knn_periodic_cuboid) Finishing rho-estimate for ', npar1, '  particles. ', n1, '(',&
                real(n1)/real(npar1)*100, '%) have   max_dist > margin;  ', n2, '(', real(n2)/real(npar2)*100, '%) have  possible neighbors lying out!'
        close(981)



        write(*,'(A)') ' (rho_knn_periodic_cuboid) Done.'


end program
