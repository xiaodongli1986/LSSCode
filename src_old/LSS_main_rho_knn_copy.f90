
program main_rho_knn

use LSS_constants_types
use LSS_chisq

implicit none

	character(len=char_len) :: inputfile, outputfile, printstr, optionstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5
	integer :: i,j,k,di, xcol, ycol, zcol, wcol, numNB, maxcol, npar1,  n1,n2
        integer, allocatable :: selected_list(:)
        real(dl) :: tmpA(1000), x,y,z, ps1(3,26), ps2(3,26), rho,drhox,drhoy,drhoz,max_dist
        real(dl), allocatable :: par2s(:,:)
        logical :: printinfo =.True., erflag, hasweight = .false.
        logical, allocatable :: compute_list(:)
	
	printstr = " Computing density around each particle. Using k nearest neighbours. Periodical cuboid. \\n"//&
                "      Usage::: EXE -inputfile inputfile -xcol 1 -ycol 2 -zcol 3 -numNB 30 -outputfile yourfilename  -hasweight T/F  \\n"//&
                " numNB: number of NNBs; outputfile name can be automatically generated (if not given)."
	if(iargc().le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif

	outputfile = ""
        xcol = 1; ycol = 2; zcol = 3; wcol = 4; numNB = 30;
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
		elseif(trim(adjustl(tmpstr1)).eq."-wcol") then
			read(tmpstr2,*) wcol
		elseif(trim(adjustl(tmpstr1)).eq."-numNB") then
			read(tmpstr2,*) numNB
		elseif(trim(adjustl(tmpstr1)).eq."-hasweight") then
			read(tmpstr2,*) hasweight
		else
			print *, "Unkown argument: ", trim(adjustl(tmpstr1))
			write(*,"(A)") trim(adjustl(printstr))
			stop
		endif
	enddo

        write(tmpstr1, *) numNB
        optionstr = ".numNB"//trim(adjustl(tmpstr1))
        if(hasweight) optionstr = trim(adjustl(optionstr))//'.weightdensity'
        write(*,'(A,A)') ' (rho_knn) options: ', trim(adjustl(optionstr))
	if(trim(adjustl(outputfile)).eq."") then
                outputfile = trim(adjustl(inputfile))//trim(adjustl(optionstr))
	endif
        write(*,'(A,A)') ' (rho_knn) inputfile: ', trim(adjustl(inputfile))
        write(*,'(A,A)') ' (rho_knn) outputfile: ', trim(adjustl(outputfile))

        maxcol = xcol; 
        if(hasweight) then
                maxcol = max(maxcol, ycol, zcol, wcol);
        else
                maxcol = max(maxcol, ycol, zcol);
        endif

        ! read in the file 
        npar1 = 0; 
        open(file=inputfile,unit=90812743)
        open(file=outputfile,unit=981)
        do while(.true.)
                read(90812743, *, end=100) tmpA(1:maxcol)
                x = tmpA(xcol); y = tmpA(ycol); z= tmpA(zcol); npar1 = npar1+1;
                if(hasweight) then
                        write(981,'(4e15.7)') x,y,z, tmpA(wcol)
                else
                        write(981,'(3e15.7)') x,y,z
                endif

                cycle
100             exit                
        enddo
        close(90812743); close(981)
        write(*,'(A,i10,A)') ' (rho_knn) processing ', npar1, ' particles; '




        ! build the grid for NNB search...
        write(*,'(A,i10,A,i10,A)') ' (rho_knn) start to building cells...'
        gb_datafile = outputfile
        !gb_i_datatype = gb_dt_xyz
        if(hasweight) then
                gb_i_datatype = gb_dt_xyzw
        else
                gb_i_datatype = gb_dt_xyz
        endif

        gb_i_rantype = gb_noran
        if(hasweight) then
                gb_usenumdensity = .false.
        else
                gb_usenumdensity = .true.
        endif

        gb_dodensitynorm = .false.
        call cosmo_funs_init(printinfo)
        call readin_dataran(printinfo)
        call init_cosmo(om_dft,w_dft,h_dft,printinfo)
        call init_mult_lists(printinfo)
        call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)

        ! compute rho near each particle
        open(file=outputfile,unit=981)
        npar1 = 0;
        allocate(selected_list(numNB))
        di = gb_numdata / 20
        do i = 1, gb_numdata
                if(mod(i,di).eq.di-1) then
                        write(*,'(i2,A,$)') int(real(i)/real(gb_numdata)*100), '%==>'
                endif
                x = gb_xyz_list(1,i); y = gb_xyz_list(2,i); z = gb_xyz_list(3,i)
                !if(0.<x .and. x<xsize .and. 0.<y .and. y<ysize .and. 0.<z .and. z<zsize) then
                if(.true.) then
                        !call NNBExactSearch_data(x,y,z, numNB,selected_list,erflag)
                        call nb_list0(x,y,z,numNB,rho,drhox,drhoy,drhoz,max_dist,erflag)
                        !print *, rho, drhox, drhoy, drhoz, max_dist
                        npar1 = npar1 + 1
                        write(981,'(8e15.7)') x,y,z,rho,drhox,drhoy,drhoz,max_dist
                endif
        enddo
        print *
        write(*,'(A,i10,A)') ' (rho_knn) Finishing rho-estimate for ', npar1, '  particles. '
        close(981)



        write(*,'(A)') ' (rho_knn) Done.'


end program
