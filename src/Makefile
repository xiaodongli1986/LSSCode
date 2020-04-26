EXE0 =  ../bin/LSS_test
EXE1 =  ../bin/LSS_calc_grad_chisqs 
EXE2 =  ../bin/LSS_BSK_calc
EXE3 =  ../bin/LSS_BSK_Scan
EXE4 =  ../bin/LSS_rmass_stat
EXE5 =  ../bin/LSS_SnapshotRSDmock
EXE6 =  ../bin/LSS_SnapshotGF
EXE7 =  ../bin/LSS_MakeLightCone
EXE8 =  ../bin/LSS_MakeCompactSample
EXE9 =  ../bin/LSS_CosmoConvert
EXE10 = ../bin/LSS_Make-1eighth-Shell-Random
EXE11 = ../bin/LSS_Rbinned-Split-Data
EXE12 = ../bin/LSS_mpi_sh
EXE13 = ../bin/LSS_TakePartofFile
EXE14 = ../bin/LSS_RaDecZToXYZ
EXE15 = ../bin/LSS_Make-allsky-Random
EXE16 = ../bin/LSS_constant-nbar-down-sampling
EXE17 = ../bin/LSS_DoCut-wrt-SingleColumn
EXE18 = ../bin/LSS_RMSv-with-masscut
EXE19 = ../bin/LSS_HR4LC
EXE20 = ../bin/LSS_HR4LC_stat
EXE21 = ../bin/LSS_HR4LC_masscut
EXE22 = ../bin/LSS_scpfromkiascluster
EXE23 = ../bin/LSS_scptokiascluster
EXE24 = ../bin/LSS_make-allsky-random
EXE25 = ../bin/LSS_add-0vxvyvz-randmass
EXE26 = ../bin/LSS_HR4LC_stat_massoutlier
EXE27 = ../bin/LSS_grid_test
EXE28 = ../bin/LSS_grid_list_fixr_rhobar
EXE29 = ../bin/LSS_boss2pcf_fmtxyzwfile
EXE30 = ../bin/LSS_HR4LC_SliceSelection
EXE31 = ../bin/LSS_HR4LC_masscut_smallbox
EXE32 = ../bin/LSS_2pcf__numline_sumweight
EXE33 = ../bin/LSS_HR4LC_SliceSelection_ThinSlices
EXE34 = ../bin/LSS_2pcf
EXE35 = ../bin/LSS_3pcf
EXE36 = ../bin/LSS_XYZToRaDecZ
EXE37 = ../bin/LSS_converttoRhct
EXE38 = ../bin/LSS_LSS_ximu
EXE39 = ../bin/LSS_CosmoConvert_CPL
EXE40 = ../bin/LSS_XYZWToRaDecZW
EXE41 = ../bin/LSS_rho_knn_periodic_box
EXE42 = ../bin/LSS_mpi_random_select
EXE43 = ../bin/LSS_add_rsd_lightcone
EXE44 = ../bin/LSS_mpi_add1
EXE45 = ../bin/LSS_mpi_lightcone_boxsplit
EXE46 = ../bin/LSS_rho_knn_rectangular_box
EXE47 = ../bin/LSS_cic_b
EXE49 = ../bin/LSS_binary_to_ascii
EXE50 = ../bin/LSS_cic
EXE51 = ../bin/LSS_gadget_to_fmt3
EXE52 = ../bin/LSS_lpicola_to_fmt3
EXE53 = ../bin/LSS_lpicola_to_gadget
EXE54 = ../bin/LSS_lpicola_lightcone_boxsplit
EXE55 = ../bin/LSS_gadget_print_head
EXE56 = ../bin/LSS_rockstar_select_xyzvxvyvz_mvir_vmax
EXE57 = ../bin/LSS_gadget_change_head
EXE58 = ../bin/LSS_cic_calcv
OBJS = LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_settings_init.o LSS_cell.o LSS_smooth.o LSS_grad_fields.o LSS_chisq.o LSS_grid_tools.o LSS_BSK.o LSS_2pcf.o LSS_3pcf.o
#F90C = ifort -heap-arrays -diag-disable 6009 #gfortran #mpif90 -heap-arrays
#OBJS = LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_settings_init.o LSS_cell.o LSS_smooth.o LSS_grad_fields.o LSS_chisq.o LSS_grid_tools.o LSS_BSK.o LSS_2pcf.o LSS_3pcf.o #LSS_kdtree3pcf.o kdtree2.o
#F90C = mpif90 -heap-arrays #gfortran #mpif90 -heap-arrays
F90C = mpiifort -heap-arrays #gfortran #mpif90 -heap-arrays
#F90C = ifort -heap-arrays #gfortran #mpif90 -heap-arrays
F90Cnompi = ifort -heap-arrays

#F90C = ifort -heap-arrays
#F90C = gfortran
#F90C = /home/xiaodongli/software/anaconda/bin/mpif90
F90FLAGS = 

default: $(EXE2)
nompi: $(EXE0)         $(EXE2)         $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11)          $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32)  $(EXE33) $(EXE34) $(EXE35) $(EXE36) $(EXE37) $(EXE38) $(EXE39) $(EXE40) $(EXE41)          $(EXE43)                   $(EXE46) $(EXE47)  $(EXE49) $(EXE50) $(EXE51) $(EXE52) $(EXE53) $(EXE54) $(EXE55) $(EXE56) $(EXE57) $(EXE58)
all:   $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32) $(EXE33) $(EXE34) $(EXE35) $(EXE36) $(EXE37) $(EXE38) $(EXE39) $(EXE40) $(EXE41) $(EXE42) $(EXE43) $(EXE44) $(EXE45) $(EXE46) $(EXE47)  $(EXE49) $(EXE50) $(EXE51) $(EXE52) $(EXE53) $(EXE54) $(EXE55) $(EXE56) $(EXE57) $(EXE58)

test: $(EXE0)
bsk: $(EXE2)
bskscan:  $(EXE3)
rmass: $(EXE4)
SnapshotRSDmock: $(EXE5)
SnapshotGF: $(EXE6)
MakeLightCone: $(EXE7)
MakeCompactSample: $(EXE8)
CosmoConvert: $(EXE9)
MakeRandom: $(EXE10)
Split-Data: $(EXE11)
mpi_sh: $(EXE12)
TakePartofFile: $(EXE13)
RaDecZToXYZ: $(EXE14)
constant-nbar-down-sampling: $(EXE16)
DoCut-wrt-SingleColumn: $(EXE17)
RMSv-with-masscut: $(EXE18)
HR4LC: $(EXE19)
HR4LC_stat: $(EXE20)
HR4LC_masscut: $(EXE21)
scpfromkiascluster: $(EXE22)
scptokiascluster: $(EXE23)
make-allsky-random: $(EXE24)
add-0vxvyvz-randmass: $(EXE25)
HR4LC_stat_massoutlier: $(EXE26)
grid_test: $(EXE27)
grid_list_fixr_rhobar: $(EXE28)
boss2pcf_fmtxyzwfile: $(EXE29)
HR4LC_SliceSelection: $(EXE30)
HR4LC_masscut_smallbox: $(EXE31)
2pcf__numline_sumweight: $(EXE32)
HR4LC_SliceSelection_ThinSlices: $(EXE33)
2pcf: $(EXE34)
3pcf: $(EXE35)
XYZToRaDecZ: $(EXE36)
converttoRhct: $(EXE37)
LSS_ximu: $(EXE38)
CosmoConvert_CPL: $(EXE39)
XYZWToRaDecZW: $(EXE40)
rho_knn_periodic_box: $(EXE41)
mpi_random_select: $(EXE42)
add_rsd_lightcone: $(EXE43)
mpi_add1: $(EXE44)
mpi_lightcone_boxsplit: $(EXE45)
rho_knn_rectangular_box: $(EXE46)
cic_b: $(EXE47)
binary_to_ascii: $(EXE49)
cic: $(EXE50)
gadget_to_fmt3: $(EXE51)
lpicola_to_fmt3: $(EXE52)
lpicola_to_gadget: $(EXE53)
lpicola_lightcone_boxsplit: $(EXE54)
gadget_print_head: $(EXE55)
rockstar_select_xyzvxvyvz_mvir_vmax: $(EXE56)
gadget_change_head: $(EXE57)
cic_calcv: $(EXE58)
AMTB: AMTB

$(EXE0): $(OBJS) LSS_main_test.o
	$(F90C) -o $(EXE0) $(OBJS) LSS_main_test.o
$(EXE1): $(OBJS) LSS_main_calc_grad_chisqs.o
	$(F90C) -o $(EXE1) $(OBJS) LSS_main_calc_grad_chisqs.o
$(EXE2): $(OBJS) LSS_main_BSK_calc.o
	$(F90C) -o $(EXE2) $(OBJS) LSS_main_BSK_calc.o
$(EXE3): $(OBJS) LSS_main_BSK_Scan.o
	$(F90C) -o $(EXE3) $(OBJS) LSS_main_BSK_Scan.o
$(EXE4): $(OBJS) LSS_main_r_mass_stat.o
	$(F90C) -o $(EXE4) $(OBJS) LSS_main_r_mass_stat.o
$(EXE5): $(OBJS) LSS_main_SnapshotRSDmock.o
	$(F90C) -o $(EXE5) $(OBJS) LSS_main_SnapshotRSDmock.o
$(EXE6): $(OBJS) LSS_main_SnapshotGF.o
	$(F90C) -o $(EXE6) $(OBJS) LSS_main_SnapshotGF.o
$(EXE7): $(OBJS) LSS_main_MakeLightCone.o
	$(F90C) -o $(EXE7) $(OBJS) LSS_main_MakeLightCone.o
$(EXE8): $(OBJS) LSS_main_MakeCompactSample.o
	$(F90C) -o $(EXE8) $(OBJS) LSS_main_MakeCompactSample.o
$(EXE9): LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_CosmoConvert.o
	$(F90C) -o $(EXE9) LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_CosmoConvert.o
$(EXE10): LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_make-1eighth-shell-random.o
	$(F90C) -o $(EXE10) LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o  LSS_main_make-1eighth-shell-random.o
$(EXE11): LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_rbinned-split-data.o
	$(F90C) -o $(EXE11) LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_rbinned-split-data.o
$(EXE12): LSS_tools.o LSS_main_mpi_sh.o
	$(F90C) -o $(EXE12) LSS_tools.o  LSS_main_mpi_sh.o
$(EXE13): LSS_tools.o LSS_main_TakePartofFile.o
	$(F90C) -o $(EXE13) LSS_tools.o  LSS_main_TakePartofFile.o
$(EXE14): LSS_tools.o LSS_cosmo_funs.o LSS_main_RaDecZToXYZ.o
	$(F90C) -o $(EXE14) LSS_tools.o  LSS_cosmo_funs.o LSS_main_RaDecZToXYZ.o
$(EXE15): LSS_tools.o LSS_main_make-allsky-random.o
	$(F90C) -o $(EXE15) LSS_tools.o  LSS_main_make-allsky-random.o
$(EXE16): LSS_tools.o LSS_main_constant-nbar-down-sampling.o
	$(F90C) -o $(EXE16) LSS_tools.o  LSS_main_constant-nbar-down-sampling.o
$(EXE17): LSS_tools.o LSS_main_DoCut-wrt-SingleColumn.o
	$(F90C) -o $(EXE17) LSS_tools.o  LSS_main_DoCut-wrt-SingleColumn.o
$(EXE18): LSS_tools.o LSS_main_RMSv-with-masscut.o
	$(F90C) -o $(EXE18) LSS_tools.o  LSS_main_RMSv-with-masscut.o
$(EXE19): LSS_tools.o LSS_main_HR4LC.o
	$(F90C) -o $(EXE19) LSS_tools.o  LSS_main_HR4LC.o
$(EXE20): LSS_cosmo_funs.o LSS_main_HR4LC_stat.o
	$(F90C) -o $(EXE20) LSS_tools.o LSS_cosmo_funs.o  LSS_main_HR4LC_stat.o
$(EXE21): LSS_cosmo_funs.o LSS_main_HR4LC_masscut.o
	$(F90C) -o $(EXE21) LSS_tools.o LSS_cosmo_funs.o  LSS_main_HR4LC_masscut.o
$(EXE22): LSS_cosmo_funs.o LSS_main_scpfromkiascluster.o
	$(F90C) -o $(EXE22) LSS_tools.o LSS_cosmo_funs.o  LSS_main_scpfromkiascluster.o
$(EXE23): LSS_cosmo_funs.o LSS_main_scptokiascluster.o
	$(F90C) -o $(EXE23) LSS_tools.o LSS_cosmo_funs.o  LSS_main_scptokiascluster.o
$(EXE24): LSS_cosmo_funs.o LSS_main_make-allsky-random.o
	$(F90C) -o $(EXE24) LSS_tools.o LSS_cosmo_funs.o  LSS_main_make-allsky-random.o
$(EXE25): LSS_cosmo_funs.o LSS_main_add-0vxvyvz-randmass.o
	$(F90C) -o $(EXE25) LSS_tools.o LSS_cosmo_funs.o  LSS_main_add-0vxvyvz-randmass.o
$(EXE26): LSS_cosmo_funs.o LSS_main_HR4LC_stat_massoutlier.o
	$(F90C) -o $(EXE26) LSS_tools.o LSS_cosmo_funs.o  LSS_main_HR4LC_stat_massoutlier.o
$(EXE27): LSS_cosmo_funs.o LSS_grid_tools.o LSS_main_grid_test.o
	$(F90C) -o $(EXE27) LSS_tools.o LSS_cosmo_funs.o  LSS_grid_tools.o LSS_main_grid_test.o
$(EXE28): LSS_cosmo_funs.o LSS_grid_tools.o LSS_main_grid_list_fixr_rhobar.o
	$(F90C) -o $(EXE28) LSS_tools.o LSS_cosmo_funs.o  LSS_grid_tools.o LSS_main_grid_list_fixr_rhobar.o
$(EXE29): LSS_cosmo_funs.o LSS_main_boss2pcf_fmtxyzwfile.o
	$(F90C) -o $(EXE29) LSS_tools.o LSS_cosmo_funs.o  LSS_main_boss2pcf_fmtxyzwfile.o
$(EXE30): LSS_cosmo_funs.o LSS_main_HR4LC_SliceSelection.o
	$(F90C) -o $(EXE30) LSS_tools.o LSS_cosmo_funs.o  LSS_main_HR4LC_SliceSelection.o
$(EXE31): LSS_cosmo_funs.o LSS_main_HR4LC_masscut_smallbox.o
	$(F90C) -o $(EXE31) LSS_tools.o LSS_cosmo_funs.o  LSS_main_HR4LC_masscut_smallbox.o
$(EXE32): LSS_cosmo_funs.o LSS_main_2pcf__numline_sumweight.o
	$(F90C) -o $(EXE32) LSS_tools.o LSS_cosmo_funs.o  LSS_main_2pcf__numline_sumweight.o
$(EXE33): LSS_cosmo_funs.o LSS_main_HR4LC_SliceSelection_ThinSlices.o
	$(F90C) -o $(EXE33) LSS_tools.o LSS_cosmo_funs.o  LSS_main_HR4LC_SliceSelection_ThinSlices.o
$(EXE34): $(OBJS) LSS_main_2pcf.o
	$(F90C) -o $(EXE34) $(OBJS) LSS_main_2pcf.o
$(EXE35): $(OBJS) LSS_main_3pcf.o
	$(F90C) -o $(EXE35) $(OBJS)  LSS_main_3pcf.o
$(EXE36): LSS_cosmo_funs.o LSS_main_XYZToRaDecZ.o
	$(F90C) -o $(EXE36) LSS_tools.o LSS_cosmo_funs.o  LSS_main_XYZToRaDecZ.o
$(EXE37): LSS_cosmo_funs.o LSS_main_converttoRhct.o
	$(F90C) -o $(EXE37) LSS_tools.o LSS_cosmo_funs.o  LSS_main_converttoRhct.o
$(EXE38): LSS_cosmo_funs.o LSS_main_LSS_ximu.o
	$(F90C) -o $(EXE38) LSS_tools.o LSS_cosmo_funs.o  LSS_main_LSS_ximu.o
$(EXE39): LSS_cosmo_funs.o LSS_main_CosmoConvert_CPL.o
	$(F90C) -o $(EXE39) LSS_tools.o LSS_cosmo_funs.o  LSS_main_CosmoConvert_CPL.o
$(EXE40): LSS_cosmo_funs.o LSS_main_XYZWToRaDecZW.o
	$(F90C) -o $(EXE40) LSS_tools.o LSS_cosmo_funs.o  LSS_main_XYZWToRaDecZW.o
$(EXE41):  $(OBJS) LSS_main_rho_knn_periodic_box.o
	$(F90C) -o $(EXE41) $(OBJS)  LSS_main_rho_knn_periodic_box.o
$(EXE42): LSS_cosmo_funs.o LSS_main_mpi_random_select.o
	$(F90C) -o $(EXE42) LSS_tools.o LSS_cosmo_funs.o  LSS_main_mpi_random_select.o
$(EXE43): LSS_cosmo_funs.o LSS_main_add_rsd_lightcone.o
	$(F90C) -o $(EXE43) LSS_tools.o LSS_cosmo_funs.o  LSS_main_add_rsd_lightcone.o
$(EXE44): LSS_cosmo_funs.o LSS_main_mpi_add1.o
	$(F90C) -o $(EXE44) LSS_tools.o LSS_cosmo_funs.o  LSS_main_mpi_add1.o
$(EXE45): LSS_cosmo_funs.o LSS_main_mpi_lightcone_boxsplit.o
	$(F90C) -o $(EXE45) LSS_tools.o LSS_cosmo_funs.o  LSS_main_mpi_lightcone_boxsplit.o
$(EXE46): LSS_cosmo_funs.o LSS_main_rho_knn_rectangular_box.o
	$(F90C) -o $(EXE46) $(OBJS)  LSS_main_rho_knn_rectangular_box.o
$(EXE47): LSS_cosmo_funs.o LSS_main_cic_b.o
	$(F90C) -o $(EXE47) LSS_tools.o LSS_cosmo_funs.o  LSS_main_cic_b.o
$(EXE49): LSS_cosmo_funs.o LSS_main_binary_to_ascii.o
	$(F90C) -o $(EXE49) LSS_tools.o LSS_cosmo_funs.o  LSS_main_binary_to_ascii.o
$(EXE50): LSS_cosmo_funs.o LSS_main_cic.o
	$(F90C) -o $(EXE50) LSS_tools.o LSS_cosmo_funs.o  LSS_main_cic.o
$(EXE51): LSS_cosmo_funs.o LSS_main_gadget_to_fmt3.o
	$(F90C) -o $(EXE51) LSS_tools.o LSS_cosmo_funs.o  LSS_main_gadget_to_fmt3.o
$(EXE52): LSS_cosmo_funs.o LSS_main_lpicola_to_fmt3.o
	$(F90C) -o $(EXE52) LSS_tools.o LSS_cosmo_funs.o  LSS_main_lpicola_to_fmt3.o
$(EXE53): LSS_cosmo_funs.o LSS_main_lpicola_to_gadget.o
	$(F90C) -o $(EXE53) LSS_tools.o LSS_cosmo_funs.o  LSS_main_lpicola_to_gadget.o
$(EXE54): LSS_cosmo_funs.o LSS_main_lpicola_lightcone_boxsplit.o
	$(F90C) -o $(EXE54) LSS_tools.o LSS_cosmo_funs.o  LSS_main_lpicola_lightcone_boxsplit.o
$(EXE55): LSS_cosmo_funs.o LSS_main_gadget_print_head.o
	$(F90C) -o $(EXE55) LSS_tools.o LSS_cosmo_funs.o  LSS_main_gadget_print_head.o
$(EXE56): LSS_cosmo_funs.o LSS_main_rockstar_select_xyzvxvyvz_mvir_vmax.o
	$(F90C) -o $(EXE56) LSS_tools.o LSS_cosmo_funs.o  LSS_main_rockstar_select_xyzvxvyvz_mvir_vmax.o
$(EXE57): LSS_cosmo_funs.o LSS_main_gadget_change_head.o
	$(F90C) -o $(EXE57) LSS_tools.o LSS_cosmo_funs.o  LSS_main_gadget_change_head.o
$(EXE58): LSS_cosmo_funs.o LSS_main_cic_calcv.o
	$(F90C) -o $(EXE58) LSS_tools.o LSS_cosmo_funs.o  LSS_main_cic_calcv.o
amtb: amtb


LSS_cosmo_funs.o: LSS_tools.o
LSS_settings_init.o: LSS_cosmo_funs.o 
LSS_cell.o: LSS_settings_init.o
LSS_smooth.o: LSS_cell.o
LSS_grad_fields.o: LSS_smooth.o LSS_mu_tools.o
LSS_chisq.o: LSS_grad_fields.o 
LSS_2pcf.o: LSS_chisq.o
LSS_main_test.o:  LSS_chisq.o 
LSS_main_calc_grad_chisqs.o:  LSS_chisq.o 
LSS_BSK.o : LSS_chisq.o
#LSS_kdtree3pcf.o : kdtree2.o
LSS_main_BSK_calc.o: LSS_BSK.o 
LSS_main_test.o: LSS_settings_init.o
LSS_main_BSK_Scan.o: LSS_BSK.o
LSS_main_r_mass_stat.o: LSS_tools.o
LSS_main_SnapshotRSDmock.o: LSS_cosmo_funs.o
LSS_main_SnapshotGF.o: LSS_chisq.o
LSS_main_MakeLightCone.o: LSS_settings_init.o
LSS_main_CosmoConvert.o: LSS_cosmo_funs.o
LSS_main_make-1eighth-shell-random.o: LSS_cosmo_funs.o
LSS_main_rbinned-split-data.o: LSS_cosmo_funs.o
LSS_main_mpi_sh.o: LSS_tools.o
LSS_main_TakePartofFile.o: LSS_tools.o
LSS_main_RaDecZToXYZ.o: LSS_cosmo_funs.o
LSS_main_make-allsky-random.o: LSS_tools.o
LSS_main_constant-nbar-down-sampling: LSS_tools.o
LSS_main_DoCut-wrt-SingleColumn: LSS_tools.o
LSS_main_RMSv-with-masscut: LSS_tools.o
LSS_main_HR4LC: LSS_tools.o
LSS_main_HR4LC_stat: LSS_cosmo_funs.o
LSS_main_HR4LC_masscut: LSS_cosmo_funs.o
LSS_main_scpfromkiascluster: LSS_cosmo_funs.o
LSS_main_scptokiascluster: LSS_cosmo_funs.o
LSS_main_make-allsky-random: LSS_cosmo_funs.o
LSS_main_add-0vxvyvz-randmass: LSS_cosmo_funs.o
LSS_main_HR4LC_stat_massoutlier: LSS_cosmo_funs.o
LSS_main_grid_test: LSS_cosmo_funs.o LSS_grid_tools.o
LSS_main_grid_list_fixr_rhobar: LSS_cosmo_funs.o
LSS_main_boss2pcf_fmtxyzwfile: LSS_cosmo_funs.o
LSS_main_HR4LC_SliceSelection: LSS_cosmo_funs.o
LSS_main_HR4LC_masscut_smallbox: LSS_cosmo_funs.o
LSS_main_2pcf__numline_sumweight: LSS_cosmo_funs.o
LSS_main_HR4LC_SliceSelection_ThinSlices: LSS_cosmo_funs.o
LSS_main_2pcf: LSS_cosmo_funs.o LSS_2pcf.o
LSS_main_3pcf: LSS_cosmo_funs.o LSS_3pcf.o
LSS_main_XYZToRaDecZ: LSS_cosmo_funs.o
LSS_main_converttoRhct: LSS_cosmo_funs.o
LSS_main_LSS_ximu: LSS_cosmo_funs.o
LSS_main_CosmoConvert_CPL: LSS_cosmo_funs.o
LSS_main_XYZWToRaDecZW: LSS_cosmo_funs.o
LSS_main_rho_knn_periodic_box: LSS_cosmo_funs.o
LSS_main_mpi_random_select: LSS_cosmo_funs.o
LSS_main_add_rsd_lightcone: LSS_cosmo_funs.o
LSS_main_mpi_add1: LSS_cosmo_funs.o
LSS_main_mpi_lightcone_boxsplit: LSS_cosmo_funs.o
LSS_main_rho_knn_rectangular_box: LSS_cosmo_funs.o
LSS_main_cic_b: LSS_cosmo_funs.o
LSS_main_binary_to_ascii: LSS_cosmo_funs.o
LSS_main_binary_to_ascii: LSS_cosmo_funs.o
LSS_main_cic: LSS_cosmo_funs.o
LSS_main_gadget_to_fmt3: LSS_cosmo_funs.o
LSS_main_lpicola_to_fmt3: LSS_cosmo_funs.o
LSS_main_lpicola_to_gadget: LSS_cosmo_funs.o
LSS_main_lpicola_lightcone_boxsplit: LSS_cosmo_funs.o
LSS_main_gadget_print_head: LSS_cosmo_funs.o
LSS_main_rockstar_select_xyzvxvyvz_mvir_vmax: LSS_cosmo_funs.o
LSS_main_gadget_change_head: LSS_cosmo_funs.o
LSS_main_cic_calcv: LSS_cosmo_funs.o
GYPS: GYPS
#LSS_main.o: LSS_settings_init.o
#LSS_main.o: LSS_tools.o
#LSS_main.o: LSS_cosmo_funs.o

%.o: %.f90
	$(F90C) -c $*.f90 $(F90FLAGS)
%.o: %.F90
	$(F90C) -c $< -o $@

clean :
	rm -rf *.o *.mod $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(OBJS) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32)  $(EXE33) $(EXE34) $(EXE35) $(EXE36) $(EXE37) $(EXE38) $(EXE39) $(EXE40) $(EXE41) $(EXE42) $(EXE43) $(EXE44) $(EXE45) $(EXE46) $(EXE47) $(EXE49) $(EXE50) $(EXE51) $(EXE52) $(EXE53) $(EXE54) $(EXE55) $(EXE56) $(EXE57) $(EXE58)
