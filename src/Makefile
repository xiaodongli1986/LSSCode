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
OBJS = LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_settings_init.o LSS_cell.o LSS_smooth.o LSS_grad_fields.o LSS_chisq.o LSS_grid_tools.o LSS_BSK.o
F90C = mpif90 -heap-arrays #gfortran #mpif90 -heap-arrays
#F90C = ifort -heap-arrays
#F90C = gfortran
#F90C = /home/xiaodongli/software/anaconda/bin/mpif90
F90FLAGS = 

default: $(EXE2)
nompi: $(EXE0)         $(EXE2)         $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11)          $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32)  $(EXE33)
all:   $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32) $(EXE33)

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
amtb: amtb


LSS_cosmo_funs.o: LSS_tools.o
LSS_settings_init.o: LSS_cosmo_funs.o 
LSS_cell.o: LSS_settings_init.o
LSS_smooth.o: LSS_cell.o
LSS_grad_fields.o: LSS_smooth.o LSS_mu_tools.o
LSS_chisq.o: LSS_grad_fields.o 
LSS_main_test.o:  LSS_chisq.o 
LSS_main_calc_grad_chisqs.o:  LSS_chisq.o 
LSS_BSK.o : LSS_chisq.o
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
GYPS: GYPS
#LSS_main.o: LSS_settings_init.o
#LSS_main.o: LSS_tools.o
#LSS_main.o: LSS_cosmo_funs.o

%.o: %.f90
	$(F90C) -c $*.f90 $(F90FLAGS)

clean :
	rm -rf *.o *.mod $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(OBJS) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32)  $(EXE33)
