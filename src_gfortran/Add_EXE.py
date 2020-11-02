#!/usr/bin/python
import sys
import commands
import os

def str_to_numbers(str1,exitcode='#', do_float_conver=True):
    floatlist = []
    str2 = ''
    numberunderconstruct = False
    for i in range(len(str1)):
        duru = str1[i]
        if duru == exitcode:
            break
        elif duru == '\n' or duru == ' ' or duru == '\t' or duru == ',':
            if numberunderconstruct == True:
                floatlist.append(str2); str2 = ''
                numberunderconstruct = False
        else:
            if numberunderconstruct == False:
                numberunderconstruct = True
                str2 = duru
            else:
                str2 += duru
    if str2 != '':
        floatlist.append(str2);
    if do_float_conver:
        floatlist = [float(floatlist[row]) for row in range(len(floatlist))]
    return floatlist

### Declare the function of this code
print 
print 'Add a new main programe to the LSS codes...'
print 'Usage: EXE Name-of-the-main-program mpi/nompi'
print

cmdargs = sys.argv
if len(cmdargs) != 3:
	print 'Error: len(cmdargs) = ', len(cmdargs)
	sys.exit()
nowEXEname = cmdargs[1]
mpi = cmdargs[2]
if mpi != 'mpi' and mpi != 'nompi':
	print 'Error: mpi must be mpi or nompi!: ', mpi
	sys.exit()


print
print 'We will add a new EXE named as  ', nowEXEname

### Make a save of the Makefile
output = commands.getoutput('ls')
output = str_to_numbers(output, do_float_conver = False)
MFs = []
MFIs = [-1]
for MF in output:
	if MF[0:13] =='Makefile.SAVE':
		MFs.append(MF)
		MFIs.append(int(MF[13:len(MF)]))
saveindex = max(MFIs) + 1
tmpstr = 'Makefile.SAVE'+str(saveindex)
print 'Copy of old Makefile saved to :', tmpstr
os.system('cp Makefile '+tmpstr)
file1 = tmpstr
os.system('rm Makefile')
file2 = 'Makefile'

### Generate new Makefile...
nowf1 = open(file1, 'r')
nowf2 = open(file2, 'w')

prevstr = ''
while True:
	nowstr = nowf1.readline()
	nowstr0 = nowstr[0:len(nowstr)-1]
	if nowstr == '':
		break
	now_str_array = str_to_numbers(nowstr, do_float_conver = False)
	prev_str_array = str_to_numbers(prevstr, do_float_conver = False)
	if len(now_str_array) == 0:
		nowf2.write(nowstr)
		prevstr = nowstr
		continue
	elif now_str_array[0] == 'OBJS':
		lastEXE = prev_str_array[0]
		print lastEXE
		lastEXEi = int(lastEXE[3:len(lastEXE)])
		nowEXEi = lastEXEi + 1
		nowEXE = 'EXE'+str(nowEXEi)
		nowf2.write(nowEXE + ' = ../bin/LSS_'+nowEXEname+'\n')
		nowf2.write(nowstr)
		prevstr = nowstr; continue
	elif now_str_array[0] == 'all:':
		nowf2.write(nowstr0 + ' $('+nowEXE+')\n')
		prevstr = nowstr; continue
	elif now_str_array[0] == 'nompi:':
		if mpi == 'nompi':
			nowf2.write(nowstr0 + ' $('+nowEXE+')\n')
			prevstr = nowstr; continue
		else:
			nowf2.write(nowstr0 + '         \n')
			prevstr = nowstr; continue
	elif now_str_array[0] == 'AMTB:':
		nowf2.write(nowEXEname+': $('+nowEXE+')\n')
		nowf2.write(nowstr)
		prevstr = nowstr; continue
	elif now_str_array[0] == 'amtb:':
		nowf2.write('$('+nowEXE+'): LSS_cosmo_funs.o LSS_main_'+nowEXEname+'.o\n')
		nowf2.write('\t$(F90C) -o $('+nowEXE+') LSS_tools.o LSS_cosmo_funs.o  LSS_main_'+nowEXEname+'.o\n')
		nowf2.write(nowstr)
		prevstr = nowstr; continue
	elif now_str_array[0] == 'GYPS:':
		nowf2.write('LSS_main_'+nowEXEname+': LSS_cosmo_funs.o\n')
		nowf2.write(nowstr)
		prevstr = nowstr; continue
	elif now_str_array[0] == 'rm':
		nowf2.write(nowstr0+ ' $('+nowEXE+')\n')
		prevstr = nowstr; continue
	else:
		prevstr = nowstr
		nowf2.write(nowstr)
nowf1.close()
nowf2.close()


nowf3 = open('LSS_main_'+nowEXEname+'.f90', 'w')
nowf3.write('program main_'+nowEXEname)
if mpi == 'mpi':
	nowf3.write('\n\nuse mpi')
nowf3.write('\n\nuse LSS_cosmo_funs\n\nimplicit none\n\n\tcharacter(len=char_len) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr\n\tinteger :: i\n\t')
nowf3.write('\n\tprintstr = "Now it is empty!"\n\tif(iargc().le.1) then\n\t\tprint *, printstr\n\t\tstop\n\tendif\n')
nowf3.write('\n\toutputfile = ""\n\tdo i = 1, iargc()\n\t\tif(mod(i,2).eq.0) cycle\n\t\tcall getarg(i,tmpstr1)\n\t\tcall getarg(i+1,tmpstr2)')
nowf3.write('\n\t\tif(trim(adjustl(tmpstr1)).eq."-inputfile") then\n\t\t\tread(tmpstr2,"(A)") inputfile\n\t\telseif(trim(adjustl(tmpstr1)).eq."-outputfile") then')
nowf3.write('\n\t\t\tread(tmpstr2,"(A)") outputfile\n\t\telse\n\t\t\tprint *, "Unkown argument: ", trim(adjustl(tmpstr1))\n\t\t\twrite(*,"(A)") trim(adjustl(printstr))\n\t\t\tstop\n\t\tendif\n\tenddo')
nowf3.write('\n\n\tif(trim(adjustl(outputfile)).eq."") then\n\t\t\n\tendif\n\n\tprint *, \'This is an empty program '+nowEXEname+'!!\'\n\nend program')
nowf3.close()

print '\nNew file created: ', 'LSS_main_'+nowEXEname+'.f90\n\n'

