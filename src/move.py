import commands;
for nowstr in commands.getoutput('ls *').split():
	if nowstr[0:3] == 'ap_':
		print nowstr
		newstr = 'LSS_'+nowstr[3:len(nowstr)]
		print newstr
		commands.getoutput('mv '+nowstr+' '+newstr)
