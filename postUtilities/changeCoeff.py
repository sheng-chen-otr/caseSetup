import os
import sys
import glob


dirlist = glob.glob(os.path.join('postProcessing','*'))
for dir in dirlist:
	if os.path.exists(os.path.join(dir,'0','coefficient_0.dat')):
		print(dir)
		os.system('mv %s %s' % (os.path.join(dir,'0','coefficient_0.dat'),os.path.join(dir,'0','coefficient.dat')))
