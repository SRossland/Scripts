#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
import sys

if sys.argv[1] == '':
	print('sys.exit()')
	sys.exit()
#	print('continue')
#	continue
try:
	print('Hello, world!'+str(sys.argv[1][0:11])+' the length of sys.argv[1] is: '+str(len(sys.argv[1])))
	if sys.argv[1] == '50311002006.obs':
		open('thisfiledoesnotexist.dat','r+')

except Exception as e:
	print('File does not exist')
	print('This was a test of exception handling')
	print(str(e))
	sys.exit
