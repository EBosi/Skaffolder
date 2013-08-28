import os,sys,optparse

################################

__doc__='''

This module use mcl (Markov CLustering), a program of S. Van Dongen 
which takes as input a network in abc format, and outputs clusters.
abc format is:
LABEL1 TAB LABEL2 TAB WEIGHT

Usage is:

python clustering.py [FILE] [OPTIONS] -> executes mcl
python clustering.py 				  -> shows this documentation
python clustering.py -h				  -> shows mcl options
'''

################################

"""if len(sys.argv)==1:
	print __doc__
	sys.exit()

for arg in sys.argv[1:]:
	to_print=''
	if arg[0]=='-':
		print to_print
		to_print=''
		to_print+=arg
	else: to_print+=arg
print to_print

line='mcl '+ ' '.join(sys.argv[1:])
print line
os.system(line)"""

def run(inp,I=2.0,out=None):
	I=float(I)
	if out==None: line='mcl %s --abc -I %s --d' %(inp,I)
	else: line='mcl %s --abc -I %s -o %s' %(inp,I,out)
	print line
	os.system(line)
