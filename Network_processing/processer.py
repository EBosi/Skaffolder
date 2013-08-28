from Noise_removal import noiser
from Network_clustering.abc_maker import abc as abcm
import Network_clustering.clustering as clustering
import os,sys,optparse
import networkx as nx

####################################################

__doc__='''
this module(s) performs computation on scaffolding network, removing 
connector nodes, then clustering the network using mcl. obviously mcl is
supposed to be installed on the machine you are working on, and to be 
in /usr/bin/. Please be sure that this happen otherwise clustering wont 
succeed.
You can download mcl (or read about it) at http://www.micans.org/mcl/. 
'''


####################################################

def mcl(inp,I=2.0,out=None):
	I=float(I)
	if out==None: line='mcl %s --abc -I %s --d' %(inp,I)
	else: line='mcl %s --abc -I %s -o %s' %(inp,I,out)
	print line 
	os.system(line)
	return

def path(path_,file_):
	check_path(path_)
	return path_+file_

def batcher(inp):
	if not os.path.isdir(opts.inp):
		print 'you selected batch option... please enter a directory as input'
		sys.exit()
	else:
		inps=[file_ for file_ in os.listdir(opts.inp) if file_.endswith(opts.format)]
		for inp in inps:
			inp=opts.inp+inp
			try: out=opts.out+inp
			except: out=opts.inp+inp
			inp=Input_file(inp,opts.abc,opts.format,opts.weight,out,opts.N,opts.clean)
			inp.run()

def helper(usage):
	if __name__=='__main__':
		print usage
		sys.exit()


####################################################


usage='''python processer.py [OPTIONS]
options are:
	-i (required) -> input: the input network (or folder)
	-w (required) -> weight: the attribute name for clustering 
	(can also not be required! read below)
	-o (optional) -> output: the output file (or folder) 
	(DEFAULT=overwrite input)
	-c (optional) -> only cleaning: skip the clustering phase
	-f (optional) -> format: the network file format (DEFAULT=.gexf)
	-b (optional) -> batch: use this option if you have multiple files
	-a (optional) -> abc: us this option if you want to provide  abc files
	(skip cleaning)
	-h (optional) -> help: print this usage and exit
	-n (optional) -> N: how conservative is network cleaning 
	(low N = more cleaning) (DEFAULT=1)
	
NOTE that:
	- if you set -a flag, you must use as input abc file(s). The file
	  format will be set to ".abc"
	  The network cleaning part will be skipped and -w argument will be not needed
	- if you set -b flag, you must use a folder as -i and -o (if used)
	  argument(s)
'''
parser=optparse.OptionParser(usage)
parser.add_option('-i','--inp',help='Input file',dest='inp',default=None)
parser.add_option('-w','--weight',help='Attribute name',dest='weight',default=None)
parser.add_option('-o','--out',help='Output file',dest='out',default=None)
parser.add_option('-f','--format',help='Input file format',dest='format',default='.gexf')
parser.add_option('-a','--abc',help='abc file(s) as input',action='store_true',dest='abc',default=False)
parser.add_option('-b','--batch',help='folder as input file',action='store_true',dest='batch',default=False)
parser.add_option('-n',help='cleaning stringency',dest='N',default=1)
parser.add_option('-c',help='only cleaning',action='store_true',dest='clean',default=False)
opts,args=parser.parse_args()

if not opts.weight:
	print '''
	''', opts.weight
	if not opts.abc: opts.inp=False
	
if not opts.inp:
	print 'GIF OUTPUTZ'
	print 'Mandatory options are found to be missing... try again please '
	helper(usage)
	
if not os.path.exists(opts.inp):
	print 'input not found'
	print 'Input file was not found... check your files or permissions...'
	helper(usage)

if opts.abc: opts.format='.abc'
	
#if not os.path.exists(opts.out):
#	opts.help_=True
#	print 'Output file was not found... check your files or permissions...'



####################################################

class Input_file(object):
	''' This class models input files of processer.py, which can be:
	- (cleaned) networks
	- abc files (see mcl documentation)
	
Basically the init funct discriminates the type of input file and its 
arguments are user-defined in the script options.

inp = input file path
  
	'''
	def __init__(self,inp,abc,tag,weight,out,N,condition):
		self.inp=inp
		self.abc=abc
		self.tag=tag
		self.weight=weight
		self.out=out
		self.inps=None
		self.cond=condition
		self.N=N
		
	def check_stop(self,condition):
		if condition:
			print 'exiting the script ("-c" option selected)...'
			return True
		else: return False

	def abc_make(self,abc,name_abc):
		if name_abc==None:
			name_abc=self.inp+'.abc'
		if abc: return self.inp
		print 'making abc file...'
		print 'network file %s will be used to make abc file %s using this weight: %s' %(self.inp,name_abc,self.weight)
		abc=abcm(self.inp)
		abc_out=name_abc+'.abc'
		abc_file=abc.write_in(name_abc,self.weight)
		return abc_file
	#def clean_batch(self,inp,tag,N,out):

	def clean(self,inp,tag,N,out):
		if self.abc: return inp
		#print 'asda', inp
		return noiser.run_single(inp,tag,N,out)
	
	def launch_mcl(self,abc_file,I=2.0,out=None):
		clustering.run(abc_file,I,out)
		return out
	
	def run(self):
		self.inp=self.clean(self.inp,self.tag,self.N,self.out)
		if self.check_stop(self.cond): return
		self.inp=self.abc_make(self.abc,None)
		self.abc=True
		self.launch_mcl(self.inp,30.0,None)		

if __name__=='__main__':

	if opts.batch:
		print '\nbatch mode selected...'
		batcher(opts.inp)
		print 'clustering ended'
		print 'Output can be found in: '
		if opts.out!=None: print opts.out
		else: print opts.inp
		sys.exit()

	inp=Input_file(opts.inp,opts.abc,opts.format,opts.weight,opts.out,opts.N,opts.clean)
	inp.run()
	
