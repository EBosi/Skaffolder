import networkx as nx
import os,sys,optparse

###########################

__doc__='''
Takes as input a directory with .gexf (default) files, and removes noisy contigs from the network.
Other format can be specified
Noisy contigs are (more formally) connector nodes, 
which are nodes with extremely high betweenness-centrality value.
The presence of such nodes may disturb community analysis, connecting distinct communities,
producing community merging, leading to poor accuracy in scaffold predictions.
This module compute a betweenness-centrality (bc) distribution for a network,
and a mean and a standard deviation values associated with the distribution,
then removes those nodes with bc value > mean + N*st-deviation.
N is an integer that determines nodes removal: lower values of N will 
lead to more nodes being cut from the network.'''

###########################

def get_graph(tag,file_):
	diz={'.gexf':nx.read_gexf,'.gml':nx.read_gml,'.yaml':nx.read_yaml}
	if tag not in diz.keys():
		print '''
		ERROR! network extension not recognized!
		please retry with .gexf, .gml or .yaml'''
		sys.exit()
	#print 'reading file', file_, 'with tag', tag
	G=diz[tag](file_)
	return G

def mean_stdev(vector):
	from math import sqrt
	if vector.__class__==dict: #dict like TAG:VALUE
		vector=[vector[k] for k in vector] #transform it into vec type
	sum_,sum2,N=0,0,len(vector)
	for num in vector:
		sum_, sum2 = sum_+num, sum2+(num**2)
	mean= sum_/float(N)
	sd= sqrt(N*sum2 - (sum_**2))/float(N)
	return mean,sd

def graph_stats(graph):
	print graph
	v=nx.betweenness_centrality(graph,normalized=True)
	m,sd=mean_stdev(v)
	return m,sd
	
def find_outliers(diz,m,s,fold=1): # find elements w/ values < m+fold*s
	bad=[node for node,val in diz.iteritems()
					if diz[node] > m + (s*fold)]
	return bad

def run_batch(net_folder,tag,fold_times,out):
	if out==None: out=net_folder
	output_list=[]
	for file_ in os.listdir(net_folder):
		if file_.endswith(tag):
			out_file=out+file_[:-len(tag)]+'_cleaned'+tag
			net=net_folder+file_
			run_single(net,tag,fold_times,out_file)
			output_list.append(out_file)
	return output_list 

def run_single(net,tag,fold_times,out_graph):
	if out_graph==None:
		#print 'it is none'
		out_graph=net+'_cleaned'
	graph=get_graph(tag,net)
	mean,stdev=graph_stats(graph)
	outliers=find_outliers(nx.betweenness_centrality(graph,normalized=True),mean,stdev,fold_times)
	for n in outliers:
		graph.remove_node(n)
	nx.write_gexf(graph,out_graph)
	return out_graph
	
def helper(usage):
	print usage
	sys.exit()

###########################
"""
usage='''python noiser.py [OPTIONS]
options are:
	-i (required): the input network (or network folder if "-b" is used)
	-o (optional): tho output network (or network folder if "-b " is used) (default=NET_NAME_cleaned.TAG)
	-f (optional): the network file format (DEFAULT=.gexf)
	-n (optional): the stringency of the removal (DEFAULT=1)
	-b (optional): flag to set batch analysis (folder of networks as input file) (DEFAULT=FALSE)

NOTE THAT if you set the flag "-b", -i and -o (if used) must be folders,
otherwise the script will FAIL!
'''
parser=optparse.OptionParser(usage)
parser.add_option('-i','--inp',help='give the Input file',dest='inp',default=None)
parser.add_option('-o','--out',help='give the Output file',dest='out',default=None)
parser.add_option('-f','--format',help='set the format of input file',dest='format',default='.gexf')
parser.add_option('-n','--N',help='set the value of N',dest='fold_times',default=1)
parser.add_option('-b','--batch',help='folder as input file',action='store_true',dest='batch',default=False)
#parser.add_option('-h','--help',help='ask for help',dest='help',default=True)
opts,args=parser.parse_args()
if not opts.inp:
	print 'CAZZONE!!!! mandatory options are found to be missing... try again please or STUDY PYTHON PROPERLY'
	helper(usage)
if not os.path.exists(opts.inp):
	print 'Input file was not found... check your files or permissions...'
	helper(usage)
"""
###########################

