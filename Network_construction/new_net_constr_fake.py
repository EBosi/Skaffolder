import optparse,sys,os,sets,numpy
import cPickle as pick
import networkx as nx
from Input_processing import inp_proc
from Utils import Blast_cl as bcl
from genome.pangenome import PanGenomer
from Synteny import synteny_simpler as synteny
from shutil import rmtree
import time


##########################
# Functions
##########################

from netfun import *
						
##########################
# Classes
##########################

class Draft(object):

	def __init__(self,fasta):
		self.fasta=fasta
		self.name=self.fasta.split('/')[-1]
		self.path='/'.join(self.fasta.split('/')[:-1])+'/'
		self.graph_name=self.name+'.gexf'
		self.contigs={}
		self.adjacency_list={}

	
	def add_contig(self,contig):
		self.contigs[contig.name[1:]]=contig

	def build(self):
		for line in open(path(self.path,self.name)).readlines():
			line=line.strip()
			if len(line)==0: continue
			if line[0]=='>':
				try: self.add_contig(cont)
				except: pass
				cont=Contig(line)
			else:
				try: cont.seq+=line
				except: print line
		self.add_contig(cont)


class Contig(object):

	def __init__(self,name,seq='',tag=''):
		self.name=name.split()[0]
		self.draft=self.name.split('_!')[0]+tag
		self.seq=seq

	def coverage(self,aligned):
		return float(aligned)/float(self.length)

	def __len__(self):
		return len(self.seq)

	def samedraft(self,cntg):
		return self.draft==cntg.draft


##########################
# Options
##########################

if __name__=='__main__':

	usage='''python new_net_constr_fake [OPTIONS]

	some common usage are: 
		python new_net_constr_fake.py -t 10 -m 50 -p 80
		python new_net_constr_fake.py -y -l TARGET
		python new_net_constr_fake.py -z -b ORTH_FILES 
	'''
	parser=optparse.OptionParser(usage)
	parser.add_option('-i','--inp',help='give the Input Sequences folder path',dest='inp',default='../Input/')
	parser.add_option('-o','--out',help='give the Output folder path',dest='out',default='../Output/')
	parser.add_option('-n','--ncpu',help='set the number of cpus to use',type="int",dest='ncpu',default=1)
	parser.add_option('-t','--bound',help='# of nucleotide threshold for contig boundaries during BLAST parsing',type='int',dest='bord_thr',default=5)
	parser.add_option('-m','--hit',help='minimum hit length',dest='hit_thr',type='int',default=100)
	parser.add_option('-p','--idperc',help='id-% threshold',dest='id_perc_thr',type='float',default=90.0)
	parser.add_option('-r','--removetemp', help='auto remove temporary folders',dest='auto_remove',default=True)
	parser.add_option('-s','--serialize', help='serialize orthologs',dest='serialize',default=None)
	parser.add_option('-l','--limitedto',help='scaffold only the following genomes',dest='target', default='all')
	parser.add_option('-b','--skipbbh',help='skip bbh by presenting orth serialized obj',dest='orth', default=None)
	parser.add_option('-y','--onlyhomology',help='only_homology',action='store_true',dest='homo', default=False)
	parser.add_option('-z','--skiphomology',help='skip homology part: give as input networks',action='store_true',dest='no_homo', default=False)	
	#parser.add_option('-h','--help',help='ask for help',dest='help',default=True)
	opts,args=parser.parse_args()
	opts.help=False
	if (not opts.inp) or (not opts.out):
		opts.help=True
		print 'mandatory options are found to be missing... try again please'
	if opts.help:
		print usage
		sys.exit()


##########################
# Main
##########################

if __name__=='__main__':

	# OPTIONS ASSIGNMENT
	input_folder=opts.inp
	output_folder=opts.out
	bord_thr=opts.bord_thr
	hit_thr=opts.hit_thr
	id_perc_thr=opts.id_perc_thr
	auto_remove=opts.auto_remove
	cpus=opts.ncpu
	out_serialized=opts.serialize
	arget_genomes=opts.target
	inp_orth=opts.orth
	only_homology=opts.homo
	only_synt=opts.no_homo
	############################



	###########################
	# Only synteny
	##

	if only_synt==True:
		### Input processing!
		#
		input_folder,output_folder=check_path(input_folder),check_path(output_folder)
		networks={}
		nets=[f for f in os.listdir(input_folder) if f.endswith('.gexf')]
		for net in nets:
			k=net[:-5]
			networks[k]=input_folder + net

		### Synteny part
		#	
		orthologs=save_orth(inp_orth=inp_orth) 
		relations=find_synteny_relations(orthologs)
		#aggiungi relazioni synteny in network e salva1e
		for network in networks:
			G=nx.read_gexf(networks[network])
			print network
			for lista in relations[network]:
				lst=relations[network][lista]
				lst.sort(key= lambda x: x[1])
				for i in range(len(lst)-1):
					contig1,contig2=lst[i][0],lst[i+1][0]
					add_synteny_to_G(contig1,contig2,G)
			graph_name=networks[network][:-5]+'_with_synteny.gexf'
			nx.write_gexf(G,graph_name)

		### Finish
		#
		print 'task completed!'
		sys.exit()
	
	### Input processing!
	#
	input_folder,output_folder=check_path(input_folder),check_path(output_folder)
	# clean input/output folders
	inp_proc.clean_folder(input_folder,'') # keep only file ending with TAG from input_folder
	inp_proc.clean_folder(output_folder,'') # keep only file ending with TAG from output_folder
	inp_proc.tag_all(input_folder, '.inp') # add .inp tag to all files in input_folder
	inp_proc.inp_rename_len(input_folder,'.inp') # rename all .inp files in input
	# make all_fasta file from .inp files
	all_fasta,fastas=make_all(input_folder,input_folder+'all_fasta','.inp') 
	print 'creating allfasta db file...', all_fasta
	#time.sleep(5)
	bcl.format_db(all_fasta) # fai il db di allfasta per il BLAST 
	networks={}
	temps=[]
	
	### Homology part
	# 
	for fasta in fastas:
		# Instanciate Draft object + BLASTall
		out=fasta[:-6]+'_blast_out'
		draft=Draft(fasta)
		draft.build()
		bcl.launch_BLAST(fasta,all_fasta,out)
		#
		# Parse BLASTall and graph construction
		# parsed=parse_blast(out)
		network,temp=parse_blast_tree(out,hit_thr,bord_thr,id_perc_thr) # create a graph
		temps.append(temp)
		add_singletons(network,draft.contigs)
		# write gexf file and append to network list
		graph_name=output_folder+draft.graph_name
		nx.write_gexf(network,graph_name)
		networks[draft.name]=graph_name
	# remove tmp folders 
	if auto_remove:
		for temp in temps: rmtree(temp)
	if only_homology == True:
		print 'only homology flag set True (-y option)...'
		print 'task completed!'
		sys.exit()
	# protein prediction + synteny
	prots=None
	if inp_orth == None:
		prots=synteny.predictor(fastas)

	### Synteny part
	#	
	orthologs=save_orth(prots=prots,ncpu=cpus,inp_orth=inp_orth) 
	relations=find_synteny_relations(orthologs)
	# add synteny relations in network and save it
	for network in networks:
		G=nx.read_gexf(networks[network])
		print network
		for lista in relations[network]:
			lst=relations[network][lista]
			lst.sort(key= lambda x: x[1])
			for i in range(len(lst)-1):
				contig1,contig2=lst[i][0],lst[i+1][0]
				add_synteny_to_G(contig1,contig2,G)
		graph_name=networks[network][:-5]+'_with_synteny.gexf'
		nx.write_gexf(G,graph_name)

	### Finish
	#
	print 'task completed!'
	sys.exit()
