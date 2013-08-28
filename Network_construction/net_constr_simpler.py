import optparse,sys,os,sets,numpy
import networkx as nx
from Input_processing import inp_proc
from Utils import Blast_cl as bcl
from genome.pangenome import PanGenomer
from Synteny import synteny_simpler as synteny
import time


##########################
# Functions
##########################

def same_border(hit1,hit2,subj,thr):
	''' checks if two contigs are lying on the
	same border of a contig or not.
	'''
	subj_len=len(subj)
	s1,e1=hit1[4],hit1[5]
	s2,e2=hit2[4],hit2[5]
	if s1 <= thr:
		if e2 >= subj_len - thr: return False
		if s2 <=  thr: return True
		return False
	if e1 >= subj_len - thr:
		if e2 >= subj_len - thr: return True
		if s2 <=  thr: return False
		return False

def length(contig):
	return int(contig.split('_&')[-1])

def border_condition(query,subj,qstart,qend,sstart,send,thr,al_length):
	qlen,slen=length(query),length(subj)
	al_perc=al_length/float(qlen)
	if al_perc >= 0.9:
		print '%s maps with al perc %s on %s' %(query,al_perc,subj)
		return True
	if qstart <= thr:
		if ((send >= slen - thr) or (send <= thr)):
			return True
		else: return False
	if qend >= qlen - thr:
		if ((sstart <= thr) or (sstart >= slen - thr)):
			return True
	return False
	
def path(path_,file_):
	path_=check_path(path_)
	return path_+file_

def make_all(inp_dir,out_name,tag=''):
	fastas=[path(inp_dir,file_) for file_ in os.listdir(inp_dir) if file_.endswith(tag)]
	print 'input files are:',' '.join(fastas)
	if len(fastas)==0: raise Exception('no %s file in directory selected') %tag
	out=open(out_name,'w')
	for i in fastas:
		for line in open(i).readlines(): out.write(line.strip()+'\n')
	return out_name,fastas

def do_overlap(hit1,hit2,threshold=0):
	s1start,s1end,s2start,s2end=hit1[4],hit1[5],hit2[4],hit2[5]
	subj1_set=sets.Set(xrange(min(s1start,s1end),max(s1start,s1end)))
	subj2_set=sets.Set(xrange(min(s2start,s2end),max(s2start,s2end)))
	return len(subj1_set.intersection(subj2_set))>threshold

def distance_between(hit1,hit2):
	if do_overlap(hit1,hit2): return 0
	all_distances=[abs(x1 - x2) for x1 in [hit1[4],hit1[5]] for x2 in [hit2[4],hit2[5]]]
	return min(all_distances)
	
def samedraft(hit1,hit2):
	return hit1.split('.fasta')[0]==hit2.split('.fasta')[0]

def check_path(path):
	if path[-1]=='/': return path
	return path+'/'
	
def average(lista):
	return sum(lista)/float(len(lista))

def weighted_average(tuple_list):
	return sum([tuple_[0]*tuple_[1] for tuple_ in tuple_list])/float(sum([tuple_[1] for tuple_ in tuple_list]))

def parse_blast(b_out,length_thresh=200):
	inv_dict={}
	print 'start parsing', b_out
	for line in open(b_out).xreadlines():
		line=line.strip()
		if line[0]=='#': continue
		line=line.split()
		# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
		query,subj,id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score=line[0],line[1],float(line[2]),int(line[3]),int(line[6]),int(line[7]),int(line[8]),float(line[9]),float(line[10]),float(line[11])
		if samedraft(query,subj):continue
		if al_length<length_thresh:continue
		if id_perc<70: continue
		inv_dict[subj]=inv_dict.get(subj,{})
		try: inv_dict[subj][query]
		except: inv_dict[subj][query]=[id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score]
		else:
			if inv_dict[subj][query][1]<al_length: inv_dict[subj][query]=[id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score]
	return inv_dict

def parse_blast_bord(b_out,length_thresh=200,bord_thr=5,id_perc_thr=90):
	''' parse function considering only hits overlapping contig borders.
	Outputs an 'inverse' dictionary: diz[subj][query]=[hit]'''
	inv_dict={}
	print 'start parsing', b_out
	for line in open(b_out).xreadlines():
		if line[0]=='#': continue
		line=line.strip()
		line=line.split()
		# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
		query,subj,id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score=line[0],line[1],float(line[2]),int(line[3]),int(line[6]),int(line[7]),int(line[8]),float(line[9]),float(line[10]),float(line[11])
		if samedraft(query,subj):
			#print 1
			continue
		if al_length<length_thresh:
			#print 2
			print al_length, length_thresh
			continue
		if id_perc<id_perc_thr:
			#print 3
			continue
		#if not border_condition(query,subj,q_start,q_end,s_start,s_end,bord_thr,al_length): continue
		if not border_condition(query,subj,q_start,q_end,s_start,s_end,bord_thr,al_length):
			#print '!!!!!!!!!!!1'
			#print query.subj
			continue
		inv_dict[subj]=inv_dict.get(subj,{})
		try: inv_dict[subj][query]
		except: inv_dict[subj][query]=[id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score]
		else:
			if inv_dict[subj][query][1]<al_length: inv_dict[subj][query]=[id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score]
	return inv_dict

def build_adj_list(diz,draft,thr):
	''' from a inverse dictionary,
	outputs a dictionary diz[contig1]:[contigui connessi].
	Contigs are connected if a same key is shared ''' 
	contigs=draft.contigs
	adj_list={}
	print 'start building graph...'
	for k in diz:
		for v1 in diz[k]:
			for v2 in diz[k]:
				if v1==v2:continue
				hit1=diz[k][v1]
				hit2=diz[k][v2]
				subj=k
				if same_border(hit1,hit2,subj,thr): continue
				adj_list[v1]=adj_list.get(v1,{})
				adj_list[v1][v2]=adj_list[v1].get(v2,{'homology':0,'distance':[]})
				adj_list[v1][v2]['homology']+=1
				adj_list[v1][v2]['distance'].append(distance_between(diz[k][v1],diz[k][v2]))
	return adj_list
	
def build_network(adj_list):
	G=nx.Graph()
	for node1 in adj_list:
		for node2 in adj_list[node1]:
			adj_list[node1][node2]['distance']=average(adj_list[node1][node2]['distance'])
			for kw in adj_list[node1][node2]:
				tupla=(node1,node2,adj_list[node1][node2][kw])
				G.add_weighted_edges_from([tupla],kw)
	return G

def add_to_network(Graph,node1,node2,keyw,value):
	G=Graph
	try: G[node1][node2]
	except:
		print 'no previous links between',node1,node2
		diz={'homology':0,'distance':float('inf'),'synteny':0}
		for k in diz:
			tupla=(node1,node2,diz[k])
			G.add_weighted_edges_from([tupla],k)
		G[node1][node2][keyw]=value
	else:
		tupla=(node1,node2,value)
		G.add_weighted_edges_from([tupla],keyw)
		
def from_gene_to_contig(gene):
	draft=gene.split('_!')[0]
	contig='_'.join(gene.split('_')[:-1])
	return draft,contig
	
def get_correlation(G,label1,label2):
	values=[]
	for edge in G.edges():
		val1,val2=G[edge[0]][edge[1]].get(label1,0),G[edge[0]][edge[1]].get(label2,0)
		values.append( (val1,val2) )
	''' FAI CORRELAZIONE SU VALUES'''
	return values
						
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
			if line[0]=='>':
				try: self.add_contig(cont)
				except: cont=Contig(line)
				else: cont=Contig(line)
			else: cont.seq+=line
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

	usage='''python net_constr_simpler [OPTIONS]
	options are:
		-i (required): the path of the folder containing input draft genomes
		-o (required): the path of the output files
		-n (optional): the number of cpus for analyzing (default=1)
		-t (optional): the threshold of bases for defining border overlap
					(default=5)
	'''
	parser=optparse.OptionParser(usage)
	parser.add_option('-i','--inp',help='give the Input Sequences folder path',dest='inp',default='../Input/')
	parser.add_option('-o','--out',help='give the Output folder path',dest='out',default='../Output/')
	parser.add_option('-n','--ncpu',help='set the number of cpus to use',dest='ncpu',default=1)
	parser.add_option('-t','--bound',help='# of nucleotide threshold for contig boundaries during BLAST parsing',dest='bord_thr',default=5)
	parser.add_option('-m','--hit',help='minimum hit length',dest='hit_thr',default=100)
	parser.add_option('-p','--idperc',help='id-% threshold',dest='id_perc_thr',default=90)
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

	# assign options to variables
	input_folder=opts.inp
	output_folder=opts.out
	bord_thr=opts.bord_thr
	hit_thr=opts.hit_thr
	id_perc_thr=opts.id_perc_thr
	############################
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
	for fasta in fastas: # per ogni file di input...
		#Instazia oggetto Draft + BLASTall
		out=fasta[:-6]+'_blast_out'
		draft=Draft(fasta)
		draft.build()
		bcl.launch_BLAST(fasta,all_fasta,out)
		#
		#Parse di BLASTall e costruzione grafo
		#parsed=parse_blast(out)
		parsed=parse_blast_bord(out,hit_thr,bord_thr,id_perc_thr)
		adj_list=build_adj_list(parsed,draft,bord_thr)
		network=build_network(adj_list)
		#
		#Scrivi file gexf e appendi alla lista dei network
		graph_name=output_folder+draft.graph_name
		nx.write_gexf(network,graph_name)
		networks[draft.name]=graph_name
	#predizione proteine e BBH seriale
	prots=synteny.predictor(fastas)
	pg=PanGenomer(prots,opts.ncpu)
	print 'now finding orthologs... it may take some time..'
	pg.run()
	orthologs=pg.orthologs
	relations={}
	#organizzazione delle relazioni: contigui di uno stesso genoma che condividono ortologhi con contigui di altri genomi sono raggruppati
	for k in orthologs:
		for v1 in orthologs[k]:
			for v2 in orthologs[k]:
				if v1==v2:continue
				draft1,contig1=from_gene_to_contig(v1)
				draft2,contig2=from_gene_to_contig(v2)
				relations[draft1],relations[draft2]=relations.get(draft1,{}),relations.get(draft2,{})
				relations[draft1][contig2]=relations[draft1].get(contig2,[])
				relations[draft2][contig1]=relations[draft2].get(contig1,[])
				if contig1 not in relations[draft1][contig2]:
					relations[draft1][contig2].append(contig1)
				if contig2 not in relations[draft2][contig1]:
					relations[draft2][contig1].append(contig2)

	#aggiungi relazioni synteny in network e salva1
	for network in networks:
		G=nx.read_gexf(networks[network])
		for lista in relations[network]:
			for contig1 in relations[network][lista]:
				for contig2 in relations[network][lista]:
					if contig1==contig2: continue
					try: G[contig1][contig2]					
					except: 
						diz={'homology':0,'distance':float('inf'),'synteny':1}
						for k in diz:
							tupla=(contig1,contig2,diz[k])
							G.add_weighted_edges_from([tupla],k)
					
					else: G[contig1][contig2]['synteny']=G[contig1][contig2].get('synteny',0) + 1
		nx.write_gexf(G,networks[network])
	#print all_fasta
	#print fastas
