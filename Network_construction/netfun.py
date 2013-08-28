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


def save_orth(prots=None,ncpu=1,out_file=None,inp_orth=None):
	if inp_orth != None:
		print 'using following file as orthologs', inp_orth
		inp_orth=pick.load(open(inp_orth))
		return inp_orth
	pg=PanGenomer(prots,ncpu)
	print 'now finding orthologs... it may take some time..'
	pg.run()
	orthologs=pg.orthologs
	if out_file != None:
		dumper=open(out_file,'w')
		pick.dump(orthologs, dumper)
	return orthologs

def find_synteny_relations(ort):
	relations={}
	#organizzazione delle relazioni: contigui di uno stesso genoma che condividono ortologhi con contigui di altri genomi sono raggruppati
	for k in ort: # gruppi di ortologhi
		for v1 in ort[k]: # proteina 1 del gruppo k
			for v2 in ort[k]: # proteina 2 del gruppo k
				if v1==v2: continue # skip se prot1 = prot2
				draft1,contig1,order1=from_gene_to_contig(v1) # trova genoma, contiguo, e posizione sul contiguo  
				draft2,contig2,order2=from_gene_to_contig(v2) # trova genoma ecc ecc
				relations[draft1],relations[draft2]=relations.get(draft1,{}),relations.get(draft2,{}) # inizializza dizionario per ogni draft
				relations[draft1][contig2]=relations[draft1].get(contig2,[]) # ogni draft ha una chiave per ogni contiguo dell' altro draft
				relations[draft2][contig1]=relations[draft2].get(contig1,[]) # idd
				if contig1 not in relations[draft1][contig2]: 
					relations[draft1][contig2].append(tuple((contig1,order2))) # metti in una lista tutti i contigui che mappano su un contiguo dell'altro genoma
				if contig2 not in relations[draft2][contig1]:
					relations[draft2][contig1].append(tuple((contig2,order1))) # idd
	return relations
	
def add_synteny_to_G(contig1,contig2,G):
	if contig1==contig2: return
	try: G[contig1][contig2]					
	except: G.add_edge(contig1,contig2,homology=0,synteny=1,overlap=0,contain=0)
	else: G[contig1][contig2]['synteny'] += 1
	

def add_singletons(G,contigs):
	for node in contigs:
		if node in G.nodes():continue
		G.add_node(node)
	
def hit_to_node(hit):
	name=hit[0]
	start,end=min(hit[6:8]),max(hit[6:8])
	if hit[6]==start: strand=1
	else: strand=-1
	return name,start,end,strand

def add_to_T(T,hit):
	'''
	add hits to a list-of-tuples pickled object
	the hits will be eventually be sorted using:
	list.sort(key=lambda x: x[1])
	'''
	if os.path.exists(T):
		tree=pick.load(open(T))
		node=hit_to_node(hit)
		tree.append(node)
		out=open(T,'w')
		pick.dump(tree,out)
		out.close()		
		return 		
	else:
		#name,start,end=hit_to_node(hit)
		node=hit_to_node(hit)
		tree=[node]
		out=open(T,'w')
		pick.dump(tree,out)
		out.close()
		return 

def build_Tree(tree_file):
	'''
	Create a Graph object (a tree) from a list.
	The list is sorted and an edge between the 
	i-th and i-th+1 nodes is added to T.
	'''
	T=nx.Graph()
	tree=pick.load(open(tree_file))
	tree.sort(key=lambda x: x[1])
	for i in range(len(tree)-1):
		node1,node2=tree[i][0],tree[i+1][0]
		orient=tuple((tree[i][3],tree[i+1][3]))
		T.add_edge(node1,node2,overlap=False,orientation=orient,contain=False)
		if tree[i+1][1] < tree[i][2]:
			T[node1][node2]['overlap']=True
		if tree[i+1][2] < tree[i][2]:
			T[node1][node2]['contain']=True
	return T

def add_T_to_G(G,T):
	'''
	Add a tree to a Graph
	'''
	T=build_Tree(T)
	for edge in T.edges():
		n1,n2=edge[0],edge[1]
		if is_edge_in(G,edge):
			G[n1][n2]['homology']+=1
			if T[n1][n2]['overlap']: G[n1][n2]['overlap']+=1
			if T[n1][n2]['contain']: G[n1][n2]['contain']+=1
		else:
			G.add_edge(n1,n2,homology=1,synteny=0,overlap=0,contain=0)
			if T[n1][n2]['overlap']: G[n1][n2]['overlap']+=1
			if T[n1][n2]['contain']: G[n1][n2]['contain']+=1
			
			
def inverse_edge(edge):
	n1,n2=edge[0],edge[1]
	return tuple((n2,n1))
	
def is_edge_in(G,edge):
	return ( (edge in G.edges()) or 
			 (inverse_edge(edge) in G.edges()) )

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
	
def distance_between_nodes(node1,hit2,G):
	hit1=[0,0,0,0,G.node[node1]['start'],G.node[node1]['end']]
	return distance_between(hit1,hit2)
	
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
	

def parse_blast_tree(b_out,length_thresh=200,bord_thr=5,id_perc_thr=90):
	'''
	function that parse a blast output, applies threshold, yields only
	border hits, builds intermediate trees and build a final Graph
	'''
	temp_dir='temp_dir_1' # create a temp directory
	while os.path.exists(temp_dir): # if it exists, create a different one 
		temp_dir= temp_dir[:-1]+str(( int(temp_dir[-1]) +1 ))
	temp_dir+='/'
	os.makedirs(temp_dir)
	print 'created temp directory', temp_dir
	print 'start parsing', b_out
	for line in open(b_out).xreadlines():
		if line[0]=='#': continue
		line=line.strip()
		line=line.split()
		query,subj,id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score=line[0],line[1],float(line[2]),int(line[3]),int(line[6]),int(line[7]),int(line[8]),float(line[9]),float(line[10]),float(line[11])
		# blast thresholds applied...
		if samedraft(query,subj): continue
		if al_length < length_thresh: continue
		if id_perc<id_perc_thr: continue
		if not border_condition(query,subj,q_start,q_end,s_start,s_end,bord_thr,al_length): continue
		# if is ok...
		T=temp_dir+subj # create pickle obj
		hit=[query,subj,id_perc,al_length,q_start,q_end,s_start,s_end,evalue,score] # keep hit
		add_to_T(T,hit) # pickle hit into T
	### Finished blast_out parsing, now build a Graph
	G=nx.Graph() # create a Graph
	for file_ in os.listdir(temp_dir): # for every pickled T:
		T=temp_dir+file_ 
		add_T_to_G(G,T) # open T and add it to G
	return G,temp_dir # return Graph

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
	order=gene.split('_')[-1]
	return draft,contig,order
	
def get_correlation(G,label1,label2):
	values=[]
	for edge in G.edges():
		val1,val2=G[edge[0]][edge[1]].get(label1,0),G[edge[0]][edge[1]].get(label2,0)
		values.append( (val1,val2) )
	''' FAI CORRELAZIONE SU VALUES'''
	return values

