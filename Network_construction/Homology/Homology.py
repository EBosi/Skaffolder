import networkx as nx
import BLAST_parser

########################

def w_average(tuple_list):
	summ=0.0
	weights=0.0
	for i in tuple_list:
		summ+=i[0]*i[1]
		weights+=i[1]
	return summ/weights


def do_overlap(hit1,hit2,threshold=0):
	return len(hit1.subj_set.intersection(hit2.subj_set))>threshold


def measure_distance(hit1,hit2):
	if do_overlap(hit1,hit2): return 0
	return min( abs(hit1.sstart - hit2.sstart),
					abs(hit1.sstart - hit2.send), 
					abs(hit1.send - hit2.sstart),
					abs(hit1.send - hit2.send) )


def samedraft(hit1,hit2):
	return hit1.split('.fasta')[0]==hit2.split('.fasta')[0]
########################


class Homology_parsing(BLAST_parser.parser):
	
	def net_construction(self,Nets):
		self.parse(Nets)
		print self.inverse_hits
		for inv in self.inverse_hits:
			print self.inverse_hits[inv]
			if len(self.inverse_hits[inv])<2: continue
			for hit1 in self.inverse_hits[inv]:
				hit1=self.inverse_hits[inv][hit1]
				print hit1
				for hit2 in self.inverse_hits[inv]:
					hit2=self.inverse_hits[inv][hit2]
					print hit2
					if hit1==hit2: continue
					if samedraft(hit1.query,hit2.query):
						print hit1, hit1.query,'!!!!'
						coverage=min(hit1.coverage,hit2.coverage)
						distance=(measure_distance(hit1,hit2),coverage)
						Nets[hit1.draft].add(hit1.query,hit2.query,coverage,distance)
						

	def net_final(self,Nets):
		for net in Nets:
			for i in Nets[net].adj_list:
				for j in Nets[net].adj_list[i]:
					Nets[net].adj_list[i][j]['distance']=w_average(Nets[net][i][j]['distance'])
					for kw in Nets[net].adj_list[i][j]:
						tupla=(k,w,Nets[net].adj_list[i][j][kw])
						Nets[net].G.add_weighted_edges_from([tupla],kw)

	def run(self,Nets):
		self.Homology_parsing.net_construction(Nets)						
		self.Homology_parsing.net_final(Nets)
		
	def parse_only(self,name):
		self.parse


class Net(object):

	def __init__(self,draft_obj):
		self.draft=draft_obj
		self.fasta=self.draft.path+self.draft.name
		self.name=self.draft.name.split('.fasta')[0]
		self.G=nx.Graph()
		self.adj_list=self.draft.adjacency_list
		self.blast_out=self.name+'blast_out'
		for i in self.adj_list:
			for j in self.adj_list:
				if i==j:continue
				self.adj_list[i][j]={'homology':0,'coverage':0.0,'distance':[],'synteny':0}
		
	def add(self,hit1,hit2,coverage,distance):
		self.adj_list[hit1][hit2]['homology']+=1
		self.adj_list[hit1][hit2]['coverage']=max(coverage,self.adj_list[hit1][hit2]['coverage'])
		self.adj_list[hit1][hit2]['distance'].append(distance)

	def add_synteny(self,hit1,hit2):
		self.adj_list[hit1][hit2]['synteny']+=1

	def add_links(self):
		diz=self.adj_list
		for k in diz:
			for w in diz[k]:
				for kw in diz[k][w]:
					tupla=(k,w,diz[k][w][kw])
					self.G.add_weighted_edges_from([tupla],kw)
		
		
	def export(self,path):
		nx.write_gexf(self.G,path)

#### NET STATISTICS ####

	@property
	def singletons(self):
		return [node for node in self.G.nodes() if len(self.G[node])==0]

	@property
	def not_singletons(self):
		return [node for node in self.G.nodes() if len(self.G[node])>0]

	def statistics(self):
		self.n_nodes=len(self.G.nodes())
		self.n_edges=len(self.G.edges())
		print 'nodes: %s edges: %s' %(len(self.not_singletons),float(self.n_edges + 1))
		self.noe=len(self.not_singletons)/float(self.n_edges + 1)

	def degrees(self):
		self.degrees={}
		for k,v in nx.degree(self.G).iteritems(): self.degrees[v]=self.degrees.get(v,0)+1
		self.K_mode=mode(degrees)
		self.K_av=average(degrees)

	def conn_comps(self):
		self.components={}
		for n in nx.connected_components(self.G):
			components[len(n)]=components.get(len(n),0)+1
		self.n_comp=len(nx.connected_components(self.G))
		self.comp_mode=mode(components)
		self.comp_average=average(components)

	def cyclez(self):
		self.cycles={}
		for k,v in nx.cycle_basis(self.G).iteritems(): self.cycles[v]=self.cycles.get(v,0)+1
		self.n_cycles=len(nx.cycle_basis(self.G))
		self.cycles_mode=mode(cycles)
		self.cycles_average=average(cycles)

