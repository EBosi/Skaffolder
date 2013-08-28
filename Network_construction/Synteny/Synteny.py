
#Gestire package in maniera =\= da sys.path.append(package)
from Utils import Blast_cl as bcl
from Homology import BLAST_parser as blp

##############################


def from_gene_to_contig(gene):
	draft=gene.split('_!')[0]
	contig='_'.join(gene.split('_')[:-1])
	return draft,contig
			

##############################


class synteny(object):

	def __init__(self,inp_folder):
		self.inpf=inp_folder

	def find_orthologs(self,nets):
		self.Blast=bcl.BLAST_make(self.inpf)
		self.Blast.BLAST_genes()
		self.parser=blp.parser(self.Blast.to_parse())
		return self.parser.BBH(nets)

	def shared_orthologs(self,nets):
		inv_bidirect=self.find_orthologs(nets)
		for k in inv_bidirect:
			for v1 in inv_bidirect[k]:
				for v2 in inv_bidirect[k]:
					
					if v1==v2: continue
					if v1.split('_!')[0]==v2.split('_!')[0]:
						nets[v1.split('_!')[0]].add_synteny(v1,v2)
						

	def add_edges_to_net(self,nets):
		for net in nets:
			nets[net].add_links()
		

	def run(self,nets):
		self.shared_orthologs(nets)
		self.add_edges_to_net(nets)
