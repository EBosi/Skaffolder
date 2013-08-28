import os


###################

def from_gene_to_contig(gene):
	draft=gene.split('_!')[0]
	contig='_'.join(gene.split('_')[:-1])
	return draft,contig

def predictor(draft_list):
	protein_list=[]
	for draft in draft_list:
		protein_list.append(predict_prots(draft))
	return protein_list
			
def predict_prots(draft):
	line ='prodigal -i %s -o %s -a %s -q' %(
					draft,draft+'_prodigal',draft+'_prodigal.prots')
	print 'running prodigal with following command line...'
	print 'cm used:',line
	os.system(line)
	return draft+'_prodigal.prots'


###################

class BBH(object):
	
	def __init__(self,drafts):
		self.drafts=drafts
		self.orthologs={}
		self.prots=()
		
	def predict_prots(self):
		self.prots = predictor(self.drafts)

	def join_prots(self):
		all_prot_path=os.path.split(self.prots[0])+'/all_prots'
		all_prot=open(all_prot_path)
		for prot in self.prots:
			for line in open(prot):
				all_prot.write(line)
	
	#def DB_creation(self):
		
