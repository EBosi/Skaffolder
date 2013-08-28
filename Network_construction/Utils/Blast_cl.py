import os
from Bio.Blast.Applications import NcbiblastnCommandline


##############################

def path(path_,file_):
	if path_[-1]!='/': path+='/'
	return path_+file_

def format_db(path, p='F'):
	os.system('formatdb -i%s -p%s -oT ' %(path,p))


def launch_BLAST(qry,rhs,out,fmt=7,ev=0.001):
	#qry,rhs,out,fmt,ev='prova'
	line='blastn -query %s -db %s -out %s -evalue %s -outfmt %s' %(qry,rhs,out,ev,fmt)
	os.system(line)
	return out

def make_all_fasta(inp_dir,out_name,tag='.fasta'):
	fastas=[path(inp_dir,file_) for file_ in os.listdir(inp_dir) if file_[-len(tag):]==tag]
	if len(fastas)==0: raise Exception('no %s file in directory selected') %tag
	out=open(out_name,'w')
	for i in fastas:
		for line in open(i).readlines(): out.write(line.strip()+'\n')
	return out_name

def make_all_genes(inp_dir,out_name):
	try: make_all_fasta(inp_dir,out_name,tag='.genes')
	except:
		for draft in [seq for seq in os.listdir(inp_dir) if seq[-6:]=='.fasta']: gene_prediction(path(inp_dir,draft))
		return make_all_fasta(inp_dir,out_name,tag='.genes')

def gene_prediction(draft_path):
	line ='prodigal -i %s -o %s -d %s -a %s -q' %(
					draft_path,draft_path+'_prodigal',draft_path+'_prodigal.genes',draft_path+'_prodigal.prots')
	print 'running prodigal with following command line...'
	print 'cm used:',line
	os.system(line)
	return draft_path+'_prodigal.genes'

##############################

class BLAST_make(object):

	def __init__(self,input_folder):
		self.inp_fold=input_folder
		self.B_out=''

	def make(self,db,query,out):
		launch_BLAST(query,db,out)
		self.B_out=out

	def all_fasta(self,all_fasta_name='all_fasta'):
		if all_fasta_name in os.listdir(self.inp_fold): return all_fasta_name
		self.all_fasta=make_all_fasta(self.inp_fold,all_fasta_name)
		format_db(self.all_fasta)
		return all_fasta_name

	def BLAST_genes(self,all_genes_name='all_genes',B_out='blastall_genes_out'):
		self.all_genes=make_all_genes(self.inp_fold,all_genes_name)
		format_db(self.all_genes)
		launch_BLAST(self.all_genes,self.all_genes,B_out)
		self.B_out=B_out

	def to_parse(self):
		return self.B_out 

