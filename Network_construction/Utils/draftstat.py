#!/usr/bin/env python

'''
DraftStat

Class to generate statistical indices and figures for draft genomes

'''

__author__='Bozzo'
#####################
#Imports

import os,sys
import matplotlib.pyplot as plt

####################
#Functions

def binner(value,bin):
	if value%bin==0:
		return bin*(value/bin)
	else:
		return bin*((value/bin)+1)

def gc(seq):
	diz={'A':0,'G':0,'C':0,'T':0}
	for nu in seq:
		diz[nu]=diz[nu]+1.0
	return diz['G']+diz['C']/diz['A']+diz['T']+diz['G']+diz['C']

def parse_prodigal(handler):
	for line in handler:
		line=line.strip()
		if line.startswith('>'):
			try: yield cntg,head,start,end,strand,seq
			except: pass
			line=line.split('#')
			cntg,head,start,end,strand,seq='_'.join(line[0].split('_')[:-2]),line[0],line[1],line[2],line[3],''
		else: seq+=line

def give_path(path):
	path=path.strip('./')
	if '/' not in path:
		return './'
	else:
		return '/'.join(path.split('/')[:-1]) +'/'

#####################
#Classes

class Contig(object):

	def __init__(self,name,seq):
		self.name=name.split()[0]
		self.seq=seq
		self.length=len(seq)
		self.gc=gc(seq)
		self.genes={}
		self.prots={}

	def num_genes(self):
		return len(self.genes)
	
	def num_prots(self):
		return len(self.prots)

	def add_gene(self,harvest):
		gene=Gene(harvest[0],harvest[1],harvest[2],harvest[3],harvest[4],harvest[5])
		self.genes[gene.name]=gene
		try: self.prots[gene.name]
		except KeyError: pass
		else: 
			self.prots[gene.name].gene=gene
			self.genes[gene.name].prot=self.prots[gene.name].gene

	def add_prot(self,harvest):
		prot=Prot(harvest[0],harvest[1],harvest[2],harvest[3],harvest[4],harvest[5])
		self.prots[prot.name]=prot
		try: self.genes[prot.name]
		except KeyError: pass
		else:
			self.genes[prot.name].prot=prot
			self.prots[prot.name].gene=self.genes[prot.name]


class Gene(object):

	def __init__(self,cntg,name,start,end,strand,seq):
		self.name=name
		self.cntg=cntg
		self.start=start
		self.end=end
		self.strand=strand
		self.seq=seq
		self.prot=None
		self.length=len(self.seq)

	def gc(self):
		return gc(self.seq)
	

class Prot(object):

	def __init__(self,cntg,name,start,end,strand,seq):
		self.name=name
		self.cntg=cntg
		self.start=start
		self.end=end
		self.strand=strand
		self.seq=seq
		self.gene=None
		self.length=len(self.seq)	

class Histo():
	def __init__(self,name,data,**kwargs):
		self.name=name
		self.data=data
		self.kw=kwargs		
	def print_(self):
		plt.hist(self.data,50)
		if self.kw.get('xlab'): plt.xlabel(self.kw['xlab'])
		if self.kw.get('ylab'): plt.ylabel(self.kw['ylab'])
		plt.title=self.name
		plt.show()
	def export(self,name,format='pdf'):
		plt.hist(self.data,50)
		if self.kw.get('xlab'): plt.xlabel(self.kw['xlab'])
		if self.kw.get('ylab'): plt.ylabel(self.kw['ylab'])
		plt.title=self.name
		plt.savefig(name,format=format)
		plt.close()

class Scatterplot():

	def __init__(self,name,data):
		self.name=name
		self.data=data
		self.x=[]
		self.y=[]
		for couple in self.data:
			self.x.append(couple[0])
			self.y.append(couple[1])
			
	def print_(self):
		plt.plot(self.x,self.y,'ro')
		plt.ylabel('number of genes')
		plt.xlabel('contig length')
		plt.title=self.name
		plt.show()

	def export(self,name,format='pdf'):
		plt.plot(self.x,self.y,'ro')
		plt.ylabel('number of genes')
		plt.xlabel('contig length')
		plt.title=self.name
		plt.savefig(name,format=format)
		plt.close()

class Draft(object):

	def __init__(self,name,path):
		self.name=name
		self.path=path
		self.contigs={}
		for line in open(self.path):
			line=line.strip()
			if line.startswith('>'):
				try:
					cntg=Contig(head,seq)
					self.contigs[cntg.name]=cntg
				except: pass
				head,seq=line,''
			else: seq+=line
		cntg=Contig(head,seq)
		self.contigs[cntg.name]=cntg
			
	def num_contigs(self):
		return len(self.contigs)

	def length(self):
		length=0.0
		for cntg in self.contigs:
			length+=self.contigs[cntg].length
		return length

	def num_genes(self):
		genes=0.0
		for cntg in self.contigs:
			genes+=len(cntg.genes)
		return genes

			
	def prodigal(self,prodigal_path='/home/bozzo/Python/prodigal.v2_50.linux'):
		try: open(prodigal_path)
		except IOError: 
			try:prodigal_path=raw_input(' please input your prodigal path... ')
			except IOError:
				print 'prodigal not found... did you install it?'
				sys.exit()
		k=self.path
		print 'launching prodigal...'
		line ='%s -i %s -d %s -a %s' %(
				prodigal_path,k,k+'_prodigal_genes',k+'_prodigal_prots')
		os.system(line)

	def get_genes(self):
		try: reader=open(self.path+'_prodigal_genes')
		except IOError: 
			self.prodigal()
			reader=open(self.path+'_prodigal_genes')
		for harvest in parse_prodigal(reader):
			self.contigs[harvest[0]].add_gene(harvest)
	
	def get_prots(self):
		try: reader=open(self.path+'_prodigal_prots')
		except IOError: 
			self.prodigal()
			reader=open(self.path+'_prodigal_prots')
		for harvest in parse_prodigal(reader):
			self.contigs[harvest[0]].add_gene(harvest)

	def codon_usage(self):
		return

	def gc(self):
		return


##############
#Statistics


	## Histograms ##

	def length_distr(self,out=[]):
		for cntg in self.contigs:
			out.append(self.contigs[cntg].length)
		plt.hist(out,50)
		title='Contig length distribution'
		xlabel='Contig length'
		return Histo(title,out,xlab=xlabel)
		

	'''def gene_cntg_length_distr(self,bins=500):
		out=[]
		diz,diz2={},{}
		for cntg in self.contigs:
			diz[binner(self.contigs[cntg].length,bins)]=diz.get(binner(self.contigs[cntg].length,bins),0)+1
			for gene in self.contigs[cntg].genes:
				diz2[binner(self.contigs[cntg].length,bins)]=diz2.get(binner(self.contigs[cntg].length,bins),0)+1
		f=diz1.keys()
		f.sort()
		out=map(lambda x: diz2[x]/float(diz[x]) , f)
		for cntg in self.contigs:
			for gene in self.contigs[cntg].genes:
				out.append(self.contigs[cntg].length)
		bins=f
		plt.hist(out,50)
		title='Genes per Contig_length distribution'
		xlabel='Contig length'
		ylabel='Number of genes'
		return Histo(title,out,xlab=xlabel,ylab=ylabel)'''

	def gene_length_distr(self,out=[]):
		for cntg in self.contigs:
			for gene in self.contigs[cntg].genes:
				out.append(self.contigs[cntg].genes[gene].length)
		title='Gene length distribution'
		xlabel='Gene length'
		return Histo(title,out,xlab=xlabel) 


	## Scatterplots ##
	
	def genes_vs_cont_length(self):
		vect=[]
		for cntg in self.contigs:
			couple=(self.contigs[cntg].length,len(self.contigs[cntg].genes))
			vect.append(couple)
		return Scatterplot('Genes vs contig length',vect)


	## Indices ##

	def draft_over_genes(self):
		return self.length()/self.num_genes()

	def cont_vs_draft(self):
		return len(self.contigs)/self.length()

	def cont_vs_mean_lenght(self):
		return (len(self.contigs)*len(self.contigs))/self.length()

	def N50(self):
		sky,ground,last=self.length()/2.0,0,0
		L=[ self.contigs[cntg].length for cntg in self.contigs ]
		for i in L:
			if ground + i < sky:
				last=i
				ground+=i
			else:
				break
		return last

	## Export ##

	def export_figs(self,path,format='pdf'):
		self.length_distr().export(path+'contig_length_distribution',format=format)
		self.gene_length_distr().export(path+'gene_length_distribution',format=format)
		self.genes_vs_cont_length().export(path+'gene_vs_cont_length',format=format)
		#self.gene_cntg_length_distr().export(path+'length_distribution',format=format)
	
