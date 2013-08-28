import sets


##############################


def get_draft(gene):
	draft=gene.split('_!')[0]
	return draft


def list_getter(lista,prefix):
	getter=[x for x in lista if get_draft(x[0])==prefix]
	if len(getter)==0: getter=True
	else: getter=getter[0]
	return getter


def get_contig(gene):
	contig='_'.join(gene.split('_')[:-1])
	return contig


def best_hits(diz):
	best_diz={}
	for k in diz:
		for v in diz[k]:
			if k==v: continue
			try: best_diz[k]
			except: best_diz[k]=[(v,diz[k][v].evalue)]
			else:
				getter=list_getter(best_diz[k],get_draft(v))
				if getter==True: best_diz[k].append(v,diz[k][v].evalue)
				else:
					if getter[1]>diz[k][v].evalue:
						best_diz[k].remove(getter)
						best_diz[k].append(v,diz[k][v].evalue)
	return best_diz


def bidirectional(best_hit_diz):
	bhd=best_hit_diz
	bidirectional={}
	for k in bhd:
		for best_hit in bhd[k]:
			if k in [x[0] for x in bhd[best_hit[0]]]:
				try: bidirectional[k]
				except: bidirectional[k]=[best_hit[0]]
				else: bidirectional[k].append(best_hit[0])
	return bidirectional


def inverse_bidirectional(diz):
	inverse_diz={}
	for k in diz:
		for j in diz[k]:
			try: inverse_diz[j].append(k)
			except: inverse_diz[j]=[k]
	return inverse_diz


##############################


class hit(object): 

	def __init__(self,query,subj,hit,nets,length=0):
		# query, subject, % id, align length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
		self.query=query
		self.subj=subj
		self.hit=hit
		self.draft=self.query.split('_!')[0]
		if length!=0: self.query_length=length
		else: self.query_length=nets[self.draft].draft.contigs['>'+self.query].length
		self.qstart=int(hit[4])
		self.qend=int(hit[5])
		self.sstart=int(hit[6])
		self.send=int(hit[7])
		self.evalue=float(hit[8])
		self.bscore=float(hit[9])
		self.id=float(self.hit[0])
		self.al_length=float(self.hit[1])
		self.coverage= (abs(self.qstart-self.qend)+1) / float(self.query_length)
		self.subj_set=sets.Set(xrange(min(self.sstart,self.send),max(self.sstart,self.send)))
		self.query_set=sets.Set(xrange(min(self.qstart,self.qend),max(self.qstart,self.qend)))


class parser(object):

	def __init__(self,blast_out):
		self.blast_out=open(blast_out)
		self.hits={}
		self.inverse_hits={}
		
	
	def parse(self,nets):
		for line in self.blast_out.readlines():
			line=line.strip()
			print line
			if line[0]!='#':
				print'!!! line'
				hit_=hit(line.split()[0],line.split()[1],line.split()[2:],nets)
				self.hits[hit_.query],self.inverse_hits[hit_.subj]=self.hits.get(hit_.query,{}),self.inverse_hits.get(hit_.subj,{})
				self.hits[hit_.query][hit_.subj],self.inverse_hits[hit_.subj][hit_.query]=hit_,hit_
				
	def parse_genes(self,nets):
		for line in self.blast_out.readlines():
			line=line.strip()
			if line[:7]=='# Query':
				l=line.split('#')
				length=(int(l[3])-int(l[2])+1)
		if line[0]!='#':
				hit_=hit(line.split()[0],line.split()[1],line.split()[2:],nets,length)
				self.hits[hit_.query],self.inverse_hits[hit_.subj]={},{}
				self.hits[hit_.query][hit_.subj],self.inverse_hits[hit_.subj][hit_.query]=hit_,hit_
			

	def BBH(self,nets):
		self.parse_genes(nets)
		return inverse_bidirectional(bidirectional(best_hits(self.hits)))

	'''def parse_single(self,draft):
		for line in self.blast_out.readlines():
			line=line.strip()
			if line[0]!='#':'''
				


