import os


##########################################

def get_seq(f):
	''' generator from fasta file It yields a list of lines
	(the first is the seq id, the others the seq)
	'''
	for line in open(f):
		line=line.strip()
		if line.startswith('>'):
			try: yield seq
			except: pass
			seq=[line]
		else:
			try: seq.append(line)
			except: pass
	try: yield seq
	except: print f,'may be not good as input...'

def tag(file_,tag):
	os.rename(file_,file_+tag)
	
def tag_all(dir_,tag_):
	for file_ in os.listdir(dir_):
		file_=dir_+file_
		tag(file_,tag_)

def path(path_,file_):
	if path_[-1]!='/': path_+='/'
	return path_+file_


def clean_folder(folder,tag='.fasta'):
	for i in os.listdir(folder):
		if not i.endswith(tag): os.remove(path(folder,i))
		
def get_draft(folder,tag='.fasta'):
	out=[]
	tag_len=len(tag)
	for i in os.listdir(folder):
		if i[-tag_len:]==tag: out.append(i)
	return out
	

def inp_rename(folder,tag='.fasta'):
	fastas=[file_ for file_ in os.listdir(folder) if file_.endswith(tag)]
	if len(fastas)==0: raise Exception('no fastas in %s' %folder) 
	for i in fastas:
		count,diz=1,{}
		text=[]
		for line in open(path(folder,i)).readlines():
			line=line.strip()
			if len(line)==0: continue
			if line[0]=='>':
				new_line='>%s_!contig_%s' %(i,count)
				count+=1
				diz[line]=new_line
				text.append(new_line+'\n')
			else: text.append(line+'\n')
		conv_table=open(path(folder,i)+'_conversion_table','w')
		out=open(path(folder,i),'w')
		for line in text:
			out.write(line)
		out.close()
		for k,v in diz.iteritems():
			conv_table.write('%s\t%s\n' %(k,v))
		conv_table.close()

def inp_rename_len(folder,tag='.fasta'):
	fastas=[file_ for file_ in os.listdir(folder) if file_.endswith(tag)]
	if len(fastas)==0: raise Exception('no fastas in %s' %folder) # raise error if no input file
	#
	for fasta in fastas: # for each input file
		count,diz=1,{}# set a dict for conversion table, and a counter for new contig name
		out=open(path(folder,'overwrite_handler'),'w') # use a temporary filehandler to overwrite stuff
		for contig in get_seq(path(folder,fasta)): # for each contig
			old_name=contig[0]
			c_length=len(''.join(contig[1:]))
			new_name='>%s_!contig_%s_&%s' %(fasta,count,c_length)
			count+=1
			diz[old_name]=new_name
			contig[0]=new_name
			to_write='\n'.join(contig)+'\n'
			out.write(to_write)
		os.rename(path(folder,'overwrite_handler'),path(folder,fasta)) # rewrite filehandler into original file
		conv_table=open(path(folder,fasta)+'_conversion_table','w') # make conversion table file and write into it the dict
		for k,v in diz.iteritems():
			conv_table.write('%s\t%s\n' %(k,v))
		conv_table.close()		
			
def instanciate_drafts(folder):
	drafts={}
	for draft in get_draft(folder):
		print 'building %s' %draft
		draft=Draft(draft,folder)
		draft.build()
		drafts[draft.name]=draft
	return drafts
	
def instanciate_single_draft(path):
	name,folder=path.split('/')[-1],'/'.join(path.split('/')[:-1])
	return Draft(name,folder)


##########################################


class Draft(object):

	def __init__(self,fasta):
		self.fasta=fasta
		self.name=self.fasta.split('/')[-1]
		self.path='/'.join(self.fasta.split('/')[:-1])+'/'
		self.contigs={}
		self.adjacency_list={}

	def add_contig(self,contig):
		self.contigs[contig.name]=contig

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

	@property
	def length(self):
		return len(self.seq)

	def samedraft(self,cntg):
		return self.draft==cntg.draft


#dir_='../../Prove/'
#inp_rename(dir_)
