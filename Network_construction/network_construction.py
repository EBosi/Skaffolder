from Homology import BLAST_parser as blp
from Homology import Homology
from Input_processing import inp_proc
from Synteny import Synteny
from Utils import Blast_cl as bcl
import optparse,sys,os
import networkx as nx


##########################################


usage="Usage: sticazzi?!?!?!?"
parser=optparse.OptionParser(usage)
parser.add_option('-i','--inp',help='give the Input Sequences folder path',dest='inp',default=None)
parser.add_option('-o','--out',help='give the Output folder path',dest='out',default=None)
#parser.add_option('-h','--help',help='ask for help',dest='help',default=True)
opts,args=parser.parse_args()
opts.help=False
if (not opts.inp) or (not opts.out):
	opts.help=True
	print 'mandatory options are found to be missing... try again please'
if opts.help:
	print usage
	sys.exit()


def check_path(path):
	if path[-1]=='/': return path
	return path+'/'

##########################################


#process folders and instanciate networks
input_folder=opts.inp
output_folder=opts.out
input_folder,output_folder=check_path(input_folder),check_path(output_folder)
to_parse='blastall_out'
inp_proc.clean_folder(input_folder)
inp_proc.clean_folder(output_folder,'')
inp_proc.inp_rename(input_folder)
print "\ninput processing finshed... now instanciating network objects"
nets={}
for name,draft in inp_proc.instanciate_drafts(input_folder).iteritems():
	nets[name]=Homology.Net(draft)

#Homology part
print "\nnow starting Homology part"
for net in nets:
	blast_all=bcl.BLAST_make(input_folder)
	db=blast_all.all_fasta()
	print 'blasting %s on db %s.... producing %s as output' %(nets[net].fasta,db,nets[net].blast_out)
	blast_all.make(db,nets[net].fasta,nets[net].blast_out)

	print "\nBLASTed',nets[net].name,'file.. now connecting contigs"
#blast_all_parsed=Homology.Homology_parsing(blast_all.to_parse())
	blast_out_parsed=Homology.Homology_parsing(nets[net].blast_out)
	print 'constructing network...'
	blast_out_parsed.net_construction(nets)
print "\nBLASTing all this shit done... now starting Synteny part"

#Synteny part
all_genes=Synteny.synteny(input_folder)
all_genes.run(nets)
print "\nBLASTing genes done... now exporting nets"

#Export nets in gexf
for net in nets:
	path=output_folder+net
	nets[net].export(path)
	
print "\nNetworks succesfully built!"

