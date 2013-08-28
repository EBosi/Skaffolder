from Homology import BLAST_parser as blp
from Homology import Homology
from Input_processing import inp_proc
from Synteny import Synteny
from Utils import Blast_cl as bcl
import optparse,sys,os
import networkx as nx



inpt='../Input/44TAC125_vecchio_ray31.fasta'
blastall='blastall_out'
blastall_genes='blastall_genes_out'
draft=inp_proc.instanciate_single_draft(inpt)
parser=Homology.Homology_parsing(blastall)


