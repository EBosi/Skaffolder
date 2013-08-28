from new_net_constr_fake import *
import sys

########################################
# Args
########################################

args=sys.argv
hit_thr=int(args[1])
bord_thr=int(args[2])
try: id_perc_thr=float(args[3])
except: id_perc_thr=0.9

########################################
# Main
########################################

fasta='../Input/44TAC125_vecchio_ray31.fasta'
all_fasta='../Input/all_fasta'
output_folder='../Output/'
out=fasta[:-6]+'_blast_out'
draft=Draft(fasta)
draft.build()
bcl.launch_BLAST(fasta,all_fasta,out)
network=parse_blast_tree(out,hit_thr,bord_thr,id_perc_thr)
add_singletons(network,draft.contigs)
graph_name='%sTAC125_%s_%s_tree_method.gexf' %(output_folder,hit_thr,bord_thr)
print graph_name
nx.write_gexf(network,graph_name)


