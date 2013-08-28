import os,sys
import networkx as nx

#############################

class abc():
	def __init__(self,graph_file):
		self.graph=graph_file
	def write_in(self,out,weight):
		self.G=nx.read_gexf(self.graph)
		self.out=open(out,'w')
		for edge in self.G.edges():
			node1,node2=edge[0],edge[1]
			try: self.G[node1][node2][weight]
			except: self.out.write('%s\t%s\t%s\n' %(node1,node2,0))
			else: self.out.write('%s\t%s\t%s\n' %(node1,node2,self.G[node1][node2][weight]))
		self.out.close()
		return out
	
#############################

if __name__=='__main__':
	path=sys.argv[1]
	weight=sys.argv[2]
	out=path+'.abc' 
	file_=abc(path)
	file_.write_in(out,weight)

