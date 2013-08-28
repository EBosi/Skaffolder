#!/usr/bin/env python
"""
Net

Kegg Library

Networks made using Kegg data
"""
import logging
import networkx as nx

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.net')

################################################################################
# Classes

class MetabolicNet(object):
    '''
    A metabolic network as a networkX graph
    
    Nodes: kegg compounds (w co_id, name; w or w/o weight)
    Edges: kegg reactions (w co1, co2, re_id, name; w or w/o weight)
    
    nodes weight indicate the activity index
    edges weight indicate the copy number 
    '''
    def __init__(self, nodes=None, edges=None, name='MetNet'):
        self.name = name
        
        self.net = nx.Graph()
        
        if nodes:
            for n in nodes:
                self.net.add_node(n.co_id, name=n.name)
                if hasattr(n, 'weight'):
                    self.net.node[n.co_id]['weight'] = n.weight
        
        if edges:
            for e in edges:
                self.net.add_edge(e.co1, e.co2, re_id=e.re_id, name=e.name)
                if hasattr(e, 'weight'):
                    self.net.edge[e.co1][e.co2]['weight'] = e.weight
                    
    def setNet(self, net):
        '''
        Use an external networkx graph
        '''
        self.net = net
                    
    def addNodes(self, nodes):
        '''
        Takes a compounds iterable and adds (or updates) the nodes
        w co_id, name; w or w/o weight attributes
        nodes weight indicate the activity index
        '''
        for n in nodes:
            self.net.add_node(n.co_id, name=n.name)
            if hasattr(n, 'weight'):
                self.net.node[n.co_id]['weight'] = n.weight
