# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt

import copy
import numpy as np


from firedrake import *
from firedrake_adjoint import *


class Clarke_Wright(object):
    
    def __init__(self, turbine_locations, substation_location, plot_output = False):
        self.turbine_locations = turbine_locations
        self.substation_location = substation_location
        self.plot_output = plot_output
      
    def routing_coordinates(self, vertices, rting):
        '''Convert a routing expressed in indexes to a routing expressed in coordinates'''        
        for i in range(len(rting)):
            rting[i] = vertices[rting[i]]
        return rting
    

    def routing_indexes(self, R, vertices):
        '''Convert a routing expressed in tuples to a routing expressed in indexes'''
        chromosome = [0]
        while len(chromosome) <= len(R):
            for i in range(len(R)):
                if R[i][0] == chromosome[-1]:
                    chromosome.append(R[i][1])
                if R[i][1] == chromosome[-1]:
                    chromosome.append(R[i][0])
        return chromosome
    
    
    # def produce_plot(self, R, vertices):
    #     '''Display with matplotlib'''
    #     R = self.routing_indexes(R, vertices)
    #     R = self.routing_coordinates(vertices, R)
    #     R = np.array(([R])) 
    #     V = np.array((vertices))
    #     plt.plot(V[:,0],V[:,1],'o')
    #     plt.plot(R[0][:,0], R[0][:,1], '-')
    #     plt.text(V[0][0], V[0][1], '%s' % (str(0)))
    #     plt.axis('equal')
    #     plt.savefig('cable_result.jpg',dpi=300) 
        

    def construct_cost_matrix(self, vertices):
        '''Constructs a matrix of costs for every potential edge - connections of vertices to themselves are set to infinity'''
        distances = []
        grouped_distances = []
        for i in range(len(vertices)):
            for j in range(len(vertices)):
                if i == j:
                    dist = 99
                else:
                    dist = (sqrt((vertices[i][0] - vertices[j][0])**2 + (vertices[i][1] - vertices[j][1])**2))
                distances.append(dist)
        for i in range(0, len(distances), len(vertices)):
            grouped_distances.append(tuple(distances[i:i+len(vertices)]))
        C = np.array(grouped_distances)
        return C
        

    def make_graph(self, Rgraph, no_depot = True):
        '''Produces a dictionary representing the graph of the form (vertex: [connected vertices], ...) can remove the depot so as to seperate the routes from one another'''
        if no_depot:
            remove_list=[]
            for edge in Rgraph:
                if edge[0] == 0 or edge[1] == 0:
                    remove_list.append(edge)
            for edge in remove_list: 
                Rgraph.remove(edge)
        vertices = []
        for a_tuple in Rgraph:
            vertices.extend(list(a_tuple))
        vertices = list(set(vertices))
        nVertices = len(vertices)
        G = {}
        for i in range(0,nVertices):
            G[vertices[i]]=[]
        for edge in Rgraph:
            G[edge[1]].append(edge[0])
            G[edge[0]].append(edge[1])
        for vertex in G:
            G[vertex] = list(set(G[vertex]))
        return G
        

    def initialise_R(self, n):
        '''produce the initial edge set - all clients connected to their nearest depot'''
        R = []
        # For each client vertex, identify the nearest depot and add (v,d) to initial edge list R
        for index in range(1, n):
            R.append((index, 0))
        return R
        

    def initialise_G(self, R, n):
        '''produce the graph of the initial edge set'''
        G = []
        Rgraph = copy.deepcopy(R)
        G = self.make_graph(Rgraph)
        return G
        

    def construct_savings_list(self, C, n):
        '''Constructs a list of the savings in making every connection and orders them'''
        S = []    
        # For each possible pair of clients in the graph, calculate the cost saving from joining with edge, and add to list S
        R = self.initialise_R(n)        
        for vertex1 in range(1, n):
            for vertex2 in range(1, n):
                if not vertex1 == vertex2:
                    depot = [j[1] for i, j in enumerate(R) if j[0] == vertex1]
                    saving = C[vertex1, depot] - C[vertex1, vertex2]
                    S.append((vertex1, vertex2, saving))
        # Sort the list S of cost savings in decreasing order
        S.sort(key = lambda s:s[2], reverse = True)
        return S


    def find_path_length(self, graph, start, end, path=[]):
        '''Returns a list of vertices on a path in order that they are connected'''        
        path = path + [start]
        if start == end:
            return path
        if not start in graph:
            return None
        for node in graph[start]:
            if node not in path:
                newpath = self.find_path_length(graph, node, end, path)
                if newpath: return newpath
        return None


    def neighbour_depot(self, edge, R, Vd):
        '''Ensures that k is neighbouring a depot - ensures no short-circuiting'''      
        nhbrDepot = False
        if (edge[0], 0) or (0, edge[0]) in R:
            nhbrDepot = True
        return nhbrDepot

        
    def one_neighbour(self, edge, R):
        '''Prevent route branching by checking that u has only one neighbour'''
        oneNhbr = True
        for edge_temp in R:
            if edge_temp[1] != 0 and edge_temp[0] == edge[0]:
                oneNhbr = False
            if edge_temp[1] == edge [1]:
                oneNhbr = False
        return oneNhbr

    
    def perform_checks(self, S, G, R, n, edge):
        '''Performs viability checks on proposed saving'''
        if not (self.find_path_length(G, edge[0], edge[1])==None):           
            return False
        if not self.neighbour_depot(edge, R, range(1, n)):
            return False
        if not self.one_neighbour(edge, R):
            return False
        return True

        
    def k_d(self, edge, R):
        '''Defines the edge whose removal is proposed'''
        if (edge[0],0) in R:
            k_d = (edge[0],0)
        elif (0,edge[0]) in R:
            k_d = (0,edge[0])
        return k_d

        
    def perform_merge(self, edge, R, k_d):
        '''Removes the superfluous old edge and inserts the new edge & updates the routing graph'''
        R.remove(k_d)
        R.append((edge[0], edge[1]))
        Rgraph = copy.deepcopy(R)
        G = self.make_graph(Rgraph)
        return G

        
    def run(self):
        vertices = self.substation_location + self.turbine_locations
        n = len(self.turbine_locations)+1
        C = self.construct_cost_matrix(vertices)
        S = self.construct_savings_list(C, n)
        R = self.initialise_R(n)
        G = self.initialise_G(R, n)
        while (not S == []) and (S[0][2][0] > 0):
            edge = (S[0][0],S[0][1])
            allowable = self.perform_checks(S, G, R, n, edge)
            if allowable:
                G = self.perform_merge(edge, R, self.k_d(edge, R))
            else:                 
                del S[0]
        if self.plot_output:
            self.produce_plot(R, vertices)
        R = self.routing_indexes(R, vertices)
        R = self.routing_coordinates(vertices, R)
        return R
