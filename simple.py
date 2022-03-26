import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from spp.convex_sets import Singleton, Polyhedron, Ellipsoid
from spp.convex_functions import TwoNorm, SquaredTwoNorm
from spp.graph import GraphOfConvexSets
from spp.shortest_path import ShortestPathProblem
import sys
import graphviz
import itertools
from itertools import tee

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


#https://en.wikipedia.org/wiki/Line_graph
def SaveGraph(G, filename):
    plt.figure(figsize=(4,5))
    G.draw_sets()
    G.draw_edges()
    plt.savefig(filename, bbox_inches='tight')


def SaveAsGraph(G, filename):
    L = graphviz.Digraph()
    for vertex, set in G.sets.items():
        L.node(f'vertex')

    for i, e in enumerate(G.edges):
        print(e[0], e[1])
        L.edge(e[0], e[1], f'{i}')

    L.render(filename)


def MakeLineGraph(G):
    L = graphviz.Graph()
    for i, e in enumerate(G.edges):
        L.node(f'{i}')

    for vertex, set in G.sets.items():
        edges_in = G.incoming_edges(vertex)[1]
        edges_out = G.outgoing_edges(vertex)[1]
        print( (vertex, edges_out))
        #for ei, ej in pairwise(edges_in):
        for ei, ej in itertools.combinations(edges_in, 2):
                L.edge(f'{ei}', f'{ej}', f'{vertex}')

        for ei, ej in itertools.combinations(edges_out, 2):
                print( (ei, ej))
                L.edge(f'{ei}', f'{ej}', f'{vertex}')

    L.render('line_graph')
    SaveAsGraph(G, 'digraph')


#def MakeLineGraph2(G):
#    L = GraphOfConvexSets()
#    sets = []
#    vertices = []
#    row = 0
#    spacing = 2
#    for i, e in enumerate(G.edges):
#        sets.append(Singleton((i*spacing, row)))
#        row =  i % 3
#        vertices.append(f'{i}')
#    print(sets)
#    print(vertices)
#    L.add_sets(sets, vertices)
#
#    l = TwoNorm(H)
#    for vertex, set in G.sets.items():
#        edges_in = G.incoming_edges(vertex)[1]
#        edges_out = G.outgoing_edges(vertex)[1]
#        for ei in edges_in:
#            for ej in edges_in:
#                print(ei, ej)
#                if (ei == ej):
#                    continue
#                L.add_edge(f'{ei}', f'{ej}', l)
#
#    SaveGraph(L, 'dual.pdf')

# convex sets
singletons = (
    Singleton((0, 0)),
    Singleton((9, 0)),
)
polyhedra = (
    Polyhedron.from_vertices(([1, 0], [1, -2], [3, -2], [3, -1])),
    Polyhedron.from_vertices(([4, -2], [5, -4], [3, -4], [2, -3])),
    Polyhedron.from_vertices(([2, 2], [1, 3], [2, 4], [4, 4], [4, 3])),
)
ellipsoids = (
    Ellipsoid((4, 1), ([1, 0], [0, 1])),
)
sets = singletons + polyhedra + ellipsoids

# label for the vertices
vertices = ['s', 't']
vertices += [f'p{i}' for i in range(len(polyhedra))]
vertices += [f'e{i}' for i in range(len(ellipsoids))]

# add convex sets to the graph
G = GraphOfConvexSets()
G.add_sets(sets, vertices)
G.set_source('s')
G.set_target('t')

# edges
H = np.hstack((np.eye(2), -np.eye(2)))
l = TwoNorm(H)
edges = {
    's': ('p0', 'p1', 'p2'),
    'p0': ('e0',),
    'p1': ('p2', 'e0'),
    'p2': ('t', 'e0'),
    'e0': ('t' ),
}
for u, vs in edges.items():
    for v in vs:
        G.add_edge(u, v, l)
        
# draw convex sets and edges
plt.figure()
G.draw_sets()
G.draw_edges()

MakeLineGraph(G)
SaveGraph(G, 'graph.pdf')
sys.exit()
#G.label_sets()

#spp = ShortestPathProblem(G, relaxation=0)
spp = ShortestPathProblem(G, relaxation=1)
sol = spp.solve()

print('Cost:', sol.cost)
print('\nFlows:')
for k, edge in enumerate(G.edges):
    flow = round(abs(sol.primal.phi[k]), 4)
    print(edge, flow)

# edge lenghts
l2 = SquaredTwoNorm(H)
G2 = deepcopy(G)
for e in G2.edges:
    G2.set_edge_length(e, l2)

spp2 = ShortestPathProblem(G2, relaxation=0)
sol2 = spp2.solve()




print('Cost:', sol2.cost)
print('\nFlows:')
for k, edge in enumerate(G2.edges):
    flow = round(abs(sol2.primal.phi[k]), 4)
    print(edge, flow)

plt.figure(figsize=(4,5))
G.draw_sets()
G.draw_edges()

offset = np.array([0, -.25])
plt.text(*(G.source_set.center + offset), r'$X_s$', ha='center', va='top')
plt.text(*(G.target_set.center + offset), r'$X_t$', ha='center', va='top')

plt.plot([np.nan] * 2, 'b--', label='Euclidean distance')
plt.plot([np.nan] * 2, 'r-.', label='Euclidean distance squared')
G.draw_path(sol.primal.phi, sol.primal.x, color='b', linestyle='--')
G.draw_path(sol2.primal.phi, sol2.primal.x, color='r', linestyle='-.')

plt.xticks(range(10))
plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0))
plt.grid()
plt.savefig('2d_setup.pdf', bbox_inches='tight')


sys.exit()

scales = np.logspace(-2, 2, 200)
c_micp = []
c_relaxation = []
c_micp2 = []
c_relaxation2 = []
for s in scales:
    
    G_scaled = deepcopy(G)
    G_scaled.scale(s)
    spp = ShortestPathProblem(G_scaled, relaxation=0)
    c_micp.append(spp.solve().cost)
    spp = ShortestPathProblem(G_scaled, relaxation=1)
    c_relaxation.append(spp.solve().cost)
    
    G2_scaled = deepcopy(G2)
    G2_scaled.scale(s)
    spp = ShortestPathProblem(G2_scaled, relaxation=0)
    c_micp2.append(spp.solve().cost)
    spp = ShortestPathProblem(G2_scaled, relaxation=1)
    c_relaxation2.append(spp.solve().cost)

def micp_vs_relaxation(c_micp, c_relaxation):
    plt.plot(scales, c_micp, label='MICP', linestyle='-', linewidth=2)
    plt.plot(scales, c_relaxation, label='Convex relaxation', linestyle='--', linewidth=2)
    plt.xlabel(r'Scale factor $r$')
    plt.ylabel('Cost')
    plt.xlim([scales[0], scales[-1]])
    plt.xscale('log')
    plt.grid(1)
    plt.legend(loc=1)

plt.figure(figsize=(3.5, 3))
micp_vs_relaxation(c_micp, c_relaxation)

plt.figure(figsize=(3.5, 3))
micp_vs_relaxation(c_micp2, c_relaxation2)
