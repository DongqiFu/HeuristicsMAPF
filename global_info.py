import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib as mpl
import numpy.random as rnd
import numpy.linalg as la
import math as ma

def loadData(fileName):
    dataSet = []
    fr = open(fileName)
    for line in fr.readlines():
        curLine=line.strip().split(',')
        dataSet.append([float(curLine[0]),float(curLine[1])])
    return dataSet

def cal_centroid(dataSet):
    for i in range(len(dataSet)):
        sum=np.sum(dataSet,axis=0)
        centroid=sum/len(dataSet)
    return centroid

def cal_eigen(dataSet):
    dataSet = np.array(dataSet)
    sum = np.sum(dataSet, axis=0)
    mean=sum/len(dataSet)
    newdataSet = dataSet - mean
    covariance = np.dot(newdataSet.T, newdataSet) / len(newdataSet)
    eigenvalue = np.linalg.eig(covariance)[0]
    eigenvector=np.linalg.eig(covariance)[1]
    return eigenvector,eigenvalue

from math import acos
from math import sqrt
from math import pi
def length(v):
    return sqrt(v[0]**2+v[1]**2)
def dot_product(v,w):
   return v[0]*w[0]+v[1]*w[1]
def determinant(v,w):
   return v[0]*w[1]-v[1]*w[0]
def inner_angle(v,w):
   cosx=dot_product(v,w)/(length(v)*length(w))
   rad=acos(cosx) # in radians
   return rad*180/pi # returns degrees
def angle_clockwise(A, B):
    inner=inner_angle(A,B)
    det = determinant(A,B)
    if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
        return inner
    else: # if the det > 0 then A is immediately clockwise of B
        return 360-inner

initial_dataSet=loadData('remote_square_data')
initial_centroid=cal_centroid(initial_dataSet)
initial_eigenvector,initial_eigenvalue=cal_eigen(initial_dataSet)

target_dataSet=loadData('circle_data')
target_centroid=cal_centroid(target_dataSet)
target_eigenvector,target_eigenvalue=cal_eigen(target_dataSet)

move_step=target_centroid-initial_centroid

for i in range(len(initial_dataSet)):
    initial_dataSet[i]=initial_dataSet[i]+move_step

initial_centroid=cal_centroid(initial_dataSet)

move_degree=angle_clockwise(initial_eigenvector[0],target_eigenvector[0])
#print move_degree

import math
rotation_matrix=np.array([[math.cos(math.radians(-move_degree)),-math.sin(math.radians(-move_degree))],[math.sin(math.radians(-move_degree)),math.cos(math.radians(-move_degree))]])

for i in range(len(initial_dataSet)):
    initial_dataSet[i]=np.dot((initial_dataSet[i]-initial_centroid),rotation_matrix)+initial_centroid

#print initial_dataSet
initial_eigenvector,initial_eigenvalue=cal_eigen(initial_dataSet)
#print angle_clockwise(initial_eigenvector[0],target_eigenvector[0])

def cal_Elipse(dataSet):
    dataSet = np.array(dataSet)
    eigenvector=[]
    eigenvalue=[]
    for i in range(len(dataSet)):
        newdataSet=np.array(dataSet-dataSet[i])
        sum=np.sum(newdataSet,axis=0)
        #mean=sum/len(dataSet)
        mean=0.0
        newdataSet=newdataSet-mean
        covariance=np.dot(newdataSet.T, newdataSet)/len(newdataSet)
        eigen=np.linalg.eig(covariance)
        #print np.matrix.round(eigen[0],3) #eigenvalue
        #print np.matrix.round(eigen[1],3) #eigenvector
        eigenvalue.append(np.matrix.round(eigen[0],3))  # the first item in eigen represents the eigenvalue
        eigenvector.append(np.matrix.round(eigen[1],3)) # the second item in eigen represents the eigenvector
    return eigenvector,eigenvalue

def degrees_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2' """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    arctan= np.arctan2(sinang, cosang)
    #if (v1[0] > 0 and v1[1] > 0) or (v1[0] < 0 and v1[1] > 0):
    return ma.degrees(arctan)
    #if (v1[0] < 0 and v1[1] < 0) or (v1[0] > 0 and v1[1] < 0):
        #return 360-ma.degrees(arctan)


initial_eigenvector,initial_eigenvalue=cal_Elipse(initial_dataSet)
print 'The coordinate of first crowd is: ',initial_dataSet
print 'The eigenvector of first crowd is: ',initial_eigenvector
print 'The eigenvalue of first crowd is: ',initial_eigenvalue

ells=[mpl.patches.Ellipse(xy=(initial_dataSet[i][0], initial_dataSet[i][1]),
                          width=initial_eigenvalue[i][0], height=initial_eigenvalue[i][1],
                          angle=degrees_between(initial_eigenvector[i][1], [0,initial_eigenvector[i][1][0]])) for i in range(len(initial_dataSet))]

fig,ax=plt.subplots()
for e in ells:
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_alpha(rnd.rand())
    e.set_facecolor(rnd.rand(3))

for i in range(len(initial_dataSet)):
    plt.plot(initial_dataSet[i][0],initial_dataSet[i][1],'bs')
    plt.axis([-10, 10, -10, 10])
    plt.annotate((initial_eigenvalue[i],initial_eigenvector[i]), xy=(initial_dataSet[i][0], initial_dataSet[i][1]))

#check degrees and print
#for i in range(len(dataSet)):
    #print degrees_between(eigenvector[i][1], [0,eigenvector[i][1][0]])

plt.xlabel('X-Coordinate')
plt.ylabel('Y-Coordinate')
plt.show()

target_eigenvector,target_eigenvalue=cal_Elipse(target_dataSet)
print 'The coordinate of second crowd is: ',target_dataSet
print 'The eigenvector of second crowd is: ',target_eigenvector
print 'The eigenvalue of second crowd is: ',target_eigenvalue
ells = [mpl.patches.Ellipse(xy=(target_dataSet[i][0], target_dataSet[i][1]), width=target_eigenvalue[i][0], height=target_eigenvalue[i][1], angle=degrees_between(target_eigenvector[i][1], [0, target_eigenvector[i][1][0]])) for i in range(len(target_dataSet))]
fig, ax = plt.subplots()
for e in ells:
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_alpha(rnd.rand())
    e.set_facecolor(rnd.rand(3))

for i in range(len(target_dataSet)):
    plt.plot(target_dataSet[i][0], target_dataSet[i][1], 'ro')
    plt.axis([-10, 10, -10, 10])
    plt.annotate((target_eigenvalue[i], target_eigenvector[i]), xy=(target_dataSet[i][0], target_dataSet[i][1]))

# check degrees and print
# for i in range(len(dataSet)):
# print degrees_between(eigenvector[i][1], [0,eigenvector[i][1][0]])

plt.xlabel('X-Coordinate')
plt.ylabel('Y-Coordinate')
plt.show()

cost=np.zeros((len(initial_dataSet),len(target_dataSet)))
for i in range(len(initial_dataSet)):
    for j in range(len(target_dataSet)):
        cost[i][j]=(sum(np.abs(np.array(initial_dataSet[i])-np.array(target_dataSet[j])))+np.sum(np.abs(initial_eigenvalue[i]-target_eigenvalue[j]))+np.sum(np.abs(initial_eigenvector[i]-target_eigenvector[j])))
print 'The cost matrix is: ',cost # the row represents the index of the first graph, and the column second

import networkx as nx
from networkx.algorithms import bipartite
B=nx.Graph()

#add nodes
nodes_in_shape1=[]
for i in range(len(initial_dataSet)):
    nodes_in_shape1.append(('S',i))
B.add_nodes_from(nodes_in_shape1, bipartite=0)  # Add the node attribute "bipartite"

nodes_in_shape2=[]
for i in range(len(target_dataSet)):
    nodes_in_shape2.append(('C',i))
B.add_nodes_from(nodes_in_shape2, bipartite=1)

B.add_nodes_from(['source','sink'])

#add edges
mid_edges=[]
for i in range(len(nodes_in_shape1)):
    for j in range(len(nodes_in_shape2)):
        mid_edges.append([nodes_in_shape1[i],nodes_in_shape2[j],{'capacity': 1000, 'weight': float(cost[i][j])}])   #changable in future iteration
B.add_edges_from(mid_edges)

source_mid_edges=[]
for i in range(len(nodes_in_shape1)):
    source_mid_edges.append(['source',nodes_in_shape1[i],{'capacity': 1, 'weight': 0}])
B.add_edges_from(source_mid_edges)

mid_sink_edges=[]
for j in range(len(nodes_in_shape2)):
    mid_sink_edges.append([nodes_in_shape2[j],'sink',{'capacity': 1, 'weight': 0}])
B.add_edges_from(mid_sink_edges)

# verification of the bipartite graph
#bottom_nodes, top_nodes = bipartite.sets(B)
#print nx.is_bipartite(B)
#print nx.is_connected(B)
#print nx.is_biconnected(B)
#print bottom_nodes, top_nodes

# verification of the maxflow
mincostFlow = nx.max_flow_min_cost(B,'source', 'sink',capacity='capacity',weight='weight')
#print mincostFlow
def mappping_and_mincost(d):
    valid_edges=[]
    for k, v in d.iteritems():
        if isinstance(v, dict):
            for kk,vv in v.iteritems():
                if vv==1:
                    valid_edges.append([k,kk,vv])
    index_for_valid_edges=[]
    for e in valid_edges:
        if isinstance(e[0],tuple) and isinstance(e[1],tuple) and e[0][0]=='S':
            index_for_valid_edges.append([e[0][1],e[1][1]])
    total_cost=0.0
    for e in index_for_valid_edges:
        total_cost+=cost[e[0]][e[1]]
    return index_for_valid_edges,total_cost

mapping, mincost=mappping_and_mincost(mincostFlow)
print 'The valid mappingss are: ',mapping
print 'The total minimal cost is: ',mincost

for e in mapping:
    plt.plot(initial_dataSet[e[0]][0], initial_dataSet[e[0]][1],'bs')
    plt.annotate(initial_dataSet[e[0]],xy=(initial_dataSet[e[0]][0], initial_dataSet[e[0]][1]))
    plt.plot(target_dataSet[e[1]][0], target_dataSet[e[1]][1], 'ro')
    plt.annotate(target_dataSet[e[1]],xy=(target_dataSet[e[1]][0], target_dataSet[e[1]][1]))
    plt.axis([-10, 10, -10, 10])
    A=initial_dataSet[e[0]]
    B=target_dataSet[e[1]]
    plt.arrow(A[0], A[1], (B[0]-A[0]), (B[1]-A[1]),head_width=0.07)
plt.show()
