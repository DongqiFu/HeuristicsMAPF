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

#initial_eigenvector,initial_eigenvalue=cal_eigen(initial_dataSet)
#print angle_clockwise(initial_eigenvector[0],target_eigenvector[0])
print 'after rotation: ',initial_dataSet
#now, data is ready in the same frame

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
        eigenvalue.append(np.matrix.round(eigen[0],3))  # the first item in eigen represents the eigenvalue
        eigenvector.append(np.matrix.round(eigen[1],3)) # the second item in eigen represents the eigenvector
    return eigenvector,eigenvalue

initial_eigenvector,initial_eigenvalue=cal_Elipse(initial_dataSet)
target_eigenvector,target_eigenvalue=cal_Elipse(target_dataSet)

def cal_cost(initial_dataSet,target_dataSet,initial_eigenvalue,target_eigenvalue,initial_eigenvector,target_eigenvector):
    return (sum(np.abs(initial_dataSet-target_dataSet))+np.sum(np.abs(initial_eigenvalue-target_eigenvalue))+np.sum(np.abs(initial_eigenvector-target_eigenvector)))

flag=False
while(flag==False):
    min_cost=10000000
    for i in range (len(initial_dataSet)):
        #pick the most similar node
        for j in range(len(target_dataSet)):
            cost=cal_cost(initial_dataSet[i],target_dataSet[j],initial_eigenvalue[i],target_eigenvalue[j],initial_eigenvector[i],target_eigenvector[j])
            if cost< min_cost:
                min_cost=cost
                min_cost_index=j
        #pick the best direction to move unit length --dictionary store, find min cost, extract key word, settle coordinate based on key word
        #up
        initial_dataSet[i]+=[0, 0.01]
        initial_eigenvector, initial_eigenvalue = cal_Elipse(initial_dataSet)
        N_cost = cal_cost(initial_dataSet[i], target_dataSet[min_cost_index], initial_eigenvalue[i], target_eigenvalue[min_cost_index],initial_eigenvector[i], target_eigenvector[min_cost_index])
        #down
        initial_dataSet[i]-=[0, 0.01] #reset from the north operation
        initial_dataSet[i]+=[0,-0.01]
        initial_eigenvector, initial_eigenvalue = cal_Elipse(initial_dataSet)
        S_cost = cal_cost(initial_dataSet[i], target_dataSet[min_cost_index], initial_eigenvalue[i], target_eigenvalue[min_cost_index],initial_eigenvector[i], target_eigenvector[min_cost_index])
        #right
        initial_dataSet[i]-=[0,-0.01] #reset from the south operation
        initial_dataSet[i]+=[0.01, 0]
        initial_eigenvector, initial_eigenvalue = cal_Elipse(initial_dataSet)
        E_cost = cal_cost(initial_dataSet[i], target_dataSet[min_cost_index], initial_eigenvalue[i], target_eigenvalue[min_cost_index],initial_eigenvector[i], target_eigenvector[min_cost_index])
        #left
        initial_dataSet[i]-=[0.01, 0] #reset from the east operation
        initial_dataSet[i]+=[-0.01,0]
        initial_eigenvector, initial_eigenvalue = cal_Elipse(initial_dataSet)
        W_cost = cal_cost(initial_dataSet[i], target_dataSet[min_cost_index], initial_eigenvalue[i], target_eigenvalue[min_cost_index],initial_eigenvector[i], target_eigenvector[min_cost_index])
        initial_dataSet[i]-=[-0.01, 0] #reset from the west operation
        cost_dict={'N':N_cost,'S':S_cost,'E':E_cost,'W':W_cost}
        best_direction=min(cost_dict, key=cost_dict.get)
        if best_direction=='N':
            initial_dataSet[i]+=[0,0.01]
        if best_direction=='S':
            initial_dataSet[i]+=[0,-0.01]
        if best_direction=='E':
            initial_dataSet[i]+=[0.01,0]
        if best_direction=='W':
            initial_dataSet[i]+=[-0.01,0]
    #after all nodes moved, calculate eigenvalue and eigenvector of initial_dataset, and report coordinate to check whether to stop
    print initial_dataSet
    initial_eigenvector,initial_eigenvalue=cal_Elipse(initial_dataSet)
    if np.array_equal(initial_eigenvector,target_eigenvector) and np.array_equal(initial_eigenvalue,target_eigenvalue):
        flag=True

print 'after global termination: ', initial_dataSet
