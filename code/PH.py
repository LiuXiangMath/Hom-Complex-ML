# -*- coding: utf-8 -*-

import numpy as np


def get_point_index(point,simplices):
    for i in range(len(simplices)):
        if simplices[i][2]==0:
            if simplices[i][3]==point:
                return i

def get_edge_index(edge,simplices):
    for i in range(len(simplices)):
        if simplices[i][2]==1:
            if (edge[0]==simplices[i][3]) and (edge[1]==simplices[i][4]):
                return i

def get_triangle_index(triangle,simplices):
    for i in range(len(simplices)):
        if simplices[i][2]==2:
            if (triangle[0]==simplices[i][3]) and (triangle[1]==simplices[i][4]) and (triangle[2]==simplices[i][5]):
                return i

def get_boundary(simplex,simplices):
    
    if simplex[2]==0:
        return []
    elif simplex[2]==1:
        one = simplex[3]
        two = simplex[4]
        index1 = get_point_index(one,simplices)
        index2 = get_point_index(two,simplices)
        if index1<index2:
            return [ [simplex[3]], [simplex[4]] ]
        else:
            return [ [simplex[4]], [simplex[3]] ]
    elif simplex[2]==2:
        one_two = [ simplex[3], simplex[4] ]
        one_three = [simplex[3], simplex[5] ]
        two_three = [simplex[4], simplex[5] ]
        index1 = get_edge_index(one_two,simplices)
        index2 = get_edge_index(one_three,simplices)
        index3 = get_edge_index(two_three,simplices)
        
        temp = [ [one_two, index1], [one_three,index2], [two_three,index3] ]
        temp2 = sorted(temp,key=lambda x:(x[1]))
        res = [temp2[0][0], temp2[1][0], temp2[2][0]]
        return res
    elif simplex[2]==3:
        one_two_three = [ simplex[3], simplex[4], simplex[5] ]
        one_two_four = [simplex[3], simplex[4], simplex[6]]
        one_three_four = [ simplex[3], simplex[5], simplex[6] ]
        two_three_four = [ simplex[4], simplex[5], simplex[6] ]
        index1 = get_triangle_index(one_two_three,simplices)
        index2 = get_triangle_index(one_two_four,simplices)
        index3 = get_triangle_index(one_three_four,simplices)
        index4 = get_triangle_index(two_three_four,simplices)
        
        temp = [ [one_two_three, index1], [one_two_four, index2], [one_three_four, index3], [two_three_four, index4] ]
        
        temp2 = sorted(temp,key=lambda x:(x[1]))
        res = [temp2[0][0], temp2[1][0], temp2[2][0], temp2[3][0]]
        return res


def get_index(simplex,simplices,m):
    t = len(simplex)
    for i in range(m,0,-1):
        if simplices[i][2]==(t-1):
            if simplices[i][3::]==simplex:
                return i


def get_largest_positive(boundary,simplices,boundary_ls,m):
    for i in range(len(boundary)-1,-1,-1):
        index = get_index(boundary[i],simplices,m)
        if boundary_ls[index]!=[-1]:
            return index
        

def add_zero_boundary(b1,b2,simplices):
    t1 = len(b1)
    t2 = len(b2)
    res = []
    for i in range(t1):
        res.append(b1[i][0])
    for i in range(t2):
        have = 0
        for j in range(t1):
            if b2[i][0]==b1[j][0]:
                res.remove(b2[i][0])
                have = 1
                break
        if have==0:
            res.append(b2[i][0])
    temp = []
    for i in range(len(res)):
        index = get_point_index(res[i],simplices)
        temp.append( [ res[i], index ] )
    temp2 = sorted(temp,key=lambda x:(x[1]))
    temp3 = []
    for i in range(len(temp2)):
        temp3.append( [temp2[i][0]] )
    return temp3


def add_one_boundary(b1,b2,simplices):
    t1 = len(b1)
    t2 = len(b2)
    res = []
    for i in range(t1):
        res.append( [ b1[i][0], b1[i][1] ] )
    for i in range(t2):
        have = 0
        for j in range(t1):
            if (b2[i][0]==b1[j][0]) and (b2[i][1]==b1[j][1]):
                res.remove( [ b2[i][0], b2[i][1] ] )
                have = 1
                break
        if have==0:
            res.append( [ b2[i][0], b2[i][1] ] )
    temp = []
    for i in range(len(res)):
        index = get_edge_index(res[i],simplices)
        temp.append( [ res[i], index ] )
    temp2 = sorted(temp,key=lambda x:(x[1]))
    temp3 = []
    for i in range(len(temp2)):
        temp3.append( temp2[i][0] )
    return temp3




def add_two_boundary(b1,b2,simplices):
    t1 = len(b1)
    t2 = len(b2)
    res = []
    for i in range(t1):
        res.append( [ b1[i][0], b1[i][1], b1[i][2] ] )
    for i in range(t2):
        have = 0
        for j in range(t1):
            if (b2[i][0]==b1[j][0]) and (b2[i][1]==b1[j][1]) and (b2[i][2]==b1[j][2]):
                res.remove( [ b2[i][0], b2[i][1], b2[i][2] ] )
                have = 1
                break
        if have==0:
            res.append( [ b2[i][0], b2[i][1], b2[i][2] ] )
    temp = []
    for i in range(len(res)):
        index = get_triangle_index(res[i],simplices)
        temp.append( [ res[i], index ] )
    temp2 = sorted(temp,key=lambda x:(x[1]))
    temp3 = []
    for i in range(len(temp2)):
        temp3.append( temp2[i][0] )
    return temp3


    
    
    

def add_boundary(b1,b2,simplices):
    t = len(b1[0])
    res = []
    if t==1:
        res = add_zero_boundary(b1,b2,simplices)
    elif t==2:
        res = add_one_boundary(b1,b2,simplices)
    elif t==3:
        res = add_two_boundary(b1,b2,simplices)
    return res
    

def get_persistence(simplices):
    # coefficient Z/2
    number = len(simplices)
    if number==0:
        return []
    boundary_ls = [] 
    I = [] # store the forever-persistent bar
    P = [] # store the pair [birth, death]
    for m in range(number):
        #print(m,number)
        boundary = get_boundary(simplices[m],simplices)
        if len(boundary)==0:
            boundary_ls.append([])
            I.append(m)
            continue
        largest_index = get_largest_positive(boundary,simplices,boundary_ls,m-1)
        
        while (boundary!=[]) and (boundary_ls[largest_index]!=[]):
            boundary = add_boundary(boundary,boundary_ls[largest_index],simplices)
            if boundary==[]:
                break
            largest_index = get_largest_positive(boundary,simplices,boundary_ls,largest_index)
            
        if boundary==[]: # positive
            boundary_ls.append([])
            I.append(m)
        else: # negative
            boundary_ls[largest_index] = boundary
            boundary_ls.append([-1])
            P.append( [largest_index,m] )
            for tt in range(len(I)):
                if I[tt] == largest_index:
                    del I[tt]
                    break
    
    zero_bar = []
    one_bar = []
    two_bar = []
    
    
    for i in I:
        if simplices[i][2]==0:   
            f = simplices[i][1] + 0.0
            zero_bar.append( [f,-1] )
        elif simplices[i][2]==1: 
            f = simplices[i][1] + 0.0
            one_bar.append( [f,-1] )
        elif simplices[i][2]==2:
            f = simplices[i][1] + 0.0
            two_bar.append( [f,-1] )
            
        # can add higher dimensional information if you need
    for p in P:
        birth = p[0]
        death = p[1]
        if simplices[death][1]>simplices[birth][1]:
            bf = simplices[birth][1] + 0.0
            df = simplices[death][1]
            if simplices[birth][2]==0:
                zero_bar.append( [bf,df] )
            elif simplices[birth][2]==1:
                one_bar.append( [bf,df] )
            elif simplices[birth][2]==2:
                two_bar.append( [bf,df] )
        
            # can add higher dimensional information if you need
   
    result = {'diagrams':[zero_bar,one_bar,two_bar]}
    return result


    
  