import pandas as pd
import numpy as np
import PH
import gudhi


Pre = './s645/'


def get_info(mut):
    chainid = mut[0]
    wildtype = mut[2]
    mutanttype = mut[-1]
    residueid = mut[3:-1]
    return chainid, wildtype, mutanttype, residueid



def distance_of_two_points(p1,p2):
    temp = pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2) + pow(p1[3]-p2[3],2)
    res = pow(temp,0.5)
    return res


def get_alpha_one_skeleton(point,filtration):
    t = point.shape
    if len(t)==1:
        return [0],[]
    alpha_complex = gudhi.AlphaComplex(points=point)
    simplex_tree = alpha_complex.create_simplex_tree()
    V = []
    E = []
    c = 0
    
    for item in simplex_tree.get_filtration():
        if item[1]<=filtration:
            if len(item[0])==2:
                m = min(item[0][0],item[0][1])
                M = max(item[0][0],item[0][1])
                temp = [ m, M, item[1] ]
                E.append(temp)
                
        else:
            break
    for i in range(t[0]):
        V.append(i)
    return V,E
    
    

def intersect2(a,b):
    return list(set(a)&set(b))  
def intersect3(a,b,c):
    return list(set(a)&set(b)&set(c))

def get_hom_complex_from_graph_new(V,E):
    res = []
    index = 0
    for i in V:
        temp = [index, 0, 0, i]
        index = index + 1
        res.append(temp)
    if len(E)<=1:
        return res
    
    
    # create filtered hom complex
    t = len(V)
    neighbor = []
    
    
    # set neighbor as empty
    for i in range(t):
        neighbor.append([])
    
    
    # add higher dimensional simplices
    for i in range(len(E)):
        one = E[i][0]
        two = E[i][1]
        dis = E[i][2]
        
        # one
        
        # edge
        for ii in range(len(neighbor[one])):
                one2 = neighbor[one][ii]
                m = min(two,one2)
                M = max(two,one2)
                Inter = intersect2(neighbor[m],neighbor[M])
                if Inter==[]:
                    temp = [index, dis, 1, m, M]
                    res.append(temp)
                    index = index + 1
            
        # triangle
        for ii in range(len(neighbor[one])):
                one2 = neighbor[one][ii]
                for jj in range(ii+1,len(neighbor[one])):
                    one3 =  neighbor[one][jj]
                    m = min(two,one2,one3)
                    M = max(two,one2,one3)
                    mid = two + one2 + one3 - m - M
                    Inter = intersect3(neighbor[m],neighbor[mid],neighbor[M])
                    if Inter==[]:
                        temp = [index, dis, 2, m, mid, M]
                        res.append(temp)
                        index = index + 1
            
            #neighbor[one].append(two)
            
        
        # two
        
            
        # edge
        for ii in range(len(neighbor[two])):
                one2 = neighbor[two][ii]
                m = min(one,one2)
                M = max(one,one2)
                Inter = intersect2(neighbor[m],neighbor[M])
                if Inter==[]:
                    temp = [index, dis, 1, m, M]
                    res.append(temp)
                    index = index + 1
            
        # triangle
        for ii in range(len(neighbor[two])):
                one2 = neighbor[two][ii]
                for jj in range(ii+1,len(neighbor[two])):
                    one3 =  neighbor[two][jj]
                    m = min(one,one2,one3)
                    M = max(one,one2,one3)
                    mid = one + one2 + one3 - m - M
                    Inter = intersect3(neighbor[m],neighbor[mid],neighbor[M])
                    if Inter==[]:
                        temp = [index, dis, 2, m, mid, M]
                        res.append(temp)
                        index = index + 1
        neighbor[one].append(two)    
        neighbor[two].append(one)
    return res
        
        

def get_one_skeleton_hom_complex_from_graph(V,E):
    res = []
    index = 0
    for i in V:
        temp = [index, 0, 0, i]
        index = index + 1
        res.append(temp)
    if len(E)<=1:
        return res
    
    
    # create filtered hom complex
    t = len(V)
    neighbor = []
    is_edge_matrix = np.zeros((t,t))
    #is_triangle_matrix = np.zeros((t,t,t))
    #triangles = []
    
    # set neighbor as empty
    for i in range(t):
        neighbor.append([])
    
    
    # add higher dimensional simplices
    for i in range(len(E)):
        one = E[i][0]
        two = E[i][1]
        dis = E[i][2]
        
        # one
        number_had = len(neighbor[one])
        if number_had==0:
            neighbor[one].append(two)
        else:
            
            # edge
            for ii in range(len(neighbor[one])):
                one2 = neighbor[one][ii]
                m = min(two,one2)
                M = max(two,one2)
                if is_edge_matrix[m][M]==0:
                    is_edge_matrix[m][M] = 1
                    temp = [index, dis, 1, m, M]
                    res.append(temp)
                    index = index + 1
            '''
            # triangle
            for ii in range(len(neighbor[one])):
                one2 = neighbor[one][ii]
                for jj in range(ii+1,len(neighbor[one])):
                    one3 =  neighbor[one][jj]
                    m = min(two,one2,one3)
                    M = max(two,one2,one3)
                    mid = two + one2 + one3 - m - M
                    if is_triangle_matrix[m][mid][M]==0:
                    #if [ m,mid,M ] not in triangles:
                        #triangles.append([ m,mid,M ])
                        is_triangle_matrix[m][mid][M] = 1
                        temp = [index, dis, 2, m, mid, M]
                        res.append(temp)
                        index = index + 1
            '''
            neighbor[one].append(two)
            
        
        # two
        number_had = len(neighbor[two])
        if number_had==0:
            neighbor[two].append(one)
        else:
            
            # edge
            for ii in range(len(neighbor[two])):
                one2 = neighbor[two][ii]
                m = min(one,one2)
                M = max(one,one2)
                if is_edge_matrix[m][M]==0:
                    is_edge_matrix[m][M] = 1
                    temp = [index, dis, 1, m, M]
                    res.append(temp)
                    index = index + 1
            '''
            # triangle
            for ii in range(len(neighbor[two])):
                one2 = neighbor[two][ii]
                for jj in range(ii+1,len(neighbor[two])):
                    one3 =  neighbor[two][jj]
                    m = min(one,one2,one3)
                    M = max(one,one2,one3)
                    mid = one + one2 + one3 - m - M
                    if is_triangle_matrix[m][mid][M]==0:
                    #if [ m,mid,M ] not in triangles:
                        #triangles.append( [m,mid,M] )
                        is_triangle_matrix[m][mid][M] = 1
                        temp = [index, dis, 2, m, mid, M]
                        res.append(temp)
                        index = index + 1
            '''
            neighbor[two].append(one)
    return res
    
    
    

def get_hom_complex_2D_from_bipartite_new(point_cloud,filtration):
    t = point_cloud.shape
    if t[0]<=1:
        return [],[]
    p1_number = 1
    p2_number = 0
    for i in range(t[0]-1):
        if point_cloud[i][0] == point_cloud[i+1][0]:
            p1_number = p1_number + 1
        else:
            p2_number = t[0] - i - 1
            break
    #print(p1_number,p2_number)
    if p1_number+p2_number!=t[0]:
        print('error99999999999999999999')
        return
    if p1_number==t[0]:
        return [],[]
    dis_list = []
    for i in range(p1_number):
        for j in range(p1_number,t[0]):
            dis = distance_of_two_points(point_cloud[i], point_cloud[j])
            if dis<=filtration:
                dis_list.append([ i,j,dis ])
    dis_list = sorted(dis_list,key=lambda x:(x[2]))
    #print(len(dis_list))
    
    
    
    simplices_p = []
    count_p = 0
    simplices_l = []
    count_l = 0
    neighbor = []
    
    for i in range(t[0]):
        neighbor.append([])
    
    # vertices
    for i in range(p1_number):
        temp = [ count_p, 0, 0, i ]
        simplices_p.append(temp)
        count_p = count_p + 1
    
    for i in range(p2_number):
        temp = [ count_l, 0, 0, i+p1_number ]
        simplices_l.append(temp)
        count_l = count_l + 1
    
    for i in range(len(dis_list)):
        p = dis_list[i][0]
        l = dis_list[i][1]
        fil = dis_list[i][2]
        
        
        # first one
        # edges
        t = len(neighbor[l])
        for tt in range(t):
                one = neighbor[l][tt]
                m = min(one,p)
                M = max(one,p)
                Inter = intersect2(neighbor[m],neighbor[M])
                if Inter==[]:
                    temp = [ count_p, fil, 1, m, M ]
                    simplices_p.append(temp)
                    count_p = count_p + 1
        # triangles
        for ii in range(t):
            one = neighbor[l][ii]
            for jj in range(ii+1,t):
                    two = neighbor[l][jj]
                    m = min(one,two,p)
                    M = max(one,two,p)
                    mid = one + two + p - m - M
                    Inter = intersect3(neighbor[m],neighbor[mid],neighbor[M])
                    if Inter==[]:
                        temp = [ count_p, fil, 2, m, mid, M ]
                        simplices_p.append(temp)
                        count_p = count_p + 1
            
            #L_neighbour[l_index].append(p)
   
    
    
        # second one
        # edges
        t = len(neighbor[p])
        for tt in range(t):
                one = neighbor[p][tt]
                m = min(one,l)
                M = max(one,l)
                Inter =  intersect2(neighbor[m],neighbor[M])
                if Inter==[]:
                    temp = [ count_l, fil, 1, m, M ]
                    simplices_l.append(temp)
                    count_l = count_l + 1
        # triangles
        for ii in range(t):
            one = neighbor[p][ii]
            for jj in range(ii+1,t):
                    two = neighbor[p][jj]
                    m = min(one,two,l)
                    M = max(one,two,l)
                    mid = one + two + l - m - M
                    Inter = intersect3(neighbor[m],neighbor[mid],neighbor[M])
                    if Inter==[]:
                        temp = [ count_l, fil, 2, m, mid, M ]
                        simplices_l.append(temp)
                        count_l = count_l + 1
        neighbor[p].append(l)
        neighbor[l].append(p)
    
    return simplices_p,simplices_l

def zero_homology_of_a_complex_to_file(complex1,filename):
    res = PH.get_persistence(complex1)
    if len(res)==0:
        zero_bar = np.array([])
        np.savetxt(filename,zero_bar,delimiter=',')
        return
    diag = res['diagrams']
    zero_bar = np.array(diag[0])
    zero = zero_bar.shape
    zero_number = zero[0]
    c = 0
    for j in range(zero_number):
        if zero_bar[j][1]==-1:
            zero_bar[j][1] = float('inf')
            c = c + 1
    #print(c,zero,zero_bar)
    np.savetxt(filename,zero_bar,delimiter=',')

def one_homology_of_a_complex_to_file(complex1,filename):
    res = PH.get_persistence(complex1)
    if len(res)==0:
        one_bar = np.array([])
        np.savetxt(filename,one_bar,delimiter=',')
        return
    diag = res['diagrams']
    one_bar = np.array(diag[1])
    one = one_bar.shape
    one_number = one[0]
    c = 0
    for j in range(one_number):
        if one_bar[j][1]==-1:
            one_bar[j][1] = float('inf')
            c = c + 1
    #print(c,zero,zero_bar)
    np.savetxt(filename,one_bar,delimiter=',')
    

def alpha_zero_homology_to_file(start,end,filtration):
    filename = Pre + 'AB-Bind_S645.xlsx'
    df1 = pd.read_excel(filename)
    t1 = df1.shape
    
    for i in range(start,end):
        print(i)
        pdbid = df1.iloc[i,0]
        chainid, wildtype, mutanttype, residueid = get_info(df1.iloc[i,1])
        folder = pdbid + '_' + chainid + '_' + wildtype + '_' + residueid + '_' + mutanttype
        print(folder)
        for typ in [ 'binding_mut', 'binding_wild', 'mutation_mut', 'mutation_wild' ]:
        #for typ in ['']
            for atom in ['C','N','O','CN','CO','NO']:
                filename1 = Pre + 's645_point_10/' + folder + '_' + typ + '_' + atom + '1.txt'
                filename2 = Pre + 's645_point_10/' + folder + '_' + typ + '_' + atom + '2.txt'
                point1 = np.loadtxt(filename1,delimiter=',')
                point2 = np.loadtxt(filename2,delimiter=',')
                #print(point1.shape,point2.shape)
                
                V1,E1 = get_alpha_one_skeleton(point1,filtration)
                #complex1 = get_hom_complex_from_graph(V1,E1)
                complex1 = get_one_skeleton_hom_complex_from_graph(V1,E1)
                #print(atom,'number1:',len(V1),len(E1),len(complex1))
                filepath1 =Pre + 's645_alpha_zero_homology_' + str(filtration) + '/' + folder + '_' + typ + '_' + atom + '1_zero_homology_' + str(filtration) + '.csv'
                zero_homology_of_a_complex_to_file(complex1,filepath1)
                
                
                
                V2,E2 = get_alpha_one_skeleton(point2,filtration)
                #complex2 = get_hom_complex_from_graph(V2,E2)
                complex2 = get_one_skeleton_hom_complex_from_graph(V2,E2)
                #print(atom,'number2:',len(V2),len(E2),len(complex2))
                filepath2 =Pre + 's645_alpha_zero_homology_' + str(filtration) + '/' + folder + '_' + typ + '_' + atom + '2_zero_homology_' + str(filtration) + '.csv'
                zero_homology_of_a_complex_to_file(complex2,filepath2)
                
        

                
def get_number(bar,left,right):
    t = bar.shape
    if (len(t)==1):
        return 0
    res = 0
    for i in range(t[0]):
        if (bar[i][0]<=left)&(bar[i][1]>=right):
            res = res + 1
    return res
                
def label_to_file():
    filename = Pre + 'AB-Bind_S645.xlsx'
    df1 = pd.read_excel(filename)
    t1 = df1.shape
    row = 645
    column = 1
    feature_matrix = np.zeros((row,column))
    
    for i in range(row):
        #print(i)
        count = 0
        feature_matrix[i,0] = df1.iloc[i,3]
    filename = Pre + 'feature/label.csv'
    np.savetxt(filename,feature_matrix,delimiter=',')
    
        
       
def alpha_h0_feature_to_file(start,end,filtration,grid_size):
    filename = Pre + 'AB-Bind_S645.xlsx'
    df1 = pd.read_excel(filename)
    t1 = df1.shape
    row = 645
    grid_number = int(filtration/grid_size)
    column = 6 * 4 * 2 * grid_number
    feature_matrix = np.zeros((row,column))
    
    for i in range(start,end):
        print(i)
        count = 0
        pdbid = df1.iloc[i,0]
        chainid, wildtype, mutanttype, residueid = get_info(df1.iloc[i,1])
        folder = pdbid + '_' + chainid + '_' + wildtype + '_' + residueid + '_' + mutanttype
        for j in range(grid_number):
            left = j * grid_size
            right = (j+1) * grid_size
            for atom in ['C','N','O','CN','CO','NO']:
                # binding_mutate
                filename1 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_binding_mut_' + atom + '1_zero_homology_' + str(filtration) + '.csv'
                filename2 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_binding_mut_' + atom + '2_zero_homology_' + str(filtration) + '.csv'
                # binding_wild
                filename3 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_binding_wild_' + atom + '1_zero_homology_' + str(filtration) + '.csv'
                filename4 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_binding_wild_' + atom + '2_zero_homology_' + str(filtration) + '.csv'
                # mutation_mut
                filename5 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_mutation_mut_' + atom + '1_zero_homology_' + str(filtration) + '.csv'
                filename6 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_mutation_mut_' + atom + '2_zero_homology_' + str(filtration) + '.csv'
                # mutation_wild
                filename7 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_mutation_wild_' + atom + '1_zero_homology_' + str(filtration) + '.csv'
                filename8 = Pre + 'S645_alpha_zero_homology_' + str(filtration) + '/' + folder \
                    + '_mutation_wild_' + atom + '2_zero_homology_' + str(filtration) + '.csv'
                
                bar1 = np.loadtxt(filename1,delimiter=',')
                bar2 = np.loadtxt(filename2,delimiter=',')
                bar3 = np.loadtxt(filename3,delimiter=',')
                bar4 = np.loadtxt(filename4,delimiter=',')
                bar5 = np.loadtxt(filename5,delimiter=',')
                bar6 = np.loadtxt(filename6,delimiter=',')
                bar7 = np.loadtxt(filename7,delimiter=',')
                bar8 = np.loadtxt(filename8,delimiter=',')
                
                
                feature_matrix[i,count] = get_number(bar1,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar2,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar3,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar4,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar5,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar6,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar7,left,right)
                count = count + 1
                feature_matrix[i,count] = get_number(bar8,left,right)
                count = count + 1
                    
                    
    #filename = Pre + 'feature/alpha_zero_homology_feature' + str(grid_size) + '.csv'
    #np.savetxt(filename,feature_matrix,delimiter=',')
    return feature_matrix
    
    

def write_edge_triangle_number_of_a_complex_to_file(complex1,filename,filtration,grid_size):
    # output format [edge_number,triangle_number,euler_number]
    complex_length = len(complex1)
    item_number = int(filtration/grid_size)
    res_equ = []
    res_leqs = []
    for i in range(item_number):
        res_equ.append( [0,0,0] )
        res_leqs.append( [0,0,0] )
    for simplex in complex1:
        fil = simplex[1]
        dim = simplex[2]
        index = int(fil/grid_size)
        res_equ[index][dim] = res_equ[index][dim] + 1
        res_leqs[index][dim] = res_leqs[index][dim] + 1
    for i in range(1,len(res_leqs)):
        res_leqs[i][0] = res_leqs[i][0] + res_leqs[i-1][0]
        res_leqs[i][1] = res_leqs[i][1] + res_leqs[i-1][1]
        res_leqs[i][2] = res_leqs[i][2] + res_leqs[i-1][2]
    res = [ res_equ,res_leqs ]
    f = open(filename,'w')
    f.write(str(res))
    f.close()
    
    


def alpha_edge_triangle_number_to_file(start,end,filtration,grid_size):
    
    filename = Pre + 'AB-Bind_S645.xlsx'
    df1 = pd.read_excel(filename)
    t1 = df1.shape
    
    for i in range(start,end):
        print(i)
        pdbid = df1.iloc[i,0]
        chainid, wildtype, mutanttype, residueid = get_info(df1.iloc[i,1])
        folder = pdbid + '_' + chainid + '_' + wildtype + '_' + residueid + '_' + mutanttype
        for typ in [ 'binding_mut', 'binding_wild', 'mutation_mut', 'mutation_wild' ]:
        
            for atom in ['C','N','O','CN','CO','NO']:
                filename1 = Pre + 's645_point_10/' + folder + '_' + typ + '_' + atom + '1.txt'
                filename2 = Pre + 's645_point_10/' + folder + '_' + typ + '_' + atom + '2.txt'
                point1 = np.loadtxt(filename1,delimiter=',')
                point2 = np.loadtxt(filename2,delimiter=',')
                #print(point1.shape,point2.shape)
                
                V1,E1 = get_alpha_one_skeleton(point1,filtration)
                complex1 = get_hom_complex_from_graph_new(V1,E1)
                #complex1 = get_one_skeleton_hom_complex_from_graph(V1,E1)
                #print(atom,'number1:',len(V1),len(E1),len(complex1))
                filepath1 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder + '_' + typ + '_' + atom + '1.txt'
                write_edge_triangle_number_of_a_complex_to_file(complex1,filepath1,filtration,grid_size)
                
                
                V2,E2 = get_alpha_one_skeleton(point2,filtration)
                complex2 = get_hom_complex_from_graph_new(V2,E2)
                #complex2 = get_one_skeleton_hom_complex_from_graph(V2,E2)
                #print(atom,'number2:',len(V2),len(E2),len(complex2))
                filepath2 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder + '_' + typ + '_' + atom + '2.txt'
                write_edge_triangle_number_of_a_complex_to_file(complex2,filepath2,filtration,grid_size)
                
                

def get_simplex_number(number,j,typ):
    if typ=='edge':
        return number[j][1]
    elif typ=='triangle':
        return number[j][2]
    elif typ=='euler':
        return number[j][0] - number[j][1] + number[j][2]

       
def alpha_simplex_feature_to_file(typ,start,end,filtration,grid_size):
    filename = Pre + 'AB-Bind_S645.xlsx'
    df1 = pd.read_excel(filename)
    t1 = df1.shape
    row = 645
    grid_number = int(filtration/grid_size)
    column = 6 * 4 * 2 * grid_number
    feature_matrix = np.zeros((row,column))
    
    for i in range(start,end):
        print(i)
        count = 0
        pdbid = df1.iloc[i,0]
        chainid, wildtype, mutanttype, residueid = get_info(df1.iloc[i,1])
        folder = pdbid + '_' + chainid + '_' + wildtype + '_' + residueid + '_' + mutanttype
        for j in range(grid_number):
            #left = j * grid_size
            #right = (j+1) * grid_size
            for atom in ['C','N','O','CN','CO','NO']:
                # binding_mutate
                filename1 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_binding_mut_' + atom + '1.txt'
                filename2 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_binding_mut_' + atom + '2.txt'
                # binding_wild
                filename3 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_binding_wild_' + atom + '1.txt'
                filename4 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_binding_wild_' + atom + '2.txt'
                # mutation_mut
                filename5 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_mutation_mut_' + atom + '1.txt'
                filename6 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_mutation_mut_' + atom + '2.txt'
                # mutation_wild
                filename7 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_mutation_wild_' + atom + '1.txt'
                filename8 = Pre + 'alpha_complex_number_' + str(filtration) + '/' + folder \
                    + '_mutation_wild_' + atom + '2.txt'
                
                f = open(filename1)
                pre_number1 = f.read()
                number1 = eval(pre_number1)
                f.close()
                f = open(filename2)
                pre_number2 = f.read()
                number2 = eval(pre_number2)
                f.close()
                f = open(filename3)
                pre_number3 = f.read()
                number3 = eval(pre_number3)
                f.close()
                f = open(filename4)
                pre_number4 = f.read()
                number4 = eval(pre_number4)
                f.close()
                f = open(filename5)
                pre_number5 = f.read()
                number5 = eval(pre_number5)
                f.close()
                f = open(filename6)
                pre_number6 = f.read()
                number6 = eval(pre_number6)
                f.close()
                f = open(filename7)
                pre_number7 = f.read()
                number7 = eval(pre_number7)
                f.close()
                f = open(filename8)
                pre_number8 = f.read()
                number8 = eval(pre_number8)
                f.close()
                
                
                
                
                feature_matrix[i,count] = get_simplex_number(number1[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number2[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number3[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number4[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number5[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number6[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number7[1],j,typ)
                count = count + 1
                feature_matrix[i,count] = get_simplex_number(number8[1],j,typ)
                count = count + 1
                    
                    
    #filename = Pre + 'feature/alpha_' + str(typ) + '_feature' + str(grid_size) + '.csv'
    #np.savetxt(filename,feature_matrix,delimiter=',')
    return feature_matrix
    
def get_topological_feature():
    d1 = alpha_h0_feature_to_file(0,645,5,0.1)
    d2 = alpha_simplex_feature_to_file('euler',0,645,5,0.1)
    d = np.hstack((d1,d2))
    filename = Pre + 'feature/alpha_h0_euler.csv'
    np.savetxt(filename,d,delimiter=',')
    
    
    
    
    
def main():
    alpha_zero_homology_to_file(0,645,5)
    alpha_edge_triangle_number_to_file(0,645,5,0.1)
    get_topological_feature()
    label_to_file()
    
    

    
    
main()
    

