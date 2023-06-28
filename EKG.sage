from time import *

##############
### CASE 1 ###
##############

"""
n=16
m = ceil(n^(1/3)*4*(8/3)^(n/3))
l=4
Poly = [0,11,13,14]
p = 2^16-15
a=5
b=5

R=GF(p)
E=EllipticCurve(R,[a,b])     #y^2 = x^3+5x+5 (65111 points) a4=5 a6=5
q = E.order()
S=GF(q)
P0 = E(30507,21693,1)

"""

"""
n=16
m = ceil(n^(1/3)*4*(8/3)^(n/3))
l=4
Poly = [0,11,13,14]
p=2^40+15
a = 1
b = 14
R=GF(p)
E=EllipticCurve(R,[a,b])
q = E.order()
S=GF(q)
"""


def output(Q,l):
    xQ = ZZ(Q[0])
    return(xQ >> l)

def knapsack(n,E,Poly,m,l):
    #init
    u=[] #we keep the whole u even if do not need it
    w=[]
    for i in range(n):
        u.append(ZZ.random_element(0,2))
        w.append(E.random_element())
    #outputs
    Y = []
    for j in range(m):
        Q = sum([u[j+i]*w[i] for i in range(n)])
        Y.append(output(Q,l))
        u.append(sum([u[j+k]for k in Poly])%2)
    return u,w,Y
    #this generator is the same wheither we use w or -w

def check_knapsack(n,Poly,m,l,u,w):
    Y = []
    for j in range(m):
        Q = sum([u[j+i]*w[i] for i in range(n)])
        Y.append(output(Q,l))
        u.append(sum([u[j+k]for k in Poly])%2)
    return Y    

def find_triplet(n,u):
    e = 1/6
    m = len(u) - n
    list_triplet = []
    weight = n*(1/3+e)
    dicosmall = {}
    contredico = {}
    for i in range(m-n):
        ui = vector(ZZ,u[i:i+n])
        if sum(ui) <= weight:
            dicosmall[i] = ui
        contredico[tuple(ui)] = i
    for i in dicosmall:
        ui = dicosmall[i]
        for j in dicosmall:
            if j > i:
                uj = dicosmall[j]
                uk = ui+uj
                k = contredico.get(tuple(uk))
                if k!= None and j<k:
                    list_triplet.append([i,j,k])              
    return(list_triplet)

def find_systeme(n,E,a,b,u,Y,l,triplets):
    m = len(u)-n
    list_index=[]
    list_point=[]
    for [i,j,k]  in triplets:
        if i in list_index:
            continue
        if j in list_index:
            continue
        list_Qi = []
        list_Qj = []

        for rQi in range(2**l):
            xQi = Y[i]*2^l+rQi
            yQi = R(xQi^3+a*xQi+b)
            if yQi.is_square():
                yQi = sqrt(yQi) #we will choose this one an not another
                Qi=E(xQi,yQi,1)
                for rQj in range(2^l):
                    xQj = Y[j]*2^l +rQj
                    yQj = R(xQj^3+a*xQj+b)
                    if yQj.is_square():
                        yQj = sqrt(yQj)
                        Qj = E(xQj,yQj,1)

                        if output(Qi+Qj,l) == Y[k]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)

                        Qj = -Qj
                        if output(Qi+Qj,l) == Y[k]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)


        if len(list_Qi)==1 :
            list_index.append(i)
            list_point.append(list_Qi[0])
        #if len(list_Qi) >1 :
            #print('l too big')           
        if len(list_Qj)==1 :
            list_index.append(j)
            list_point.append(list_Qj[0]) 

    return(list_index,list_point) 

def attack(n,E,a,b,u,Y,l):
    t = n//2
    m=len(u)-n
    triplets = find_triplet(n,u)
    list_index,list_point = find_systeme(n,E,a,b,u,Y,l,triplets)
    if len(list_index) < n:
        return 

    M = matrix(S,n,n)
    for i in range(n):
        ii = list_index[i]
        for j in range(n):
            M[i,j] = u[ii+j]
    
    if M.rank() !=n :
        return

    while(M.rank() != n):
        list_index = list_index[1:]
        if len(list_index) < n:
            return 
        list_point = list_point[1:]
        M = matrix(S,n,n)
        for i in range(n):
            ii = list_index[i]
            for j in range(n):
                M[i,j] = u[ii+j]

    antiM = M^(-1)
    antiM = matrix(ZZ,antiM)
    
    for temp in range(2^(t-1)):
        sign = bin(temp)[2:]
        sign = '0'*(t-1-len(sign)) + sign
        
        B = [] #list points with signs 

        Qi = list_point[0]
        Qj = list_point[1]

        B.append(Qi)
        B.append(Qj)

        for i in range(1,t):
            exp = int(sign[i-1])
            Qi = (-1)^exp * list_point[2*i]
            Qj = (-1)^exp * list_point[2*i+1]
            B.append(Qi)
            B.append(Qj)

        #if the sign is correct we have M*w = B hence w = antiM*B
        
        P=[]
        for i in range(n):
            Pi = sum([antiM[i,j]*B[j] for j in range(n)])
            P.append(Pi)

        ubis = copy(u)[:n]
        Ybis = check_knapsack(n,Poly,100,l,ubis,P)
        if Ybis == Y[:100]:
            return P
            #sometime we return -P but it is also correct

def test_triplets(n,E,Poly,m,l,run):
    instances = []
    for i in range(run) :
        u,w,Y = knapsack(n,E,Poly,m,l)
        instances.append(u)
    T=time()
    res = 0
    for i in range(run):
        u=instances.pop()
        liste = find_triplet(n,u)
        res+= len(liste)
    T=time()-T 
    return(T/run,res/run)

def test_attack(n,E,a,b,Poly,m,l,run):
    instances = []
    for i in range(run) :
        u,w,Y = knapsack(n,E,Poly,m,l)
        instances.append([u,Y])
    T=time()
    res=0
    for i in range(run):
        [u,Y]=instances.pop()
        A = attack(n,E,a,b,u,Y,l)
        if A!= None:
            res+=1

    T=time()-T 
    return(T/run,res/run) 

def find_systeme_tricky(n,E,a,b,u,Y,l,r):
    m=len(u)-n
    list_index = []
    list_Q=[]

    #we need to force one of the Qi
    [i,j,k] = r[0]
    list_Qi = []
    list_Qj = []
    list_Qk = []
    for rQi in range(2**l):
        xQi = Y[i]*2^l+rQi
        yQi = R(xQi^3+a*xQi+b)
        if yQi.is_square():
            yQi = sqrt(yQi) #we will choose this one an not another
            Qi=E(xQi,yQi,1)
            for rQj in range(2^l):
                xQj = Y[j]*2^l +rQj
                yQj = R(xQj^3+a*xQj+b)
                if yQj.is_square():
                    yQj = sqrt(yQj)
                    Qj = E(xQj,yQj,1)
                    Qk = Qi+Qj
                    if output(Qk,l) == Y[k]:
                        if Qi not in list_Qi:
                            list_Qi.append(Qi)
                        if Qj not in list_Qj:
                            list_Qj.append(Qj)
                        if Qk not in list_Qk:
                            list_Qk.append(Qk)


                    Qj = -Qj
                    Qk = Qi+Qj
                    if output(Qk,l) == Y[k]:
                        if Qi not in list_Qi:
                            list_Qi.append(Qi)
                        if Qj not in list_Qj:
                            list_Qj.append(Qj)
                        if Qk not in list_Qk:
                            list_Qk.append(Qk)                        

    if len(list_Qi)==1:
        list_index.append(i)
        list_Q.append(list_Qi[0])                   
    if len(list_Qj)==1:
        list_index.append(j)
        list_Q.append(list_Qj[0])  
    if len(list_Qk)==1:
        list_index.append(k)
        list_Q.append(list_Qk[0])

    r.remove(r[0])
    t=0
    c=0
    while(len(r)!=0 and c<10000):
        c += 1
        if t==len(r) :
            t=0
        [i,j,k] = r[t]

        #all points in list_Q
        if i in list_index and j in list_index and k in list_index:
            r.remove(r[t])
        
        #two points are in list_Q
        elif i in list_index and j in list_index and k not in list_index:
            ii = list_index.index(i)
            Qi = list_Q[ii] 
            jj = list_index.index(j)
            Qj = list_Q[jj]
            list_index.append(k)
            list_Q.append(Qi+Qj)
            r.remove(r[t])

        elif i in list_index and j not in list_index and k  in list_index:
            ii = list_index.index(i)
            Qi = list_Q[ii] 
            kk = list_index.index(k)
            Qk = list_Q[kk]
            list_index.append(j)
            list_Q.append(Qk-Qi)
            r.remove(r[t])

        elif i not in list_index and j in list_index and k  in list_index:
            kk = list_index.index(k)
            Qk = list_Q[kk] 
            jj = list_index.index(j)
            Qj = list_Q[jj]
            list_index.append(i)
            list_Q.append(Qk-Qj)
            r.remove(r[t])

        #one point is in list_Q
        elif i in list_index and j not in list_index and k not in list_index:
            ii = list_index.index(i)
            Qi = list_Q[ii]
            list_Qj = []
            list_Qk = []
            for rQj in range(2^l):
                xQj = Y[j]*2^l +rQj
                yQj = R(xQj^3+a*xQj+b)
                if yQj.is_square():
                    yQj = sqrt(yQj)
                    Qj = E(xQj,yQj,1)
                    Qk = Qi+Qj
                    if output(Qk,l) == Y[k]:
                        if Qj not in list_Qj:
                            list_Qj.append(Qj)
                        if Qk not in list_Qk:
                            list_Qk.append(Qk)
                    Qj = -Qj
                    Qk = Qi+Qj
                    if output(Qk,l) == Y[k]:
                        if Qj not in list_Qj:
                            list_Qj.append(Qj)
                        if Qk not in list_Qk:
                            list_Qk.append(Qk)
            if len(list_Qj) == 1:
                list_index.append(j)
                list_Q.append(list_Qj[0])
                list_index.append(k)
                list_Q.append(list_Qk[0])
            r.remove(r[t])


        elif i not in list_index and j in list_index and k not in list_index:
            jj = list_index.index(j)
            Qj = list_Q[jj]
            list_Qi = []
            list_Qk = []
            for rQi in range(2^l):
                    xQi = Y[i]*2^l +rQi
                    yQi = R(xQi^3+a*xQi+b)
                    if yQi.is_square():
                        yQi = sqrt(yQi)
                        Qi = E(xQi,yQi,1)
                        Qk = Qi+Qj
                        if output(Qk,l) == Y[k]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qk not in list_Qk:
                                list_Qk.append(Qk)
                        Qi=-Qi
                        Qk = Qi+Qj
                        if output(Qk,l) == Y[k]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qk not in list_Qk:
                                list_Qk.append(Qk)
            if len(list_Qi) == 1:
                list_index.append(i)
                list_Q.append(list_Qi[0])
                list_index.append(k)
                list_Q.append(list_Qk[0]) 
            r.remove(r[t])

        elif i not in list_index and j not in list_index and k in list_index:
            kk = list_index.index(k)
            Qk = list_Q[kk]
            list_Qi = []
            list_Qj = []
            for rQi in range(2^l):
                    xQi = Y[i]*2^l +rQi
                    yQi = R(xQi^3+a*xQi+b)
                    if yQi.is_square():
                        yQi = sqrt(yQi)
                        Qi = E(xQi,yQi,1)
                        Qj = Qk-Qi
                        if output(Qj,l) == Y[j]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)
                        Qi=-Qi
                        Qj = Qk-Qi
                        if output(Qj,l) == Y[j]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)
            if len(list_Qi) == 1:
                list_index.append(i)
                list_Q.append(list_Qi[0])
                list_index.append(j)
                list_Q.append(list_Qj[0]) 
            r.remove(r[t])

        #no point are in the list
        elif i not in list_index and j not in list_index and k not in list_index:
            t=t+1

    return list_index,list_Q

def attack_tricky(n,E,a,b,u,Y,l):
    m=len(u)-n
    triplets = find_triplet(n,u)
    #print("nombre de triplets = ",len(triplets)) 
    list_index=[]
    while len(list_index) < 3*n:
        try : 
            list_index,list_point = find_systeme_tricky(n,E,a,b,u,Y,l,triplets)
        except :
            return
    #print("nombre de points =",len(list_index))
    M = matrix(S,n,n)
    B=[0]*n
    for i in range(n):
        while M.rank() != i+1 :
            ii = list_index[0]
            Qi = list_point[0]
            list_index.remove(ii)
            if len(list_index) == 0:
                return
            list_point.remove(Qi)
            B[i] = Qi
            for j in range(n):
                M[i,j] = u[ii+j]
    #print("rang de M = ",M.rank())

    antiM = M^(-1)
    antiM = matrix(ZZ,antiM)

    P=[]
    for i in range(n):
        Pi = sum([antiM[i,j]*B[j] for j in range(n)])
        P.append(Pi)

    ubis = copy(u)[:n]
    Ybis = check_knapsack(n,Poly,100,l,ubis,P)
    if Ybis == Y[:100]:
        return P
        


def test_attack_tricky(n,E,a,b,Poly,m,l,run):
    instances = []
    for i in range(run) :
        u,w,Y = knapsack(n,E,Poly,m,l)
        instances.append([u,Y])
    T=time()
    res=0
    for i in range(run):
        [u,Y]=instances.pop()
        A = attack_tricky(n,E,a,b,u,Y,l)
        if A!= None:
            res+=1

    T=time()-T 
    return(T/run,res/run) 



def refine_triplet(list_triplet):
    t=len(list_triplet)
    res=[]
    list_nb_index=[]
    list_index=[]
    for compteur in range(t):
        [i,j,k] = list_triplet[compteur]
        if i not in list_index:
            list_index.append(i)
            list_nb_index.append(1)
        else:
            ii = list_index.index(i)
            list_nb_index[ii] += 1
        if j not in list_index:
            list_index.append(j)
            list_nb_index.append(1)
        else:
            jj = list_index.index(j)
            list_nb_index[jj] += 1
        if k not in list_index:
            list_index.append(k)
            list_nb_index.append(1)
        else:
            kk = list_index.index(k)
            list_nb_index[kk] += 1



    for [i,j,k] in list_triplet:
        ii = list_index.index(i)
        jj = list_index.index(j)
        kk = list_index.index(k)
        if list_nb_index[ii] > 1 or list_nb_index[jj] > 1 or list_nb_index[kk] > 1:
            res.append([i,j,k])

    return res  

        


def attack1(u,Y,m):
    list_index = []
    list_Q = []
    lt1 = find_triplet_exo(u,m)
    r = refine_triplet(lt1)
    list_index,list_Q = find_systeme(u,Y,m,r)
    t = len(list_index)
    M = matrix(S,t,n)
    b = []
    for i in range(t):
        ii = list_index[i]
        Qi = list_Q[i]
        b.append(list_log.index(Qi))
        for j in range(n):
            M[i,j] = u[ii+j]
    b=vector(S,b)
    a = M.solve_right(b)
    return a

def attack2(u,Y,m,e):
    list_index = []
    list_Q = []
    lt2 = find_triplet_clever(u,m,e)
    r = refine_triplet(lt2)
    list_index,list_Q = find_systeme(u,Y,m,r)
    t = len(list_index)
    M = matrix(S,t,n)
    b = []
    for i in range(t):
        ii = list_index[i]
        Qi = list_Q[i]
        b.append(list_log.index(Qi))
        for j in range(n):
            M[i,j] = u[ii+j]
    b=vector(S,b)
    a = M.solve_right(b)
    return a



##############
### CASE 2 ###
##############

"""
n=16
l=4
Poly = [0,11,13,14]
p = 2**224*(2**32-1) + 2**192 + 2**96 - 1
a = p - 3
b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
R=GF(p)
E=EllipticCurve(R,[a,b])     
q = E.order()
S=GF(q)
P0 = E(0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296, 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5,1)

"""

def coppersmith(G,monomials,N,ell,m):
    """ retrieve a little common root of polynomials over a certain modulo
    
    Input:  G a list of polynomials
                monomials the list of monomials appearing in G (in a certain order)
                N the moduli
                ell the size of the unknown root
                m the multiplicity of the root in the polynomials of G (in our case it is the same for all polynomials)"""

    nb_poly = len(G)
    nb_mono = len(monomials)

    # Construction of the matrix described in section 2
    
    M = matrix(QQ,nb_mono+nb_poly)
    for i in range(nb_mono):
        if monomials[i] == 1: 
            M[i,i] = 1
        else:
            M[i,i] = (2^ell)^(-1*monomials[i].degree())

    for j in range(nb_poly):
        M[nb_mono+j,nb_mono+j] = N^m
        for i in range(nb_mono):
            M[i,nb_mono+j] = G[j].monomial_coefficient(monomials[i])
              
    # LLL
    M=M.LLL()
    
    # recuperation of the smallest vector with the correct form
    
    i=0
    
    while i < nb_mono+nb_poly and  abs(M[i,0]) == 0:
        i=i+1
    if i != nb_mono+nb_poly : 
        v = M[i,0]^-1*M[i]
        v_final = 2^ell*v

        return v_final
        
    return 0 #In the case the algorithm fails

def find_systeme_coppersmith(u,Y,m,l,r): 
    Z1 = var('Z1')
    Z2 = var('Z2')
    Z3 = var('Z3')
    Poly(Z1,Z2,Z3) = (Z1-Z2)^2*Z3^2 - 2*((Z1+Z2)*(Z1*Z2 + a) + 2*b)*Z3 + (Z1*Z2 - a)^2 - 4*b*(Z1+Z2)
    PP.<d1,d2,d3> = ZZ[]

    #r=list_triplet_refined
    list_index = []
    list_Q=[]

    #Finding the first point to force its sign
    [i,j,k] = r[0]

    list_Qi = []
    list_Qj = []
    list_Qk = []

    X1 = Y[i]*2**l + d1
    X2 = Y[j]*2**l + d2
    X3 = Y[k]*2**l + d3
    G = Poly(X1,X2,X3) 
    G=PP(G)%p
    monomials = G.monomials()
    monomials.reverse()
    # the orders of monomials in sage seems to be [1, d3, d2, d1,...]
    res = (coppersmith([G],monomials,p,l,1))[1:4]
    xQi = R(Y[i]*2**l+res[2])
    xQj = R(Y[j]*2**l+res[1])
    xQk = R(Y[k]*2**l+res[0])
    yQi = R(xQi^3+a*xQi+b)
    if yQi.is_square():
        yQi = sqrt(yQi) #we will choose this one an not yQi = -sqrt(yQi)
        Qi=E(xQi,yQi,1)
        yQj = R(xQj^3+a*xQj+b)
        if yQj.is_square():
            yQj = sqrt(yQj)
            Qj = E(xQj,yQj,1)
            Qk = Qi+Qj

            if output(Qk,l) == Y[k]:
                if Qi not in list_Qi:
                    list_Qi.append(Qi)
                if Qj not in list_Qj:
                    list_Qj.append(Qj)
                if Qk not in list_Qk :
                    list_Qk.append(Qk)

            Qj=E(xQj,-yQj,1)
            Qk = Qi+Qj
            if output(Qk,l) == Y[k]:
                if Qi not in list_Qi:
                    list_Qi.append(Qi)
                if Qj not in list_Qj:
                    list_Qj.append(Qj)
                if Qk not in list_Qk :
                    list_Qk.append(Qk)

    if len(list_Qi)==1:
        list_index.append(i)
        list_Q.append(list_Qi[0])                   
    if len(list_Qj)==1:
        list_index.append(j)
        list_Q.append(list_Qj[0])   
    if len(list_Qk)==1:
        list_index.append(k)
        list_Q.append(list_Qk[0])

    r.remove(r[0])
    t=0
    last_t=0
    while(len(r)!=0 and last_t<2):
        if t==len(r) :
            t=0
            last_t+=1
        [i,j,k] = r[t]
        if i in list_index:
            ii = list_index.index(i)
            Qi = list_Q[ii]
            if j in list_index:
                jj = list_index.index(j)
                Qj = list_Q[jj]
                if k in list_index:
                    kk = list_index.index(k)
                    Qk = list_Q[kk]
                    assert(Qi+Qj==Qk)
                else :
                    Qk = Qi+Qj
                    list_index.append(k)
                    list_Q.append(Qk)
                    r.remove(r[t])
            else : #j is not in list_index
                if k in list_index:
                    kk = list_index.index(k)
                    Qk = list_Q[kk]
                    Qj = Qk-Qi
                    list_index.append(j)
                    list_Q.append(Qj)
                    r.remove(r[t])
                else : 
                    list_Qj = []
                    list_Qk = []
                    xQi = ZZ(Qi[0])
                    X1 = xQi
                    X2 = Y[j]*2**l + d2
                    X3 = Y[k]*2**l + d3
                    G = Poly(X1,X2,X3)
                    G=PP(G)%p
                    monomials = G.monomials()
                    monomials.reverse()
                    #[1,d3,d2, ... ]
                    res = (coppersmith([G],monomials,p,l,1))[1:3]
                    xQj = Y[j]*2^l+R(res[1])
                    xQk = Y[k]*2^l+R(res[0])
                    yQj = R(xQj^3+a*xQj+b)
                    if yQj.is_square():
                        yQj = sqrt(yQj)
                        Qj = E(xQj,yQj,1)
                        Qk = Qi+Qj
                        if output(Qk,l) == Y[k]:
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)
                            if Qk not in list_Qk:
                                list_Qk.append(Qk)
                        yQj = -yQj
                        Qj=E(xQj,yQj,1)
                        Qk = Qi+Qj
                        if output(Qk,l) == Y[k]:
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)
                            if Qk not in list_Qk:
                                list_Qk.append(Qk)
                    if len(list_Qj) == 1:
                        list_index.append(j)
                        list_Q.append(list_Qj[0])
                        list_index.append(k)
                        list_Q.append(list_Qk[0])  
                        r.remove(r[t])                                   
        else : #i is not in list_index
            if j in list_index:
                jj = list_index.index(j)
                Qj = list_Q[jj]
                if k in list_index:
                    kk = list_index.index(k)
                    Qk = list_Q[kk]
                    Qi = Qk - Qj
                    list_index.append(i)
                    list_Q.append(Qi) 
                    r.remove(r[t])                
                else :
                    list_Qi = []
                    list_Qk = []
                    xQj = ZZ(Qj[0])
                    X1 = Y[i]*2**l + d1
                    X2 = xQj
                    X3 = Y[k]*2**l + d3
                    G = Poly(X1,X2,X3)
                    G=PP(G)%p
                    monomials = G.monomials()
                    monomials.reverse()
                    #[1,d3,d1,...]
                    res = (coppersmith([G],monomials,p,l,1))[1:3]
                    xQi = Y[i]*2^l+R(res[1])
                    xQk = Y[k]*2^l+R(res[0])
                    yQi = R(xQi^3+a*xQi+b)
                    if yQi.is_square():
                        yQi = sqrt(yQi)
                        Qi = E(xQi,yQi,1)
                        Qk = Qi+Qj
                        if output(Qk,l) == Y[k]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qk not in list_Qk:
                                list_Qk.append(Qk)
                        yQi = -yQi
                        Qi=E(xQi,yQi,1)
                        Qk = Qi+Qj
                        if output(Qk,l) == Y[k]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qk not in list_Qk:
                                list_Qk.append(Qk)

                    if len(list_Qi) == 1:
                        list_index.append(i)
                        list_Q.append(list_Qi[0])
                        list_index.append(k)
                        list_Q.append(list_Qk[0])  
                        r.remove(r[t])                   
            else: #nor i neither j are in list_index
                if k in list_index:
                    kk = list_index.index(k)
                    Qk=list_Q[kk]
                    list_Qi = []
                    list_Qj = []
                    xQk = ZZ(Qk[0])
                    X1 = Y[i]*2**l + d1
                    X2 = Y[j]*2**l + d2
                    X3 = xQk
                    G = Poly(X1,X2,X3)
                    G=PP(G)%p
                    monomials = G.monomials()
                    monomials.reverse()
                    #[1,d2,d1,...]
                    res = (coppersmith([G],monomials,p,l,1))[1:3]
                    xQi = Y[i]*2^l+R(res[1])
                    xQj = Y[j]*2^l+R(res[0])
                    yQi = R(xQi^3+a*xQi+b)
                    if yQi.is_square():
                        yQi = sqrt(yQi)
                        Qi = E(xQi,yQi,1)
                        Qj = Qi-Qk
                        if output(Qj,l) == Y[j]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)
                        yQi = -yQi
                        Qi=E(xQi,yQi,1)
                        Qj = Qi-Qk
                        if output(Qj,l) == Y[j]:
                            if Qi not in list_Qi:
                                list_Qi.append(Qi)
                            if Qj not in list_Qj:
                                list_Qj.append(Qj)

                    if len(list_Qi) == 1:
                        list_index.append(i)
                        list_Q.append(list_Qi[0])
                        list_index.append(j)
                        list_Q.append(list_Qj[0])  
                        r.remove(r[t])                                      
                else:
                    t=t+1
    return list_index,list_Q        
                                    

