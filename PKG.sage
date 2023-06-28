from time import *

"""
Case 1

n=16
l=4
Poly = [0,11,13,14]
q = 99839
p=2*q+1
R = Integers(p)
"""

def knapsack(n,R,q,Poly,m,l):
    #init
    u=[] #we keep the whole u even if do not need it
    w=[]
    for i in range(n):
        u.append(ZZ.random_element(0,2))
        wi = R.random_element()
        while wi.multiplicative_order() != q:
            wi = R.random_element()
        w.append(wi)
    #outputs
    Y = []
    V = []
    for j in range(m):
        v=product([w[i]^u[i+j] for i in range(n)]) 
        V.append(v)
        Y.append(ZZ(v)>>l)
        u.append(sum([u[j+k]for k in Poly])%2)
    return u,w,Y

def check_knapsack(n,Poly,m,l,u,w):
    Y = []
    for j in range(m):
        Q = prod([w[i]^u[i+j] for i in range(n)])
        Y.append(ZZ(Q)>>l)
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


def find_systeme(n,R,u,Y,l,triplets):
    m = len(u)-n
    list_index=[]
    list_point=[]
    for [i,j,k]  in triplets:
        if i in list_index:
            continue
        if j in list_index:
            continue
        list_hi = []
        list_hj = []

        for ri in range(2**l):
            hi = R(Y[i]*2^l+ri)
            for rj in range(2^l):
                hj = R(Y[j]*2^l +rj)
                if ZZ(hi*hj)>> l == Y[k]:
                    if hi not in list_hi:
                        list_hi.append(hi)
                    if hj not in list_hj:
                        list_hj.append(hj)

        if len(list_hi)==1 :
            list_index.append(i)
            list_point.append(list_hi[0])
        #if len(list_Qi) >1 :
            #print('l too big')           
        if len(list_hj)==1 :
            list_index.append(j)
            list_point.append(list_hj[0]) 

    return(list_index,list_point) 

def attack(n,q,R,u,Y,l):
    t = n//2
    S=Integers(q)
    m=len(u)-n
    triplets = find_triplet(n,u)
    list_index,list_point = find_systeme(n,R,u,Y,l,triplets)
    if len(list_index) < n:
        return 
 
    M = matrix(S,n,n)
    for i in range(n):
        ii = list_index[i]
        for j in range(n):
            M[i,j] = u[ii+j]
    
    if M.rank() !=n :
        return
    antiM = M^(-1)
    antiM = matrix(ZZ,antiM)
        
    P=[]
    for i in range(n):
        Pi = prod([list_point[j]^antiM[i,j] for j in range(n)])
        P.append(Pi)
    ubis = copy(u)[:n]
    Ybis = check_knapsack(n,Poly,100,l,ubis,P)
    if Ybis == Y[:100]:
        return P
            
def test_attack(n,R,q,Poly,m,l,run):
    #warning initialization phase is long!
    instances = []
    for i in range(run) :
        u,w,Y = knapsack(n,R,q,Poly,m,l)
        instances.append([u,Y])
    T=time()
    res=0
    print("end of initialization phase")
    for i in range(run):
        [u,Y]=instances.pop()
        A = attack(n,q,R,u,Y,l)
        if A!= None:
            res+=1

    T=time()-T 
    return(T/run,res/run) 



"""
Case 2

n=16
l=4
Poly = [0,11,13,14]
q = 72536599031050480402372360602698911648481683373808860129469667649180998227293
p=2*q+1
R = Integers(p)
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

def find_systeme_cop(n,R,u,Y,l,triplets):
    PP.<d1,d2,d3> = ZZ[]

    m = len(u)-n
    list_index=[]
    list_point=[]
    for [i,j,k]  in triplets:
        if i in list_index:
            continue
        if j in list_index:
            continue
        si = Y[i]*2^l
        sj = Y[j]*2^l
        sk = Y[k]*2^l

        G = (si+d1)*(sj+d2)-(sk+d3)
        G = PP(G)%p 
        monomials = [d1^0,d1,d2,d3,d1*d2]
        res = (coppersmith([G],monomials,p,l,1))[1:4]
        try :
            hi = R(si + res[0])
            hj = R(sj + res[1])
        except :
            continue
        list_index.append(i)
        list_index.append(j)
        list_point.append(hi)
        list_point.append(hj)

    return(list_index,list_point) 

def attack_cop(n,q,R,u,Y,l):
    t = n//2
    S=Integers(q)
    m=len(u)-n
    triplets = find_triplet(n,u)
    list_index,list_point = find_systeme_cop(n,R,u,Y,l,triplets)
    if len(list_index) < n:
        return 
 
    M = matrix(S,n,n)
    for i in range(n):
        ii = list_index[i]
        for j in range(n):
            M[i,j] = u[ii+j]
    
    if M.rank() !=n :
        return
    antiM = M^(-1)
    antiM = matrix(ZZ,antiM)
        
    P=[]
    for i in range(n):
        Pi = prod([list_point[j]^antiM[i,j] for j in range(n)])
        P.append(Pi)
    ubis = copy(u)[:n]
    Ybis = check_knapsack(n,Poly,100,l,ubis,P)
    if Ybis == Y[:100]:
        return P
            
def test_attack_cop(n,R,q,Poly,m,l,run):
    #warning initialization phase is long!
    instances = []
    for i in range(run) :
        u,w,Y = knapsack(n,R,q,Poly,m,l)
        instances.append([u,Y])
    T=time()
    res=0
    print("end of initialization phase")
    for i in range(run):
        [u,Y]=instances.pop()
        A = attack_cop(n,q,R,u,Y,l)
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

def find_systeme_tricky(u,Y,m,r):
    #r=list_triplet_refined
    list_index = []
    list_v=[]
    #we need to force one of the Qi
    [i,j,k] = r[0]
    print(r[0])
    list_vi = []
    list_vj = []
    list_vk = []

    for ri in range(2**l):
        vi = R(Y[i]*2^l+ri)
        for rj in range(2^l):
            vj = R(Y[j]*2^l +rj)
            vk = vi*vj
            if (vk >> l) == Y[k]:
                print(1)
                if vi not in list_vi:
                    list_vi.append(vi)
                if vj not in list_vj:
                    list_vj.append(vj)
                if vk not in list_vk :
                    list_vk.append(vk)

    if len(list_vi)==1:
        list_index.append(i)
        list_v.append(list_vi[0])                   
    if len(list_vj)==1:
        list_index.append(j)
        list_v.append(list_vj[0])   
    if len(list_vk)==1:
        list_index.append(k)
        list_v.append(list_vk[0])

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
            vi = list_v[ii]
            if j in list_index:
                jj = list_index.index(j)
                vj = list_v[jj]
                if k in list_index:
                    kk = list_index.index(k)
                    vk = list_v[kk]
                    assert(vi*vj==vk)
                else :
                    vk = vi*vj
                    list_index.append(k)
                    list_v.append(vk)
                    r.remove(r[t])
            else : #j is not in list_index
                if k in list_index:
                    kk = list_index.index(k)
                    vk = list_v[kk]
                    vj = vk*(vi^(-1))
                    list_index.append(j)
                    list_v.append(vj)
                    r.remove(r[t])
                else : 
                    list_vj = []
                    list_vk = []
                    for rj in range(2^l):
                        vj = R(Y[j]*2^l +rj)
                        vk = vi*vj
                        if (vk >> l) == Y[k]:
                            if vj not in list_vj:
                                list_vj.append(vj)
                            if vk not in list_vk:
                                list_vk.append(vk)
                    if len(list_vj) == 1:
                        list_index.append(j)
                        list_v.append(list_vj[0])
                        list_index.append(k)
                        list_v.append(list_vk[0])  
                        r.remove(r[t])                                   
        else : #i is not in list_index
            if j in list_index:
                jj = list_index.index(j)
                vj = list_v[jj]
                if k in list_index:
                    kk = list_index.index(k)
                    vk = list_v[kk]
                    vi = vk*(vj^(-1))
                    list_index.append(i)
                    list_v.append(vi) 
                    r.remove(r[t])                
                else :
                    list_vi = []
                    list_vk = []
                    for ri in range(2^l):
                        vi = R(Y[i]*2^l +ri)
                        vk = vi*vj
                        if( vk >> l) == Y[k]:
                            if vi not in list_vi:
                                list_vi.append(vi)
                            if vk not in list_vk:
                                list_vk.append(vk)
                    if len(list_vi) == 1:
                        list_index.append(i)
                        list_v.append(list_vi[0])
                        list_index.append(k)
                        list_v.append(list_vk[0])  
                        r.remove(r[t])                   
            else: #nor i neither j are in list_index
                if k in list_index:
                    kk = list_index.index(k)
                    vk=list_v[kk]
                    list_vi = []
                    list_vj = []
                    for ri in range(2^l):
                        vi = R(Y[i]*2^l +ri)
                        vj = vk*(vi^(-1))
                        if (vj >> l) == Y[j]:
                            if vi not in list_vi:
                                list_vi.append(vi)
                            if vj not in list_vj:
                                list_vj.append(vj)
                    if len(list_vi) == 1:
                        list_index.append(i)
                        list_v.append(list_vi[0])
                        list_index.append(j)
                        list_v.append(list_vj[0])  
                        r.remove(r[t])                                      
                else:
                    t=t+1
    return list_index,list_v