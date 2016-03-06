import numpy
import itertools, random
variables=['x0','x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12']

pvarsC3 = variables[1:7]
prelationsC3 = [[['x1','x2','x3'],['x4','x6'],[1,1]],
              [['x4','x5','x6'],[],[]]]

pvarsD3 = variables[1:8]
prelationsD3 = [[['x1','x2','x3'],['x4','x6'],[1,1]],
                [['x4','x5'],['x6'],[1]],
                [['x6','x7'],[],[]]]

pvarsD9 = variables[1:8]
prelationsD9 = [[['x1','x2','x3'],['x4','x6'],[1,1]],
                [['x4','x5'],[],[]],
                [['x6','x7'],[],[]]]

pvarsD16 = variables[1:8]
prelationsD16 = [[['x1','x2','x3'],['x4','x7'],[1,1]],
                 [['x4','x5'],['x6'],[1]],
                 [['x6','x7'],[],[]]]

pvarsG1 = variables[1:8]
prelationsG1 = [[['x1','x7'],[],[]],[['x2','x3','x4'],['x1'],[1]],
                [['x4','x5','x6'],['x1'],[2]],[['x5','x6','x7'],['x2','x3'],[1,1]],
                [['x1','x2','x3'],['x5','x6'],[1,1]]]

pvarsG3 = variables[1:8]
prelationsG3 = [[['x1','x7'],[],[]],[['x2','x3','x4'],[],[]],
                [['x4','x5','x6'],['x1'],[1]],[['x5','x6','x7'],['x2','x3'],[1,1]],
                [['x1','x2','x3'],['x5','x6'],[1,1]]]

pvarsG4 = variables[1:8]
prelationsG4 = [[['x1','x7'],['x4'],[1]],[['x2','x3','x4'],['x7'],[1]],
                [['x4','x5','x6'],['x1','x2'],[1,1]],[['x5','x6','x7'],['x2'],[1]],
                [['x1','x2','x3'],[],[]]]

pvarsG5 = variables[1:8]
prelationsG5 = [[['x1','x7'],['x4'],[1]],[['x2','x3','x4'],['x7'],[1]],
                [['x4','x5','x6'],[],[]],[['x5','x6','x7'],['x2','x3'],[1,1]],
                [['x1','x2','x3'],[],[]]]

pvarsH2 = variables[0:8]
prelationsH2 = [[['x1','x2'],['x0'],[1]],[['x0','x7'],['x1'],[1]],
                [['x1','x6'],['x7'],[1]],[['x2','x7'],[],[]],
                [['x0','x6'],[],[]],[['x3','x4','x5'],['x0','x1'],[1,1]]]

pvarsH6 = variables[0:8]
prelationsH6 = [[['x1','x2'],['x0'],[1]],[['x0','x7'],['x1'],[1]],
                [['x1','x6'],['x7'],[1]],[['x2','x7'],[],[]],
                [['x0','x6'],[],[]],[['x3','x4','x5'],['x0','x2'],[1,1]]]

pvarsH10 = variables[0:8]
prelationsH10 = [[['x1','x2'],['x0'],[1]],[['x0','x7'],['x1'],[1]],
                [['x1','x6'],['x7'],[1]],[['x2','x7'],[],[]],
                [['x0','x6'],[],[]],[['x3','x4','x5'],['x2','x6'],[1,1]]]

pvarsI3 = variables[0:8]
prelationsI3 = [[['x0','x7'],['x3'],[1]],[['x3','x6'],['x7'],[1]],
                [['x0','x6'],[],[]],[['x1','x2'],['x0'],[1]],
                [['x3','x4','x5'],['x0','x1'],[1,1]],
                [['x4','x5','x7'],['x1'],[1]]]

pvarsI8 = variables[0:8]
prelationsI8 = [[['x0','x7'],['x3'],[1]],[['x3','x6'],['x7'],[1]],
                [['x0','x6'],[],[]],[['x1','x2'],['x6'],[1]],
                [['x3','x4','x5'],['x0','x1'],[1,1]],
                [['x4','x5','x7'],['x1'],[1]]]

pvarsI11 = variables[0:8]
prelationsI11 = [[['x0','x7'],['x3'],[1]],[['x3','x6'],['x7'],[1]],
                [['x0','x6'],[],[]],[['x1','x2'],[],[]],
                [['x3','x4','x5'],['x0','x1'],[1,1]],
                [['x4','x5','x7'],['x1'],[1]]]

pvarsJ1 = variables[0:8]
prelationsJ1 = [[['x3','x6'],['x7'],[1]], [['x0','x1','x2'],['x4','x5'],[1,1]],
                [['x4','x5','x6'],['x1','x2'],[1,1]],[['x0','x7'],['x3'],[1]],
                [['x0','x6'],[],[]],[['x3','x4','x5'],[],[]],[['x4','x5','x7'],['x6'],[1]],
                [['x1','x2','x3'],['x6'],[1]],[['x1','x2','x7'],['x6'],[2]]]

pvarsJ2 = variables[0:8]
prelationsJ2 = [[['x3','x6'],['x7'],[1]], [['x0','x1','x2'],['x4','x5'],[1,1]],
                [['x4','x5','x6'],['x1','x2'],[1,1]],[['x0','x7'],['x3'],[1]],
                [['x0','x6'],[],[]],[['x3','x4','x5'],['x0'],[1]],[['x4','x5','x7'],[],[]],
                [['x1','x2','x3'],[],[]],[['x1','x2','x7'],['x6'],[1]]]

pvarsK2 = variables[0:9]
prelationsK2 = [[['x0','x7'],['x1'],[1]],[['x1','x8'],['x7'],[1]],
                [['x0','x8'],[],[]],[['x8','x2'],['x6'],[1]],
                [['x7','x6'],['x8'],[1]],[['x1','x6'],[],[]],
                [['x0','x6'],['x2'],[1]],[['x1','x2'],['x0'],[1]],
                [['x7','x2'],[],[]],[['x3','x4','x5'],['x0','x1'],[1,1]]]

pvarsZ1 = variables[1:9]
prelationsZ1 = [[['x1','x2','x5'],[],[]],[['x1','x2','x6'],['x7'],[1]],
                [['x2','x4','x5'],['x8'],[1]],
                [['x2','x4','x6'],['x7','x8'],[1,1]],
                [['x3','x8','x7'],[],[]],[['x3','x4','x6'],['x1','x5'],[1,1]],
                [['x3','x4','x7'],['x1'],[1]],[['x3','x6','x8'],['x5'],[1]],
                [['x1','x8'],['x4'],[1]],[['x5','x7'],['x6'],[1]]]

pvarsZ2 = variables[1:9]
prelationsZ2 = [[['x1','x2','x5'],[],[]],[['x1','x2','x6'],['x7'],[1]],
                [['x2','x4','x5'],['x8'],[1]],
                [['x2','x4','x6'],['x7','x8'],[1,1]],
                [['x3','x8','x7'],['x2'],[1]],[['x3','x4','x6'],[],[]],
                [['x3','x4','x7'],['x1','x2'],[1,1]],[['x3','x6','x8'],['x2','x5'],[1,1]],
                [['x1','x8'],['x4'],[1]],[['x5','x7'],['x6'],[1]]]

#note: 3rd primitive relation for M polytopes given in paper has typo
#(paper says x4+x6+x8 = x1+x2, not x4+x6+x8 = x2+x3)
pvarsM1 = variables[1:9]
prelationsM1 = [[['x1','x8'],[],[]],[['x1','x2','x3'],['x4','x6'],[1,1]],
                [['x4','x6','x8'],['x2','x3'],[1,1]],
                [['x4','x5'],[],[]],[['x6','x7'],[],[]],
                [['x2','x3','x5'],['x6','x8'],[1,1]],
                [['x2','x3','x7'],['x4','x8'],[1,1]]]

pvarsM4 = variables[1:9]
prelationsM4 = [[['x1','x8'],[],[]],[['x1','x2','x3'],['x4','x6'],[1,1]],
                [['x4','x6','x8'],['x2','x3'],[1,1]],
                [['x4','x5'],['x1'],[1]],[['x6','x7'],[],[]],
                [['x2','x3','x5'],['x6'],[1]],
                [['x2','x3','x7'],['x4','x8'],[1,1]]]

pvarsM = variables[1:9]
prelationsM = [[['x1','x8'],['x5'],[1]],[['x4','x5'],['x7'],[1]],
               [['x6','x7'],['x1'],[1]],[['x1','x2','x3'],['x6'],[1]],
               [['x2','x3','x5'],['x6','x8'],[1,1]],[['x2','x3','x7'],[],[]],
               [['x4','x8','x6'],[],[]]]

#Last primitive relation is x1+x2+x8=0 in paper, should be x2+x4+x8=0
pvarsP341 = variables[0:9]
prelationsP341 = [[['x0','x7'],[],[]],[['x0','x8'],['x1'],[1]],[['x3','x5'],['x4'],[1]],
                  [['x4','x6'],['x5'],[1]],[['x1','x7'],['x8'],[1]],[['x3','x6'],[],[]],
                  [['x1','x2','x5'],['x0','x6'],[1,1]],[['x1','x2','x4'],['x0'],[1]],
                  [['x2','x5','x8'],['x6'],[1]],[['x2','x4','x8'],[],[]]]

#For type R, need to include x4+x6=x0
pvarsR1 = variables[0:9]
prelationsR1 = [[['x1','x2','x7'],['x3','x8'],[1,1]],[['x1','x2','x4'],['x3'],[1]],
                [['x3','x5'],['x4'],[1]],[['x1','x2','x5'],[],[]],
                [['x0','x1','x2'],['x3','x6'],[1,1]],[['x3','x6','x8'],['x1','x2'],[1,1]],
                [['x0','x7'],['x4'],[1]],[['x4','x8'],['x7'],[1]],[['x0','x8'],[],[]],
                [['x6','x7'],[],[]],[['x4','x6'],['x0'],[1]]]

pvarsR2 = variables[0:9]
prelationsR2 = [[['x1','x2','x7'],['x3','x8'],[1,1]],[['x1','x2','x4'],['x3'],[1]],
                [['x3','x5'],['x0'],[1]],[['x1','x2','x5'],['x6'],[1]],
                [['x0','x1','x2'],['x3','x6'],[1,1]],[['x3','x6','x8'],['x1','x2'],[1,1]],
                [['x0','x7'],['x4'],[1]],[['x4','x8'],['x7'],[1]],[['x0','x8'],[],[]],
                [['x6','x7'],[],[]],[['x4','x6'],['x0'],[1]]]

pvarsR3 = variables[0:9]
prelationsR3 = [[['x1','x2','x7'],['x3','x8'],[1,1]],[['x1','x2','x4'],['x3'],[1]],
                [['x3','x5'],[],[]],[['x1','x2','x5'],['x6','x8'],[1,1]],
                [['x0','x1','x2'],['x3','x6'],[1,1]],[['x3','x6','x8'],['x1','x2'],[1,1]],
                [['x0','x7'],['x4'],[1]],[['x4','x8'],['x7'],[1]],[['x0','x8'],[],[]],
                [['x6','x7'],[],[]],[['x4','x6'],['x0'],[1]]]

pvarsP358ii = variables[0:10]
prelationsP358ii = [[['x0','x4'],[],[]],[['x1','x5'],[],[]],[['x2','x6'],[],[]],
                    [['x3','x7'],[],[]],[['x8','x9'],[],[]],
                    [['x0','x1','x2'],['x7','x8'],[1,1]],
                    [['x0','x1','x3'],['x6','x8'],[1,1]],
                    [['x0','x2','x3'],['x5','x8'],[1,1]],
                    [['x1','x2','x3'],['x4','x8'],[1,1]],
                    [['x0','x1','x9'],['x6','x7'],[1,1]],
                    [['x0','x2','x9'],['x5','x7'],[1,1]],
                    [['x0','x3','x9'],['x5','x6'],[1,1]],
                    [['x1','x2','x9'],['x4','x7'],[1,1]],
                    [['x1','x3','x9'],['x4','x6'],[1,1]],
                    [['x2','x3','x9'],['x4','x5'],[1,1]],
                    [['x4','x5','x6'],['x3','x9'],[1,1]],
                    [['x4','x5','x7'],['x2','x9'],[1,1]],
                    [['x4','x6','x7'],['x1','x9'],[1,1]],
                    [['x5','x6','x7'],['x0','x9'],[1,1]],
                    [['x4','x5','x8'],['x2','x3'],[1,1]],
                    [['x4','x6','x8'],['x1','x3'],[1,1]],
                    [['x4','x7','x8'],['x1','x2'],[1,1]],
                    [['x5','x6','x8'],['x0','x3'],[1,1]],
                    [['x5','x7','x8'],['x0','x2'],[1,1]],
                    [['x6','x7','x8'],['x0','x1'],[1,1]]]

pvarsP358iii = variables[0:9]
prelationsP358iii = [[['x0','x4'],[],[]],[['x1','x5'],[],[]],[['x2','x6'],[],[]],
                    [['x3','x7'],[],[]],
                    [['x0','x1','x2'],['x7','x8'],[1,1]],
                    [['x0','x1','x3'],['x6','x8'],[1,1]],
                    [['x0','x2','x3'],['x5','x8'],[1,1]],
                    [['x1','x2','x3'],['x4','x8'],[1,1]],
                    [['x4','x5','x8'],['x2','x3'],[1,1]],
                    [['x4','x6','x8'],['x1','x3'],[1,1]],
                    [['x4','x7','x8'],['x1','x2'],[1,1]],
                    [['x5','x6','x8'],['x0','x3'],[1,1]],
                    [['x5','x7','x8'],['x0','x2'],[1,1]],
                    [['x6','x7','x8'],['x0','x1'],[1,1]]]

pvarsW = variables[1:10]
prelationsW = [[['x1','x4'],['x7'],[1]],[['x2','x5'],['x8'],[1]],
               [['x3','x6'],['x9'],[1]],[['x1','x2','x3'],[],[]],
               [['x4','x5','x6'],[],[]],[['x7','x8','x9'],[],[]],
               [['x1','x2','x9'],['x6'],[1]],[['x4','x5','x9'],['x3'],[1]],
               [['x1','x3','x8'],['x5'],[1]],[['x4','x6','x8'],['x2'],[1]],
               [['x2','x3','x7'],['x4'],[1]],[['x5','x6','x7'],['x1'],[1]],
               [['x1','x8','x9'],['x5','x6'],[1,1]],[['x4','x8','x9'],['x2','x3'],[1,1]],
               [['x2','x7','x9'],['x4','x6'],[1,1]],[['x5','x7','x9'],['x1','x3'],[1,1]],
               [['x3','x7','x8'],['x4','x5'],[1,1]],[['x6','x7','x8'],['x1','x2'],[1,1]]]

def rowreduce(m, col, prevpivots, pd, pv, varlist):
    dims = m.shape
    i = 0
    if col in varlist:
        if col < dims[1]-1:
            rowreduce(m, col+1, prevpivots, pd, pv, varlist)
            return
        else:
            return
    while (m[i,col]==0 or i in prevpivots):
        i = i+1
    p = m[i,col]
    pindex = i
    i = 0
    while i < dims[0]:
        if (i != pindex):
            m[i,...] = -p*m[i,...]+m[i,col]*m[pindex,...]
        i = i+1
    if (col < dims[1]-1):
        prevpivots.append(pindex)
        rowreduce(m, col+1, prevpivots, pd, pv, varlist)
    pd[pv[col]] = (-m[pindex,varlist[0]]*pd[pv[varlist[0]]]
                          -m[pindex,varlist[1]]*pd[pv[varlist[1]]]
                          -m[pindex,varlist[2]]*pd[pv[varlist[2]]]
                          -m[pindex,varlist[3]]*pd[pv[varlist[3]]])/m[pindex,col]

def getmaxfans(pv, pr, polytopename):
    #Get all primitive relations of the form v1+v2+v3=v4+v5
    #and put them into a list
    print polytopename
    relevantprs = []
    for i in pr:
        if len(i[0])==3 and len(i[1])==2 and set(i[2])==set([1]):
            relevantprs.append(i)
    #print(relevantprs)
    pdict = dict()
    idmat = numpy.identity(4)
    initialvars = list(relevantprs[0][0])
    initialvars.append(relevantprs[0][1][0])
    pdict[initialvars[0]] = idmat[0,...]
    pdict[initialvars[1]] = idmat[1,...]
    pdict[initialvars[2]] = idmat[2,...]
    pdict[initialvars[3]] = idmat[3,...]
    varlist = [pv.index(initialvars[0]),
               pv.index(initialvars[1]),
               pv.index(initialvars[2]),
               pv.index(initialvars[3])]
    m = numpy.zeros([len(pr),len(pv)])
    for i in range(len(pr)):
        for j in pr[i][0]:
            m[i,pv.index(j)] = 1
        for j in range(len(pr[i][1])):
            m[i,pv.index(pr[i][1][j])] = -pr[i][2][j]
    rowreduce(m, 0, [], pdict, pv, varlist)
    ell = []
    for i in pv:
        ell.append(pdict[i])
    a=numpy.array(ell)
    print(a)
    ell = itertools.combinations(pv,4)
    ell = [set(i) for i in ell]
    for i in pr:
        print(i)
        ell = filter(lambda k: not set(i[0]).issubset(k), ell)
    #print(ell)
    #Do a consistency check to make sure that the vertices and faces
    #really come from a smooth fano polytope, all primitive relations
    #hold, and every primitive relation is included
    '''curnonfaces = []
    for i in range(3):
        for j in itertools.combinations(pv,i+2):
            flag = 0
            for k in ell:
                if set(j).issubset(k):
                    flag = 1
            for q in curnonfaces:
                if q.issubset(set(j)):
                    flag = 1
            if flag==0:
                curnonfaces.append(set(j))
    if (len(curnonfaces)==len(pr)):
        print("All primitive relations included!")
    else:
        print("Not all primitive relations included!")
    dualvertices = []
    v = numpy.ones(4)
    for n in ell:
        q = sorted(list(n),key = lambda k: variables.index(k))
        mat = numpy.array([pdict[i] for i in q])
        dualvertices.append(numpy.dot(numpy.linalg.inv(mat),v))
    print(dualvertices)
    flag = 0
    for i in dualvertices:
        dots = [numpy.dot(i,pdict[k]) for k in pv]
        #Every dual vertex should eval to 1 on exactly four vertices
        #and < 1 on everything else
        dots = filter(lambda k: k >= 1,dots)
        #print(dots)
        if not (len(dots)==4 and sum(dots)==4):
            flag = 1
    if flag==0:
        print('Vertex check OK')
    else:
        print('Failed vertex check')
    flag = 0
    for i in pr:
        #Check RHS and LHS equal each other in all primitive relations
        lhs = sum([pdict[j] for j in i[0]])
        if (len(i[2])==0):
            rhs = numpy.zeros(4)
        else:
            rhs = sum([i[2][k]*pdict[i[1][k]] for k in range(len(i[2]))])
        if not (lhs[0]==rhs[0] and lhs[1]==rhs[1] and lhs[2]==rhs[2] and
                lhs[3]==rhs[3]):
            flag = 1
    if flag==0:
        print('Primitive relations check OK')
    else:
        print('Failed primitive relations check')'''
    dmaxfanslist = [ell]
    #Generate all the possible combinations of primitive relations that can
    #give a delta-maximal fan
    prlist = []
    #n = 0
    for i in range(len(relevantprs)):
        #Loop through all possible subsets of the set of relevant primitive relations
        #Elements of the subset will be primitive relations on which we do wall crossings
        for j in itertools.combinations(relevantprs,i+1):
            flag = 0
            #Check that the variables in every pair of primitive relations in the current subset
            #overlap on no more than 3 variables
            for x in itertools.combinations(j,2):
                s1 = set(x[0][0]).union(set(x[0][1]))
                s2 = set(x[1][0]).union(set(x[1][1]))
                if len(s1.intersection(s2)) > 3:
                    flag = 1
            if flag==0:
                prlist.append(j)
                #n = n+1
    #print(prlist)
    #print('Number of delta-maximal fans:',n)
    print('Number of delta-maximal fans:',len(prlist)+1)
    #Generate the cones in the new fans by doing wall-crossings
    for i in prlist:
        newflist = list(ell)
        for j in i:
            bs = set(j[0])
            cs = set(j[1])
            nset1 = bs.union(set([j[1][0]]))
            nset2 = bs.union(set([j[1][1]]))
            iset1 = cs.union(set([j[0][0],j[0][1]]))
            iset2 = cs.union(set([j[0][0],j[0][2]]))
            iset3 = cs.union(set([j[0][1],j[0][2]]))
            newflist.remove(iset1)
            newflist.remove(iset2)
            newflist.remove(iset3)
            newflist.append(nset1)
            newflist.append(nset2)
        dmaxfanslist.append(newflist)
    #print(dmaxfanslist)
    nonfaces = []
    for n in dmaxfanslist:
        curnonfaces = []
        for i in range(3):
            for j in itertools.combinations(pv,i+2):
                flag = 0
                for k in n:
                    if set(j).issubset(k):
                        flag = 1
                for q in curnonfaces:
                    if q.issubset(set(j)):
                        flag = 1
                if flag==0:
                    curnonfaces.append(set(j))
        nonfaces.append(curnonfaces)
    #print(nonfaces)
    coneinverses = []
    for m in dmaxfanslist:
        curinverses = []
        for n in m:
            q = sorted(list(n),key = lambda k: variables.index(k))
            mat = numpy.array([pdict[i] for i in q])
            curinverses.append(numpy.linalg.inv(mat))
        coneinverses.append(curinverses)
    prelationslist = []
    #generate the primitive relations for all the delta-maximal fans so we can check if they are projective
    numvars = len(pv)
    idmat = numpy.identity(numvars)
    stdbasisdict = dict()
    for i in range(numvars):
        stdbasisdict[pv[i]] = idmat[i,...]
    for m in range(len(nonfaces)):
        curprelations = []
        for n in range(len(nonfaces[m])):
            cr = numpy.zeros(numvars)
            curvec = numpy.zeros(4)
            for i in nonfaces[m][n]:
                curvec = curvec+pdict[i]
                cr = cr+stdbasisdict[i]
            flag = 0
            j = 0
            while (flag != 1):
                x = numpy.dot(curvec,coneinverses[m][j])
                if (x[0] >= 0 and x[1] >= 0 and x[2] >= 0 and x[3] >= 0):
                    flag = 1
                    v = sorted(list(dmaxfanslist[m][j]))
                    #print(v)
                    #print(x)
                    cr = cr-x[0]*stdbasisdict[v[0]]-x[1]*stdbasisdict[v[1]]-x[2]*stdbasisdict[v[2]]-x[3]*stdbasisdict[v[3]]
                j = j+1
            curprelations.append(cr)
        prelationslist.append(curprelations)
    #print(prelationslist)
    numprojective = 0
    numnonprojective = 0
    projrelationdata = []
    nonprojrelationdata = []
    projrelations = []
    nonprojrelations = []
    for m in prelationslist:
        #print(m)
        result = isPointedCone(m)
        if result==None:
            raise Exception
        elif (result[0]==0):
            if (result[1] < 0).any():
                raise Exception
            if (result[1]==0).all():
                raise Exception
            if not (numpy.dot(result[1],numpy.array(m))==0).all():
                raise Exception
            #print("Is not projective!",result)
            numnonprojective = numnonprojective+1
            nonprojrelationdata.append(getTuples(m))
            nonprojrelations.append(m)
        else:
            flag = 0
            for i in m:
                if numpy.dot(i,result[1]) <= 0:
                    flag = 1
            if (flag==1):
                raise Exception
            #print("Is projective!",result)
            numprojective = numprojective+1
            projrelationdata.append(getTuples(m))
            projrelations.append(m)
    #Uncomment next line to generate processing files for Macaulay2
    #m2fileoutput(pv, a, [projrelations[0],projrelations[1]], polytopename)
    print('Number of nonprojective fans:',numnonprojective)
    print('Number of projective fans:',numprojective)
    print('Projective isomorphism types:')
    if len(projrelationdata) <= 1:
        counter = len(projrelationdata)
    else:
        previndex = 0
        i = 1
        counter = 0
        while i < len(projrelationdata):
            if projrelationdata[i] != projrelationdata[i-1]:
                k = checkAreIsomorphicTranspose(projrelations[previndex:i])
                #k = checkAreIsomorphic(projrelations[previndex:i])
                if (k >1):
                    print k
                previndex = i
                counter = counter+k
            i = i+1
        k = checkAreIsomorphicTranspose(projrelations[previndex:i])
        #k = checkAreIsomorphic(projrelations[previndex:i])
        if (k >1):
            print k
        counter = counter+k
    print(counter)
    print('Nonprojective isomorphism types:')
    #print(nonprojrelations[0])
    if len(nonprojrelationdata) <= 1:
        counter = len(nonprojrelationdata)
    else:
        previndex = 0
        i = 1
        counter = 0
        while i < len(nonprojrelationdata):
            if nonprojrelationdata[i] != nonprojrelationdata[i-1]:
                k = checkAreIsomorphicTranspose(nonprojrelations[previndex:i])
                if (k >1):
                    print k
                previndex = i
                counter = counter+k
            i = i+1
        k = checkAreIsomorphicTranspose(nonprojrelations[previndex:i])
        if (k > 1):
            print k
        counter = counter+k
    print(counter)
    '''print('Checking if the projective types have isolated rays:')
    for i in projrelations:
        flag = 0
        q = numpy.ones([len(pv),len(pv)])
        q = q-numpy.identity(len(pv))
        for m in i:
            maxes = map(lambda k: max(k,0), m)
            #print(maxes)
            if (sum(maxes)==2):
                firstone = maxes.index(1)
                maxes[firstone] = 0
                secondone = maxes.index(1)
                maxes[secondone] = 0
                q[firstone, secondone] = 0
                q[secondone, firstone] = 0
        #print(q)
        for j in range(len(pv)):
            if (sum(q[j,...])<= 4):
                mat = []
                for k in range(len(pv)):
                    if q[j,k]==1:
                        mat.append(pdict[pv[k]])
                if numpy.linalg.det(mat) != 0:
                    flag = 1
        if (flag==1):
            print('Has isolated ray!')
        else:
            print('No isolated ray!')'''
            
#Given rectangular matrix m, get the two matrices
#formed by sorting each row entry-wise, then sorting the
#sorted rows, and likewise sorting each column entry-wise
#and then sorting the columns.  If two matrices m1 and
#m2 are related by some permutation of rows and columns,
#then getTuples(m1)==getTuples(m2).  This is a necessary
#but not sufficient condition.  A simple counterexample
#to see that it is not sufficient is:
#m1 = [[1,0,1,0],
#      [0,1,0,0],
#      [0,0,1,1],
#      [0,0,0,1]]
#m2 = [[1,0,1,0],
#      [0,1,0,1],
#      [0,0,1,0],
#      [0,0,0,1]]
def getTuples(m):
    rowdata = []
    for i in m:
        k = list(i)
        rowdata.append(sorted(k))
    j = numpy.array(m)
    j = numpy.transpose(j)
    rowdata = sorted(rowdata)
    coldata = []
    for i in j:
        k = list(i)
        coldata.append(sorted(k))
    coldata = sorted(coldata)
    return [rowdata,coldata]

#Check if list of rectangular matrices are equivalent
#up to arbitrary permutation of rows and columns (in other words,
#arbitrary left and right multiplication by permutation matrices of the
#appropriate dimensions).
#The algorithm works by first identifying the "isomorphism class"
#of each row in the first matrix, where two rows are considered "isomorphic" if they
#are the same after sorting element-wise.  Then we sort the rows according
#to their isomorphism class.  (This sorting is not unique since there is no
#definite order for isomorphic rows.)

#We also sort all of the subsequent matrices in the list with the same method.

#Then we consider all permutations of the sorted first matrix obtained by
#permuting isomorphic rows.
#If the permuted first matrix is equivalent to any of the other matrices in the list,
#then after sorting both matrices by columns, both matrices must be equal.
def checkAreIsomorphic(relations):
    if (len(relations)==1):
        return 1
    k = list(relations[0])
    k = map(sorted, k)
    s = sorted(range(len(k)),key = lambda y: k[y])
    k.sort()
    groups = []
    cgroup = [s[0]]
    for i in range(1,len(k)):
        if k[i] != k[i-1]:
            groups.append(cgroup)
            cgroup = [s[i]]
        else:
            cgroup.append(s[i])
    groups.append(cgroup)
    comps = []
    for i in range(1,len(relations)):
        k = list(relations[i])
        k = map(sorted, k)
        s = sorted(range(len(k)),key = lambda y: k[y])
        ncomps = numpy.transpose(numpy.array([relations[i][j] for j in s]))
        ncomps = [list(j) for j in ncomps]
        ncomps.sort()
        comps.append(ncomps)
    iterobj = map(itertools.permutations, groups)
    resultvec = numpy.zeros(len(relations)-1)
    for x in itertools.product(*iterobj):
        curpermutation = []
        for i in x:
            curpermutation.extend(i)
        reorder = numpy.zeros([len(relations[0]),len(relations[0][0])])
        for i in range(len(relations[0])):
            reorder[i] = relations[0][curpermutation[i]]
        reorder = numpy.transpose(reorder)
        reorder = [list(i) for i in reorder]
        reorder.sort()
        for i in range(len(comps)):
            if reorder==comps[i]:
                resultvec[i]=1
        if (sum(resultvec)==len(relations)-1):
            return 1
    if (sum(resultvec)==len(relations)-1):
        return 1
    else:
        newrelations = []
        for i in range(len(resultvec)):
            if resultvec[i]==0:
                newrelations.append(relations[i+1])
        return 1+checkAreIsomorphic(newrelations)

def checkAreIsomorphicTranspose(relations):
    relations = [numpy.array([[r[j][k] for j in range(len(relations[0]))] for k in range(len(relations[0][0]))]) for r in relations]
    return checkAreIsomorphic(relations)

#Determine if list of vectors generates a pointed cone
#Pointed means that the cone contains no nontrivial linear subspaces of the
#ambient vector space
#The algorithm is as follows: we look for an "extremal ray" of the cone.
#Extremal ray means that if v1 is a point in the ray r, and
#v2, v3 are points in the cone such that v2+v3=v1, then
#v2, v3 are also in r.  It is easy to show that a cone is pointed
#if and only if it has at least one extremal ray.
#Assume the cone is pointed.
#Given any coordinate in the ambient vector space, if S is the subset
#of the list of vectors with positive (resp. negative) values of this coordinate,
#then at least one of the vectors in the list must generate an extremal ray.
#If a vector v does generate an extremal ray, then we can quotient out by the
#extremal ray to get a pointed cone in the quotient space.
#Conversely, if the quotient cone is pointed, then the original cone
#is pointed provided that it doesn't contain the one-dimensional
#subspace spanned by v, which is equivalent in this case to containing
#positive scalar multiples of both v and -v.
#The algorithm works recursively based on this process.
#If the cone is pointed, it returns 1, along with (as a certificate) a vector in the ambient space
#whose dot product is positive with all nonzero vectors in the list.
#If it isn't pointed, it returns 0, along with a positive linear combination of the vectors that
#equals the zero vector.            
def isPointedCone(veclist):
    #print('length of veclist',len(veclist))
    ambientdimension = len(veclist[0])
    curlowest = len(veclist)+1
    curp = -1
    curlist = []
    if len(veclist)==1:
        if (veclist[0]==0).all():
            return [1,numpy.ones(ambientdimension)]
        else:
            return [1, veclist[0]]
    #To minimize number of recursive function calls, look for the coordinate with smallest possible nonzero number of pos or neg entries
    for i in range(ambientdimension):
        numneg = 0
        numpos = 0
        neglist = []
        poslist = []
        for j in range(len(veclist)):
            if veclist[j][i] < 0:
                numneg = numneg+1
                neglist.append(j)
            elif veclist[j][i] > 0:
                numpos = numpos+1
                poslist.append(j)
        if numneg > 0 and numneg < curlowest:
            curlowest = numneg
            curp = i
            curlist = neglist
        if numpos > 0 and numpos < curlowest:
            curlowest = numpos
            curp = i
            curlist = poslist
    #print('curlist',curlist)
    #If curp is still -1 then all of the vectors are zero, so cone is pointed
    if (curp == -1):
        return([1,numpy.ones(ambientdimension)])
    for j in curlist:
        curvec = veclist[j]
        pivot = curvec[curp]
        if pivot < 0:
            sgn = -1
        else:
            sgn = 1
        quotientvecs = []
        for k in range(len(veclist)):
            if k != j:
                #print(pivot,veclist[k],curvec)
                r = abs(pivot)*veclist[k]-sgn*veclist[k][curp]*curvec
                #print(r)
                #If r is the zero vector, then one of the vectors is either zero or a scalar multiple of the pivot vector
                if (r==0).all() and veclist[k][curp]*pivot < 0:
                        a = numpy.zeros(len(veclist))
                        a[k] = abs(pivot)
                        a[j] = abs(veclist[k][curp])
                        return [0,a]
                else:    
                    r = list(r)
                    del r[curp]
                    quotientvecs.append(numpy.array(r))
        #print('quotientvecs',quotientvecs)
        result = isPointedCone(quotientvecs)
        #print(result)
        #If we find that the quotient cone is not pointed, then ray j is not extremal, so try to find a linear dependence, or delete the ray if it's
        #a positive linear combination of the other vectors
        #If the quotient cone is pointed, use the result to find a certificate vector
        if (result[0]==0):
            m = result[1]
            l = list(m[0:j])
            l.append(0)
            l.extend(list(m[j:]))
            r = 0
            for k in range(len(veclist)):
                r = r+l[k]*veclist[k][curp]
            #print('r',r,'pivot',pivot,'l',l)
            if (r*pivot <= 0):
                dep = numpy.array(l)
                dep = abs(pivot)*dep+abs(r)*numpy.identity(len(veclist))[j,...]
                return [0,dep]
            else:
                vcopy = list(veclist)
                del vcopy[j]
                result = isPointedCone(vcopy)
                if (result[0]==1):
                    return result
                else:
                    m = result[1]
                    l = list(m[0:j])
                    l.append(0)
                    l.extend(list(m[j:]))
                    return [0,numpy.array(l)]
        else:
            m = result[1]
            #print(m,curp)
            l = list(m[0:curp])
            l.append(0)
            l.extend(list(m[curp:]))
            mlift = numpy.array(l)
            #print(mlift, curvec)
            dotp = numpy.dot(mlift,curvec)
            mlift = abs(pivot)*mlift
            mlift[curp] = -sgn*dotp
            #print('mlift',mlift)
            maxabs = 0
            for k in range(len(veclist)):
                if k != j:
                    maxabs = max(abs(veclist[k][curp]),maxabs)
            #print('maxabs',maxabs)
            if pivot > 0:
                mlift = mlift+maxabs*mlift+numpy.identity(ambientdimension)[curp,...]
            else:
                mlift = mlift+maxabs*mlift-numpy.identity(ambientdimension)[curp,...]
            return [1,mlift]
        

def macaulay2format(m):
    output = "M = transpose matrix {{"
    for i in range(len(m)):
        for j in range(len(m[i])-1):
            output = output+str(int(m[i][j]))+','
        if i < len(m)-1:
            output = output+str(int(m[i][-1]))+'},{'
        else:
            output = output+str(int(m[i][-1]))+'}}'
    return output

def m2fileoutput(pv, a, relations, polytopename):
    filestr1 = polytopename+'a.m2'
    filestr2 = polytopename+'b.m2'
    pname1 = polytopename+'a'
    pname2 = polytopename+'b'
    commands = '''ff = (v, M, sr) ->
    (
    R = QQ[v]; 
    I = ideal(M*(transpose(vars(R))));
    sr = value(sr);
    sr = apply(sr, product);
    J = ideal sr;
    c = I+J;
    QR = R/c;
    QQR = QR/(annihilator((sum(gens(R)))_QR));
    print(numgens(source((basis(QQR)))));
    if (numgens(source(basis(2,QQR)))==numgens(source(basis(1,QQR)))) then print "ONTO" else print "NOT ONTO";
    ell = basis(1,QQR);
    QQRc = QQR[vars(0..(numgens(source(ell))-1))];
    x = vars(QQRc)*(transpose(ell));
    tf = (x_(0,0))^3;
    T = factor tf;
    print(T);
    L = {};
    for i from 0 to (#T)-1 do (
    L = append(L, (T#i#0)^(T#i#1));
    );
    L = select(L, x->((degree(x))#0 > 0));
    filename = concatenate(pname,"milnor.m2");
    filename << "clearAll;" << endl;
    filename << "v=" << gens(QQRc) << ";" << endl;
    filename << "R=QQ[v];" << endl;
    filename << "p=" << toString(product(L)) << ";" << endl;
    filename << "I = ideal(p); J = ideal(jacobian(I)); hilbertSeries J" << endl;
    filename << close;
    );'''
    with open(filestr1, 'w') as f:
        f.write('clearAll;');
        f.write('pname = "'+pname1+'"; \n');
        f.write(commands+"\n")
        f.write(macaulay2format(a)+"; \n");
        f.write("v = {");
        f.write(pv[0])
        for i in range(1, len(pv)):
            f.write(","+pv[i])
        f.write("};\n")
        f.write('sr="{');
        for i in range(len(relations[0])):
            if i>0:
                f.write(',')
            f.write('{')
            for j in range(0,len(relations[0][i])):
                if relations[0][i][j]==1:
                    if j > (list(relations[0][i]).index(1)):
                        f.write(',')
                    f.write(pv[j])
            f.write('}')
        f.write('}"; \n')
        f.write('ff(v, M, sr);')
    with open(filestr2, 'w') as f:
        f.write('clearAll;')
        f.write('pname = "'+pname2+'"; \n');
        f.write(commands+"\n")
        f.write(macaulay2format(a)+"; \n");
        f.write("v = {");
        f.write(pv[0])
        for i in range(1, len(pv)):
            f.write(","+pv[i])
        f.write("};\n")
        f.write('sr="{');
        for i in range(len(relations[1])):
            if i>0:
                f.write(',')
            f.write('{')
            for j in range(0,len(relations[1][i])):
                if relations[1][i][j]==1:
                    if j > (list(relations[1][i]).index(1)):
                        f.write(',')
                    f.write(pv[j])
            f.write('}')
        f.write('}"; \n')
        f.write('ff(v, M, sr);')

getmaxfans(pvarsC3, prelationsC3, "C3")
getmaxfans(pvarsD3, prelationsD3, "D3")
getmaxfans(pvarsD9, prelationsD9, "D9")
getmaxfans(pvarsD16, prelationsD16, "D16")
getmaxfans(pvarsG1, prelationsG1, "G1")
getmaxfans(pvarsG3, prelationsG3, "G3")
getmaxfans(pvarsG4, prelationsG4, "G4")
getmaxfans(pvarsG5, prelationsG5, "G5")
getmaxfans(pvarsH2, prelationsH2, "H2")
getmaxfans(pvarsH6, prelationsH6, "H6")
getmaxfans(pvarsH10, prelationsH10, "H10")
getmaxfans(pvarsI3, prelationsI3, "I3")
getmaxfans(pvarsI8, prelationsI8, "I8")
getmaxfans(pvarsI11, prelationsI11, "I11")
getmaxfans(pvarsJ1, prelationsJ1, "J1")
getmaxfans(pvarsJ2, prelationsJ2, "J2")
getmaxfans(pvarsK2, prelationsK2, "K2")
getmaxfans(pvarsZ1, prelationsZ1, "Z1")
getmaxfans(pvarsZ2, prelationsZ2, "Z2")
getmaxfans(pvarsM1, prelationsM1, "M1")
getmaxfans(pvarsM4, prelationsM4, "M4")
getmaxfans(pvarsM, prelationsM, "M5")
getmaxfans(pvarsP341, prelationsP341, "P341")
getmaxfans(pvarsR1, prelationsR1, "R1")
getmaxfans(pvarsR2, prelationsR2, "R2")
getmaxfans(pvarsR3, prelationsR3, "R3")
#getmaxfans(pvarsP358ii, prelationsP358ii, "P358ii")
getmaxfans(pvarsP358iii, prelationsP358iii, "P358iii")
getmaxfans(pvarsW, prelationsW, "W")
