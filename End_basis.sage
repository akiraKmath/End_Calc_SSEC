def first_cycle_conversion_prime(E, p, cycle, deg, is_elkies):
    try:
        cycle_a = trace_cal(E, cycle, deg, is_elkies, False)
        R = cycle[0][0].parent()
        x,y = R.gens()
        cycle_c = trace_cal(E, cycle, deg, is_elkies, True)
        cycle_c = -cycle_c/p
        if deg - cycle_a^2-p*(cycle_c^2) == 0:
            return False
        Q = DiagonalQuadraticForm(QQ, [1,p])
        print("tenary_test : ", tenary_quadratic_solve(1, p, -(deg - cycle_a^2-p*(cycle_c^2))))
        try:
            sov = Q.solve(deg-cycle_a^2-p*(cycle_c^2))
            a = cycle_a; b = sov[0]; c = cycle_c; d = sov[1]
            if (a^2+b^2+p*(c^2+d^2))==deg:
                    return (a,b,c,d)
            else:
                return False;
        except ArithmeticError:
            return False

    except ZeroDivisionError:
        return False

def first_cycle_conversion_prime_v2(E, p, cycle, deg, prime_point, is_elkies):
    #print("prime_point in fccp: ", prime_point)
    try:
        cycle_a = trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        R = cycle[0][0].parent()
        x,y = R.gens()
        cycle_c = trace_cal(E, cycle, deg, prime_point, is_elkies, True)
        cycle_c = -cycle_c/p
        if deg - cycle_a^2-p*(cycle_c^2) == 0:
            return False
        
        sov = tenary_quadratic_solve(1, p, -(deg - cycle_a^2-p*(cycle_c^2)))
        if sov != False:
            a = cycle_a; b = sov[0]; c = cycle_c; d = sov[1]
            if (a^2+b^2+p*(c^2+d^2))==deg:
                print("deg = ", deg.nbits())
                return (a,b,c,d)
            else:
                return False;
        else: return False

    except ZeroDivisionError:
        return False

def first_cycle_conversion_ext(E, p, cycle, deg, is_elkies):
    try:
        cy_a =  trace_cal(E, cycle, deg, is_elkies, False)
        cy_a = QQ(cy_a)
        if cy_a == None:
            return False
        if deg-cy_a^2 == 0:
            return False
        print("tenary_test : ", tenary_quadratic_solve(1, p, -(deg - cy_a^2)))
        Q = DiagonalQuadraticForm(QQ,  [1,p])
        try:
            sov = Q.solve(deg-cy_a^2)
            a = cy_a; b = sov[0]; c = 0; d = sov[1]
            if (a^2+b^2+p*(c^2+d^2))==deg:
                return (a,b,c,d)
            else:
                return False;
        except ArithmeticError:
            # print("here")
            return False;
    except ZeroDivisionError:
        return False;

def first_cycle_conversion_ext_v2(E, p, cycle, deg, prime_point, is_elkies):
    #print("prime_point in fcce: ", prime_point)
    try:
        cy_a =  trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        cy_a = QQ(cy_a)
        if cy_a == None:
            return False
        if deg-cy_a^2 == 0:
            return False
        sov = tenary_quadratic_solve(1, p, -(deg - cy_a^2))
        if sov != False:
            a = cy_a; b = sov[0]; c = 0; d = sov[1]
            if (a^2+b^2+p*(c^2+d^2))==deg:
                print("deg = ", deg.nbits())
                return (a,b,c,d)
            else:
                return False;
        else: return False
    except ZeroDivisionError:
        return False;

def second_cycle_conversion_prime(E, p, cycle, deg, cycle0, deg0, coefs0, prime_point, is_elkies):
    a0,b0,c0,d0 = coefs0
    try:
        cycle_a = trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        tr = trace_cal(E, cycle0+cycle, deg*deg0, prime_point, is_elkies, False)
        R = cycle[0][0].parent()
        x,y = R.gens()
        cycle_c = trace_cal(E, cycle, deg, prime_point, is_elkies, True)
        trc = trace_cal(E, cycle0+cycle, deg0*deg, prime_point, is_elkies, True) 
        cycle_c = -cycle_c/p; trc = -trc/p
        if deg - cycle_a^2-p*(cycle_c^2) == 0:
            return False
        A = tr-a0*cycle_a+p*c0*cycle_c
        C = trc-a0*cycle_c-cycle_a*cycle_c
        nrm0 = deg0-a0^2-p*c0^2
        nrm1 = deg-cycle_a^2-p*(cycle_c^2)
        var('xx yy')
        f = d0*xx-b0*yy-C;
        g = -b0*xx-p*d0*yy-A
        ans = solve([f,g], xx, yy, solution_dict = True)
        try:
            a1 = cycle_a; b1 = ans[0][xx]; c1 = cycle_c; d1 = ans[0][yy]
            if (b1^2+p*c1^2+p*d1^2)==nrm1 and (-b0*b1-p*c0*c1-p*d0*d1) == A:
                return (a1,b1,c1,d1);
        except ArithmeticError:
            # print("here")
            return False

    except ZeroDivisionError:
        return False

def second_cycle_conversion_ext(E, p, cycle, deg, cycle0, deg0, coefs0, is_elkies):
    a0,b0,c0,d0 = coefs0
    try:
        a1 = trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        tr = trace_cal(E, cycle0+cycle, deg*deg0, prime_point, is_elkies, False)
        if deg - a1^2 == 0 or tr == None or a1 == None:
            return False
        A = tr-a0*a1
        nrm0 = deg0-a0^2
        nrm1 = deg-a1^2
        alpha = (A*c0)/(nrm0)
        beta = (A*d0)/(nrm0)
        eq = nrm1*b0^2+alpha^2*(b0^2+c0^2*p)*p+beta^2*(b0^2+d0^2*p)*p+2*c0*d0*p^2*alpha*beta-A^2
        rtcd = prod([f[0]^(f[1]/2) for f in (c0^2+d0^2).factor()])
        print("tenary_test : ", tenary_quadratic_solve(b0^2*p, p*nrm0, -(eq)))
        Q2 = DiagonalQuadraticForm(QQ,  [b0^2*p, p*nrm0])
        try:
            sov = Q2.solve(eq)
            prexy = ((1/rtcd)*(sov[0]*d0+sov[1]*c0)-alpha, (1/rtcd)*(-sov[0]*c0+sov[1]*d0)-beta)
            b1 = -((A+c0*p*prexy[0]+d0*p*prexy[1])/b0); c1 = prexy[0]; d1=prexy[1]
            if (b1^2+p*c1^2+p*d1^2)==nrm1 and (-b0*b1-p*c0*c1-p*d0*d1) == A:
                return (a1,b1,c1,d1);
        except ArithmeticError:
            return False;

    except ZeroDivisionError:
        return False

def second_cycle_conversion_ext_v2(E, p, cycle, deg, cycle0, deg0, coefs0, prime_point, is_elkies):
    a0,b0,c0,d0 = coefs0
    try:
        a1 = trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        tr = trace_cal(E, cycle0+cycle, deg*deg0, prime_point, is_elkies, False)
        if deg - a1^2 == 0 or tr == None or a1 == None:
            return False
        A = tr-a0*a1
        nrm0 = deg0-a0^2
        nrm1 = deg-a1^2
        alpha = (A*c0)/(nrm0)
        beta = (A*d0)/(nrm0)
        eq = nrm1*b0^2+alpha^2*(b0^2+c0^2*p)*p+beta^2*(b0^2+d0^2*p)*p+2*c0*d0*p^2*alpha*beta-A^2
        rtcd = prod([f[0]^(f[1]/2) for f in (c0^2+d0^2).factor()])
        sov = tenary_quadratic_solve(b0^2*p, p*nrm0, -(eq))
        if sov != False:
            prexy = ((1/rtcd)*(sov[0]*d0+sov[1]*c0)-alpha, (1/rtcd)*(-sov[0]*c0+sov[1]*d0)-beta)
            b1 = -((A+c0*p*prexy[0]+d0*p*prexy[1])/b0); c1 = prexy[0]; d1=prexy[1]
            if (b1^2+p*c1^2+p*d1^2)==nrm1 and (-b0*b1-p*c0*c1-p*d0*d1) == A:
                print("deg = ", deg.nbits())
                return (a1,b1,c1,d1);
            else: return False
        else: return False

    except ZeroDivisionError:
        return False

def third_cycle_conversion_ext(E, p, cycle, deg, cycle0, deg0, coefs0, cycle1, deg1, coefs1, prime_point, is_elkies):
    a0,b0,c0,d0 = coefs0; a1,b1,c1,d1 = coefs1; 
    try:
        a2 = trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        tr02 = trace_cal(E, cycle0+cycle, deg*deg0, prime_point, is_elkies, False)
        tr12 = trace_cal(E, cycle1+cycle, deg*deg1, prime_point, is_elkies, False)
        if deg - a2^2 == 0 or tr02 == None or tr12 == None:
            return False
        var('x y z')
        f = a0*a2 - b0*x - p*c0*y - p*d0*z-tr02
        g = a1*a2 - b1*x - p*c1*y - p*d1*z-tr12
        h_sub = a2^2 + x^2 + p*y^2 +p*z^2 - deg
        ans = solve([f,g,h_sub], x, y, z, solution_dict = True)
        if Matrix(QQ, [(ans[0][x], ans[0][y], ans[0][z]),(-b0, -c0, -d0),(-b1, -c1, -d1)]).rank()<3:
            return False
        return (a2, ans[0][x], ans[0][y], ans[0][z])

    except TypeError:
        return False

def other_cycle_conversion_ext(E, p, cycle, deg, cycle0, deg0, coefs0, cycle1, deg1, coefs1, cycle2, deg2, coefs2, prime_point, is_elkies):
    a0, b0, c0, d0 = coefs0; a1, b1, c1, d1 = coefs1; a2, b2, c2, d2 = coefs2;
    try:
        a3 = trace_cal(E, cycle, deg, prime_point, is_elkies, False)
        tr03 = trace_cal(E, cycle0+cycle, deg0*deg, prime_point, is_elkies, False)
        tr13 = trace_cal(E, cycle1+cycle, deg1*deg, prime_point, is_elkies, False)
        tr23 = trace_cal(E, cycle2+cycle, deg2*deg, prime_point, is_elkies, False)
        if a3 == None or tr03 == None or tr13 == None or tr23 == None or deg == a3^2:
            return False;
        tM = Matrix(QQ,[[b0,-p*c0,-p*d0,a0*a3-tr03],[b1,-p*c1,-p*d1,a1*a3-tr13],[b2,-p*c2,-p*d2,a2*a3-tr23]])
        if tM.rank() < 3:
            print("here: line84")
            return False
        var('x y z')
        f = a0*a3 - b0*x - p*c0*y - p*d0*z-tr03
        g = a1*a3 - b1*x - p*c1*y - p*d1*z-tr13
        h = a2*a3 - b2*x - p*c2*y - p*d2*z-tr23
        ans = solve([f,g,h], x, y, z, solution_dict = True)
        if ans[0][z] not in QQ:
            ans[0].substitute(r1 = 1)
        return (a3, ans[0][x], ans[0][y], ans[0][z])
    except ZeroDivisionError:
        return False;

#======================================================================
def QArep(cofs, p):
    QA.<i,j,k> = QuaternionAlgebra(QQ, -1,-p)
    a = cofs[0];
    b = cofs[1];
    c = cofs[2];
    d = cofs[3];
    x = QQ(a) + QQ(b)*i + QQ(c)*j + QQ(d)*k
    return x;

#======================================================================

#======================================================================

def product_phi_mat_ver_qua(baslis):
    prelist = baslis
    last = []
    for x in prelist:
        inlist = []
        for y in prelist:
            inlist.append(2*(x*(y.conjugate()))[0])

        last.append(inlist)

    return matrix(last)

#======================================================================

#======================================================================

def is_maximal_prime(totos, k, p):
    korpp = totos
    for i in range(k):
        for j in range(i+1,k):
            karix = QArep(totos[i], p)
            kariy = QArep(totos[j], p)
            kari = karix*kariy
            korpp.append([kari[0],kari[1],kari[2],kari[3]])
    korpp.append([1,0,0,0])
    korpp.append([0,0,1,0])  

    M1 = Matrix(QQ, korpp)
    M2 = M1.LLL()
    baslis = []
    for cof in M2:
        x = QArep(cof, p)
        if x != 0:
            baslis.append(x)
    N = product_phi_mat_ver_qua(baslis)
    if N.det() != 0:
        print("discriminant of subring of End(E): ", factor(N.det()))
        return N.det(), baslis

def is_maximal_ext(totos, k, p):
    korpp = totos
    for i in range(k):
        for j in range(i+1,k):
            karix = QArep(totos[i], p)
            kariy = QArep(totos[j], p)
            kari = karix*kariy
            korpp.append([kari[0],kari[1],kari[2],kari[3]])
    korpp.append([1,0,0,0]) 

    M1 = Matrix(QQ, korpp)
    M2 = M1.LLL()
    baslis = []
    for cof in M2:
        x = QArep(cof, p)
        if x != 0:
            baslis.append(x)
    N = product_phi_mat_ver_qua(baslis)
    if N.det() != 0:
        print("discriminant of subring of End(E): ", factor(N.det()))
        return N.det(), baslis



def End_basis_prime(E, p, ell_set, bound, collect_nums, is_elkies):
    j0 = E.j_invariant()
    check = 0
    cycles = []; coeflist = []
    found_cycles = []; done_cycles = []
    rejects = [[],[]]
    prime_point = prime_point_collect(E, bound)
    while (1):
        if check == 0:
            cycle_set = cycles_pair(E, ell_set, bound/3, collect_nums,rejects)
            prime_cycle_set = cycle_set[0] + cycle_set[1]
            for cycle in prime_cycle_set:
                rats = rats_isogeny(E, cycle)
                if rats != False:
                    found_cycles.append([rats, prod(cycle[1])])
            print(len(found_cycles))
            if len(found_cycles) != 0:
                check = 1
        else:
            cycle_set = cycles_pair(E, ell_set, bound/3, collect_nums, rejects)
            prime_cycle_set = cycle_set[0] + cycle_set[1]
            for cycle in prime_cycle_set:
                rats = rats_isogeny(E, cycle)
                if rats != False and [rats, prod(cycle[1])] not in done_cycles:
                    found_cycles.append([rats, prod(cycle[1])])
        for cycle in found_cycles:
            if len(coeflist) == 0:
                if  cycle not in done_cycles:
                    cycle0 = cycle[0]; deg0 = cycle[1]
                    coefs0 = first_cycle_conversion_prime_v2(E, p, cycle0, deg0, prime_point, is_elkies)
                    if coefs0 != False:
                        if coefs0 != None:
                            coeflist.append(coefs0); 
                        done_cycles.append(cycle)
                        
            else:
                if cycle not in done_cycles:
                    cycle1 = cycle[0]; deg1 = cycle[1]
                    coefs1 = second_cycle_conversion_prime(E, p, cycle1, deg1, cycle0, deg0, coefs0, prime_point, is_elkies)
                    if coefs1 != False:
                        if coefs1 != None:
                            coeflist.append(coefs1); 
                        done_cycles.append(cycle)
                        
                        
            #print("coeflist = ", coeflist)
            if len(coeflist) >= 2:
                det, basis = is_maximal_prime(copy(coeflist), len(coeflist), p)

                if det == p^2:
                    print()
                    print("basis : ", basis);
                    print()
                    print("#generators = ", len(coeflist))
                    return basis, len(coeflist)
            print()

##===============================================================##


def End_basis_ext(E, p, ell_set, bound, collect_nums, is_elkies):
    j = E.j_invariant(); ellset0 = copy(ell_set)
    cycles = []; coflist = []
    trys = 1; predet = 0; check = 0
    pretrys = 1; cycs_cand = []; removed_cycs = []
    cyc_find_time = 0
    rejects = [[],[]]
    prime_point = prime_point_collect(E, bound)
    #print("prime_point = ", prime_point)
    while(1):
        st = time.time()
        if check == 0:
            print()
            sets = cycles_pair(E, ell_set, bound/3, collect_nums, rejects)
            pre_cycles = cycles_by_pair(sets[0], sets[1], bound/3, p)
            print()
            for cycle in pre_cycles:
                rats = rats_isogeny(E, cycle)
                if rats != False:
                    cycs_cand.append([rats, prod(cycle[1])])
        else:
            mm = max([c[1].nbits() for c in cycles])
            print(mm)
            old_cycles = copy(cycs_cand)
            sets = cycles_append(E, ell_set, sets[0], sets[1], bound-mm,rejects)
            pre_cycles = cycles_by_pair(sets[0], sets[1], bound-mm, p)
            for cycle in pre_cycles:
                if cycle not in cycs_cand and cycle not in cycles:
                    rats = rats_isogeny(E, cycle)
                    if rats != False:
                        cycs_cand.append([rats, prod(cycle[1])])
            for c in old_cycles:
                cycs_cand.remove(c)
        et = time.time()-st
        cyc_find_time += et
        print("pair found : ", cyc_find_time)
        if len(coflist) == 0:
            for cycle in cycs_cand:
                cycle0 = cycle[0]; deg0 = cycle[1]
                print("degree of cycle : ", cycle[1].nbits())
                cofs = first_cycle_conversion_ext_v2(E, p, cycle0, deg0, prime_point, is_elkies)
                if cofs != False and cofs[3] != 0:
                    coflist.append(cofs)
                    rv_cycle = cycle
                    check = 1
                    break;
                print()
            if check == 1:
                cycles.append(rv_cycle)
                cycs_cand.remove(rv_cycle)
            #print("coflist : ", coflist)
        
        if check == 1:
            for cycle in cycs_cand:
                print()
                print("degree of cycle : ", cycle[1].nbits())
                subcheck = 0
                if cycle not in cycles:
                    if len(coflist) == 1:
                        cycle1 = cycle[0]; deg1 = cycle[1]
                        print("second cycle cof calc")
                        cofs = second_cycle_conversion_ext_v2(E, p, cycle1, deg1, cycle0, deg0, coflist[0], prime_point, is_elkies)
                        
                    elif len(coflist) == 2:
                        cycle2 = cycle[0]; deg2 = cycle[1]
                        print("third cycle cof calc")
                        cofs = third_cycle_conversion_ext(E, p, cycle2, deg2, cycle0, deg0, coflist[0], cycle1, deg1, coflist[1], prime_point, is_elkies)
                        
                    else:
                        cycle3 = cycle[0]; deg3 = cycle[1]
                        print("other cycle cof calc")
                        cofs = other_cycle_conversion_ext(E, p, cycle3, deg3, cycle0, deg0, coflist[0], cycle1, deg1, coflist[1], cycle2, deg2, coflist[2], prime_point, is_elkies)
                        
                    
                    if cofs != False:
                            coflist.append(cofs)
                            rv_cycle = cycle
                            subcheck = 1
                
                if subcheck == 1:
                    cycles.append(rv_cycle)
                    cycs_cand.remove(rv_cycle)

                #print("coflist : ", coflist)
                if len(coflist) >= 2:
                    det, baslis = is_maximal_ext(copy(coflist), len(coflist), p)
                    predet = det
                    if det == p^2:
                        print()
                        print("basis : ", baslis);
                        print()
                        print("#generators = ", len(coflist))
                        sum_of_len_of_cycles = 0
                        for res_cyc in cycles:
                            sum_of_len_of_cycles += len(res_cyc[0])
                        return baslis, len(coflist)






