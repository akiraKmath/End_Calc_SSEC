import sys
import itertools

def ell_path_to_Ep(E, ell, ell_set, max_length, reject_jinvs):
    F = E.base_field()
    p = F.characteristic()

    # polynomials on F for modular poly.
    Rxy.<X,Y> = PolynomialRing(F)
    Ry.<Y> = PolynomialRing(F)

    #It uses the database downloded from https://aur.archlinux.org/packages/sage-data-kohel
    PHI = ClassicalModularPolynomialDatabase()
    ModPoly = dict([(ell,PHI[ell](X,Y)) for ell in ell_set])
    jt = j0 = E.j_invariant()
    j0p = j0^p
    dont_reverse1 = 0
    trail1 = [j0]
    degs1 = []
    pr_degs1 = 1
    while True:

        # Update jt with avoiding 0,1728
        jt_dummy = F([1728])
        while jt_dummy == F([1728]) or jt_dummy == F([0]) or jt_dummy == F([287496]):
            jt_dummy = choice(ModPoly[ell](jt,Y).roots())[0] #It picks a random root of modular poly.
        jt = jt_dummy

        if jt == j0p and jt not in reject_jinvs:#Find a path from j0 to j0p
            #print("trail1 : ", trail1)
            if j0 == j0p and len(trail1) == 2:
                continue
            trail1.append(jt)
            rj = jt
            degs1.append(ell)
            dont_reverse1 = 1
            break
        
        #cut self-loop or self-cycle
        if jt in trail1:
            temp_ind = trail1.index(jt)
            trail1 = trail1[:temp_ind+1]
            degs1 = degs1[:temp_ind]
            continue
        if jt^p in trail1:
            temp_ind = trail1.index(jt.conjugate())
            trail1 = trail1[:temp_ind+1]
            degs1 = degs1[:temp_ind]
            continue

        trail1.append(jt)
        degs1.append(ell)
        pr_degs1 *= ell
        #print("check246 : ", len(trail1), max_length)

        ### length check
        # print("trail1 : ", trail1)
        if pr_degs1.nbits() > max_length:
            return False

        ### termination check ###
        Check = 0
        for el in ell_set:
            if ModPoly[el](jt,jt^p) == 0 and el != ell and jt not in reject_jinvs:
                #print("jt : ", jt, el)
                rj = jt
                Check = 1
                degs1.append(el)
                break
        if Check == 1:
            break

    #print("trail1 : ", trail1)
    if j0^p == j0:#In this case,it is enough to only get the first path
        if dont_reverse1 == 0: #Connect to the path of j.conjugate()

            trail1p = list(reversed(trail1))
            for j in trail1p:
                trail1.append(j.conjugate())

            degs1p = list(reversed(degs1))
            temp = len(degs1p)
            for ind in range(1,temp):
                degs1.append(degs1p[ind])

        return [trail1,degs1], rj

    if dont_reverse1 == 0:

        trail1p = list(reversed(trail1))
        for j in trail1p:
            #print("j : ", j, j.conjugate(), j.parent())
            trail1.append(j.conjugate())

        degs1p = list(reversed(degs1))
        temp = len(degs1p)
        for ind in range(1,temp):
            degs1.append(degs1p[ind])

    return [trail1,degs1], rj

def cycles_pair(E, ell_set, max_length, setcard, rejects):
    f_set = []
    chars = itertools.cycle(r'/-\|')
    print("finding path : ", end = "")
    while (1):
        path = ell_path_to_Ep(E, ell_set[0], ell_set, max_length/2, rejects[0])
        sys.stdout.write('\b'+next(chars))
        sys.stdout.flush() 
        if path != False:
            if path[0] not in f_set:
                f_set.append(path[0])
                rejects[0].append(path[1])
                sys.stdout.write('\b'+' ')
                sys.stdout.flush() 
                print(repr(ell_set[0]) + "-path found : ", len(f_set), " / " + repr(setcard))
                print("finding path : ", end = "")
        if len(f_set) >= setcard:
            break;
    mm = max([prod(pp[1]).nbits() for pp in f_set])
    #print("mm = ", mm)
    s_set = []
    while (1):
        path = ell_path_to_Ep(E, ell_set[1], ell_set, max_length-mm, rejects[1])
        sys.stdout.write('\b'+next(chars))
        sys.stdout.flush() 
        if path != False:
            sys.stdout.write('\b'+' ')
            sys.stdout.flush() 
            if path[0] not in s_set:
                s_set.append(path[0])
                rejects[1].append(path[1])
                print(repr(ell_set[1]) + "-path found : ", len(s_set), " / " + repr(setcard))
                print("finding path : ", end = "")
        if len(s_set) >= setcard:
            sys.stdout.write('all done \n')
            sys.stdout.flush() 
            break;

    return f_set, s_set

def cycles_by_pair(f_set,s_set, max_len, p):
    cycles = []
    for f in f_set:
        for s in s_set:
            cycles.append([list(f[0][:-1])+list(reversed(s[0])),f[1]+list(reversed(s[1]))])
    for i in range(len(f_set)):
        f0 = f_set[i]
        for j in range(i+1, len(f_set)):
            f1 = f_set[j]
            if f0 != f1:
                cycles.append([list(f0[0][:-1])+list(reversed(f1[0])),f0[1]+list(reversed(f1[1]))])
    for i in range(len(s_set)):
        f0 = s_set[i]
        for j in range(i+1, len(s_set)):
            f1 = s_set[j]
            if f0 != f1:
                cycles.append([list(f0[0][:-1])+list(reversed(f1[0])),f0[1]+list(reversed(f1[1]))])
    cycles.sort(key=lambda x: prod(x[1]))
    if p.nbits() > 20:
        ncycles = []
        for cycle in cycles:
            cyc_deg = prod(cycle[1]).nbits()
            if cyc_deg >= 5 and cyc_deg <= max_len:
                ncycles.append(cycle)
    
    else:
        ncycles = cycles
    print("degs of cycles : ", [prod(c[1]).nbits() for c in ncycles])
    return ncycles

def cycles_append(E, ell_set, f_set, s_set, max_length, rejects):
    chars = itertools.cycle(r'/-\|')
    print("finding extra path : ", end = "")
    while (1):
        path = ell_path_to_Ep(E, ell_set[0], ell_set, max_length/2, rejects[0])
        sys.stdout.write('\b'+next(chars))
        sys.stdout.flush()
        if path != False:
            if path[0] not in f_set:
                sys.stdout.write('\b'+' ')
                sys.stdout.flush()
                f_set.append(path[0])
                rejects[0].append(path[1])
                print("path found")
                break;
    print("finding extra path : ", end = "")
    while (1):
        path = ell_path_to_Ep(E, ell_set[1], ell_set, max_length/2, rejects[1])
        sys.stdout.write('\b'+next(chars))
        sys.stdout.flush()
        if path != False:
            if path[0] not in s_set:
                sys.stdout.write('\b'+' ')
                sys.stdout.flush()
                s_set.append(path[0])
                rejects[1].append(path[1])
                print("path found")
                break;

    return f_set, s_set