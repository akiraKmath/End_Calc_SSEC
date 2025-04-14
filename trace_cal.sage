"""
An implementation program of computing the endomorphism ring of a given supersingular elliptic curve besed on finding isogeny cycle.
(C) 2025 Mitsubisi Electric, Rikkyo University, Created by Yuta Kambe, Akira Katayama, Kazuki Komine, Yusuke Aikawa, Yuki Ishihara, Masaya Yasuda, Kazuhiro Yokoyama.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
"""

import sys
import itertools

def ext_isogeny(E, phi, extension_degree):      ###同種写像列の基礎体の拡大
    ext_phi = []
    F = E.base_field()
    exF = F.extension(extension_degree)
    R = phi[0][0].parent()
    #print("R = ", R)
    for i in range(len(phi)):
        psi = (R(phi[i][0]), R(phi[i][1])); expsi = []
        for rat in psi:
            ratn = rat.numerator().change_ring(exF)
            ratd = rat.denominator().change_ring(exF)
            exR = ratd.parent()
            exQR = FractionField(exR)
            expsi.append(exQR(ratn)/exQR(ratd))
        ext_phi.append(expsi)
    return ext_phi

def prime_point_collect(E, ext_bound):
    print("ext_bound : ", ext_bound)
    chars = itertools.cycle(r'/-\|')
    print("finding ell-torsion points : ", end = "")
    qq = E.base_field().order(); pp = (E.cardinality()-(qq+1))/2 #pp = qq.factor()[0][0]
    prime_point = dict(); d = 0
    elkies_prime = Primes()[:27]
    a = E.a4(); b = E.a6(); T = prod(elkies_prime)
    card_E_all_factors = []
    while True:
        d += 1
        card_E = (qq^d+ 1 - (-1)^d*2*pp^d) #E(F_p^2d)
        card_Et = sqrt(qq^d+ 1 - (-1)^d*2*pp^d) #E(F_p^2d)
        #print("d = ", d)
        #sys.stdout.write('\b'+repr(d))
        #sys.stdout.flush()
        
        for card_E_factor in card_E_all_factors:  #E(F_p^(2^i)) | E(F_p^2^(i+1))
            while card_E % card_E_factor == 0:
                card_E //= card_E_factor
        
        card_E_factors = [ell[0] for ell in ZZ(sqrt(card_E)).factor()]

        for factor in card_E_factors:
            card_E_all_factors.append(factor) 
        use_ells = [ell for ell in card_E_factors if (ell < 10000000)] #min{p, 10000000}   
        
        if d > 1:
            exF = E.base_field().extension(d); q = exF.order()
            EE = EllipticCurve(exF, [a,b])
        else:
            EE = E
        check = 0 
        for ell in use_ells:
            if ell not in elkies_prime:
                while(1):
                    sys.stdout.write('\b'+next(chars))
                    sys.stdout.flush()
                    P = EE.random_element()
                    #print("P = ", P, P.order().factor())
                    preP = ZZ(card_Et/ell) * P
                    #print("preP = ", preP)
                    if preP != EE(0) and ell*preP == EE(0):
                        Q = preP
                        break;
                T *= ell; prime_point[ell] = Q
                if T.nbits() > (ext_bound/2+3):
                    check = 1
                    break
        if check == 1:
            sys.stdout.write('\b'+' ')
            sys.stdout.flush()
            print("points found")
            return prime_point
            

def isogeny_point_data(E, P, phi):           ###サイクルの座標計算
    xt = P.xy()[0]; yt = P.xy()[1]
    for i in range(len(phi)):
        if i != 0:          
            psi = (R(phi[i][0]), R(phi[i][1]))
        else:
            psi = phi[0]
        
        R = psi[0].parent()
        FF = R.base_ring()
        x,y = R.gens()
        xt0 = psi[0].numerator()(x = xt, y = yt); xt1 = psi[0].denominator()(x = xt, y = yt);
        yt0 = psi[1].numerator()(x = xt, y = yt); yt1 = psi[1].denominator()(x = xt, y = yt);
        xt =xt0*(FF(xt1).inverse_of_unit()); yt = yt0*(FF(yt1).inverse_of_unit())
        #print("xt, yt = ", xt, yt)
    else:
        PP = E((xt, yt))
    return PP

def ell_tor_prime(E, deg, T, done_ells):
    qq = E.base_field().order(); pp = (E.cardinality()-(qq+1))/2
    d = 0; limit = 4*sqrt(deg)
    card_E_all_factors = []
    use_ells_list = []
    while True:
        d += 1
        card_E = (qq^d+ 1 - (-1)^d*2*pp^d) #E(F_p^2d)
        #print("d = ", d)
        
        for card_E_factor in card_E_all_factors:  #E(F_p^(2^i)) | E(F_p^2^(i+1))
            while card_E % card_E_factor == 0:
                card_E //= card_E_factor
        
        card_E_factors = [ell[0] for ell in ZZ(sqrt(card_E)).factor()]

        for factor in card_E_factors:
            card_E_all_factors.append(factor) 
        use_ells = [ell for ell in card_E_factors if (ell < 10000000) and (ell not in done_ells)] #min{p, 10000000}    
        #use_ells = [ell for ell in card_E_factors if (ell < abs(pp)) and (ell not in done_ells)] #min{p, 10000000}
        for i in range(len(use_ells)):
            done_ells.append(use_ells[i])
            T *= use_ells[i]
            if T >= limit:
                use_ells = use_ells[:i+1]
                use_ells_list.append(use_ells)
                return d, use_ells_list
        use_ells_list.append(use_ells)

def Point_data(E, phi, deg, d, ells, prime_point, is_frobenius):     ###Schoofに用いる点などのデータの計算
    data = []; qq = E.base_field().order(); pp = (E.cardinality()-(qq+1))/2 #pp = qq.factor()[0][0]
    Et_card = sqrt(qq^d+ 1 - (-1)^d*2*pp^d) 
    #print("d = ", d, "\nells = ", ells, "\nprime_point = ", prime_point)
    if d > 1:
        a = E.a4(); b = E.a6(); exF = E.base_field().extension(d); q = exF.order()
        EE = EllipticCurve(exF, [a,b])
        exphi = ext_isogeny(E, phi, d)
    else:
        EE = E; exphi = phi
    for ell in ells:
        st = time.time()
        try:
            Q = prime_point[ell]
        except KeyError:
            while(1):
                P = EE.random_element()
                preP = ZZ(Et_card/ell) * P
                    
                if preP != EE(0) and ell*preP == EE(0):
                    Q = preP
                    break;
            
        et = time.time()-st
        if is_frobenius:
            p = abs(pp)
            phiQ = isogeny_point_data(EE, Q, exphi)
            phiQ = EE(phiQ.xy()[0]^p, phiQ.xy()[1]^p)
            dphiQ= isogeny_point_data(EE, phiQ, exphi)
            dphiQ = EE(dphiQ.xy()[0]^p, dphiQ.xy()[1]^p)
        else:
            phiQ= isogeny_point_data(EE, Q, exphi)
            dphiQ= isogeny_point_data(EE, phiQ, exphi)
        data.append([dphiQ, phiQ, Q, ell, EE])
    return data

def mtor(m, deg, Points): #Points = (R, Q, P, m, ExtE, ExtCycle), isogeny_traceの一部
    l_m = deg % m;
    if l_m > (m/2):
        l_m =  l_m - m;

    def BSGS(r,S,T,E):
        E0 = E
        bsm = round(sqrt(r))+1; R = bsm*S;
        L = list(range(bsm));
        tmp=E0(0);
        for i in range(bsm):
            L[i] = tmp; tmp += S;

        flag = 0; tmp = T;
        for j in range(bsm):
            if flag != 0:
                break;
            if tmp in L:
                for i in range(len(L)):

                    if tmp == L[i]:
                        I = i; J = j; flag = 1;
                        break;
            tmp -= R
        d = (I+J*bsm)%r;
        return d;

    def standard_trace(m, S, T, E):
        Rw = S
        for i in range(m):
            if Rw == T:
                return i
            Rw += S


    E = Points[4] #Eの演算用に表示
    pi = Points[1]
    pis = Points[0]
    l_mP = l_m * Points[2]
    pi_d = pis + l_mP
    jo = BSGS(m,pi,pi_d,E)
    #print("Schoof_ check: ", pis - jo*pi + l_m * Points[2])
    if jo > m/2:
        jo = jo-m;
    return jo;

def elliptic_quotient(E, phi):
    F = E.base_field(); a = E.a4(); b = E.a6(); R = phi[0].numerator().parent()
    s, u = R.gens()
    def_eq = s^3+a*s+b
    phiss = [phi[0].numerator(), phi[0].denominator(), phi[1].numerator(), phi[1].denominator()]
    nphiss = []
    for i in range(len(phiss)):
        phi_part = phiss[i]
        lis_phi_part = list(phi_part)
        if phi_part == 1:
            nphiss.append(R(1))
        else:
            nlis = R(0)
            for cor in lis_phi_part:
                d = cor[1]; y_sq = 0
                while True:
                    if (d/(u^2)).denominator() == R(1):
                        d = d / (u^2)
                        y_sq += 1
                    else: break;
                nlis += cor[0]*d*((def_eq)^y_sq)
            nphiss.append(nlis)
    NR.<s> = F[]
    nphiss = [NR(nphi_p(x=s, y= 1)) for nphi_p in nphiss]
    return nphiss

def Schoof_Hybrid_data(E, phi, deg, prime_point):  #RandomSamplingで使う素数を事前に集める．Elkies素数はそれ以外のものを使う．
    #========== time data ===============================
    Elkies_list = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
    Elkies_time = [0.46527299284935, 0.43057015538215637, 0.5715641677379608, 0.7338555256525675, 0.6296719551086426, 0.8843600273132324, 1.0611517429351807, 1.5139443635940553, 1.6602551698684693, 2.680908226966858, 3.0691380977630613, 3.2790233850479127, 3.8630520343780517, 5.0347225904464725, 6.481257033348084, 6.9038142442703245, 8.45649094581604, 9.379428505897522, 9.890097093582153, 11.766771030426025, 13.093309259414672, 15.296850228309632, 17.797670602798462, 19.846211624145507, 21.15657832622528, 23.242745447158814, 24.819961524009706, 27.34380567073822]
    exdeg_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    exdeg_time = [0.056093454360961914, 0.3169085184733073, 0.4095543920993805, 0.5012552056993756, 0.9880895614624023, 1.5507220824559529, 3.100886901219686, 4.010195398330689, 8.352376752429539, 14.300596740510729, 38.491298834482826, 83.52607655525208]
    #=====================================================
    print("degree : ", deg.nbits())
    #print("prime_point in SHd: ", prime_point)
    F = E.base_field()
    T = 1; S = []; Z = []; R.<t> = F[]
    done_ells = [ell[0] for ell in deg.factor()]  
    mod_phi = []
    for psi in phi:
        npsi = elliptic_quotient(E, psi)
        mod_phi.append(npsi)
    exdeg, use_ells_list = ell_tor_prime(E, deg, T, done_ells) #collecting prime for RandomSampling

    count = 0
    for d in range(1, exdeg+1):
        #================== do Elkies while time of Elkies prime less than extension degree d time =======
        if d <= 12:
            #print()
            time_sum = Elkies_time[count]
            while time_sum < exdeg_time[d-1]:
                if Elkies_list[count] not in done_ells:
                    ell = Elkies_list[count]
                    check = Elkies_test(E, ell)
                    if check == True:
                        print("Elkies prime : ", ell); 
                        st = time.time()
                        tr = trace_mod_elkies(E, deg, mod_phi, ell)
                        endt = time.time()-st
                        done_ells.append(ell)
                        #print(ell, "time = " , endt)
                        T *= ell
                        a = tr
                        S.append(a); Z.append(ell)
                        count += 1
                        time_sum += Elkies_time[count]
                        if T > 4*sqrt(deg):
                            trace = crt(S, Z)
                            if trace > T/2:
                                trace = trace - T
                            return trace
                    else: count += 1
                else: count += 1
            #print()
        #===============================================================================================

        #========== Random Sampling ============================================================
        stime = time.time()
        if use_ells_list[d-1] == []: continue
        print("extension degree : ", d)
        Data = Point_data(E, phi, deg, d, use_ells_list[d-1], prime_point, False) 
        etime = time.time() - stime
        #print("time(RandomSampling) =", etime)
        if Data != None:
            for Point in Data:
                print("ell = ", Point[3])
                done_ells.append(Point[3])
                st = time.time()
                a = mtor(Point[3], deg, Point)
                et = time.time() - st
                S.append(a); Z.append(Point[3])
                T *= Point[3]
                #print("progress : ", T.nbits(), "/", round(RR(4*sqrt(deg))).nbits())
                #print("progress : ", int(T.nbits()/round(RR(4*sqrt(deg))).nbits()*100), "%")
            etime = time.time() - stime
            if T > 4*sqrt(deg): 
                trace = crt(S, Z)
                if trace > T/2:
                    trace = trace - T
                return trace
        #============================================================================================

def Schoof_data(E, phi, deg, prime_point, is_frobenius):  #RandomSamplingで使う素数を事前に集める．Elkies素数はそれ以外のものを使う．
    print("degree : ", deg.nbits())
    F = E.base_field()
    T = 1; S = []; Z = []; R.<t> = F[]
    done_ells = [ell[0] for ell in deg.factor()]  
    mod_phi = []
    exdeg, use_ells_list = ell_tor_prime(E, deg, T, done_ells) #collecting prime for RandomSampling
    count = 0
    for d in range(1, exdeg+1):
        #========== Random Sampling ============================================================
        stime = time.time()
        if use_ells_list[d-1] == []: continue
        print("extension degree : ", d)
        Data = Point_data(E, phi, deg, d, use_ells_list[d-1], prime_point, is_frobenius) 
        etime = time.time() - stime
        #print("time(RandomSampling) =", etime)
        if Data != None:
            for Point in Data:
                print("ell = ", Point[3])
                done_ells.append(Point[3])
                st = time.time()
                a = mtor(Point[3], deg, Point)
                et = time.time() - st
                S.append(a); Z.append(Point[3])
                T *= Point[3]
                #print("progress : ", T.nbits(), "/", round(RR(4*sqrt(deg))).nbits())
                #print("progress : ", int(T.nbits()/round(RR(4*sqrt(deg))).nbits()*100), "%")
            etime = time.time() - stime
            if T > 4*sqrt(deg): 
                trace = crt(S, Z)
                if trace > T/2:
                    trace = trace - T
                return trace
        #============================================================================================

def trace_cal(E, cycle, deg, prime_point, is_elkies, is_frobenius):
    p = E.base_field().characteristic()
    if is_frobenius: return Schoof_data(E, cycle, deg*p, prime_point, is_frobenius)/2
    if is_elkies:
        return Schoof_Hybrid_data(E, cycle, deg, prime_point)/2
    else:
        return Schoof_data(E, cycle, deg, prime_point, False)/2
