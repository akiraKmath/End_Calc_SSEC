#based on http://www.sfb013.uni-linz.ac.at/reports/2004/pdf-files/rep_04-32_pilnikova.pdf


#SolubilityCertificate(可解性検証): done
def SC(b1,b2,b3):
    P_set = [p[0] for p in b3.factor()]
    rp_set = []
    for p in P_set:
        F = GF(p)
        r_sq = F(-b1/b2)
        if r_sq.is_square():
            rp_set.append(ZZ(sqrt(r_sq)))
        else:
            return False
    if rp_set == [] and P_set == []:
        return 1
    else:
        return crt(rp_set, P_set)

#find_sol_mod2(a1,a2,a3): done
def find_sol_mod2(a1,a2,a3):
    set_sol = [(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)]
    for sol in set_sol:
        if a1*sol[0]^2+a2*sol[1]^2+a3*sol[2]^2 % 2 == 0:
            return vector(ZZ, sol)
    return False
    
#find_sol_mod4(a1,a2,a3): done
def find_sol_mod4(a1,a2,a3):
    set_sol = [(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)]
    for sol in set_sol:
        if a1*sol[0]^2+a2*sol[1]^2+a3*sol[2]^2 % 4 == 0:
            return vector(ZZ, sol)
    return False

def rational_squarefree(b):
    b_num = b.numerator(); b_den = b.denominator()
    if b == 0:
        return True
    if b_num.is_squarefree() and b_den.is_squarefree():
        return True
    else: return False

#make_squarefree: 
def make_squarefree(c1,c2,c3):
    M = vector(QQ, [1,1,1])
    a1,a2,a3 = c1,c2,c3
    #while !((a1*a2*a3).is_squarefree()):
    if not (rational_squarefree(a1)):
        a1_square_primes = [(p[0],p[1]//2) for p in a1.factor() if p[1] % 2 == 0]
        a1 = a1//(prod([p[0]^(2*p[1]) for p in a1_square_primes]))
        M[0] = 1 / prod([p[0]^(p[1]) for p in a1_square_primes])
    if not (rational_squarefree(a2)):
        a2_square_primes = [(p[0],p[1]//2) for p in a2.factor() if p[1] % 2 == 0]
        a2 = a2//(prod([p[0]^(2*p[1]) for p in a2_square_primes]))
        M[1] = 1 / prod([p[0]^(p[1]) for p in a2_square_primes])
    if not (rational_squarefree(a3)):
        a3_square_primes = [(p[0],p[1]//2) for p in a3.factor() if p[1] % 2 == 0]
        a3 = a3//(prod([p[0]^(2*p[1]) for p in a3_square_primes]))
        M[2] = 1 / prod([p[0]^(p[1]) for p in a3_square_primes])
    if gcd(a1,a2) != 1:
        g3 = gcd(a1,a2)
        a1, a2, a3 = a1/g3, a2/g3, a3/g3
    if gcd(a2,a3) != 1:
        g1 = gcd(a2,a3)
        a1, a2, a3 = a1/g1, a2/g1, a3/g1
    if gcd(a3,a1) != 1:
        g2 = gcd(a3,a1)
        a1, a2, a3 = a1/g2, a2/g2, a3/g2
    return a1,a2,a3,M
            


#IsIsotropic3(3変数二次形式の有理数解の有無判定)
def find_sols_mod(c1,c2,c3,m):
    a0 = vector(ZZ,[-c3-c2,c1,c1])%abs(m)
    a1 = vector(ZZ,[c2,-c3-c1,c2])%abs(m)
    a2 = vector(ZZ,[c3,c3,-c1-c2])%abs(m)
    M = Matrix([a0,a1,a2])
    if M.rank() == 3:
        return  [a0,a1,a2]
    else:
        return [a0,a1]
    
def find_sols_mod4(cofs):
    sols0 = []
    sols1 = []
    sols = []
    flag = 0 
    check = [cofs[1]+cofs[2], cofs[0]+cofs[2], cofs[0]+cofs[1]]
    for i in range(3):
        if check[i] % 4 == 0:
            for j in range(4):
                for k in range(2):
                    pre_sol = vector(ZZ, [j,j,j])
                    pre_sol[i] = ZZ(2)*k
                    mm = [ZZ(cofs[t]*pre_sol[t]^2) for t in range(3)]
                    
                    if (mm[0]+mm[1]+mm[2]) % 4 == 0 and pre_sol != 0:
                        sols.append(pre_sol)
#             break;
            if flag == 0:
                sols0 = sols
                sols = []
                flag = 1
            else:
                sols1 = sols
    return sols0, sols1

def find_sol_gen(sol_a,sol_4, m):
    sols = []
    for sol0 in sol_a:
        for sol1 in sol_4:
            v = vector(ZZ, 3)
            for i in range(3):
                v[i] = crt([sol0[i],sol1[i]],[abs(m),4])
            sols.append(v%(4*abs(m)))
#     print(Matrix(sols).LLL())
    msl = Matrix(sols).LLL()
    return msl[-1], msl[-2], msl[-3], sols

def GSO_expand(B, n, cofs):
    GS = Matrix(QQ, n, 3)
    mu = Matrix(QQ, n, n)
#     print("n = ", n)
    for i in range(n):
#         print("i = ", i)
        GS[i] = B[i]
        mu[i,i] = 1
        for j in range(i):
            try:
                mu[i,j] = (cofs[0]*B[i][0]*GS[j][0]+cofs[1]*B[i][1]*GS[j][1]+cofs[2]*B[i][2]*GS[j][2])/(cofs[0]*(GS[j][0]^2)+cofs[1]*(GS[j][1]^2)+cofs[2]*(GS[j][2]^2))
                GS[i] -= mu[i,j]*GS[j]
            except ZeroDivisionError:
                #print("done")
                return 0, GS[j]
#             print("mu = ", mu)
#             print("mu[i,j] = ", mu[i,j])
            
    return GS, mu

def LLL_expand(B, delta, n, cofs):
    GS, mu = GSO_expand(B, n, cofs)
    if GS == 0:
        return mu
    BB = vector(QQ, n)
    for i in range(n):
        BB[i] = abs(cofs[0]*GS[i][0]^2+cofs[1]*GS[i][1]^2+cofs[2]*GS[i][2]^2)
    k = 1
    while k <= n-1:
        for j in range(k)[::-1]:
            if abs(mu[k,j]) > 0.50:
                q = round(mu[k,j])
                B[k] -= q*B[j]
                for l in range(j+1):
                    mu[k,l] -= q*mu[j,l]
        if BB[k] >= (delta-mu[k,k-1]^2)*BB[k-1]:
            k += 1
        else:
            v = B[k-1]; B[k-1]=B[k]; B[k] = v;
            GS, mu = GSO_expand(B, n, cofs)
            if GS == 0:
                return mu
            for i in range(n):
                BB[i] = abs(cofs[0]*GS[i][0]^2+cofs[1]*GS[i][1]^2+cofs[2]*GS[i][2]^2)
                k = max(k-1, 1)
    #print("GS : ", GS)
    return True

def ENUM_expand(B, n, cofs, R):
    GS, mu = GSO_expand(B, n, cofs)
    BB = vector(QQ, n)
    for i in range(n):
        BB[i] = abs(cofs[0]*GS[i][0]^2+cofs[1]*GS[i][1]^2+cofs[2]*GS[i][2]^2)
    sigma = Matrix(QQ, n+1, n)
    r = vector(ZZ, n+1)
    rho = vector(QQ, n+1)
    v = vector(ZZ, n)
    c = vector(QQ, n)
    w = vector(ZZ, n)
    for i in range(n+1):
        r[i] = i
    v[0] = 1
    last_nonzero = 1
    k = 1
    while(1):
        rho[k-1] = rho[k] + (v[k-1] - c[k-1])^2*BB[k-1]
        #print("k = ", k)
        if RR(rho[k-1]) <= RR(R):
            if k == 1:
                return v
            k = k-1
            r[k-1] = max(r[k-1], r[k])
            for i in range(k+1, r[k]+1)[::-1]:
                sigma[i-1, k-1] = sigma[i,k-1] + mu[i-1, k-1]*v[i-1]
            c[k-1] = -sigma[k,k-1]
            v[k-1] = round(c[k-1])
            w[k-1] = 1
        else:
            k = k + 1
            if k == n+1:
                return False
            r[k-1] = k
            if k >= last_nonzero:
                last_nonzero = k
                v[k-1] = v[k-1] + 1
            else:
                if RR(v[k-1]) > RR(c[k-1]):
                    v[k-1] = v[k-1] - w[k-1]
                else:
                    v[k-1] = v[k-1] + w[k-1]
                w[k-1] = w[k-1] + 1
    

def is_isotropic3_2(c1,c2,c3):
    #part3: done
    a1, a2, a3, MF = make_squarefree(ZZ(c1),ZZ(c2),ZZ(c3))
    #part4: done
    ta1, ta2, ta3 = a1, a2, a3
    #print("ta1, ta2, ta3 = ", ta1, ta2, ta3)
    if a1 % 2 == 0:
        ta1 /= 2
    if a2 % 2 == 0:
        ta2 /= 2
    if a3 % 2 == 0:
        ta3 /= 2
    r3 = SC(a1,a2,ta3); r1 = SC(a2,a3,ta1)
    r2 = SC(a3,a1,ta2)
    ff3 = vector(ZZ, [r3,-1,0]); ff1 = vector(ZZ, [0,r1,-1])
    ff2 = vector(ZZ, [-1,0,r2])
    #print(ff1,ff2,ff3)
    gg0 = vector(ZZ, 3)
    for i in range(3):
        gg0[i] = crt([ff3[i],ff1[i],ff2[i]],[ta3,ta1,ta2])
    #print("gg0 = ", gg0%abs(ta1*ta2*ta3))
    sol_a = find_sols_mod(gg0[0],gg0[1],gg0[2],abs(ta1*ta2*ta3))
    sol_4 = find_sols_mod4((ta1, ta2, ta3))
    for i in range(2):
        if sol_4[i] != []:
            sols = find_sol_gen(sol_a,sol_4[i], abs(ta1*ta2*ta3))
            M = Matrix([sols[0],sols[1],sols[2]])
            #print("M = ", M)
            for j in range(3):
                if ta1*(sols[j][0]^2)+ta2*(sols[j][1]^2)+ta3*(sols[j][2]^2) == 0:
                    return MF[0]*sols[j][0], MF[1]*sols[j][1], MF[2]*sols[j][2]
            LC = LLL_expand(M, 0.99, 3, (ta1, ta2, ta3))
            if LC != True:
                return MF[0]*LC[0], MF[1]*LC[1], MF[2]*LC[2]
            #print("M = ", M)
            for j in range(3):
                if ta1*(M[j][0]^2)+ta2*(M[j][1]^2)+ta3*(M[j][2]^2) == 0:
                    return MF[0]*sols[j][0], MF[1]*sols[j][1], MF[2]*sols[j][2]
            if LC:
                R = 0.99*abs(a1*(M[0][0]^2)+a2*(M[0][1]^2)+a3*(M[0][2]^2))
                while(1):
                    v = vector(ZZ, 3)
                    print()
                    #print("ENUM")
                    v = ENUM_expand(M, 3, (ta1, ta2, ta3), R)
                    
                    
                    if v != False:
                        if v[0] == 1 and v[1] == 0:
                            return False
                        print(v)
                        vec = v[0]*M[0]
                        for i in range(1,3):
                            vec += v[i]*M[i]
                        R = 0.99*abs(a1*(vec[0]^2)+a2*(vec[1]^2)+a3*(vec[2]^2))
                        if R == 0:
                            return MF[0]*vec[0], MF[1]*vec[1], MF[2]*vec[2]
                    else:

                        break
    return False

def vol_beta(tp, p):
    #print(tp.factor())
    alpha = [ell[1] for ell in tp.factor() if ell[0] == p]
    if len(alpha) != 0: 
        return alpha % 2
    else: return 0

def make_squarefree4(c1,c2,c3,c4):
    M = vector(QQ, [1,1,1,1])
    a1,a2,a3,a4 = c1,c2,c3,c4
    #while !((a1*a2*a3).is_squarefree()):
    if not (a1.is_squarefree()):
        a1_square_primes = [(p[0],p[1]//2) for p in a1.factor() if p[1] % 2 == 0]
        a1 = a1//(prod([p[0]^(2*p[1]) for p in a1_square_primes]))
        M[0] = 1 / prod([p[0]^(p[1]) for p in a1_square_primes])
    if not (a2.is_squarefree()):
        a2_square_primes = [(p[0],p[1]//2) for p in a2.factor() if p[1] >= 2]
        a2 = a2//(prod([p[0]^(2*p[1]) for p in a2_square_primes]))
        M[1] = 1 / prod([p[0]^(p[1]) for p in a2_square_primes])
    if not (a3.is_squarefree()):
        a3_square_primes = [(p[0],p[1]//2) for p in a3.factor() if p[1] >= 2]
        a3 = a3//(prod([p[0]^(2*p[1]) for p in a3_square_primes]))
        M[2] = 1 / prod([p[0]^(p[1]) for p in a3_square_primes])
    if not (a4.is_squarefree()):
        a4_square_primes = [(p[0],p[1]//2) for p in a4.factor() if p[1] >= 2]
        a4 = a4//(prod([p[0]^(2*p[1]) for p in a4_square_primes]))
        M[3] = 1 / prod([p[0]^(p[1]) for p in a4_square_primes])
    return a1,a2,a3,a4, M 

#IsIsotropic4(4変数二次形式の有理数解の有無判定)
def is_isotropic4(a1,a2,a3,a4):
    a1,a2,a3,a4, M = make_squarefree4(ZZ(a1),ZZ(a2),ZZ(a3),ZZ(a4))
    #print((2*a1*a2*a3*a4).factor())
    P_set = [-1] + [p[0] for p in (2*a1*a2*a3*a4).factor() if p[1] > 0]
    tp_set = []
    
    q = 1
    for p in P_set:
        check = 0
       #print("p = ", p)
        
        h12 = hilbert_symbol(a1, a2, p); h34 = hilbert_symbol(-a3, -a4, p);
        if p == -1:
            if hilbert_symbol(1,-a1*a2, p) == hilbert_symbol(a1, a2, p):
                tp = 1
                check = 1
            else:
                tp = -1; check = 1
        else:
            Q = Qp(p)
            if h12 == h34 or not((Q(a1*a2)/Q(a3*a4)).is_square()):
                
                for i in range(1,p+1):
                    #print(hilbert_symbol(i,-a1*a2, p), hilbert_symbol(i,-a3*a4, p))
                    if hilbert_symbol(i,-a1*a2, p) == h12 and hilbert_symbol(i,-a3*a4, p) == h34:
                        check = 1
                        tp = i
                        break
                        
                if check == 1:
                    tp_set.append((tp, p))
                    if p == -1:
                        q *= tp
                    else:
                        q *= p^(vol_beta(ZZ(tp), p))
                else: 
                    print("line: 174")
                    return False
            else:
                print("line: 177")
                return False
    #part5: done
    ell = 2
    
    while(1):
        check = 0
        t = ell*q
        for p in P_set:
            if p != -1:
                h12 = hilbert_symbol(a1, a2, p); h34 = hilbert_symbol(-a3, -a4, p);
                if hilbert_symbol(t, -a1*a2, p) != h12 or hilbert_symbol(t, -a3*a4, p) != h34:
                    check = 1
                    break
        if check == 1: ell = next_prime(ell)
        
        if check == 0:
            break;
    print("t = ", t)
    A = is_isotropic3_2(a1,a2,-t)
    if A != False:
        y1, y2, u = A
        print("A_test : ", a1*(y1^2)+ a2*(y2^2)+ (-t)*(u^2))
    else: 
        return False
    B = is_isotropic3_2(a3,a4,t)
    if B != False:
        y3, y4, v = B
        print("B_test : ", a3*(y3^2)+ a4*(y4^2)+ (t)*(v^2))
    else: 
        return False
    print("y1, y2, y3, y4 = ", y1, y2, y3, y4)
    print("u, v = ", u, v)
    #part7: done
    g = gcd(u, v)
    x1, x2 = y1/u, y2/u
    x3, x4 = y3/v, y4/v
    print("a1,a2,a3,a4 = ", a1,a2,a3,a4)
    print("F_test : ", a1*(x1^2)+ a2*(x2^2)+ a3*(x3^2)+ a4*(x4^2))
    return x1*M[0], x2*M[1], x3*M[2], x4*M[3]
    
                
        
        
def tenary_quadratic_solve(a1, a2, a3):
    #print("a1,a2,a3 = ", a1,a2,a3)
    if a1*a2*a3 == 0: return False
    a_d = a1.denominator()*a2.denominator()*a3.denominator()
    result = is_isotropic3_2(a_d*a1,a_d*a2,a_d*a3)
    if result != False:
        return (result[0]/result[2], result[1]/result[2])
    else: return False
    
