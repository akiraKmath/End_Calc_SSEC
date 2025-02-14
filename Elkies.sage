def sum(P, Q, A, f):
    global div
    if not P:
        return Q;
    if not Q:
        return P;
    x1 = P[0]; y1 = P[1]
    x2 = Q[0]; y2 = Q[1]
    if x1 == x2:
        if y1 == y2:
            return twi(P,A,f)
        else:
            return ();
    try:
        m = (y2-y1)/(x2-x1);
    except ZeroDivisionError:
        div = x2-x1
        raise
    x3 = f * m^2 - x1 - x2;
    y3 = m * (x1 - x3) - y1;
    return (x3, y3);

def twi(P, A, f):
    global div
    if not P: return P
    x1 = P[0]; y1 = P[1]
    try:
        m = (3 * x1^2 + A) / (2 * y1 * f)
    except ZeroDivisionError:
        div = 2 * y1 * f
        raise
    x3 = f*m^2 - 2*x1
    y3 = m*(x1-x3) - y1
    return (x3,y3)

def scl(P, n, A, f):
    m = abs(n);
    x1 = P[0]; y1 = P[1]
    if n == 0:
        return ()
    if m != 1:
        i = 2;
        while i <= m:
            P = sum(P, (x1, y1), A, f);
            i += 1;
    if n >= 0:
        return P;
    else:
        P = (P[0], -P[1])
        return P;


def Elkies_test(E, ell):   #入力ellがElkies素数かどうかの判定
    FF = E.base_field()
    Ea = E.a4(); Eb = E.a6()
    j = E.j_invariant()

    db = ClassicalModularPolynomialDatabase()
    P.<x, y> = FF[]
    Phi = db[ell].subs(j0 = x, j1 = y)

    R.<z> = FF[]
    tmp = Phi.subs(x = z, y = E.j_invariant())
    tmp2 = tmp.roots()
    if len(tmp2) != 0:
        return True
    else: return False

def Factor_Poly(E, ell):
    FF = E.base_ring(); p = FF.cardinality()
    a = E.a4(); b = E.a6()

    j = E.j_invariant()
    A = -48*a; B = 864*b
    
    try: jj = -j*B/A
    except ZeroDivisionError: return False

    db = ClassicalModularPolynomialDatabase()
    P.<x, y> = FF[]
    Phi = db[ell].subs(j0 = x, j1 = y)

    R.<kk> = FF[]
    tmp = Phi.subs(x = kk, y = FF(j))
    tmp2 = tmp.roots()
    if len(tmp2) == 0:
        return False
    
    tilde_j = tmp2[0][0] #A root of Phi(x, j)

    Tmp = Phi.gradient()
    Phi_x = Tmp[0]; Phi_y = Tmp[1]
    tmp = FF(Phi_x(x = j, y = tilde_j)); tmp1 = FF(Phi_y(x = j, y = tilde_j))
    try: 
        tilde_jj = -jj*tmp 
        tilde_jj *= FF(ell*tmp1).inverse_of_unit()
    except: return False
    try: tilde_a = -(tilde_jj)^2/(FF(48)*tilde_j*(tilde_j - FF(1728)))
    except: return False
    tilde_b = -(tilde_jj)^3/(FF(864)*tilde_j^2*(tilde_j - FF(1728)))
    tilde_A = -48*tilde_a; tilde_B = 864*tilde_b

    Tmp = Phi_x.gradient()
    Tmp1 = Phi_y.gradient()
    Phi_xx = Tmp[0]; Phi_xy = Tmp[1]; Phi_yy = Tmp1[1]
    tmp = FF((jj)^2*Phi_xx(x = j, y = tilde_j) + 2*ell*jj*tilde_jj*Phi_xy(x = j, y = tilde_j) + ell^2*(tilde_jj)^2*Phi_yy(x = j, y = tilde_j))
    try:
        J = -tmp/(jj*FF(Phi_x(x = j, y = tilde_j)))
        p1 = ell*J/FF(2) + ell/FF(4)*(A^2/B - ell*tilde_A^2/tilde_B) + ell/FF(3)*(B/A - ell*tilde_B/tilde_A)
    except ZeroDivisionError: return False
    d = ZZ((ell-1)/2)
    c = vector(FF, d+2); hat_c = vector(FF, d+2)
    c[1] = -a/FF(5); c[2] = -b/FF(7)
    hat_c[1] = -ell^4*tilde_a/FF(5); hat_c[2] = -ell^6*tilde_b/FF(7)
    for k in range(3, d+1):
        tmp = FF(3)/FF((k-2)*(2*k+3)); tmp1 = tmp2 = FF(0)
        for h in range(1, k-1):
            tmp1 += c[h]*c[k-1-h]
            tmp2 += hat_c[h]*hat_c[k-1-h]
        c[k] = tmp*tmp1; hat_c[k] = tmp*tmp2

    Q.<t> = PolynomialRing(FF)
    R.<w> = Q.quotient(ideal(t^(d+1)))
    F = -R(p1/FF(2))*w
    for k in range(1, d+1): F -= R((hat_c[k] - ell*c[k])/FF((2*k+1)*(2*k+2)))*w^(k+1)
    A = R(1)
    for k in range(1, d+1): A += F^k/R(factorial(k))
    C = R(0)
    for k in range(1, d+1): C += R(c[k])*w^k

    vec_A = lift(A).coefficients()
    list_C = list(range(d+1)); list_C[0] = Q(0); tmp = C
    for i in range(1, d+1):
        list_C[i] = lift(tmp); tmp *= C

    vec_F = vector(FF, d+1) #Coefficient vector of a factor of the ell-th division polynomial
    vec_F[d] = FF(1)
    for i in range(1, d+1):
        tmp = FF(0)
        for k in range(1, i+1):
            tmp2 = FF(0)
            for h in range(k+1):
                tmp1 = list_C[k-h]
                tmp2 += FF(binomial(d-i+k, k-h))*tmp1.monomial_coefficient(t^h)
            tmp += tmp2*vec_F[d-i+k]
        vec_F[d-i] = lift(A).monomial_coefficient(t^i) - tmp

    tmp = Q(vec_F[0])
    for i in range(1, d+1): tmp += vec_F[i]*t^i
    return tmp #a factor of the ell-th division polynomial

def BSGS_ell(r,S,T, a, ff):         ###④計算時間を測る際はこの中の処理藻お願いします。
        bsm = round(sqrt(r))+1; 
        L = list(range(bsm));
        tmp=()
        for i in range(bsm):
            L[i] = tmp; tmp = sum(tmp, S, a, ff);
        R = tmp;
        flag = 0; tmp = T;
        for j in range(bsm):
            if flag != 0:
                break;
            #print(tmp, L[i])
            if tmp in L:
                for i in range(len(L)):
                    
                    if tmp == L[i]:
                        I = i; J = j; flag = 1;
                        break;
            tmp = sum(tmp, (R[0], -R[1]), a, ff)
        if flag == 0: return False
        d = (I+J*bsm)%r;
        return d;

def rat_isogeny_ell(E, ell, mod_phi, FC):
    fcts = FC
    FF = E.base_field(); R.<s> = FF[]; RQ.<x> = R.quotient(ideal(fcts(t=s)))
    R0 = mod_phi[0][0].numerator().parent()
    for i in range(len(mod_phi)):
        psi = mod_phi[i]
        #print("line 157: psi = ", psi[1].parent())
        if i == 0:
            pphis = [RQ(f(s = x)) for f in psi]
            nphi = (pphis[0]*pphis[1].inverse_of_unit(), pphis[2]*pphis[3].inverse_of_unit())
        else:
            prepsi = [RQ(f(s = nphi[0])) for f in psi]
            prepsi[2] *= nphi[1] 
            nphi = (prepsi[0]*(prepsi[1].inverse_of_unit()), prepsi[2]*(prepsi[3].inverse_of_unit()))
    for i in range(len(mod_phi)):
        psi = mod_phi[i]
        if i == 0:
            pphis = [RQ(f(s = nphi[0])) for f in psi]
            dphi = (pphis[0]*pphis[1].inverse_of_unit(), nphi[1]*pphis[2]*(pphis[3].inverse_of_unit()))
        else:
            prepsi = [RQ(f(s = dphi[0])) for f in psi]
            prepsi[2] *= dphi[1] 
            dphi = (prepsi[0]*(prepsi[1].inverse_of_unit()), prepsi[2]*(prepsi[3].inverse_of_unit()))
    
    return ((nphi[0].lift(), nphi[1].lift()), (dphi[0].lift(), dphi[1].lift()))#ここ直す#

def trace_mod_elkies(E, deg, phi, ell):      ### ③ell > 50以上の実験はこの関数を用いてください
    a = E.a4(); b = E.a6()
    F = E.base_field()
    qq = F.order(); pp = qq.factor()[0][0]
    exF = F.extension(ZZ((ell-1)/2),'tt')
    exE = E.change_ring(exF)
    FC = Factor_Poly(E, ell)
    
    stime = time.time()
    nphi, dphi = rat_isogeny_ell(E, ell, phi, FC)
    etime = time.time()-stime
    MR = nphi[0].parent()
    x = MR.gen()

    h = FC(t=x)
    deg_ell = deg % ell
    if deg_ell >= ell/2:
        deg_ell -= ell
    while true:
        try:
            MMR.<x> = MR.quotient(ideal(h))
            ff = x^3+a*x+b
            MP2 = scl((x, 1), (deg_ell), a, ff)
            S = sum(dphi, MP2, a, ff)
            phi_im = (nphi[0](s = x), nphi[1](s = x))
            
            tst = BSGS_ell(ell, phi_im, S, a, ff)
            if tst == False:
                print("line 221")
                return False
            if tst > ell/2:
                tst -= ell
            return tst
                    
            assert false
        except ZeroDivisionError:
            h = gcd(h, div.lift())