def Factor_Poly_dst(E, ell, dst):
    FF = E.base_ring();
    a = E.a4(); b = E.a6(); ell = FF(ell)

    j = E.j_invariant()
    A = -48*a; B = 864*b
    
    try: jj = -j*(B/A)
    except ZeroDivisionError: return False

    db = ClassicalModularPolynomialDatabase()
    P.<x, y> = FF[]
    Phi = db[ell].subs(j0 = x, j1 = y)
    tilde_j = dst #A root of Phi(x, j)
    Tmp = Phi.gradient()
    Phi_x = Tmp[0]; Phi_y = Tmp[1]
    tmp = FF(Phi_x(x = j, y = tilde_j)); tmp1 = FF(Phi_y(x = j, y = tilde_j))
    try: 
        tilde_jj = -jj*tmp 
        tilde_jj *= (FF(ell*tmp1).inverse_of_unit())
    except: return False
    try: tilde_a = -(tilde_jj)^2/(FF(48)*tilde_j*(tilde_j - FF(1728)))
    except: return False
    tilde_b = -(tilde_jj)^3/(FF(864)*tilde_j^2*(tilde_j - FF(1728)))
    tilde_A = -48*tilde_a; tilde_B = 864*tilde_b
    Tmp = Phi_x.gradient()
    Tmp1 = Phi_y.gradient()
    Phi_xx = Tmp[0]; Phi_xy = Tmp[1]; Phi_yy = Tmp1[1]
    tmp = (jj)^2*Phi_xx(x = j, y = tilde_j) + 2*ell*jj*tilde_jj*Phi_xy(x = j, y = tilde_j) + ell^2*(tilde_jj)^2*Phi_yy(x = j, y = tilde_j)
    try:
        J = -FF(tmp)/(jj*FF(Phi_x(x = j, y = tilde_j)))
        p1 = FF(ell*(J/FF(2))) + FF((ell/FF(4))*(A^2/B - ell*tilde_A^2/tilde_B)) + FF((ell/FF(3))*(B/A - ell*tilde_B/tilde_A))
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

    Q.<t> = FF[]
    R.<w> = Q.quotient(ideal(t^(d+1)))

    F = -FF(p1/FF(2))*w
    for k in range(1, d+1): 
        F -= FF((hat_c[k] - ell*c[k])/FF((2*k+1)*(2*k+2)))*w^(k+1)
    A = R(1)
    for k in range(1, d+1): 
        A += F^k/FF(factorial(k))
    C = R(0)
    for k in range(1, d+1): C += FF(c[k])*w^k

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

def ell_kernel_poly(E, ell, dst):
    div_ell = E.division_polynomial(ell).monic()
    fcts = div_ell.factor()
    fcts = [f[0].monic() for f in fcts]
    for f in fcts:
        try:
            psi = E.isogeny(f)
            if psi.codomain().j_invariant() == dst:
                return psi
        except ValueError:
            #print("dst : ", ell, dst)
            return False
    return False

def rats_isogeny(E, cycle):
    j_invs = cycle[0]; degs = cycle[1]
    F = E.base_field()
    Et = copy(E)
    rats = []
    for i in range(len(degs)):
        ell = degs[i]
        if ell >= 5:
            KP = Factor_Poly_dst(Et, ell, j_invs[i+1])
            if KP == False:
                return False
            else:
                try:
                    KP = KP.monic()
                    prepsi = Et.isogeny(KP)
                    if prepsi.codomain().j_invariant() == j_invs[i+1]:
                        rats.append(prepsi.rational_maps())
                        Et = prepsi.codomain()
                except ValueError:
                    print("dst : ", ell, j_invs[i+1])
                    return False
        else:
            KP = ell_kernel_poly(Et, ell, j_invs[i+1])
            if KP == False:
                if ell == 3:
                    KP = Factor_Poly_dst(Et, ell, j_invs[i+1])
                    if KP == False:
                        return False
                    else:
                        try:
                            KP = KP.monic()
                            prepsi = Et.isogeny(KP)
                            if prepsi.codomain().j_invariant() == j_invs[i+1]:
                                rats.append(prepsi.rational_maps())
                                Et = prepsi.codomain()
                        except ValueError:
                            print("dst : ", ell, j_invs[i+1])
                            return False
                else:
                    return False
            else:  
                rats.append(KP.rational_maps())
                Et = KP.codomain()
    if Et != E:
        iso = Et.isomorphism_to(E).rational_maps()
        tmpR = rats[-1][0].parent()
        rats.append([tmpR(iso[0]), tmpR(iso[1])])
    return rats