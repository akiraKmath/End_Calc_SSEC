#======================================
# Action of a point by an endomorphism
#======================================
def action(P, L, g, u, p):
    EE = P.curve(); infinity = EE(0)
    if P == infinity:
        return infinity
    T1 = EE(P[0]^p, P[1]^p)
    #print("T1 = ", T1)
    T2 = EE(-P[0], u*P[1])
    #print("T2 = ", T2)
    T3 = EE(-T1[0], u*T1[1])
    #print("T3 = ", T3)
    A = (g[0]%L)*P + (g[1]%L)*T2 + (g[2]%L)*T1 + (g[3]%L)*T3
    return A

#===========================================
# Finding a kernel point by random sampling
#===========================================
def KernelPoint(E, ell, W, z):
	FF = E.base_field(); p = ZZ(sqrt(FF.order()))
	a = E.a4(); b = E.a6();

	d = Mod(-p, ell).multiplicative_order()
	order = (p^d - (-1)^d)^2; #print("order : ", order)

	# print("degree =", d)
	FFext.<w> = FF.extension(d); u = FFext(z); #print("u^2 = ", u^2, FFext(-1))
	Eext = EllipticCurve(FFext, [a, b]); #print("card : ", Eext.cardinality())

	m = 1
	while(1):
		if order % ell^(2*m) != 0: break
		m += 1
	n = ZZ(order/ell^(m))

	while(1):
		flag = 1
		tmp = Eext.random_element(); P = n*tmp
		if P != Eext(0) and ell*P == Eext(0):
			g = W[0]; gg = vector([-g[0], g[1], g[2], g[3]])
			R = action(P, ell, gg, u, p)
			if R == Eext(0):
				# print("point1")
				Q = Eext(P[0]^p, P[1]^p)
				R = action(Q, ell, gg, u, p)
			if action(R, ell, g, u, p) == Eext(0) and R != Eext(0):
				flag = 0
				for k in range(1, 4):
					if action(R, ell, W[k], u, p) != Eext(0):
						flag = 1; break
		if flag == 0: break
		# print("point: m = ", m)
	return R

#====================================
# Computation of a kernel polynomial
#====================================
def KernelPolynomial(R, ell, E):
	FFext = (R.curve()).base_field(); L = ZZ(ell)
	EEext = EllipticCurve(FFext, [E.a4(), E.a6()])
	t = EEext.isogeny(R).kernel_polynomial(); d = ZZ((ell-1)/2)

	# RR.<X> = FFext[]
	# d = ZZ((L-1)/2); t = FFext(1); TT = R
	# for i in range(1, d+1):
		# t *= (X-TT[0]); TT += R

	FF = E.base_field(); z = FF.gen(); u = FFext(z)
	LL = t.list()
	RR.<T> = FF[]; Poly = T^d
	for i in range(d):
		tmp = LL[i].polynomial().quo_rem(u.polynomial())
		Poly += (FF(tmp[0])*z + FF(tmp[1]))*T^i
	return Poly

#====================================
# Main
#====================================

def Constructive_Deuring_ver3(Alphas,Sd,p,field,E):
    W = Alphas
    z = field.gen()
    L = list(Sd.keys()); r = len(L); R = list(range(r))
    for i in range(r):
        ell = L[i]
        start_time = time.perf_counter()
        R[i] = KernelPoint(E, ell, W, z)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        # print("ell = ", ell, ": KernelPoint = ", elapsed_time)

    print(); EE = E;
    isogenies = [EllipticCurveIsogeny(EE, EE(0), degree=1)]

    for i in range(r):
        ell = L[i]
        start_time = time.perf_counter()
        Poly = KernelPolynomial(R[i], ell, EE)
        phi = EllipticCurveIsogeny(EE, Poly, degree=ell)
        isogenies.append(phi)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        EE = phi.codomain()
        # print("l = ", ell, ": j-inv =", EE.j_invariant())
        print("ell = ", ell, ": Isogeny = ", elapsed_time)

        Phi = phi.rational_maps()
        for j in range(i+1, r):
            if R[j] != False:
                FFext = R[j].curve().base_field()
                EEext = EllipticCurve(FFext, [EE.a4(), EE.a6()])
                T = R[j]
                R[j] = EEext((Phi[0])(x=T[0], y=T[1]), (Phi[1])(x=T[0], y=T[1]))

    return [EE,isogenies]
