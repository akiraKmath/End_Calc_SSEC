"""
An implementation program of computing the endomorphism ring of a given supersingular elliptic curve besed on finding isogeny cycles.
(C) 2025 Mitsubisi Electric, Rikkyo University, Created by Yuta Kambe, Akira Katayama, Kazuki Komine, Yusuke Aikawa, Yuki Ishihara, Masaya Yasuda, Kazuhiro Yokoyama.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
"""

################ 1.1 Extension degree ################

def Ext_deg(ell,char):

    # This program returns the extension degree of the ell-torsion group for F_p of a supersingular elliptic curve

    d = 1
    while(1):
        order = (char^d - (-1)^d)^2
        if order % ell^2 == 0: break
        else: d += 1
    return 2*d

################ 1.1 Extension degree End ################

################ 1.2 Power smooth number ################

def Sgenerator(char,lower_bound,S1,Sd1,max_ext_deg, excepts):

    # Return: An positive integer S = l_1 * l_2 * ... * l_a where {l_1,...,l_a} is distinct a primes,
    #       : s.t. extension degree of E[l_i] <= max_ext_deg, S1*S satisfies the same condition

    # Input
    # lower_bound: S > lower_bound
    # S1,Sd1: given positive integer and Sd1 = dict(factor(S1))
    # max_ext_deg: upper bound of extension degrees

    Primeset = Primes() #An instance of all prime numbers set.
    S = 1
    Exp_dict = dict()
    l = 3
    while True:
        if l not in excepts:
            break;
        else :
            l = Primeset.next(l)
    while S <= lower_bound:

        if l in Sd1.keys():
            l = Primeset.next(l)
            while True:
                if l not in excepts:
                    break;
                else :
                    l = Primeset.next(l)
            continue
        #print("here : ", Ext_deg(ell=l,char=char), max_ext_deg)
        if Ext_deg(ell=l,char=char) > max_ext_deg:
            l = Primeset.next(l)
            while True:
                if l not in excepts:
                    break;
                else :
                    l = Primeset.next(l)
            continue

        Exp_dict[l] = 1
        S = S*l
        l = Primeset.next(l)
        while True:
            if l not in excepts:
                break;
            else :
                l = Primeset.next(l)

    return [S,Exp_dict]

################ 1.2 Power smooth number End ################

########### 2.1 Prime Norm algorithm ###########

def PrimeNormAlgorithm(order,ideal,char,iteration = 5000,lower_bound=100):

    # Return an ideal J s.t. Nrd(J) is prime and J and the input ideal are equivalent

    # This program finds an element D of the ideal s.t. Nrd(D) = Nrd("ideal")*prime
    # This program generates "iteration" elements of the ideal,
    # and choices an element D s.t. the norm Nrd(D) is minimum among found elements.

    NrI = ideal.norm() #The norm of I.
    # print("NrI : ", NrI)
    m = ceil(log(char))   #The boundary of random elements.
    Zbasis = ideal.basis() #The rows are basis of I over Z

    #New========================================================================
    Zbasis  = strtoSage(Magma_Minkowski(Zbasis, order, char),char)
    #===========================================================================
    D = 0
    N = Infinity

    for iter in range(iteration):
        TempD = 0
        for idx in range(0,4):
            TempD = TempD + randint(-m,m)*Zbasis[idx] #A linear form of Zbasis
        Nrd = TempD.reduced_norm()
        # print("Nrd_check",Nrd, NrI)
        # print("sub_check: ", Nrd/NrI)
        if (Nrd / NrI) in ZZ:
            Temp = ZZ(Nrd/NrI) #A candidate of prime norm
            # print("Nrd_check",iter)
            if Temp.is_prime():
                if Temp < N and lower_bound <= Temp:
                    D = TempD
                    N = Temp

    # Here D is an element s.t. Nrd(D) = Nrd("ideal")*prime

    Gamma = (1/NrI)*D.conjugate()

    basis_of_J = []

    for idx in range(0,4):
        basis_of_J.append(Zbasis[idx]*Gamma)
    return order.left_ideal(basis_of_J)

########### 2.1 Prime Norm algorithm End ###########

########### 2.2 Make Beta_1,Beta_2 ###########

def CornacchiaAlgorithm(M):
    # Return a pair [a,b] such that a^2 + b^2 = M

    if M < 0:
        return False
    if M == 0:
        return [0,0]
    F = DiagonalQuadraticForm(QQ, [1, 1])
    v = F.solve(M)
    return [v[0], v[1]]

def Find_beta_1(B_pinf,order,char,N,S):

    # This program returns an element "beta1" in the given order s.t. Nrd("beta1")= NS
    # It assumes that the norm form is a^2+b^2+p(c^2+d^2).

    m = floor(sqrt((N*S)/(2*char)))
    Zbasis_order = order.basis()
    # Zbasis_order = [O]
    # print(Zbasis_order)

    while True:
        TempD = 0
        for idx in range(0,4):
            TempD = TempD + randint(-m,m)*Zbasis_order[idx]
        # print("TempD = ", TempD)
        # TempD = O_to_ZZ(order, TempD)
        Denomi = lcm(lcm(denominator(TempD[0]),denominator(TempD[1])),lcm(denominator(TempD[2]),denominator(TempD[3])))
        Denomi = lcm(denominator(TempD[2]),denominator(TempD[3]))
        TempD = Denomi*TempD
        CD = [TempD[2],TempD[3]] # the coefficients of TempD at j,k
        # print("CD : ", CD)
        MR = QQ(N*S-char*(CD[0]^2+CD[1]^2))
        # print("MR : ", MR)
        if MR < 0:
            # print("MR<0", MR)
            continue
        try: AB = CornacchiaAlgorithm(MR)
        except: continue

        return B_pinf([AB[0],AB[1],CD[0],CD[1]])

def Find_beta_2(order, B_pinf,char,ideal,beta1,NrI):

    # This program returns an element "beta2" s.t. O(beta_1*beta_2) = I mod NO
    # With very low propability, for beta1, it returns False

    # First, find an element "alpha" in I s.t. I = O(alpha) mod NO

    m = ceil(log(char))
    Zbasis = ideal.basis() #The rows are basis of I over Z

    #New========================================================================
    Zbasis  = strtoSage(Magma_Minkowski(Zbasis, order, char),char)
    #===========================================================================

    # The following finds a random element "alpha" of I s.t. gcd(Nrd("alpha"), NrI^2) = NrI

    while True:
        TempD = 0
        for idx in range(0,4):
            TempD = TempD + randint(-m,m)*Zbasis[idx] #A random linear form of Zbasis

        if gcd(TempD.reduced_norm(),NrI^2) == NrI:
            alpha = TempD
            break
    # print("alpha : ", alpha, alpha in ideal, alpha in order)
    # Nalpha = O_to_ZZ(ideal, alpha)

    # Here "alpha" is a principal generator of I/NO
    # Next, find "beta_2" s.t. O*(beta_1*beta_2) + NO = O*(alpha) + NO for given beta_1
    X = vector([alpha[0],alpha[1],alpha[2],alpha[3]])
    # X = vector(Nalpha)
    # Nbeta1 = O_to_ZZ(ideal, beta1)
    # print("ss : ", [(ZZ(f.denominator())).factor() for f in Nbeta1])
    Y = vector([ZZ(beta1[0]),ZZ(beta1[1]),ZZ(beta1[2]),ZZ(beta1[3])])
    # Y = vector([ZZ(beta1[0].numerator()*inverse_mod(beta1[0].denominator(),NrI)),ZZ(beta1[1].numerator()*inverse_mod(beta1[1].denominator(),NrI)),(beta1[2].numerator()*inverse_mod(beta1[2].denominator(),NrI)),ZZ(beta1[3].numerator()*inverse_mod(beta1[3].denominator(),NrI))])
    ZN = IntegerModRing(NrI)
    # Y = vector(ZN, Nbeta1)
    # print("Y = ", Y)
    # Y = vector(beta1)
    # ZN = IntegerModRing(NrI)
    M = matrix(ZN,[
                    [-char*Y[2],-char*Y[3],-X[0],X[1],char*X[2],char*X[3]],
                    [-char*Y[3],char*Y[2],-X[1],-X[0],-char*X[3],char*X[2]],
                    [Y[0],-Y[1],-X[2],X[3],-X[0],-X[1]],
                    [Y[1],Y[0],-X[3],-X[2],X[1],-X[0]]
                    ])

    # print("M = ", M)

    Rank = M.rank()
    Kernel = (M.right_kernel()).basis_matrix()
    # print("Kernel : ", Kernel)
    # print("Rank = ", Rank)

    if Rank == 2:
        i = 0
        while Kernel[i][0] == 0 and Kernel[i][1] == 0:
            i = i + 1

        return B_pinf([0,0,ZZ(Kernel[i][0]),ZZ(Kernel[i][1])])
        # return B_pinf([0,0,(Kernel[i][0]),(Kernel[i][1])])

    if Rank == 3:
        i = 0
        while Kernel[i][0] == 0 and Kernel[i][1] == 0:
            i = i + 1

        if (ZZ(Kernel[i][0])^2 + ZZ(Kernel[i][1])^2) % NrI != 0:
            return B_pinf([0,0,ZZ(Kernel[i][0]),ZZ(Kernel[i][1])])
        # if ((Kernel[i][0])^2 + (Kernel[i][1])^2).numerator() % NrI != 0:
        #     return B_pinf([0,0,(Kernel[i][0]),(Kernel[i][1])])

        return False

########### 2.2 Make Beta_1,Beta_2 End ###########

########### 2.3 S2_equation ###########

def S2_eq(beta2,Sd1,S2,Sd2,N,char,max_ext_deg):

    # Input
    # beta2: an element of the form Cj+Dk such that O(Beta1*Beta2) = I_prime mod NO
    # Sd1,Sd2: the dictionary of the exponentials of the prime divisors of S1,S2

    # Output
    # (Beta2)' such that Nrd((Beta2)') = S2 and (Beta2)' = Beta2 mod N, and changed S2.

    C = beta2[2]
    D = beta2[3]
    if C != 1:
        D = ZZ(Mod((Mod(C,N)^-1)*D,N))
        C = 1

    SquareD = 1+D^2
    ModNSquareD = SquareD % N
    ConstDN = D*N

    NN = N^2

    Inv = Mod(char*ModNSquareD,N)^-1

    #The core algorithm. Until finding valid solution, this program updates S2 by multiplying a suitable prime
    Primeset = Primes() #An instance of all prime numbers set.
    TempS2 = S2
    TempSd2 = copy(Sd2)

    while true:

        try:
            Lam = ZZ(sqrt(Mod(TempS2*Inv,N)))
        except:
            # There is no the solution for lambda. Update S2
            multi_prime = list(TempSd2.keys())[-1]
            while True:
                iter_num = 0
                multi_prime = Primeset.next(multi_prime)
                if Ext_deg(multi_prime,char) <= max_ext_deg:
                    TempS2 = TempS2*multi_prime
                    TempSd2[multi_prime] = 1
                    break
                if iter_num > 1000:
                    return False
                iter_num = iter_num + 1
            continue

        # Determine the minimum value of the quad form p*{(Lam*C + c*N)^2 + (Lam*D + d*N)^2}
        e = ZZ((TempS2 - Lam^2*char*SquareD)/N)
        v = ZZ(Mod(2*char*Lam,N)^-1)

        y0 = -(Lam*SquareD+v*e*N)/NN
        if y0 - floor(y0) <= 1/2:
            y = floor(y0)
            epsi = y0 - y
            epsiepsi = epsi^2
        else:
            y = ceil(y0)
            epsi = y0 - y
            epsiepsi = epsi^2

        if TempS2 - char*NN*(SquareD+NN*epsiepsi/SquareD) < 0: #In this case, there is no solution for TempS2.
            multi_prime = list(TempSd2.keys())[-1]
            while True:
                iter_num = 0
                multi_prime = Primeset.next(multi_prime)
                if Ext_deg(multi_prime,char) <= max_ext_deg:
                    TempS2 = TempS2*multi_prime
                    TempSd2[multi_prime] = 1
                    break
                if iter_num > 1000:
                    return False
                iter_num = iter_num + 1
            continue

        #This x is the best bound for d.
        x = sqrt((TempS2/(NN*char) - NN*epsiepsi/SquareD)/SquareD)


        #This d determine the minimum value of (Lam*C + c*N)^2 + (Lam*D + d*N)^2
        d = -D*(Lam - NN*epsi/SquareD)/N
        if d - floor(d) <= 1/2:
            dm = floor(d)
        else:
            dm = ceil(d)

        MinimumVector = [Lam + (-D*dm+v*e+y*N)*N, Lam*D + dm*N]

        c0 = MinimumVector[0]
        d0 = MinimumVector[1]

        #Determin the area of d and c such that TempS2 > p*((Lam*C+c*N)^2+(Lam*D + d*N)^2)
        RatioofS2andP = TempS2/char
        for t in range(ceil((-sqrt(RatioofS2andP)-d0)/N),floor((sqrt(RatioofS2andP)-d0)/N)+1):
            TempoConst1 = (d0+t*N)^2
            TempoConst2 = c0-ConstDN*t
            for u in range(ceil((-sqrt(RatioofS2andP-TempoConst1)-TempoConst2)/NN),floor((sqrt(RatioofS2andP-TempoConst1)-TempoConst2)/NN)+1):
                #Find the solution
                TempM = TempS2 - char*(TempoConst1+(TempoConst2+u*NN)^2)
                if TempM < 0:
                    continue
                if TempM == 0:
                    AB = CornacchiaAlgorithm(ZZ(TempM/NN))
                    d = ZZ((d0+t*N-Lam*D)/N)
                    c = ZZ((TempoConst2+u*NN-Lam)/N)
                    return [[AB[0]*N,AB[1]*N,Lam+c*N,Lam*D+d*N],TempS2,TempSd2]

                try: AB = CornacchiaAlgorithm(ZZ(TempM/NN))
                except: continue
                AB = CornacchiaAlgorithm(ZZ(TempM/NN))
                d = ZZ((d0+t*N-Lam*D)/N)
                c = ZZ((TempoConst2+u*NN-Lam)/N)
                return [[AB[0]*N,AB[1]*N,Lam+c*N,Lam*D+d*N],TempS2,TempSd2]

        # There is no the solution for TempS2, Update TempS2
        multi_prime = list(TempSd2.keys())[-1]
        while True:
            multi_prime = Primeset.next(multi_prime)
            if Ext_deg(multi_prime,char) <= max_ext_deg:
                TempS2 = TempS2*multi_prime
                TempSd2[multi_prime] = 1
                break
            if iter_num > 1000:
                return False
            iter_num = iter_num + 1

########### 2.3 S2_equation End ###########

########### 2.4 KLPT algorithm ###########

def KLPTalgorithm_ver3(order,ideal,char,max_ext_deg,N_iteration = "default",N_lower=100):
    # print("test")
    # Input
    # order: a maximal order corresponding to the start curve
    # ideal: an left ideal of the order
    # max_ext_deg: upper bound of extension degree of torsion groups
    # iteration: the number of iteration of prime norm algorithm to find a prime norm as small as possible

    # Return: an ideal J equivalent to the input ideal

    B_pinf = ideal.basis()[0].parent()
    # print("B_pinf : ", B_pinf)

    if N_iteration == "default":
        N_iteration = 4*(ceil(log(char))^4)

    print("start KLPT algorithm")
    total_start_time = time.time()

    print("")
    # print("start computing minimum prime norm")
    start = time.time()
    I_prime = PrimeNormAlgorithm(order=order,ideal=ideal,char=char,iteration = N_iteration,lower_bound = N_lower)
    # print("here: line349")
    # susu = strtoSage(Magma_Minkowski(I_prime.basis_matrix(),char),char)
    # NI_prime = order.left_ideal(susu)
    # print("I_prime:",I_prime)
    # print("NI_prime : ", NI_prime)
    # print("NI_prime_basis : ", susu)
    N = I_prime.norm()
    # print("prime norm N:",N)
    # print("prime norm N:",NI_prime.norm())

    elapsed_time = time.time() - start
    # print("prime norm:{0}".format(elapsed_time) + "[sec]")
    # print("finish computing prime norm")
    exss = lcm([lcm([denominator(DD) for DD in list(Ds)]) for Ds in I_prime.basis()])
    exslis = [x[0] for x in exss.factor()]
    # print("exss : ", exss, exslis)

    print("")
    print("start computing powersmooth numbers")
    start = time.time()

    S1,Sd1 = Sgenerator(char=char,lower_bound=char*log(char),S1=1,Sd1=dict(),max_ext_deg=max_ext_deg, excepts=exslis)
    # print("S1,Sd1 = ", S1, Sd1)
    S2,Sd2 = Sgenerator(char=char,lower_bound=max(char^3*log(char),char*sqrt(N^3/pi)),S1=S1,Sd1=Sd1,max_ext_deg=max_ext_deg, excepts=exslis)
    # print("S2,Sd2 = ", S2, Sd2)

    # print("S1:",Sd1)
    # print("S2:",Sd2)

    Primeset = Primes()

    elapsed_time = time.time() - start
    print("powersmooth numbers:{0}".format(elapsed_time) + "[sec]")
    print("finish computing powersmooth numbers")
    print("")


    while True:
        # print("start computing beta1,beta2,beta'2")
        start = time.time()

        # Find beta_1,beta_2
        # beta1 is an element of the given order s.t. Nrd(beta1) = N*S1
        beta1 = Find_beta_1(B_pinf=B_pinf,order = order,char=char,N=N,S=S1)
        # print("beta1 : ", (beta1[0]^2+beta1[1]^2+char*(beta1[2]^2+beta1[3]^2)).factor(), beta1)
        # print(beta1 in B_pinf)
        # print("beta1 done!")
        # beta2 = Cj+Dk such that O*(beta1*beta2) = I_prime mod O/O*N as ideals
        beta2 = Find_beta_2(order = order, B_pinf=B_pinf,char=char,ideal=I_prime,beta1=beta1,NrI=N)
        # print("beta2 done!", beta2)
        # print("bata1*beta2 : ", beta1*beta2)
        # print(beta1*beta2 in I_prime)
        if beta2 != False:
            if beta2[2] != 0:

                Result_S2 = S2_eq(beta2=beta2,Sd1=Sd1,S2=S2,Sd2=Sd2,N=N,char=char,max_ext_deg=max_ext_deg)

                if Result_S2 != False and beta1*beta2 in I_prime:
                    beta2_prime,NewS2,NewSd2 = Result_S2
                    beta2_prime = B_pinf(beta2_prime)
                    # print("beta2_prime : ", beta2_prime, ZZ(beta2_prime*(beta2_prime.conjugate())).factor())
                    # print("bata1*beta2_prime : ", beta1*beta2_prime)
                    # print(beta1*beta2_prime in I_prime)

                    # print("NewS2:",NewSd2)

                    elapsed_time = time.time() - start
                    # print("beta1,beta2,beta'2:{0}".format(elapsed_time) + "[sec]")
                    # print("finish computing beta1,beta2,beta'2")

                    # Compute the ideal of norm S1*S2
                    delta = (1/N)*(beta1*beta2_prime).conjugate()
                    Zbasis_I_prime = I_prime.basis()

                    #New========================================================
                    Zbasis_I_prime  = strtoSage(Magma_Minkowski(Zbasis_I_prime, order, char),char)
                    #===========================================================

                    Result_ideal_basis = []
                    for idx in range(0,4):
                        Result_ideal_basis.append(Zbasis_I_prime[idx]*delta)
                    NRIB = strtoSage(Magma_Minkowski(Result_ideal_basis,order, char),char)
                    # print("here : 422")
                    # Converge Sd1 and NewSd2
                    norm_divisios = dict()
                    for l in Sd1.keys():
                        norm_divisios[l] = 1
                    for l in NewSd2.keys():
                        norm_divisios[l] = 1

                    if beta1*beta2_prime in I_prime:
                        total_finish_time = time.time() - total_start_time
                        print("")

                        ("total:{0}".format(total_finish_time) + "[sec]")
                        print("finish KLPT")

                        return [order.left_ideal(NRIB),norm_divisios]
                else:
                    print("There is no solution, retry for beta1,beta2")
                    continue

########### 2.4 KLPT algorithm End ###########

def Magma_Minkowski(W,O,B):
    # Return minkowski reduced basis by Magma.
    No = O.basis()
    return magma_free(f'p:={B};B<i,j,k>:=QuaternionAlgebra<RationalField()|-1,-p>;O:=QuaternionOrder([{No[0][0]}+{No[0][1]}*i+{No[0][2]}*j+{No[0][3]}*k,{No[1][0]}+{No[1][1]}*i+{No[1][2]}*j+{No[1][3]}*k,{No[2][0]}+{No[2][1]}*i+{No[2][2]}*j+{No[2][3]}*k,{No[3][0]}+{No[3][1]}*i+{No[3][2]}*j+{No[3][3]}*k]);I:=LeftIdeal(O,[{W[0][0]}+{W[0][1]}*i+{W[0][2]}*j+{W[0][3]}*k,{W[1][0]}+{W[1][1]}*i+{W[1][2]}*j+{W[1][3]}*k,{W[2][0]}+{W[2][1]}*i+{W[2][2]}*j+{W[2][3]}*k,{W[3][0]}+{W[3][1]}*i+{W[3][2]}*j+{W[3][3]}*k]);ReducedBasis(I)')

def strtoSage(state, p):
    B_quot.<i,j,k> = QuaternionAlgebra(QQ,-1,-p)
    preall = []
    presub = ""
    for sub in state:
        if sub != "[" and sub != "]":
            if sub != "," :
                # print("sub : ", sub)
                if sub != ' ':
                    presub += sub
                    # print(presub)
            else:
                preall.append(presub)
                presub = ""
    preall.append(presub)
    # print("preall: ", preall)
    alls = []
    # for pusu in preall:
    for tstdat in preall:
        kom = [0,0,0,0]
        nuuth = tstdat.split('-')
        # print("nuuth : ", nuuth)
        strm = []
        if nuuth[0] == '':
            nuuth.pop(0)
            for thh in nuuth:
                thh = '-'+thh
                strm.append(thh)
        else:
            strm.append(nuuth[0])
            for tu in range(1,len(nuuth)):
                rhu = nuuth[tu]
                rhu = '-' + rhu
                strm.append(rhu)
        knss = []
        for mook in strm:
            knss += mook.split('+')
                # shu.append(knss)
        # print("knss : ", knss)
        for pod in knss:
            t = 0
            coff = ""
            # print("kokossu", len(pod)-1)
            # print("pod : ", pod)
            if len(pod) != 1:
                while t < len(pod):
                    if pod[t] not in ["*", "i", "j", "k"]:
                        # print(pod[t])
                        coff += pod[t]
                    # print("t = ", t)
                    t += 1
                    # print("kof", coff)
            # print("here0")
            if pod[-1] == "i":
                if coff != "":
                    if coff == "-":
                        kom[1] == QQ(-1)
                    else:
                        kom[1] = QQ(coff);
                else:
                    kom[1] = QQ(1);
            elif pod[-1] == "j":
                if coff != "":
                    if coff == "-":
                        kom[2] == QQ(-1)
                    else:
                        kom[2] = QQ(coff);
                    # kom[2] = QQ(coff);
                else:
                    kom[2] = QQ(1);
            elif pod[-1] == "k":
                if coff != "":
                    if coff == "-":
                        kom[3] == QQ(-1)
                    else:
                        kom[3] = QQ(coff);
                    # kom[3] = QQ(coff);
                else:
                    kom[3] = QQ(1);
            else:
                # print("here")
                if coff != "":
                    if coff == "-":
                        kom[0] == QQ(-1)
                    else:
                        kom[0] = QQ(coff);
                    # kom[0] = QQ(coff)
                else:
                    kom[0] = QQ(1)
        # print("kom: ", kom)
        alls.append(B_quot(kom))
    return alls;

def O_to_ZZ(order, x):
    Ob0,Ob1,Ob2,Ob3 = order.basis()
    xlis = list(x)
    Ob0=list(Ob0); Ob1=list(Ob1); Ob2=list(Ob2); Ob3=list(Ob3);
    a,b,c,d = var('a b c d')
    f0 = a*Ob0[0]+b*Ob1[0]+c*Ob2[0]+d*Ob3[0]-xlis[0]
    f1 = a*Ob0[1]+b*Ob1[1]+c*Ob2[1]+d*Ob3[1]-xlis[1]
    f2 = a*Ob0[2]+b*Ob1[2]+c*Ob2[2]+d*Ob3[2]-xlis[2]
    f3 = a*Ob0[3]+b*Ob1[3]+c*Ob2[3]+d*Ob3[3]-xlis[3]
    solns = solve([f0,f1,f2,f3],a,b,c,d, solution_dict=True)
    # print("solution : ", solns)
    return [solns[0][a], solns[0][b], solns[0][c], solns[0][d]]
