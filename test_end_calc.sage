from sage.functions.log import logb
from sage.crypto.util import random_blum_prime
import time


load("construct_isogeny.sage")
load("Elkies.sage")
load("End_basis.sage")
load("Deuring.sage")
load("diophantus_solve.sage")
load("find_isogeny_cycle.sage")
load("KLPT.sage")
load("random_curves.sage")
load("trace_cal.sage")



def End_basis_test(k, ells_num, collect_nums, Fp_defined, is_elkies, D = 0):
    bound = 9*k                            #Schoof Algorithm 用の制限(基準)
    time_stamp = []
    print("=== Fix a prime of bit k ================================")
    p = random_blum_prime(2^(k-1),2^k)
    print("p = ", p)                                     ###k-bitのランダムな素数###
    F.<z> = GF(p^2,name="z",modulus=x^2+1)               ###有限体Fの定義###
    ell_set = Primes()[:ells_num]
    print(ell_set)
    E_1728 = EllipticCurve_from_j(F(1728))
    #exp1 = floor(logb(p/36,2)) 
    exp1 = floor(logb(p/24,2))  
    print("=========================================================")
    print()
    #=======================================================
    if D == 0:
        D = RR((12.0*log(p))/(log(log(p))))
        #D = RR((8.0*log(p))/(log(log(p))))
    #======== Fp上の超特異楕円曲線の決定 =======================
    print("=== Fix a supersingular elliptic curve ==================")
    while(1) :
        phis = random_isogeny_power_of_2(exp=exp1,field=F,start=E_1728)
        ###Eの定義の表示###

        E = phis[randint(1,exp1)].codomain()  
        j = E.j_invariant()
        
        if j != j^p:
            if Fp_defined == 0 or Fp_defined == 2:
                Fp = 2
                break
            else: continue
        elif j == 0 or j == 1728 or j == 287496:
            continue
        else:
            if Fp_defined == 0 or Fp_defined == 1:
                Fp = 1
                break
            else: continue
    print("j-invariant of E = ", j)
    E_start = EllipticCurve_from_j(j)
    E = E_start
    print("E: ", E)

    print("|E| = ", E_start.cardinality().factor())
    print("=========================================================")
    print()
    print("=== Computing a basis of End(E) (Step 1)=================")
    
    start_time1 = time.time()
    if Fp == 1:
        basis, num_gen = End_basis_prime(E, p, ell_set, bound, collect_nums, is_elkies)
    else:
        basis, num_gen = End_basis_ext(E, p, ell_set, bound, collect_nums, is_elkies)
    end_time1 = time.time()-start_time1
    time_stamp.append(end_time1)

    print("=========================================================")
    print()
    print("=== Verifying via KLPT + Deuring algorithm ==============")
    print("=== KLPT algorithm (Step 2)==============================")
    start_time2 = time.time()
    result_target = basis

    B_quot = result_target[0].parent()
    Gene_of_E6 = [B_quot([1,0,0,0]),B_quot([0,1,0,0]),B_quot([1/2,0,0,1/2]),B_quot([0,1/2,1/2,0])]
    Order_0 = B_quot.quaternion_order(Gene_of_E6)
    Prod = []
    for idx in range(0,4):
        for jdx in range(0,4):
            Prod.append(Gene_of_E6[idx]*result_target[jdx])
    Prod = Order_0.left_ideal(Prod).basis()
    Order_1 = B_quot.quaternion_order(result_target)
    Intersection = Order_0.intersection(Order_1)
    U =Order_0.left_ideal(Intersection.basis()).basis_matrix()*matrix(QQ,[[1,0,0,0],[0,1,0,0],[0,-1,0,2],[-1,0,2,0]])
    NNN = (1/2)*U.det().abs()
    NNNProd = []
    for idx in range(4):
        NNNProd.append(NNN*Prod[idx])
    connecting_ide = Order_0.left_ideal(NNNProd)

    Result_KLPT = KLPTalgorithm_ver3(order=Order_0,ideal=connecting_ide,char=p,max_ext_deg=D)
    end_time2 = time.time()-start_time2
    time_stamp.append(end_time2)
    print("=========================================================")

    print("=== Deuring algorithm (Step 3)===========================")
    start_time3 = time.time()
    Result_Deuring = Constructive_Deuring_ver3(Alphas=Result_KLPT[0].basis_matrix(),Sd=Result_KLPT[1],p=p,field = F, E = EllipticCurve(F, [1, 0]))

    print("answer:", E.j_invariant())
    print("END+KLPT:", Result_Deuring[0].j_invariant(), Result_Deuring[0].j_invariant().conjugate())
    end_time3 = time.time()-start_time3
    time_stamp.append(end_time3)
    print()
    print("=========================================================")
    print()
    print("=== Result ==============================================")

    if E.j_invariant() == Result_Deuring[0].j_invariant().conjugate() or E.j_invariant() == Result_Deuring[0].j_invariant():
        total = time_stamp[0]+time_stamp[1]+time_stamp[2]
        time_stamp.append(total)
        print("Verified")
        print("Result: ", basis)
        print("num_gen : ", num_gen)
        print("=========================================================")
        return basis, time_stamp, num_gen
    else:
        print("=========================================================")
        print("Can't be verified")
        print("Result: False")
        print("=========================================================")
        return False


