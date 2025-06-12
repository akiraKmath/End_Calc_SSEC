"""
An implementation program of computing the endomorphism ring of a given supersingular elliptic curve besed on finding isogeny cycles.
(C) 2025 Mitsubisi Electric, Rikkyo University, Created by Yuta Kambe, Akira Katayama, Kazuki Komine, Yusuke Aikawa, Yuki Ishihara, Masaya Yasuda, Kazuhiro Yokoyama.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
"""

from sage.categories.morphism import SetMorphism

def random_isogeny_power_of_2(exp,field,start):
    ### exp: = log_2(degree_of_isogeny)
    ### start: the domain of the result
    ### Return: An isogeny start -> target

    ## Basis of ell-torsion points ##
    # This program is made by Yasuda Masaya
    # It returns a list [P,Q,Ext], where [P,Q] is a basis of E[ell]
    # and Ext is the curve that P,Q in Ext.
    def Basis_of_Torsions(E,ell,z):
        FF = E.base_field(); p = FF.order().factor()[0][0]
        a = E.a4(); b = E.a6(); d = 1;

        while(1):
            order = (p^d - (-1)^d)^2
            if order % ell^2 == 0: break
            else: d += 1
        
        FFext.<w> = FF.extension(d); u = FFext(z)
        Eext = EllipticCurve(FFext, [a, b])

        m = 1
        while(1):
            if order % ell^(2*m) != 0: break
            m += 1
        n = ZZ(order/ell^(m))

        while(1):
            tmp = Eext.random_element(); P = n*tmp
            if P != Eext(0) and ell*P == Eext(0):
                break
        
        while(1):
            tmp = Eext.random_element(); Q = n*tmp
            if Q != Eext(0) and ell*Q == Eext(0):
                try: v = discrete_log(Q,P,ell,operation="+")
                except ValueError: break
        return [P,Q,Eext]

    ## Make a random isogeny from start

    phi = EllipticCurveIsogeny(start,start(0)) # identity map
    Temp_E = start

    ## Compute random 2-isogenies
    FF = start.base_field()
    char = FF.characteristic()
    z = FF.gen()
    for idx in range(exp):
        ## Pick a basis of the 2-torsion group of Temp_E
        P,Q,Temp_Eext = Basis_of_Torsions(Temp_E,2,z)

        ## Determine the kernel point
        R = 0*P+0*Q
        while R == Temp_Eext(0):
            R = randint(0,2)*P+randint(0,2)*Q
                
        psi = EllipticCurveIsogeny(Temp_Eext,R) # Isogeny Temp_Eext -> Temp_Eext/<R>
        ker_poly = psi.kernel_polynomial().change_ring(FF).list()

        #Tranrate to setmorphism for composite to phi
        Temp_phi = EllipticCurveIsogeny(Temp_E,ker_poly)
        Temp_phi2 = SetMorphism(Hom(Temp_phi.domain(), Temp_phi.codomain(), Sets()),Temp_phi)
 
        phi = Temp_phi2*phi
        Temp_E = phi.codomain()

    return phi