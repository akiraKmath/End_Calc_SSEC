load("test_end_calc.sage")

k = 25                                  #実験に用いる基礎体の標数
ells_num = 6                            #同種写像の次数の素因数の個数
collect_nums = 3
Fp_defined = 0                          #Fp=1なら素体上, Fp=2なら拡大体
D = 0
is_elkies = True
test_time = 5

for i in range(test_time):
    result = End_basis_test(k, ells_num, collect_nums, Fp_defined, is_elkies, D = 0)
    if result != False:
        basis, time_stamp, num_gen = result