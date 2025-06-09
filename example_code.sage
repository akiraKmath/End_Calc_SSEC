"""
An implementation program of computing the endomorphism ring of a given supersingular elliptic curve besed on finding isogeny cycle.
(C) 2025 Mitsubisi Electric, Rikkyo University, Created by Yuta Kambe, Akira Katayama, Kazuki Komine, Yusuke Aikawa, Yuki Ishihara, Masaya Yasuda, Kazuhiro Yokoyama.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
"""

load("test_end_calc.sage")

k = 30                                  #実験に用いる基礎体の標数
ells_num = 8                            #同種写像の次数の素因数の個数
collect_nums = 3                        #2,3-isogenyの数
Fp_defined = 2                          #1なら素体上, 2なら拡大体, 0ならランダム
D = 0                                   #KLPTでのパラメータ, Defaultは0(関数内で自動で決定)
is_elkies = True                        #Elkies素数を用いるか？
test_time = 5                           #実行回数

#=== Don't change the above parameter ====================

time_calc_end = 0
time_KLPT = 0
time_Deuring = 0
time_total = 0
all_num_gen = 0

timeset_calc_end = []
timeset_KLPT = []
timeset_Deuring = []
timeset_total = []
allset_num_gen = []
count = 0
#=========================================================

while count < test_time:
    result = End_basis_test(k, ells_num, collect_nums, Fp_defined, is_elkies, D = 0)
    if result != False:
        basis, time_stamp, num_gen = result
        print("Computing End: ", time_stamp[0])
        print("KLPT: ", time_stamp[1])
        print("Deuring: ", time_stamp[2])
        print("Total: ", time_stamp[3])
        time_calc_end += time_stamp[0]; timeset_calc_end.append(time_stamp[0])
        time_KLPT += time_stamp[1]; timeset_KLPT.append(time_stamp[1])
        time_Deuring += time_stamp[2]; timeset_Deuring.append(time_stamp[2])
        time_total += time_stamp[3]; timeset_total.append(time_stamp[3])
        print()
        print("number of generator: ", num_gen)
        all_num_gen += num_gen; allset_num_gen.append(num_gen)
        print("=========================================================")
        print()
        count+=1
print()
print("=== Experimental result =================================")
print("min Computing End: ", min(timeset_calc_end))
print("max Computing End: ", max(timeset_calc_end))
print("avr Computing End: ", time_calc_end/test_time)
print()
print("min KLPT: ", min(timeset_KLPT))
print("max KLPT: ", max(timeset_KLPT))
print("avr KLPT: ", time_KLPT/test_time)
print()
print("min Deuring: ", min(timeset_Deuring))
print("max Deuring: ", max(timeset_Deuring))
print("avr Deuring: ", time_Deuring/test_time)
print()
print("min Total: ", min(timeset_total))
print("max Total: ", max(timeset_total))
print("avr Total: ", time_total/test_time)
print()
print("min num_gen: ", min(allset_num_gen))
print("max num_gen: ", max(allset_num_gen))
print("avr num_gen: ", RR(all_num_gen/test_time))
print("=========================================================")