load("test_end_calc.sage")

k = 40                                  #実験に用いる基礎体の標数
ells_num = 8                            #同種写像の次数の素因数の個数
collect_nums = 3
Fp_defined = 0                          #Fp=1なら素体上, Fp=2なら拡大体
D = 0
is_elkies = True
test_time = 5

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

#=========================================================

for i in range(test_time):
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