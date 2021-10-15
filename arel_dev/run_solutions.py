import os

solution_list = ['child_70_20','child_71_13','child_71_16','child_72_32',
                 'child_73_22','child_73_37','child_74_23','child_75_22',
                 'child_76_15','child_76_19','child_76_25','child_76_34',
                 'child_76_5','child_77_1','child_77_22','child_77_23',
                 'child_77_28','child_77_3','child_78_0','child_78_10',
                 'child_78_2','child_78_24','child_78_30','child_78_4',
                 'child_78_7','child_79_1','child_79_17','child_79_20',
                 'child_79_22','child_79_25','child_79_26','child_79_29',
                 'child_79_3','child_79_30','child_79_32','child_79_34',
                 'child_79_35','child_79_36','child_79_39','child_79_4']

for name in solution_list:
    file_ = open(f"{name}/{name}_sim.inp",'r')
    file_lines = file_.readlines()
    file_.close()

    new_file = open(f"{name}/new_{name}.inp",'w')
    for line in file_lines:
        new_file.write(line)
        if "'ITE.SRC'" in line:
            new_file.write("'PRI.STA' '2EXP' '2RPF' /\n")
            new_file.write("'PIN.EDT' 'ON' '3PIN' '3PXP' '3KWF' /\n")
            new_file.write("'PIN.ASM' 'ALL-3D'/\n")
            new_file.write("'PIN.ZED' 25*1/\n")
    new_file.close()
    command = f"cd {name} ; simulate3 new_{name}.inp ;"
    command += f" python ~/veraCS/simulate_exchanger.py --file new_{name}.out ; "
    command += f"mv new_{name}.h5 deck.h5 ; cp vera_new_{name}.inp deck.inp ; rm new_{name}.out"
    os.system(command)





