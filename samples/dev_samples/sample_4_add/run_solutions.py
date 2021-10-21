import os
import yaml
with open('optimized_solutions.yaml') as yaml_file:
    solution_list = yaml.safe_load(yaml_file)

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

os.system('rm -r initial_*')

for i in range(1000):
  for j in range(1000):
    name = f"child_{i}_{j}"
    if os.path.isdir(name):
      if name in solution_list:
        pass
      else:
        os.system(f'rm -r {name}')



