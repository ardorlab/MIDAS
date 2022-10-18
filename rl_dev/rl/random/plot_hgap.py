import csv
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

hgap=np.loadtxt("/cm/shared/databases/MC_gap_conductance/test_MARS/fort.111")
y=hgap[:,4]
x=np.arange(1,y.shape[0]+1,1)

sns.set_style('whitegrid')
plt.rcParams.update({'axes.titlesize': 18})
plt.rcParams.update({'font.size': 18})

# Read Input and Output


xlabel = r'Learning Iterations'
ylabel = r'Best Fitness'
filename = "hgap_test.png"
x=np.arange(1,y.shape[0]+1,1)
fig=plt.figure(figsize=(9, 6))
y=hgap[:,1]

plt.plot(x,y, linewidth=2)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig(filename,bbox_inches='tight',dpi=300) 
plt.close(fig)