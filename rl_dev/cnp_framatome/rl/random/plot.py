import csv
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

fitness=np.loadtxt("random.txt",skiprows=1)


sns.set_style('whitegrid')
plt.rcParams.update({'axes.titlesize': 18})
plt.rcParams.update({'font.size': 18})

# Read Input and Output

y=np.array(fitness)
x=np.arange(1,len(fitness)+1)


# Plot

xlabel = r'Learning Iterations'
ylabel = r'Fitness'
filename = "fitplot.png"

fig=plt.figure(figsize=(9, 6))
plt.plot(x,y, linewidth=2)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig(filename,bbox_inches='tight',dpi=300) 
plt.close(fig)

nep = 30
fitav = []
for i in range(len(fitness)):
    if i%nep ==0 and i!=0:
        val = np.mean(np.array(fitness[i-nep:i]))
        fitav.append(val)

y=np.array(fitav)
x=np.arange(1,len(fitav)+1)

xlabel = r'Learning Epochs'
ylabel = r'Fitness'
filename = "fitepoch.png"

fig=plt.figure(figsize=(9, 6))
plt.plot(x,y, linewidth=2)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig(filename,bbox_inches='tight',dpi=300) 
plt.close(fig)


nep = len(fitness)
bestfit = []
besti=-np.inf
for i in range(len(fitness)):
    fiti = fitness[i]
    if fiti>besti:
        bestfit.append(fiti)
        besti=fiti
    else:
        bestfit.append(besti)

y=np.array(bestfit)
x=np.arange(1,len(bestfit)+1)

xlabel = r'Learning Iterations'
ylabel = r'Best Fitness'
filename = "fitbest.png"

fig=plt.figure(figsize=(9, 6))
plt.plot(x,y, linewidth=2)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig(filename,bbox_inches='tight',dpi=300) 
plt.close(fig)