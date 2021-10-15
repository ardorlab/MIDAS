import numpy as np 
import matplotlib.pyplot as plt 

with open('rewards.txt') as f:
    lines = f.readlines()
    epochs = [line.split()[0] for line in lines]
    rewards = [line.split()[1] for line in lines]

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("plot")
ax1.set_xlabel('epochs')
ax1.set_ylabel('reward')

ax1.plot(epochs, rewards, c='r', label='data')

leg = ax1.legend()

plt.show