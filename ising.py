from lua import *
import matplotlib.pyplot as plt
from IPython import display
import time
from matplotlib.animation import FuncAnimation, PillowWriter



Nthermalization = 100
n_array = np.arange(0,Nthermalization)
beta = 2.0

frames = Nthermalization
dim = L
spin = np.zeros((L,L))

hotstart(spin,L)
plt.imshow(spin)
#plt.show()
grid = np.zeros((frames,dim, dim))
pause_time = 0.1
energy_array = []
for n in range(Nthermalization):
    for i in range(L):
        for j in range(L):
            metropolis(spin,beta,[i,j])
    grid[n] = spin 
    '''plt.imshow(spin)
    plt.xticks([])
    plt.yticks([])
    #display.display(plt.gcf())
    #display.clear_output(wait=True)
    #time.sleep(pause_time)
    energy_array.append(energy(spin))
    plt.savefig("image_%03d.png" %n)'''




def animate(i):
    #grid[i] = np.random.choice([0, 1], size = (dim, dim))

    im = ax.imshow(grid[i], cmap = "bwr")

    return [im]

fig, ax = plt.subplots()
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])

ani = FuncAnimation(fig, animate, frames = frames, interval = 200, blit = True)
ani.save("l_Ising_2.gif", writer = "pillow", fps = 5)

    
'''
print(energy_array)        
plt.plot(n_array,energy_array)
plt.show()
#print(spin)
#spin(1,:) == spin[1]
'''
