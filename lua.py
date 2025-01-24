import numpy as np
import random

L = 10
im = np.array(L)
ip = np.array(L)

ip = [(i+1)%L for i in range(L)]
im = [(i-1)%L for i in range(L)]

def hotstart(spin: int,L: int):
    for i in range(L):
        for j in range(L):
            if random.random() < 0.5: 
                spin[i,j] = 1
            else:
                spin[i,j] = -1

def DS(spin: int,x : int, beta: float) -> float:
    return 2*beta*spin[x[0],x[1]] * ( spin[ip[x[0]],x[1]] + spin[x[0],ip[x[1]]] + spin[im[x[0]],x[1]] + spin[x[0],im[x[1]]]  )
                
def metropolis(spin: int, beta: float , x: int):
    p = min(1,np.exp(-DS(spin,x,beta)))
    if random.random() <= p:
        spin[x[0],x[1]] = -spin[x[0],x[1]]

def energy(spin):
    energy = 0.0
    for i in range(L):
        for j in range(L):
            energy += spin[i,j]*(spin[i,ip[j]]+spin[ip[i],j])
    return energy/L**2
    
