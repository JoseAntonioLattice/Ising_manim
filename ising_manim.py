from manim import *
#from lua import *

import random

L = 10
im = np.array(L)
ip = np.array(L)
eps = 0.8

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
    


class grid_arrow(Scene):
    def construct(self):
        '''x_a = np.arange(-3, 4)
        for x in x_a:
            for y in x_a:
                dot = Dot(np.array([x, y, 0]) , radius = 0.08)
                self.add(dot)

        dot = Dot(point = ORIGIN, radius = 0.08, color = YELLOW)
        self.add(dot)
        self.wait(0.5)
        '''
        samples = 10
        spin = np.zeros((10,10))
        hotstart(spin,10)
        #spin[:,:] = 1
        def draw_lattice(s):
            group_arrow = VGroup()
            for i in range(10):
                for j in range(10):
                #x_r = np.random.choice(x_a)
                #y_r = np.random.choice(x_a)
                #r = np.random.rand()
                    if s[i,j] == 1:
                        color = RED
                        arr = "UP"
                        y_r_d = j + 0.5
                    else:
                        color = BLUE
                        arr = "DOWN"
                        y_r_d = j - 0.5

                    arr = Arrow(start = [i*eps, j*eps, 0], end = [i*eps, y_r_d*eps, 0],max_stroke_width_to_length_ratio = 55, buff = 3, color = color).shift(eps*LEFT*L/2+eps*DOWN*L/2)
                    group_arrow += arr
                    #self.add(arr)
            return group_arrow


        ar = draw_lattice(spin)
        self.add(ar)
        self.remove(ar)

        for n in range(20):
            for i in range(L):
                for j in range(L):
                    metropolis(spin,1.0,[i,j])

            ar = draw_lattice(spin)
            self.add(ar)
            self.wait(0.5)
            self.remove(ar)
        #self.add(xy_p, xy_arr)
            
            #self.remove(xy_p, xy_arr)
