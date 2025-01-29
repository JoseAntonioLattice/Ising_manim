from manim import *
#from lua import *

import random

L = 10
im = np.array(L)
ip = np.array(L)
eps = 0.7
beta = 0.5

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
    


class grid_arrow(ZoomedScene):
    def __init__(self, **kwargs):
        ZoomedScene.__init__(
            self,
            zoom_factor=0.3,
            zoomed_display_height=3*eps,
            zoomed_display_width=3*eps,
            image_frame_stroke_width=20,
            zoomed_camera_config={
                "default_frame_stroke_width": 3,
            },
            **kwargs
        )

    
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
        spin = np.zeros((L,L))
        hotstart(spin,L)
        #spin[:,:] = 1
        def draw_lattice(s):
            group_arrow = VGroup()
            for i in range(L):
                for j in range(L):
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

        beta_text = MathTex(f"\\beta = {beta}").shift(3.5*UP)
        self.add(beta_text)
        ar = draw_lattice(spin)
        self.add(ar)
        self.wait(0.5)
        self.remove(ar)

        for n in range(10):
            for i in range(L):
                for j in range(L):
                    metropolis(spin,0.1,[i,j])

            ar = draw_lattice(spin)
            self.add(ar)
            self.wait(0.5)
            self.remove(ar)
        #self.add(xy_p, xy_arr)
        self.add(ar)
        dot = Dot(point = [-3*eps,0,0], radius = 0.08, color = YELLOW)
        self.add(dot)
        self.wait(3)

        self.activate_zooming(animate=False)
        self.play(self.zoomed_camera.frame.animate.scale(4))
        self.play(self.zoomed_camera.frame.animate.shift(3*eps * LEFT))
            #self.remove(xy_p, xy_arr)
