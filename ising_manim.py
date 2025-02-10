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

        #self.play(self.zoomed_camera.frame.animate.)
        self.play(self.zoomed_camera.frame.animate.scale(4).shift(3*eps * LEFT))
        #self.remove(xy_p, xy_arr)


class ising_metropolis(MovingCameraScene):
    def construct(self):
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

        dot_coordinates = [-3*eps,0,0]
        dot = Dot(point = dot_coordinates, radius = 0.08, color = YELLOW)
        self.add(dot)
        self.wait(3)

        self.play(
            self.camera.frame.animate.scale(0.3).move_to(dot_coordinates)
        )

        point = 25
        
        self.wait(3)
        self.play(Indicate(ar[point - L],scale_factor=2.0),
                  Indicate(ar[point + L],scale_factor=2.0),
                  Indicate(ar[point + 1],scale_factor=2.0),
                  Indicate(ar[point - 1],scale_factor=2.0),run_time=3
                  )
        self.wait(3)

        nine_spins = VGroup(ar[point],ar[point+1],ar[point-1],ar[point+L],ar[point-L], \
                            ar[point+L+1],ar[point+L-1],ar[point-L+1],ar[point-L-1])
        #nine_spins = ar[point]

        #for x in [-1,1]:
        #    for y in [-L,L]:
        #        nine_spins += ar[x+y+point]

        #for x in [-1,1,-L,L]:
        #    nine_spins += ar[x+point]
        
 


        self.remove(ar)
        self.add(nine_spins)
        
        #self.play(FadeIn(nine_spins),FadeOut(ar))

        i = point//L
        j = point%L
        beta_text = MathTex(f"\\beta = {beta}").shift([-2*eps,eps,0]).scale(0.5)
        DS_text = MathTex("\\Delta S = 2\\beta s_{i,j}( s_{i+1,j} +  s_{i-1,j} +  s_{i,j+1} + s_{i,j-1}) ").shift([-3*eps,-0.5*eps,0]).scale(0.3)
        DS_text2 = MathTex(f"\\Delta S = 2({beta})({int(spin[i,j])}) (({int(spin[i+1,j])}) +  ({int(spin[i-1,j])}) +  ({int(spin[i,j+1])}) + ({int(spin[i,j-1])})) = {DS(spin,[i,j],beta)}").shift([-3*eps,-0.5*eps,0]).scale(0.3)

        self.play(nine_spins.animate.shift(LEFT+0.5*UP).scale(0.5),dot.animate.shift(LEFT+0.5*UP).scale(0.5))
        self.play(FadeIn(beta_text),FadeIn(DS_text))
        self.wait(2)
        self.play(TransformMatchingTex(DS_text,DS_text2))
        self.wait(3)

        prob_text = MathTex(f"p = \\exp(-\\Delta S)").shift([-3*eps,-eps,0]).scale(0.5)
        prob_text2 = MathTex(f"p = \\exp({-DS(spin,[i,j],beta)})").shift([-3*eps,-eps,0]).scale(0.5)
        prob_text3 = MathTex(f"p = {np.exp(-DS(spin,[i,j],beta)):1.4f}").shift([-3*eps,-eps,0]).scale(0.5)
        
        self.play(FadeIn(prob_text))
        self.wait(2)

        self.play(TransformMatchingTex(prob_text,prob_text2))
        self.wait(3)
        
        self.play(TransformMatchingTex(prob_text2,prob_text3))
        self.wait(3)
        r = np.random.rand(20)

        
        p = np.exp(-DS(spin,[i,j],beta))
        for l in range(20):
            text = MathTex(f"r = ",f"{r[l]:1.4f}",f"\\in [0,1)").shift([-1.8*eps,0,0]).scale(0.5)
            self.add(text)
            self.wait(0.2)
            self.remove(text)
        
        self.add(text)
        self.wait(2)

        if r[19] < p:
                   prob3 =  MathTex(f"p = {p:1.4f} \\geq {r[19]:1.4f}").shift([-3*eps,-eps,0]).scale(0.5)
                   
        else:
                   prob3 =  MathTex(f"p = {p:1.4f} < {r[19]:1.4f}").shift([-3*eps,-eps,0]).scale(0.5)

        
        self.play(Transform(prob_text3,prob3))
        self.wait(3)
