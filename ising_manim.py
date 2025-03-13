from manim import *
#from lua import *

import random

L = 10
im = np.array(L)
ip = np.array(L)
eps = 0.7
beta = 0.5

Nconfs=11

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
    return -energy/L**2
    


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


class configurations(Scene):
    def construct(self):

        spin = np.zeros((Nconfs,L,L))
        spin = np.load("configurations.npy")
        squares = VGroup()
        confs = VGroup()

        k = 0
        for j in [1,-1]:
            for i in [-1, 0, 1]:
                cir = Circle(1.5)
                sq = Square(3)
                #hotstart(spin,L)
                lat = draw_lattice2(spin[k]).scale(3/L).move_to([0,0,0])
                squares += sq
                confs += lat
                self.play(FadeIn(confs[k]),FadeIn(squares[k]))
                self.play(confs[k].animate.to_edge([i,j,0]),squares[k].animate.to_edge([i,j,0]), run_time = 2)
                k +=1
                
        self.wait(2)
        self.play(squares.animate.move_to(10*UP), confs.animate.move_to(10*UP))

        squares_2 = VGroup()
        confs_2 = VGroup()

        k_2 = 0
        for i in [-1, 0, 1]:
            sq_2 = Square(3)
            #hotstart(spin,L)
            lat_2 = draw_lattice2(spin[k]).scale(3/L).move_to([0,0,0])
            squares_2 += sq_2
            confs_2 += lat_2
            self.play(FadeIn(confs_2[k_2]),FadeIn(squares_2[k_2]))
            self.play(confs_2[k_2].animate.to_edge([i, 1, 0]),squares_2[k_2].animate.to_edge([i, 1, 0]), run_time = 2)
            k_2 +=1
            k+=1

        sq_2 = Square(3)
        hotstart(spin,L)
        lat_2 = draw_lattice2(spin[k]).scale(3/L).move_to([0,0,0])
        squares_2 += sq_2
        confs_2 += lat_2
        self.play(FadeIn(confs_2[k_2]),FadeIn(squares_2[k_2]))
        self.play(confs_2[k_2].animate.to_edge([-1, -1, 0]),squares_2[k_2].animate.to_edge([-1, -1, 0]), run_time = 2)
        k_2 +=1
        k+=1

        
        
        #self.play(FadeIn(inf), shift = DOWN)
        
        sq_2 = Square(3)
        #hotstart(spin,L)
        lat_2 = draw_lattice2(spin[k]).scale(3/L).move_to([0,0,0])
        squares_2 += sq_2
        confs_2 += lat_2
        self.play(FadeIn(confs_2[k_2]),FadeIn(squares_2[k_2]))
        self.play(confs_2[k_2].animate.to_edge([0, -1, 0]),squares_2[k_2].animate.to_edge([0, -1, 0]), run_time = 2)
        k_2 +=1
        k+=1

        cont = MathTex('\\cdots', font_size = 120)
        cont.to_edge([1,-1,0])
        #inf = MathTex('\\infty', font_size = 120)
        #inf.shift(2*DOWN)
        
        #VGroup(cont, inf)
        
        self.play(FadeIn(cont))
        self.wait()

        cf = confs_2[0].copy()
        sqf = squares_2[0].copy()
        
        self.add(cf,sqf)
        self.play(FadeOut(squares_2,confs_2,cont))
        self.play(cf.animate.move_to([0,0,0]).scale(L/4),sqf.animate.move_to([0,0,0]).scale(L/4),run_time = 3)


class configurations2(Scene):
    def construct(self):

        spin = np.zeros((Nconfs,L,L))
        spin = np.load("configurations.npy")
        

        squares = VGroup()
        confs = VGroup()

        k = 0
        for j in [1,-1]:
            for i in [-1, 0, 1]:
                sq = Square(3).move_to([0,0,0])
                #hotstart(spin,L)
                lat = draw_lattice2(spin[k]).scale(3/L).move_to([0,0,0])
                squares += sq
                confs += lat
                self.add(squares[k].to_edge([i,j,0]),confs[k].to_edge([i,j,0]))
                k += 1

        self.wait(0.5)
        for n in range(6):
            self.remove(confs[n])
            
    

        for l in range(10):
            for j in [1,-1]:
                for i in [-1, 0, 1]:    
                    #hotstart(spin,L)
                    lat = draw_lattice2(spin[k]).scale(3/L).move_to([0,0,0])
                    confs += lat
                    self.add(confs[k].to_edge([i,j,0]))
            
                    k += 1

            self.wait(0.5)
            for n in [k-6,k-5,k-4,k-3,k-2,k-1]:
                self.remove(confs[n])
    
            
        

class choose_point(Scene):
    def construct(self):

        spin = np.zeros((L,L))
        confs = np.load("configurations.npy")
        spin = confs[6]
        lat = draw_lattice2(spin).move_to([0,0,0]).scale(3/4)
        square = Square(3*L/4).move_to([0,0,0])
        

        self.add(lat,square)
        self.wait(2)
        text = Tex("Choose point lexicographically")
        self.play(Write(text))
        self.wait(2)
        self.remove(text)
        for j in np.linspace(-L/2+1,L/2,L):
            for i in np.linspace(-L/2+1,L/2,L): 
                dot = Dot(point = [3/4*(i-0.5),3/4*(j-0.5),0] , radius = 0.08 ,color = YELLOW)
                self.add(dot)
                self.wait(0.1)
                self.remove(dot)

        self.wait(2)


class choose_point_random(Scene):
    def construct(self):

        spin = np.zeros((L,L))
        confs = np.load("configurations.npy")
        spin = confs[6]
       
        lat = draw_lattice2(spin).move_to([0,0,0]).scale(3/4)
        square = Square(3*L/4).move_to([0,0,0])
        

        self.add(lat,square)
        self.wait(2)
        text = Tex("..or choose point randomly")
        self.play(Write(text))
        self.wait(2)
        self.remove(text)
    
        for x in range(L):
            for y in range(L):
                i = np.random.randint(L)-L/2+1
                j = np.random.randint(L)-L/2+1
                dot = Dot(point = [3/4*(i-0.5),3/4*(j-0.5),0] , radius = 0.08 ,color = YELLOW)
                self.add(dot)
                self.wait(0.1)
                self.remove(dot)


class create_lattice(Scene):
    def construct(self):

        spin = np.zeros((L,L))
        
        confs = np.load("configurations.npy")
        spin = confs[0]
        lat = draw_lattice2(spin).move_to([0,0,0]).scale(3/4)
        square = Square(3*L/4).move_to([0,0,0])

        self.add(square)
        self.wait(2)
        k = 0
        lattice_points = VGroup()
        for j in np.linspace(-L/2+1,L/2,L):
            for i in np.linspace(-L/2+1,L/2,L): 
                lattice_points += Dot(point = [3/4*(i-0.5),3/4*(j-0.5),0] , radius = 0.08 ,color = WHITE)

        self.play(Create(lattice_points))
        self.wait(2)
        self.play(FadeOut(lattice_points),Create(lat))
        self.wait(2)
        
        

class metropolis_step(MovingCameraScene):
    def construct(self):
        
        spin = np.zeros((L,L))
        
        confs = np.load("configurations.npy")
        spin = confs[6]
        lat = draw_lattice2(spin).move_to([0,0,0]).scale(3/4)
        square = Square(3*L/4).move_to([0,0,0])
        point = 38
        i = 3/4*(point//L - L/2+0.5)
        j = 3/4*(point%L - L/2+0.5)
        dot =  Dot(point = [i,j,0] , radius = 0.08 ,color = YELLOW)

        nine_spins = VGroup(lat[point].copy(),lat[point+1].copy(),lat[point-1].copy(),lat[point+L].copy(),lat[point-L].copy(), \
                            lat[point+L+1].copy(),lat[point+L-1].copy(),lat[point-L+1].copy(),lat[point-L-1].copy())

        

        self.add(square,nine_spins,lat)
        self.wait()

        

        self.play(Create(dot))

        self.wait(3)

        
        self.play(self.camera.frame.animate.scale(0.3).move_to([i,j,0]))

        self.wait(2)
        
        
        self.play(FadeOut(lat,square))
        self.wait()
        
        self.play(Indicate(nine_spins[1],scale_factor=2.0),
                  Indicate(nine_spins[2],scale_factor=2.0),
                  Indicate(nine_spins[3],scale_factor=2.0),
                  Indicate(nine_spins[4],scale_factor=2.0),run_time=3
                  )
        self.wait(2)

        beta = 0.5
        
        beta_text = MathTex(f"\\beta = {beta}").shift([i+0.5,j+0.5,0]).scale(0.5)
        DS_text = MathTex("\\Delta S = 2\\beta s_{i,j}( s_{i+1,j} +  s_{i-1,j} +  s_{i,j+1} + s_{i,j-1}) ").shift([i-0.1,j-0.5,0]).scale(0.3)
        
      
        i = point//L
        j = point%L

        DS_text2 = MathTex(f"\\Delta S = 2({beta})({int(spin[i,j])}) (({int(spin[i+1,j])}) +  ({int(spin[i-1,j])}) +  ({int(spin[i,j+1])}) + ({int(spin[i,j-1])})) = {DS(spin,[i,j],beta)}").shift(DS_text.get_center()).scale(0.3)

        self.play(nine_spins.animate.shift(LEFT+0.5*UP).scale(0.5),dot.animate.shift(LEFT+0.5*UP).scale(0.5))
        self.play(FadeIn(beta_text),FadeIn(DS_text))
        self.wait(2)
        self.play(TransformMatchingTex(DS_text,DS_text2))
        self.wait(3)
     
        prob_text = MathTex(f"p = \\exp(-\\Delta S)").shift(DS_text.get_center()+0.5*DOWN).scale(0.5)
        prob_text2 = MathTex(f"p = \\exp({-DS(spin,[i,j],beta)})").shift(prob_text.get_center()).scale(0.5)
        prob_text3 = MathTex(f"p = {np.exp(-DS(spin,[i,j],beta)):1.4f}").shift(prob_text.get_center()).scale(0.5)
        
        self.play(FadeIn(prob_text))
        self.wait(2)

        self.play(TransformMatchingTex(prob_text,prob_text2))
        self.wait(3)
        
        self.play(TransformMatchingTex(prob_text2,prob_text3))
        self.wait(3)
        r = np.random.rand(20)

        
        p = np.exp(-DS(spin,[i,j],beta))
        for l in range(20):
            text = MathTex(f"r = ",f"{r[l]:1.4f}",f"\\in [0,1)").shift(beta_text.get_center()+0.5*DOWN).scale(0.5)
            self.add(text)
            self.wait(0.2)
            self.remove(text)
        
        self.add(text)
        self.wait(2)
        self.remove(DS_text2,text)

        i = 3/4*(point//L - L/2+0.5)
        j = 3/4*(point%L - L/2+0.5)
        proposal = Rotate(nine_spins[0],angle=PI,about_point=dot.get_center())
        if r[19] < p:
            prob3 =  MathTex(f"p = {p:1.4f} \\geq {r[19]:1.4f}").shift(prob_text.get_center()).scale(0.5)
            text = Tex("Accept the change").shift(prob_text.get_center()+0.5*UP).scale(0.5)
            if nine_spins[0].get_color() == RED:
                t2 = Indicate(nine_spins[0])#Tex("rojo").shift([-3*eps,-0.1*eps,0]).scale(0.5)
                nine_spins[0].set_color(BLUE)
            else:
                t2 = Wiggle(nine_spins[0])#Tex("azul").shift([-3*eps,-0.1*eps,0]).scale(0.5)
                nine_spins[0].set_color(RED)        
            self.play(Transform(prob_text3,prob3),Create(text),proposal)
            self.wait(3)
        else:
            prob3 =  MathTex(f"p = {p:1.4f} < {r[19]:1.4f}").shift(prob_text.get_center()).scale(0.5)
            text = Tex("Keep configuration").shift(prob_text.get_center()+0.5*UP).scale(0.5) 
            #t2 = Tex("puta").shift([-3*eps,-0.1*eps,0]).scale(0.5)
            self.play(Transform(prob_text3,prob3),Create(text))
            self.wait(3)
       
    
        
def draw_lattice2(s):
    group_arrow = VGroup()
    for i in range(L):
        for j in range(L):
            #x_r = np.random.choice(x_a)
            #y_r = np.random.choice(x_a)
            #r = np.random.rand()
            
            if s[i,j] == 1:
                color = RED
                y_r_d = j + 0.5
            else:
                color = BLUE
                y_r_d = j - 0.5
            
            arr = Arrow(start = [i, j, 0], end = [i, y_r_d, 0],max_stroke_width_to_length_ratio = 55, buff = 1, color = color)
            group_arrow += arr
                
    return group_arrow

class pacman(Scene):
    def construct(self):

        openmouth = AnnularSector(inner_radius = 0, outer_radius = 1, angle = 315*DEGREES, start_angle = 22.5*DEGREES).set_fill(YELLOW).set_stroke(color = WHITE, width = 1.4).scale(0.5)
        closedmouth = Dot(radius = 0.5, color = YELLOW)
        dot = Dot(radius = 0.5, color = YELLOW)
        dot_left = dot.copy().shift([-8, 0, 0])
        text = Tex("Periodic boundary conditions in all directions")
        #self.play(dot.animate.move_to([8,0,0]),dot_left.animate.move_to([0,0,0]),run_time = 3 )
        text2 = Tex("Create lattice")
        square = Square(8)
        rect_right = Rectangle(height = 8, width = 4, color = BLACK, fill_opacity = 1).shift([6, 0, 0])
        rect_left = Rectangle(height = 8, width = 4, color = BLACK, fill_opacity = 1).shift([-6, 0, 0])
        rect_up = Rectangle(height = 4, width = 8, color = BLACK, fill_opacity = 1).shift([0, 6, 0])
        rect_down = Rectangle(height = 4, width = 8, color = BLACK, fill_opacity = 1).shift([0, -6, 0])
        box = Rectangle(color = BLACK,width=text.width,height=text.height,fill_opacity=1)
        box.surround(text)
        self.play(Create(text2))
        self.play(Create(square), run_time=2)
        self.wait(2)
        self.play(FadeOut(text2))
      

        lattice_points = VGroup()
        for j in np.linspace(-3.5, 3.5, 7):
            for i in np.linspace(-3.5, 3.5, 7): 
                lattice_points += Dot(point = [i, j, 0] , radius = 0.1, color = WHITE)

        self.play(Create(lattice_points))
        #self.add(rect_right, rect_left, rect_up, rect_down)
        self.wait(2)
        self.add(box)
        self.play(Create(text))
        self.wait(2)
        self.remove(text,box)
        
        self.play(FadeOut(lattice_points[24]),Create(openmouth), runt_time = 1)
        self.wait()
        self.remove(openmouth)
        
        pacman = closedmouth.copy()
        pacmanleft = pacman.copy().move_to([-8,0,0])

        dt = 0.1
        t = 0.0
        x = 0.0
        v = 2.0

        tp_ep_x = [0.6, 1.2, 1.8]
        tn_ep_x = [2.4, 3.0, 3.6]
        epp_x = 1
        epn_x = 3       
        
        for i in range(40):
            newpacman = pacman.copy().move_to([x,0,0])
            newpacmanl = pacman.copy().move_to([-8 + x,0,0])

            self.add(pacman.move_to(newpacman),pacmanleft.move_to(newpacmanl)) 
            self.add(rect_right, rect_left, rect_up, rect_down)
            self.add(square)
            self.wait(dt)
            t += dt
            x += v*dt
            self.remove(lattice_points[24 + epp_x]) if round(t, 1) in tp_ep_x else None
            epp_x += 1 if round(t, 1) in tp_ep_x else 0
            self.remove(lattice_points[24 - epn_x]) if round(t, 1) in tn_ep_x else None
            epn_x -= 1 if round(t, 1) in tn_ep_x else 0
            self.remove(pacman,pacmanleft,square,rect_right, rect_left, rect_up, rect_down)
            pacman = openmouth.copy() if i%2 == 0 else closedmouth.copy()
            pacmanleft = openmouth.copy() if i%2 == 0 else closedmouth.copy()
            pacman.move_to(newpacman)
            pacmanleft.move_to(newpacmanl)

        openmouthup = openmouth.copy()
        pacman = openmouthup
        self.add(pacman,square)
    
        
        #self.remove(pacmanleft)
        y = 0
        t = 0

        tp_ep_y = [0.6, 1.2, 1.8]
        tn_ep_y = [2.4, 3.0, 3.6]
        epp_y = 7
        epn_y = 21
        
        self.play(Rotate(pacman, angle = PI/2, rate_func = linear))
        
        
        pacmand = pacman.copy().move_to([0,-8,0])
        for i in range(40):
            newpacman = pacman.copy().move_to([0,y,0])
            newpacmand = pacman.copy().move_to([0,y-8,0])

            self.add(pacman.move_to(newpacman),pacmand.move_to(newpacmand))
            self.add(rect_right, rect_left, rect_up, rect_down)
            self.add(square)
            self.wait(dt)
            t += dt 
            y += v*dt
            self.remove(lattice_points[24 + epp_y]) if round(t, 1) in tp_ep_y else None
            epp_y += 7 if round(t, 1) in tp_ep_y else 0
            self.remove(lattice_points[24 - epn_y]) if round(t, 1) in tn_ep_y else None
            epn_y -= 7 if round(t, 1) in tn_ep_y else 0
            self.remove(pacman,pacmanleft,square,rect_right, rect_left, rect_up, rect_down)
            self.remove(pacman,pacmand,square,rect_right, rect_left, rect_up, rect_down)
            pacman = openmouthup.copy() if i%2 == 0 else closedmouth.copy()
            pacmand = openmouthup.copy() if i%2 == 0 else closedmouth.copy()
            pacman.move_to(newpacman)
            pacmand.move_to(newpacmand)

        self.add(openmouthup,square)
        self.wait(2)
        
        self.play(FadeOut(openmouthup))
        self.clear()



def draw_conf(s):
    group_sq = VGroup()
    for i in range(L):
        for j in range(L):
            #x_r = np.random.choice(x_a)
            #y_r = np.random.choice(x_a)
            #r = np.random.rand()
            
            if s[i,j] == 1:
                color = RED
                y_r_d = j + 0.5
            else:
                color = BLUE
                y_r_d = j - 0.5
            
            sq=Square(side_length=1,color = color, fill_opacity=1).move_to([i-L/2,j-L/2,0])
            group_sq += sq
    group_sq.scale(0.5)            
    return group_sq
        
class evolution(Scene):
    def construct(self):

        spin = np.zeros((Nconfs,L,L))
        confs = np.load("configurations.npy")

        spin = confs
        E = np.random.rand(11)
        coords = []
        dots = VGroup()
        
        for i in range(len(E)):
            ax = Axes(
                x_range=[0, 10, 1],
                y_range=[0, 1, 0.1],
                tips=False,
                axis_config={"include_numbers": True}
            )
            
            # x_min must be > 0 because log is undefined at 0.
            
            #self.add(ax, graph)dots = VGroup()
            #c = draw_conf(spin[i])
            #self.add(c)
            #self.wait(0.5)
            #self.remove(c)
            #coords.append(x)
            dots += Dot(radius=0.1,point=[i,E[i],0])
            self.add(ax,dots)
            self.wait(0.5)
        
        

        


class en_graph(Scene):
    def construct(self):

        confs = np.load("configurations.npy")
        coords = VGroup()
        samples = 100
        x = np.arange(1, samples, 1)
        en = [energy(confs[i]) for i in range(len(confs))] #np.random.rand(samples)
        y =  en
    
        ''' ax = Axes(
              x_range = [0, samples, 10],
              y_range = [0, 1, 1],
              tips = False,
              axis_config = {'include_numbers': True},
              )

        self.add(ax)
       '''
        plane = Axes(
            x_range = (0, len(confs),10),
            y_range = (-2,1,0.2),
            x_length = 7,
            y_length = 7,
            axis_config={"include_numbers": True},
            tips=False
        )
        self.add(plane)
        plane.center()
        for i in range(len(confs)):
            
            line_graph = plane.plot_line_graph(
                x_values = x[0:i+1],
                y_values = y[0:i+1],
                line_color=GOLD_E,
                vertex_dot_style=dict(stroke_width=0,  fill_color=PURPLE),
                #stroke_width = 4,
            )
            self.add(line_graph)
            self.wait(0.1)
            self.remove(line_graph)
        
        '''for i in range(samples):
            coord = Dot(ax.c2p(i, y[i]), color = GREEN)
            coords += coord
            self.add(coords)
            self.wait(0.1)
        '''

