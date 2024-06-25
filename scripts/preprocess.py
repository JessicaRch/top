from math import dist
import numpy as np

class Ellipse:

    def __init__(self):
        pass

    def coefficients(self, foci0: list, foci1: list, maj_axis: float) -> list:
        '''
            this function finds the ellipse's: center, minimun axis and rotation's sin, cos

            return min_axis, center, rotation

            foci0,foci1 : ellipse focis points
            maj_axis : ellipse's major axis
            
        '''
        
        f = dist(foci0,foci1)/2
    
        if f == 0: # then its a circle
            center = foci0
            min_axis = maj_axis
            sinA = 0
            cosA = 1
        else:
            # finding min axis value
            min_axis = np.sqrt(maj_axis**2 - f**2)

            # finding ellipse center
            cx = (foci1[0]+foci0[0])/2
            cy = (foci1[1]+foci0[1])/2
            center = [cx,cy]

            # finding rotated cos, and sin 
            sinA = (foci1[1]-foci0[1])/(2*f)
            cosA = (foci1[0]-foci0[0])/(2*f)
        
        return min_axis, center, [sinA, cosA]
    

    def rotation(self,point:tuple,sin0:float, cos0:float) -> tuple:
        '''
            rotate point by sin, cos of rotation's ang

            return rot_point

            point: point to be rotated
            sin0, cos0: rotation's angle sin, cos 
        '''
        rot = [[cos0,-sin0],[sin0,cos0]]
        rot_point = np.dot(rot,point)
        return rot_point
    

    def translate(self, point:tuple, ref:tuple):
        '''
            translate a point given a point reference 
          
            return trans_point

            point: point to be translated
            ref: translation's reference point
        '''
        trans_point = [pi-ri for (pi,ri) in zip(point, ref)]
        return trans_point

    def is_inside(self, point: tuple, center: tuple, a: float, b: float, sin0: float, cos0: float) -> bool:
        '''
            verify if a given point is inside a given ellipse
            True if it is in, False, otherwise

            return in

            center: ellipse center
            a: maximum ellipse axis
            b: minimum ellipse axis
            sin0, cos0:  sin, cos of angle rotation
        '''

        rot_point = self.rotation(point, sin0, cos0)
        rot_center = self.rotation(center, sin0, cos0)

        trans_point = self.translate(rot_point, rot_center)

        f = (trans_point[0]**2)/(a**2) + (trans_point[1]**2)/(b**2)
        if f <= 1:
            return True
        else:
            return False

class Preprocess:

    def __init__(self, tmax:float, n:int, m:int, data:list):
        self.tmax = tmax
        self.n = n
        self.m = m
        self.data = self.viable_points(data)
        self.A = self.create_A()

    def viable_points(self,data):
        '''
            return available points 

            this function will select all points within a range area
            that are possible to arrive given a time interval and a start and stop points
        '''
        
        start = data[0]
        stop = data[-1]

        ellipse = Ellipse()

        min_axis, center, rotation = ellipse.coefficients(start[:2], stop[:2], self.tmax)
        ava_points = [start]

        for i in range(1, self.n-1):
            
            point = data[i]

            is_in = ellipse.is_inside(point[:2], center, self.tmax, min_axis, rotation[0], rotation[1])
            if is_in:
                ava_points.append(point)

        ava_points.append(stop)

        return ava_points
    

    def create_A(self):
        '''
            return a distance matrix from each vertice
        '''
        V = set(range(self.n))
        A = np.zeros((self.n, self.n))

        for i in V[:-1]:
            for j in V[1:]:
                A[i][j] = dist(self.data[i][:2],self.data[j][:2])

        return A


    def chao_initialization(self):
        '''
            return initialized k-paths, remain vertices

            given an amount of paths and an ordered list of points based on the furthest distance
            to the start, construct the firsts k-paths
        '''
        remain = list(range(1, self.n-1))
        N = self.n-2
        kpaths = []

        for _ in range(self.m):
       
            path = [0,N+1] # inserting start and stop indexes
            time = self.A[0][self.n-1] # inicial time
            len_remain = len(remain)
            i = 0
            j = 0
            
            while i < len_remain:
                
                r = remain[i]
                prev = path[j]
                next = path[j+1]
                t = -self.A[prev][next] + self.A[next][r] + self.A[prev][r]

                if (t+time) < self.tmax:
                    j+=1
                    path.insert(-1,r)
                    remain.pop(i)
                    i-=1
                    len_remain-=1
                    time += t    
                i+=1

            kpaths.append(path)
        return kpaths, remain