from dataclasses import dataclass, field
from math import dist
import numpy as np

@dataclass
class Constructor:

    def __init__(self, tmax:float, inst:list):
        '''
        '''
        self.tmax = tmax
        self.inst = inst
        self.k = 0
        self.V = self.sort()

    def sort(self) -> list:
        '''
            return an ordered list of it's furthest points based on the start position V[0]

            return ordered

            V: matrix with coordinates of all the other points
        '''

        empty = True    
        start = self.inst[0]
        
        for i in range(1,len(self.inst)-1):
            point = self.inst[i]
            d = dist(start[:2],point[:2])
            new_point = [np.concatenate((point, [d]))]
            if empty:
                ordered=np.array(new_point)
                empty = False
            else: ordered = np.concatenate((ordered, new_point))
        
        ordered = ordered[(-ordered[:,3]).argsort()]  
        ordered = np.concatenate(([start+[0]],ordered, [self.inst[-1]+[0]]))      
        return ordered

    def initialization(self, k:int, interval:list):
        '''
            returning the k paths initialized

            return kpaths

            k: amount of paths to initialize
            interval: interval to search for points
            tmax: maximum path time allowed between start and stop time 
            V: available points to choose, sorted by the furthest distances from start point 

            if minimum value in the interval is less than 0, interval started in 1 
            is assumed, if maximum value in the intervall exceeds max lenght, max 
            lenght-1 is assumed
        '''

        start = self.V[0][:3]
        stop = self.V[-1][:3]
        n = len(self.V)

        kpaths = []
        remain = list(range(1,n-1))
     
        interval.sort()
        si = max(interval[0], 1)
        sf = min(interval[1], n-1)
        intPoints = list(range(si,sf)) # points in given interval

        counter = 0
        while intPoints != [] and counter < k:

            p = np.random.choice(intPoints, replace=False)
            t = dist(start[:2], self.V[p][:2]) + dist(self.V[p][:2], stop[:2]) - dist(start[:2], stop[:2])
            if t < self.tmax:
                kpaths.append([0,p,n-1])
                remain.remove(p)
                counter +=1     
            intPoints.remove(p)

        if counter < k:
            for _ in range(k-counter):
                kpaths.append([0,n-1])

        return kpaths, remain

    def path_len(self,path:list) -> float:
        '''
            return the path length

            path : indexes of visited vertices
            V : [[x0,y0, reward0]....[xn,yn, rewardn]]
            
        '''
        totalDist = 0
        for i in range(len(path)-1):
            vi = path[i]
            vj = path[i+1]
            totalDist += dist(self.V[vi][:2],self.V[vj][:2])
        return totalDist
    
    def cheapest_insertion(self, kpaths:list, remain: list, time_kpaths:list):
        '''
            return populated kpaths based on cheapest insertion between all available points

            kpath : k-path indexes of visited vertices
            remain: remain vertices not visited
            V : [[x0,y0, reward0]....[xn,yn, rewardn]]
        '''
        lenremain = len(remain)
        j = 0 
        while j < lenremain:

            ri = remain[j]
            best = []
            z1 = self.tmax

            for i in range(self.k):

                time = self.path_len(kpaths[i], self.V)

                for m in range(len(kpaths[i])-1):
                    p1 = kpaths[i][m]
                    p2 = kpaths[i][m+1]
                    t = dist(self.V[p1][:2], self.V[ri][:2])+dist(self.V[p2][:2], self.V[ri][:2])-dist(self.V[p2][:2], self.V[p1][:2])
                    if (t+time < self.tmax) and (z1 > t):
                        z1 = t
                        best = [i,m+1, t+time]
            
            if len(best) != 0:
                i,m,t = best
                time_kpaths[i] = t
                kpaths[i].insert(m, ri)
                lenremain -=1
                remain.remove(ri)
            else:
                j+=1

        return True

    def cost_benefit_insertion(self, kpaths: list,remain:list, time_kpaths:list) -> bool: 
        '''
            return populated kpaths based on cost_benefit of insertion (z2)

            kpath : k-path indexes of visited vertices
            remain: remain vertices not visited
            time_kpaths: list of times of each path
            V : [[x0,y0, reward0]....[xn,yn, rewardn]]
        '''           

        insert = True
        while insert:
            z2 = 0
            best = []
            for i in range(self.k):
                time = self.path_len(kpaths[i])
                for m in range(len(kpaths[i])-1):
                    p1 = kpaths[i][m]
                    p2 = kpaths[i][m+1]
                    for ri in remain:
                        t = dist(self.V[p1][:2], self.V[ri][:2])+dist(self.V[p2][:2], self.V[ri][:2])-dist(self.V[p2][:2], self.V[p1][:2])
                        if t==0: t=10e-10
                        z = self.V[ri][2]*np.sqrt(self.V[ri][2])/t
                        if (z > z2) and (t+time < self.tmax):
                            best = [i,m+1, ri, t+time]
                            z2 = z
            
            if len(best) != 0:
                i,m, ri, t = best
                time_kpaths[i] = t
                kpaths[i].insert(m, ri)
                remain.remove(ri)
                insert = True
            else:
                insert = False

        return True

    def best_economy(self, kpaths: list, gamma: float, remain: list, time_kpaths:list) -> bool:
        '''
            return populated kpaths based on best economy insertion between all available points

            kpath : k-path indexes of visited vertices
            remain: remain vertices not visited
            V : [[x0,y0, reward0]....[xn,yn, rewardn]]
        '''
        lenremain = len(remain)
        j = 0
        start = kpaths[0][0]

        while j < lenremain:

            ri = remain[j]
            best = []
            z3 = self.tmax
            
            for i in range(self.k):

                time = self.path_len(kpaths[i], self.V)

                for m in range(len(kpaths[i])-1):
                    p1 = kpaths[i][m]
                    p2 = kpaths[i][m+1]
                    t = dist(self.V[p1][:2], self.V[ri][:2])+dist(self.V[p2][:2], self.V[ri][:2])-dist(self.V[p2][:2], self.V[p1][:2])
                    z = t - gamma*2*(dist(self.V[start][:2], self.V[ri][:2]))
                    if (t+time < self.tmax) and (z3 > z):
                        z3 = z
                        best = [i,m+1,t+time]
            
            if len(best) != 0:
                i,m, t = best
                time_kpaths[i] = t
                kpaths[i].insert(m, ri)
                lenremain -=1
                remain.remove(ri)
            else:
                j+=1

        return True

    def lrc(self, k:int, interval:list, insertion=3, gamma=0.1):
        '''
            return k-paths populated 
            
            k: amount of paths to initialize
            interval: interval to search for points
            V: set of points to choose
            insertion: insertion method for initialization, 
                       1 - cheaper insertion with bonus 
                       2 - best cost-benneif of insertion 
                       3 or anything else - cheapest insertion
        '''
        self.k = k
        kpaths, remain = self.initialization(k, interval)
        time = np.zeros(k)
        if insertion==1:
            self.best_economy(kpaths, gamma, remain, time)
        elif insertion==2:
            self.cost_benefit_insertion(kpaths, remain, time)
        else:
            self.cheapest_insertion(kpaths, remain, time)
        
        return kpaths, remain, time
        
    def reward(self,kpaths:list):
        '''
            return the team reward 

            kpaths :  index of V that team member k run
            V: set of points to choose

        '''
        r = 0
        for path in kpaths:
            for pi in path:
                r += self.V[pi][2]
        return r
