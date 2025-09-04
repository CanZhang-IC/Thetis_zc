# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import math
import copy
import numpy as np
from ad import adnumber

from firedrake import *
from firedrake_adjoint import *

import random
from the_CW import Clarke_Wright


class CableCostGA(object):
    
    def __init__(self, turbine_locations, substation_location = [[0,0]], capacity = 6, pop_size = 120, num_iter = 200, convergence_definition = 20, converged = False, show_prog = False, show_result = False, figname ='fig'): 
        self.turbine_locations = []
        for i in range(int(len(turbine_locations)/2)):
            self.turbine_locations.append([turbine_locations[2*i],turbine_locations[2*i+1]])      
        self.turbine_locations_check = []
        for i in self.turbine_locations:
            if i not in self.turbine_locations_check:
                self.turbine_locations_check.append(i)
        self.substation_location = [[substation_location[0],substation_location[1]]]
        self.capacity = capacity
        self.pop_size = pop_size
        self.num_iter = num_iter
        self.convergence_definition = convergence_definition
        self.show_prog = show_prog
        self.show_result = show_result
        self.figname = figname


    def convert_to_adnumber(self, coordinate_list):
        '''Convert the location vectors from floats into adnumbers to enable differentiation'''
        adnumber_coordinate_list = []
        for i in range(len(coordinate_list)):
            adnumber_coordinate_list.append([adnumber(coordinate_list[i][0]), adnumber(coordinate_list[i][1])])
        coordinate_list = adnumber_coordinate_list
        return coordinate_list
        
    
    def differentiate(self, best_chromosome):
        '''Differentiate the length of the routing w.r.t. the position of the turbines, produces a n x 2 array, ((dC/dx1, dC/dy1), ...)'''
        print ('Determining dC/dX: Performing automatic differentiation...')
        turbine_locations = self.convert_to_adnumber(self.turbine_locations)
        substation_location = self.convert_to_adnumber(self.substation_location)
        vertices = substation_location + turbine_locations
        C = self.construct_cost_matrix(vertices)
        route = best_chromosome[0]
        breaks = best_chromosome[1]    
        rting = self.routing(route, breaks)
        total_dist = self.routing_distance(rting, C)
        n = len(turbine_locations)
        dC_dX = np.zeros((n,2)) 
        for i in range(n):
            dC_dX[i][0] = total_dist.d(turbine_locations[i][0])
            dC_dX[i][1] = total_dist.d(turbine_locations[i][1])
        print ('Automatic differentiation complete...')
        return dC_dX
    
    
    def construct_cost_matrix(self, vertices):
        '''Constructs a matrix of costs for every potential edge - connections of vertices to themselves are set to infinity'''
        distances = []
        grouped_distances = []
        for i in range(len(vertices)):
            for j in range(len(vertices)):
                if i == j:
                    dist = 99
                else:
                    #print(vertices[i].values())
                    dist = (sqrt(float((vertices[i][0] - vertices[j][0])**2 + (vertices[i][1] - vertices[j][1])**2)))
                distances.append(dist)
        for i in range(0, len(distances), len(vertices)):
            grouped_distances.append(tuple(distances[i:i+len(vertices)]))
        C = np.array(grouped_distances)
        return C
    
    
    def produce_plot(self, vertices, best_chromosome, global_min):
        '''Display with matplotlib'''
        plt.figure()
        opt_route = best_chromosome[0]
        opt_breaks = best_chromosome[1]
        R = self.routing_coordinates(vertices, self.routing(opt_route, opt_breaks))
        plt.title('Total distance='+str(global_min))
        V = np.array((vertices))
        plt.plot(V[:,0],V[:,1],'o')
        for i in range(len(R)):
            plt.plot(R[i][:,0], R[i][:,1], '-')
        for i in range(len(V)):
            plt.text(V[i][0], V[i][1], '%s' % (str(i)))
        plt.axis('equal')
        plt.savefig('./'+self.figname + '.jpg',dpi=300)
        
    
    def rand_breaks(self, n, min_route, n_routes, n_breaks):
        '''Produces a list of random, but valid, route-breaks'''
        RB = [np.random.randint(min_route, self.capacity+1)]
        for i in range(1, n_breaks):
            RB.append(np.random.randint(min_route, self.capacity+1) + RB[i-1])
        if RB[-1] < (n - self.capacity):
            short = (n - self.capacity) - RB[-1]
            add_each = int(np.ceil(0.5 + short / len(RB)))
            for i in range(len(RB)):
                RB[i] = RB[i] + add_each * (i+1)
        return RB    
    
        
    def initialise_population(self, n, min_route, n_routes, n_breaks):
        '''Randomly produces an initial population'''
        popbreaks = []
        poproutes = []
        for i in range(self.pop_size):
            popbreaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
            poproutes.append(random.sample(range(1, n+1), n))
        return poproutes, popbreaks
    
    
    def routing(self, route, breaks):
        '''Combine route and breaks into an array of the routes described as vertices in the order in which they are toured'''
        rting = [[0] + route[0:breaks[0]]]
        if len(breaks) > 1:
            for f in range(1, len(breaks)):
                rting.append([0] + route[breaks[f-1]:breaks[f]])
        rting.append([0] + route[breaks[-1]:])
        return rting
    
    
    def routing_distance(self, routing, C):
        '''Return the geometric length of the routing'''
        d = 0        
        for r in range(len(routing)):
            for v in range(len(routing[r]) - 1):
                d += C[routing[r][v]][routing[r][v+1]]
        return d
    
    
    def routing_coordinates(self, vertices, rting):
        '''Convert a routing expressed in indexes to a routing expressed in coordinates''' 
        for i in range(len(rting)):
            for j in range(len(rting[i])):
                rting[i][j] = vertices[rting[i][j]]
            rting[i] = np.array(rting[i])
        return rting
    
    
    def evaluate_population(self, pop, C, n):
        '''produce a list of the functional evaluations of each member of the population'''
        D = []
        # print(len(pop[0]),self.pop_size)
        for i in range(len(pop[0])):
            D.append(self.routing_distance(self.routing(pop[0][i], pop[1][i]), C))
        return D
        # D = []
        # for i in range(self.pop_size):
        #     D.append(self.routing_distance(self.routing(pop[0][i], pop[1][i]), C))
        # return D
        
        
    def clarke_wright(self, route, breaks, vertices):
        '''apply clarke and wright algorithm to chromosome'''
        rting = self.routing(route, breaks)
        rting_coords = self.routing_coordinates(vertices, rting)     
        rting = []
        for i in range(len(rting_coords)):
            temp_r = rting_coords[i].tolist()
            #print(temp_r)
            del temp_r[0]
            CW = Clarke_Wright(temp_r, self.substation_location)
            temp_r = CW.run()
            del temp_r[0]
            #print(temp_r)
            for j in range(len(temp_r)):
                # print(temp_r)
                temp_r[j] = self.turbine_locations.index(temp_r[j]) + 1
            rting += temp_r
        return rting
    
    
    def transformations(self, bestof8route, bestof8breaks, n, min_route, n_routes, n_breaks, vertices, iteration):
        '''Returns a copy of the original chromosome and 7 mutations'''
        n = 8
        selector = random.sample(range(1, n), n-1)
        randlist = [selector[int(n/3)], selector[int(2*n/3)]]
        I = min(randlist)
        J = max(randlist)
        temp_pop_route = []
        temp_pop_breaks = []
        ### 添加上一代最好的
        temp_pop_route.append(bestof8route)
        temp_pop_breaks.append(bestof8breaks)
        # transformation 1  
        trans_1 = copy.deepcopy(bestof8route)
        trans_1[I], trans_1[J] = trans_1[J], trans_1[I]
        temp_pop_route.append(trans_1)
        temp_pop_breaks.append(bestof8breaks)
        ## transformation 2
        trans_2 = copy.deepcopy(bestof8route)
        Temp = trans_2[I:J]; Temp = Temp[::-1]; trans_2 = trans_2[:I] + Temp + trans_2[J:]
        temp_pop_route.append(trans_2)
        temp_pop_breaks.append(bestof8breaks)
        ## transformation 3
        trans_3 = copy.deepcopy(bestof8route)
        # trans_3.remove(I); trans_3.insert(J, I)
        Temp = trans_3[I:J]; trans_3 = trans_2[:I]  + trans_2[J:] + Temp
        temp_pop_route.append(trans_3)
        temp_pop_breaks.append(bestof8breaks)
        # transformation 4
        if iteration<15 :
            trans_4 = copy.deepcopy(bestof8route)
            trans_4 = self.clarke_wright(trans_4, bestof8breaks, vertices)
            temp_pop_route.append(trans_4)
            temp_pop_breaks.append(bestof8breaks)
        else:
            temp_pop_route.append(bestof8route)
            temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))

        ## transformation 5
        temp_pop_route.append(trans_1)
        temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
        ## transformation 6        
        temp_pop_route.append(trans_2)
        temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
        ## transformation 7        
        temp_pop_route.append(trans_3)
        temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
        
        return temp_pop_route, temp_pop_breaks
        
     
    def breed_population(self, pop, D, n, min_route, n_routes, n_breaks, vertices, iteration):
        '''Create a new population based upon the best performing specimens of the old population'''
        n = len(self.turbine_locations)        
        chrome_pop = []
        for i in range(self.pop_size):
            chrome_pop.append(pop[0][i] + pop[1][i] + [D[i]])
        random.shuffle(chrome_pop) ###打乱列表中元素的排列顺序
        new_pop_route = []
        new_pop_breaks = []
        for i in range(0, self.pop_size, 8):
            Dists = [] 
            for j in range(8):
                Dists.append(chrome_pop[i+j][-1])
            MIndx = np.argmin(Dists)#min(xrange(len(Dists)),key=Dists.__getitem__)
            bestof8Chrome = chrome_pop[MIndx+i]
            bestof8route = bestof8Chrome[:n]
            bestof8breaks = bestof8Chrome[n:-1]
            new_8 = self.transformations(bestof8route, bestof8breaks, n, min_route, n_routes, n_breaks, vertices, iteration)
            new_pop_route += new_8[0]
            new_pop_breaks += new_8[1]
        return new_pop_route, new_pop_breaks
        
            
    def run_GA(self):
        '''Run the genetic algorithm'''
        ###通过给seed赋值，保证每次随机数都是固定的，保证了可重复性。
        np.random.seed(100)
        random.seed(100)        
        vertices = self.substation_location + self.turbine_locations
        n_routes = int(math.ceil(float(len(vertices)) / self.capacity))
        min_route = len(vertices) / n_routes
        n = len(self.turbine_locations)
        n_breaks = n_routes - 1
        converged = False
        convergence_counter = 0
        #print ('Running Genetic Algorithm...')
        pop = self.initialise_population(n, min_route, n_routes, n_breaks) ###随机生成染色体（点的顺序和打断点的位置）
        C = self.construct_cost_matrix(vertices)
        D = self.evaluate_population(pop, C, n)  ### 计算每一段染色体代表的距离的长度
        MIndx = np.argmin(D)#min(xrange(len(D)),key=D.__getitem__)    ###返回axis方向上最小值的index，默认axis=0即横向。
        global_min = copy.deepcopy(D[MIndx]) ###deepcopy保证列表中的子列表也完成了复制，子列表也不会再跟随原来的列表发生变化
        best_chromosome = [copy.deepcopy(pop[0][MIndx]), copy.deepcopy(pop[1][MIndx])]
        dist_history = [global_min]
    
        for iteration in range(self.num_iter):
            if not converged:
                if self.show_prog:
                    print ('Current Routing Length', global_min)
                ###返回最优的点顺序和间断点顺序，格式为([[点顺序1],[点顺序2]],[[间断点顺序1],[间断点顺序2]])
                pop = self.breed_population(pop, D, n, min_route, n_routes, n_breaks, vertices, iteration)
                ###返回每条染色体（点顺序+间断点顺序）
                D = self.evaluate_population(pop, C, n)
                ###返回最短路径对应的index
                MIndx = np.argmin(D)#min(xrange(len(D)),key=D.__getitem__)
                ###对比此时的最短路径与历史最低的最短路径，如果更小，说明找到了更好的，进行替换。
                if D[MIndx] < global_min:
                    global_min = copy.deepcopy(D[MIndx])
                    best_chromosome = [copy.deepcopy(pop[0][MIndx]), copy.deepcopy(pop[1][MIndx])]
                    dist_history.append(global_min)
                    convergence_counter = 0
                else: 
                    convergence_counter += 1
                if convergence_counter > self.convergence_definition:
                    converged = True
        if self.show_result:
            self.produce_plot(vertices, best_chromosome, global_min)
        #print ('GA Complete, Routing Length is: ', global_min)
        return best_chromosome, global_min
        
    def compute_cable_cost(self):
        if len(self.turbine_locations) != len(self.turbine_locations_check):
            cost_out = 1e10
        else:
            out = self.my_run_GA()
            cost_out = out[1]
            self.best_chromosome = out[0]
            rting = self.routing(self.best_chromosome[0],self.best_chromosome[1])
        return cost_out,self.best_chromosome[0],self.best_chromosome[1]
     
    def compute_cable_cost_derivative(self):
        dC_dX = self.differentiate(self.best_chromosome)
        return dC_dX


    def my_run_GA(self):
        '''Run the genetic algorithm'''
        ###通过给seed赋值，保证每次随机数都是固定的，保证了可重复性。
        # np.random.seed(100)
        # random.seed(100)        
        vertices = self.substation_location + self.turbine_locations
        n_routes = int(math.ceil(float(len(vertices)) / self.capacity))
        min_route = len(vertices) / n_routes
        n = len(self.turbine_locations)
        n_breaks = n_routes - 1
        converged = False
        convergence_counter = 0
        #print ('Running Genetic Algorithm...')
        pop = self.initialise_population(n, min_route, n_routes, n_breaks) ###随机生成染色体（点的顺序和打断点的位置）
        C = self.construct_cost_matrix(vertices)

        best_chromosome = []
        global_min = 10e20
        betterNumof8 = 0
        for i in range(0, self.pop_size, 8):
            popof8Chrome = [pop[0][i:i+8]] + [pop[1][i:i+8]] ###每8个一组
            # print(i,popof8Chrome)
            
            Dof8Chrome = self.evaluate_population(popof8Chrome, C, n)###算出一组8个的距离
            MIndxof8Chrome = np.argmin(Dof8Chrome) ### 找出8个中最小的index
            ###8个中最好的添加到最好的组里面，用于保存遗传到下一代
            best_chromosome.append([copy.deepcopy(popof8Chrome[0][MIndxof8Chrome]), copy.deepcopy(popof8Chrome[1][MIndxof8Chrome])]) 
            ###对比此时的最短路径与历史最低的最短路径，如果更小，说明找到了更好的，进行替换。
            if Dof8Chrome[MIndxof8Chrome] < global_min:
                global_min = copy.deepcopy(Dof8Chrome[MIndxof8Chrome])
                betterNumof8 += 1

        dist_history = [global_min]
    
        for iteration in range(self.num_iter):
            if not converged:
                if self.show_prog:
                    print ('Current Routing Length', global_min)
                # print(best_chromosome)
                ###返回最优的点顺序和间断点顺序，格式为([[点顺序1],[点顺序2]],[[间断点顺序1],[间断点顺序2]])
                pop = self.my_breed_population(best_chromosome, n, min_route, n_routes, n_breaks, vertices, iteration)
                ###返回每条染色体（点顺序+间断点顺序）
                
                best_chromosome = []
                global_min = 10e20
                betterNumof8 = 0
                for i in range(0, self.pop_size, 8):
                    popof8Chrome = [pop[0][i:i+8]] + [pop[1][i:i+8]] ###每8个一组
                    # print(popof8Chrome)
                    Dof8Chrome = self.evaluate_population(popof8Chrome, C, n)###算出一组8个的距离
                    MIndxof8Chrome = np.argmin(Dof8Chrome) ### 找出8个中最小的index
                    ###8个中最好的添加到最好的组里面，用于保存遗传到下一代
                    best_chromosome.append([copy.deepcopy(popof8Chrome[0][MIndxof8Chrome]), copy.deepcopy(popof8Chrome[1][MIndxof8Chrome])]) 
                    ###对比此时的最短路径与历史最低的最短路径，如果更小，说明找到了更好的，进行替换。
                    if Dof8Chrome[MIndxof8Chrome] < global_min:
                        global_min = copy.deepcopy(Dof8Chrome[MIndxof8Chrome])
                        global_best_index = [i,MIndxof8Chrome]
                        betterNumof8 += 1
                if betterNumof8 == 0: 
                    convergence_counter += 1
                if convergence_counter > self.convergence_definition:
                    converged = True
        print(len(pop[0]),global_best_index[0],global_best_index[1])
        global_best_chromosome = [copy.deepcopy(pop[0][global_best_index[0]+7+global_best_index[1]]), copy.deepcopy(pop[1][global_best_index[0]+7+global_best_index[1]])]
        if self.show_result:
            self.produce_plot(vertices, global_best_chromosome, global_min)
        #print ('GA Complete, Routing Length is: ', global_min)
        return global_best_chromosome, global_min

    def my_breed_population(self, best_chromosome, n, min_route, n_routes, n_breaks, vertices, iteration):
        '''Create a new population based upon the best performing specimens of the old population'''
        n = len(self.turbine_locations)        
        new_pop_route = []
        new_pop_breaks = []
        bestofeach8Chrome = best_chromosome
        # print(best_chromosome)
        for i in range(len(best_chromosome)):
            bestof8route = bestofeach8Chrome[i][0]
            bestof8breaks = bestofeach8Chrome[i][1]
            new_8 = self.transformations(bestof8route,bestof8breaks, n, min_route, n_routes, n_breaks, vertices, iteration)
            new_pop_route += new_8[0]
            new_pop_breaks += new_8[1]
        return new_pop_route, new_pop_breaks
    
    def my_transformations(self, bestofeach8Chrome, n, min_route, n_routes, n_breaks, vertices, iteration):
        '''Returns a copy of the original chromosome and 7 mutations'''
        temp_pop_route = []
        temp_pop_breaks = []
        for threetimes in range(3):
            ### 随机从最好的里面选出两个，按照scikit-opt的crossover算子进行变异
            a,b = random.sample(range(len(bestofeach8Chrome)), 2)
            # print(a,b,bestofeach8Chrome)
            Chrom1, Chrom2 = copy.deepcopy(bestofeach8Chrome[a][0]), copy.deepcopy(bestofeach8Chrome[b][0])
            break1, break2 = copy.deepcopy(bestofeach8Chrome[a][1]), copy.deepcopy(bestofeach8Chrome[b][1])
            # print(a,b,Chrom1, Chrom2)
            cxpoint1, cxpoint2 = np.random.randint(0, len(Chrom1) - 1, 2)
            if cxpoint1 >= cxpoint2:
                cxpoint1, cxpoint2 = cxpoint2, cxpoint1 + 1
            # crossover at the point cxpoint1 to cxpoint2
            pos1_recorder = {value: idx for idx, value in enumerate(Chrom1)}
            pos2_recorder = {value: idx for idx, value in enumerate(Chrom2)}
            for j in range(cxpoint1, cxpoint2):
                value1, value2 = Chrom1[j], Chrom2[j]
                pos1, pos2 = pos1_recorder[value2], pos2_recorder[value1]
                Chrom1[j], Chrom1[pos1] = Chrom1[pos1], Chrom1[j]
                Chrom2[j], Chrom2[pos2] = Chrom2[pos2], Chrom2[j]
                pos1_recorder[value1], pos1_recorder[value2] = pos1, j
                pos2_recorder[value1], pos2_recorder[value2] = j, pos2
            # 第一、三、五条变异后的染色体
            temp_pop_route.append(Chrom1)
            temp_pop_breaks.append(break2)
            # 第二、四、六条变异后的染色体
            temp_pop_route.append(Chrom2)
            temp_pop_breaks.append(break1)
        ###第七条染色体
        theone = np.random.randint(len(bestofeach8Chrome)-1,size=1)
        # print(theone)
        theChrome = bestofeach8Chrome[int(theone)][0]
        theBreak = bestofeach8Chrome[int(theone)][1]
        if iteration== 20:
            Chrom_CW = copy.deepcopy(theChrome)
            Chrom_CW = self.clarke_wright(Chrom_CW, theBreak, vertices)
            temp_pop_route.append(Chrom_CW)
            temp_pop_breaks.append(theBreak)
        else:
            Chrom_Mu = copy.deepcopy(theChrome)
            if np.random.rand() < 0.5:
                selector = random.sample(range(1, n), n-1)###返回1～n范围内的一段长度为n-1的数列
                randlist = [selector[int(n/3)], selector[int(2*n/3)]]
                I = min(randlist)
                J = max(randlist)
                Temp = Chrom_Mu[I:J]; Temp = Temp[::-1]; Chrom_Mu = Chrom_Mu[:I] + Temp + Chrom_Mu[J:]
            temp_pop_route.append(Chrom_Mu)
            temp_pop_breaks.append(self.rand_breaks(n, min_route, n_routes, n_breaks))
       
        return temp_pop_route, temp_pop_breaks
        


        
        
if __name__ == '__main__':
    xmin,ymin,xmax,ymax = 443340, 3322634, 443592, 3322848 
    import time
    for x_space in [50]:
        for pop_size in range(2000,10001,1000):
            time_s = time.time()
            turbine_location = []
            # for x in range(xmin,xmax,x_space):
            #     for y in range(ymin,ymax,x_space):
            #         turbine_location.append([x,y])
            for x in range(xmin+20,xmax-20,x_space):
                for y in range(ymin+20,ymax-20,x_space*2):
                    turbine_location.append([x,y])
            for x in range(xmin+20+int(x_space/2),xmax-20,x_space):
                for y in range(ymin+20+x_space,ymax-20,x_space*2):
                    turbine_location.append([x,y])
            turbine_locations = [x for coord in turbine_location for x in coord]
            landpointlocation = [444000,3323000]

            fig_name = str(len(turbine_location))
            CC = CableCostGA(turbine_locations, substation_location=landpointlocation,show_prog = False, pop_size=pop_size, num_iter = 1000, convergence_definition = 20,show_result = True, figname='test'+str(pop_size))
            # print ('Cable length:',CC.compute_cable_cost())
            cost = CC.compute_cable_cost()
            time_e = time.time()
            print(len(turbine_location),pop_size,cost,time_e-time_s)
    