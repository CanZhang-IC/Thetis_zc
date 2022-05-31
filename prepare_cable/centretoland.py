from matplotlib import pyplot as plt
import math
import copy
import numpy as np
from ad import adnumber



def two_length(turbine_locations):
    '''
    length1: sum of each point to their centre point.
    length2: distance between centre point to land point
    return: length1/turbine_number+length2
    '''
    turbine_number = int(len(turbine_locations)/2)
    x,y = [], []
    for i in range(turbine_number):
        x.append(turbine_locations[2*i])
        y.append(turbine_locations[2*i+1])
    x_centre = sum(x)/turbine_number
    y_centre = sum(y)/turbine_number

    length1 = 0
    for i in range(turbine_number):
        dis = float((x_centre-turbine_locations[2*i])**2+(y_centre-turbine_locations[2*i+1])**2)
        length1 += dis

    total_dist = np.sqrt(length1/turbine_number) 
    return total_dist

def convert_to_adnumber(coordinate_list):
        '''Convert the location vectors from floats into adnumbers to enable differentiation'''
        adnumber_coordinate_list = []
        for i in range(len(coordinate_list)):
            adnumber_coordinate_list.append(adnumber(float(coordinate_list[i].values())))
        coordinate_list = adnumber_coordinate_list
        return coordinate_list
        
    
def dlength_dx(turbine_locations):
    '''Differentiate the length of the routing w.r.t. the position of the turbines, produces a n x 2 array, ((dC/dx1, dC/dy1), ...)'''

    turbine_locations = convert_to_adnumber(turbine_locations)

    turbine_number = int(len(turbine_locations)/2)
    x,y = [], []
    for i in range(turbine_number):
        x.append(turbine_locations[2*i])
        y.append(turbine_locations[2*i+1])
    x_centre = sum(x)/turbine_number
    y_centre = sum(y)/turbine_number

    length1 = 0
    for i in range(turbine_number):
        dis = (x_centre-turbine_locations[2*i])**2+(y_centre-turbine_locations[2*i+1])**2
        length1 += dis

    total_dist = np.sqrt(length1/turbine_number) 
    dC_dX = np.zeros(turbine_number*2) 
    for i in range(turbine_number):
        dC_dX[2*i] = total_dist.d(turbine_locations[2*i])
        dC_dX[2*i+1] = total_dist.d(turbine_locations[2*i+1])
    return dC_dX

def convergence_rates(E_values, eps_values, show=True):
    from numpy import log
    r = []
    for i in range(1, len(eps_values)):
        r.append(log(E_values[i] / E_values[i - 1])
                 / log(eps_values[i] / eps_values[i - 1]))
    if show:
        print("Computed convergence rates: {}".format(r))
    return r

if __name__ == "__main__":
    ###turbine locations

    turbine_locations1 = [[x,y] for x in np.arange(443032+20, 443288-20, 60) for y in np.arange(3322891+20, 3323091-20, 40)]
    #substation_location is the location on the island, not a location for turbine.
    substation_location1 = [[442500,3322750]]

    turbine_locations = [i for j in turbine_locations1 for i in j]
    substation_location = [i for j in substation_location1 for i in j]


    # perturbation for taylor test
    x_add = 10
    y_add = 10
    Delta_x = []
    for x in np.arange(443032+20, 443288-20,60):
        for y in np.arange(3322891+20, 3323091-20, 40):
            Delta_x.append(x_add)
            Delta_x.append(y_add)
    Delta_x = np.array(Delta_x)
    epsilons = [0.01 / 2 ** i for i in range(4)]

    #begin taylor test
    Jm = two_length(turbine_locations,substation_location)
    dJdm = np.dot(dlength_dx(turbine_locations,substation_location),Delta_x.T)

    residuals = []
    for eps in epsilons:
        locations_plus_perturbe1 = [[x+x_add*eps,y+y_add*eps] for x in np.arange(443032+20, 443288-20, 60) for y in np.arange(3322891+20, 3323091-20, 40)]
        locations_plus_perturbe = [i for j in locations_plus_perturbe1 for i in j]

        Jp = two_length(locations_plus_perturbe,substation_location)
        res = np.abs(Jp-Jm-dJdm*eps)
        print(Jm,dJdm)
        residuals.append(res)

    np.min(convergence_rates(residuals, epsilons))
    print(residuals)




# def dlength_dx(turbine_locations,substation_location):
#     turbine_number = int(len(turbine_locations)/2)
#     x,y = [], []
#     for i in range(turbine_number):
#         x.append(turbine_locations[2*i])
#         y.append(turbine_locations[2*i+1])
#     x_centre = sum(x)/turbine_number
#     y_centre = sum(y)/turbine_number

#     grad_l = np.zeros(int(turbine_number*2))
#     # for i in range(turbine_number):
#     #     grx = 0
#     #     gry = 0
#     #     for j in range(turbine_number):
#     #         the_sqrt_1 = sqrt(float((x_centre-turbine_locations[2*j])**2+(y_centre-turbine_locations[2*j+1])**2))
#     #         the_sqrt_2 = sqrt(float((x_centre-substation_location[0])**2+(y_centre-substation_location[1])**2))
#     #         if i == j:
#     #             grx += (1-turbine_number)/turbine_number*(x_centre-turbine_locations[2*j])/the_sqrt_1#/turbine_number #+(x_centre-substation_location[0])/turbine_number/the_sqrt_2
#     #             gry += (1-turbine_number)/turbine_number*(y_centre-turbine_locations[2*j+1])/the_sqrt_1#/turbine_number #+(y_centre-substation_location[1])/turbine_number/the_sqrt_2
#     #         else:
#     #             grx += 1/turbine_number*(x_centre-turbine_locations[2*j])/the_sqrt_1#/turbine_number #+(x_centre-substation_location[0])/turbine_number/the_sqrt_2
#     #             gry += 1/turbine_number*(y_centre-turbine_locations[2*j+1])/the_sqrt_1#/turbine_number #+(y_centre-substation_location[1])/turbine_number/the_sqrt_2
#     #     grad_l[2*i] = grx
#     #     grad_l[2*i+1] = gry
#     for i in range(turbine_number):
#         grx = 0
#         gry = 0
#         for j in range(turbine_number):
#             the_sqrt_1 = sqrt(float((x_centre-turbine_locations[2*j])**2+(y_centre-turbine_locations[2*j+1])**2))
#             the_sqrt_2 = sqrt(float((x_centre-substation_location[0])**2+(y_centre-substation_location[1])**2))
#             if i == j:
#                 grx += (x_centre-substation_location[0])/turbine_number/the_sqrt_2
#                 gry += (y_centre-substation_location[1])/turbine_number/the_sqrt_2
#             else:
#                 grx += (x_centre-substation_location[0])/turbine_number/the_sqrt_2
#                 gry += (y_centre-substation_location[1])/turbine_number/the_sqrt_2
#         grad_l[2*i] = grx
#         grad_l[2*i+1] = gry
#     return grad_l


