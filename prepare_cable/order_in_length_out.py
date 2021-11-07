import numpy
from firedrake import *
from firedrake_adjoint import *



def cablelength(turbine_locations,substation_location,order_w):
    vertices = substation_location + turbine_locations
    order_pair = []
    for i in range(len(order_w)-1):
        if int(order_w[i+1].values()) == 0:
            pass
        else:
            order_pair.append([int(order_w[i].values()),int(order_w[i+1].values())]) 
    #length value matrix for every two vertices
    length = 0
    for i in range(int(len(vertices)/2)):
        for j in range(int(len(vertices)/2)):
            if [i,j] in order_pair:
                dis = sqrt(float((vertices[2*i]-vertices[2*j])**2+(vertices[2*i+1]-vertices[2*j+1])**2))
                length += dis
    return length


def dlength_dx(turbine_locations,substation_location,order_w):
    vertices = substation_location + turbine_locations
    order_pair = []
    for i in range(len(order_w)-1):
        if int(order_w[i+1].values()) == 0:
            pass
        else:
            order_pair.append([int(order_w[i].values()),int(order_w[i+1].values())]) 
    
    t_number = int(len(vertices)/2)
    grad_l = numpy.zeros((int(len(order_pair)),int(t_number*2)))
    row = 0
    for i in range(t_number):
        for j in range(t_number):
            if [i,j] in order_pair:
                the_sqrt = sqrt(float((vertices[2*i]-vertices[2*j])**2+(vertices[2*i+1]-vertices[2*j+1])**2))
                if the_sqrt == 0 :
                    the_sqrt = sqrt(1/2)
                    if i == 0 :
                        grad_l[row, 2*j] = -the_sqrt
                        grad_l[row, 2*j+1] = -the_sqrt
                        row +=1
                    elif j == 0:
                        grad_l[row, 2*i] = the_sqrt
                        grad_l[row, 2*i+1] = the_sqrt
                        row += 1
                    else:
                        grad_l[row, 2*i] = the_sqrt
                        grad_l[row, 2*j] = -the_sqrt
                        grad_l[row, 2*i+1] = the_sqrt
                        grad_l[row, 2*j+1] = -the_sqrt
                        row += 1
                else:
                    if i == 0 :
                        grad_l[row, 2*j] = -(vertices[2*i]-vertices[2*j])/the_sqrt
                        grad_l[row, 2*j+1] = -(vertices[2*i+1]-vertices[2*j+1])/the_sqrt
                        row +=1
                    elif j == 0:
                        grad_l[row, 2*i] = (vertices[2*i]-vertices[2*j])/the_sqrt
                        grad_l[row, 2*i+1] = (vertices[2*i+1]-vertices[2*j+1])/the_sqrt
                        row += 1
                    else:
                        grad_l[row, 2*i] = (vertices[2*i]-vertices[2*j])/the_sqrt
                        grad_l[row, 2*j] = -(vertices[2*i]-vertices[2*j])/the_sqrt
                        grad_l[row, 2*i+1] = (vertices[2*i+1]-vertices[2*j+1])/the_sqrt
                        grad_l[row, 2*j+1] = -(vertices[2*i+1]-vertices[2*j+1])/the_sqrt
                        row += 1

    gradrow = numpy.zeros(int(t_number*2))
    for j in range(int(t_number*2)):
        for i in range(row):
            gradrow[j] += grad_l[i][j]

    return gradrow

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
    ###link order and turbine locations
    order_w1 = [[0, 1, 2, 3, 4, 8, 12, 11,10], [0, 5, 6, 7, 9, 13, 14, 15, 16]]
    turbine_locations1 = [760.        , 160.        , 779.49681974, 182.80135892,
       791.0495674 , 210.48787958, 814.28994971, 229.45986336,
       816.45835488, 160.        , 845.72732647, 166.59408795,
       866.99395457, 189.12785533, 839.60214818, 245.56676241,
       891.72879469, 206.31966837, 911.58308915, 292.29731629,
       883.78474286, 281.00657136, 855.49174295, 271.0273383 ,
       920.17188078, 215.99847063, 950.18210255, 217.5150043 ,
       965.85573043, 243.30627277, 995.75868127, 248.28074774]
    #substation_location is the location on the island, not a location for turbine.
    substation_location1 = [[500,0]]
    order_w = [Constant(i) for j in order_w1 for i in j]
    turbine_locations = [Constant(i) for i in turbine_locations1]# for i in j]
    substation_location = [Constant(i) for j in substation_location1 for i in j]


    # # perturbation for taylor test
    # x_add = 10
    # y_add = 10
    # Delta_x = [x_add,y_add]
    # for x in numpy.arange(443032+20, 443288-20,60):
    #     for y in numpy.arange(3322891+20, 3323091-20, 40):
    #         Delta_x.append(x_add)
    #         Delta_x.append(y_add)
    # Delta_x = numpy.array(Delta_x)
    # epsilons = [0.01 / 2 ** i for i in range(4)]

    #begin taylor test
    Jm = cablelength(turbine_locations,substation_location,order_w)
    print(Jm)
    # print(dlength_dx(turbine_locations,substation_location,order_w))
    # dJdm = numpy.dot(dlength_dx(turbine_locations,substation_location,order_w),Delta_x.T)

    # residuals = []
    # for eps in epsilons:
    #     locations_plus_perturbe1 = [[x+x_add*eps,y+y_add*eps] for x in numpy.arange(443032+20, 443288-20, 60) for y in numpy.arange(3322891+20, 3323091-20, 40)]
    #     locations_plus_perturbe = [Constant(i) for j in locations_plus_perturbe1 for i in j]

    #     Jp = cablelength(locations_plus_perturbe,substation_location,order_w)
    #     res = numpy.abs(Jp-Jm-dJdm*eps)
    #     residuals.append(res)
    # numpy.min(convergence_rates(residuals, epsilons))
    # print(residuals)





