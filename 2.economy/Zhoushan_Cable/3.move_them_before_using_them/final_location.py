import h5py
import numpy as np
import matplotlib.pyplot as plt
from thetis import *
from firedrake_adjoint import *

op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare_cable.Hybrid_Code
from prepare_cable.cable_overloaded import cablelength


for BE in range(7):
    
    BE = float(BE)

    file_dir = '../../'
    mesh2d = Mesh(file_dir+'mesh/mesh.msh')
    P1 = FunctionSpace(mesh2d, "CG", 1)

    turbine_density = Function(P1).assign(0.0)
    final_turbine_density = Function(P1, name = 'optimal_density').assign(0.0)

    x = SpatialCoordinate(mesh2d)

    result_output_dir = '../../../outputs/2.economy/discrete/intermediate-1-l/yaw-cable-BE-'+str(float(BE))[:-2]
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        def_file = h5py.File(result_output_dir+'/diagnostic_'+'controls'+'.hdf5','r+')
        for name, data in def_file.items():
            all_controls = list(data[-1])
            iteration_numbers = len(data)
    else:
        all_controls = None
        iteration_numbers = None
    all_controls = comm.bcast(all_controls,root = 0)
    iteration_numbers = comm.bcast(iteration_numbers, root = 0)

    for i in range(12):
        coord = [all_controls[2*i],all_controls[2*i+1]]
        dx0 = (x[0] - coord[0])/10
        dx1 = (x[1] - coord[1])/10
        psi_x = conditional(lt(abs(dx0), 1), exp(1-1/(1-dx0**2)), 0)
        psi_y = conditional(lt(abs(dx1), 1), exp(1-1/(1-dx1**2)), 0)
        bump = psi_x * psi_y
        unit_bump_integral = 1.45661 # integral of bump function for radius=1 (copied from OpenTidalFarm who used Wolfram)

        turbine_density = turbine_density + bump/(10**2 * unit_bump_integral)

    final_turbine_density.interpolate(turbine_density)
    File('turbinedensity.pvd', mode = 'a').write(final_turbine_density)

    turbine_locations = all_controls[:24]#[x for coord in turbine_location for x in coord]
    landpointlocation = [444000,3323000]

    CC = prepare_cable.Hybrid_Code.CableCostGA(turbine_locations, substation_location=landpointlocation)
    print (BE,'Cable length:',CC.compute_cable_cost())

     
