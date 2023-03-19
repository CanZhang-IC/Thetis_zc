from firedrake import *
from thetis import *
from thetis.callback import  DiagnosticCallback
from thetis.optimisation import DiagnosticOptimisationCallback
          
class RMSECallback(DiagnosticCallback):
    """
    The RMSE between the interested velocity and simulated velocity
    """
    name = 'velocity_errors2'
    variable_names = ['RMSEall', 'RMSEaverage', 'RMSE_current']
    def __init__(self,solver_obj,original_dir,E_area_centre_point, E_area_circle, **kwargs):
        """
        :arg solver_obj: a :class:`.FlowSolver2d` object containing the tidal_turbine_farms
        :arg kwargs: see :class:`DiagnosticCallback`
        """
        
        super().__init__(solver_obj,array_dim=1,**kwargs)
        
        self.solver_obj = solver_obj
        self.dir = original_dir
        # self.dx_index = dx_index

        self.uv, elev = split(self.solver_obj.fields.solution_2d)
        self.P1 = VectorFunctionSpace(self.solver_obj.mesh2d, 'DG', 1)
        self.area = Function(FunctionSpace(self.solver_obj.mesh2d,'DG',1)).assign(Constant(0.0))

        self.coord = E_area_centre_point
        self.radius = E_area_circle

        x = SpatialCoordinate(self.solver_obj.mesh2d)
        dx0 = (x[0] - self.coord[0])/self.radius
        dx1 = (x[1] - self.coord[1])/self.radius
        psi_x = conditional(lt(abs(dx0), 1), 1, 0)
        psi_y = conditional(lt(abs(dx1), 1), 1, 0)
        bump = psi_x * psi_y
        # unit_bump_integral = 1.45661 # integral of bump function for radius=1 (copied from OpenTidalFarm who used Wolfram)
        self.area = self.area + bump

        self.RMSEall = 0
        self.RMSEaverage = 0
        self.RMSE_current = [0]
        self.t_initial = self.solver_obj.simulation_time

    def _cal_error(self):
        t = self.solver_obj.simulation_time / self.solver_obj.options.simulation_export_time 
        if t - int(t) == 0 :
            while len(str(t)) < 7:
                t = '0' + str(t)
            chk = DumbCheckpoint(self.dir+'/hdf5/Velocity2d_'+t[:5],mode=FILE_READ)
            original_velocity = Function(self.P1)
            chk.load(original_velocity, name = 'uv_2d')
            chk.close()
            # File('./original_velocity.pvd').write(original_velocity)

            # diff_vx = self.uv[0]-original_velocity[0]
            # diff_vy = self.uv[1]-original_velocity[1]
            # self.RMSE_current.append(assemble(sqrt(diff_vx**2+diff_vy**2)*dx(self.dx_index)))
            # self.RMSE_current.append(assemble(abs(dot(self.uv,self.uv)-dot(original_velocity,original_velocity))*dx(self.dx_index)))

            v_diff = assemble(abs(dot(self.uv,self.uv)-dot(original_velocity,original_velocity))*self.area*dx)

            self.RMSE_current.append(v_diff)

            self.RMSEall = sum(self.RMSE_current)
            if float(t)-self.t_initial/ self.solver_obj.options.simulation_export_time == 0:
                self.RMSEaverage = self.RMSEall
            else:
                self.RMSEaverage = self.RMSEall / (float(t)-self.t_initial/ self.solver_obj.options.simulation_export_time)  
        return self.RMSEall, self.RMSEaverage, self.RMSE_current[-1]
    
    def __call__(self):
        return self._cal_error()

    def message_str(self, RMSEall, RMSEaverage, RMSE_current):
        return 'RMSEall is: {}, RMSEaverage is: {}, RMSE_current is: {}'.format(RMSEall, RMSEaverage, RMSE_current)

class RMSEOptimisationCallback(DiagnosticOptimisationCallback):
    """
    :class:`DiagnosticOptimisationCallback` that evaluates the performance of each tidal turbine farm during an optimisation.

    See the :py:mod:`optimisation` module for more info about the use of OptimisationCallbacks."""
    name = 'RMSE_optimisation2'
    variable_names = ['RMSEall', 'RMSEaverage', 'RMSE_current']

    def __init__(self, solver_obj, RMSECallback, **kwargs):
        """
        :arg solver_obj: a :class:`.FlowSolver2d` object
        :arg turbine_functional_callback: a :class:`.TurbineFunctionalCallback` used in the forward model
        :args kwargs: see :class:`.DiagnosticOptimisationCallback`"""
        self.tfc = RMSECallback
        super().__init__(solver_obj, **kwargs)

    def compute_values(self, *args):
        RMSE_currents = [RMSE_current.block_variable.saved_output for RMSE_current in self.tfc.RMSE_current]
        RMSEall = self.tfc.RMSEall.block_variable.saved_output
        RMSEaverage = self.tfc.RMSEaverage.block_variable.saved_output
        return RMSEall, RMSEaverage, RMSE_currents[-1]

    def message_str(self, RMSEall, RMSEaverage, RMSE_current):
        return 'RMSEall is: {}, RMSEaverage is: {}, RMSE_current is: {}'.format(RMSEall, RMSEaverage, RMSE_current)