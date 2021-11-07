from firedrake import *
from thetis import *
from thetis.callback import  DiagnosticCallback
from thetis.optimisation import DiagnosticOptimisationCallback
          
class BathymetryCallback(DiagnosticCallback):
    """
    The Bathymetry changes between the initial value and affected value
    """
    name = 'bathymetry changes'
    variable_names = ['b_c']
    def __init__(self,solver_obj,initial_bathy,E_area_centre_point, E_area_circle, **kwargs):
        """
        :arg solver_obj: a :class:`.FlowSolver2d` object containing the tidal_turbine_farms
        :arg kwargs: see :class:`DiagnosticCallback`
        """
        
        super().__init__(solver_obj,array_dim=1,**kwargs)
        
        self.solver_obj = solver_obj
        self.initial_bathy = initial_bathy

        self.uv, elev = split(self.solver_obj.fields.solution_2d)
        self.bathy = self.solver_obj.fields.bathymetry_2d


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

        self.b_c_all = 0
        self.b_c_average = 0
        self.b_c_current = [0]
        self.t_initial = self.solver_obj.simulation_time

    def _cal_changes(self):
        t = self.solver_obj.simulation_time / self.solver_obj.options.simulation_export_time 
        
        b_c = assemble(abs(self.bathy - self.initial_bathy)*self.area*dx)

        self.b_c_current.append(b_c)
        self.b_c_all = sum(self.b_c_current)
        self.b_c_average = self.b_c_all / (float(t)-self.t_initial/ self.solver_obj.options.simulation_export_time) 
        return self.b_c_all, self.b_c_average, self.b_c_current[-1]
    
    def __call__(self):
        return self._cal_changes()

    def message_str(self, b_c_all, b_c_average, b_c_current):
        return 'b_c_all is: {}, b_c_average is: {}, b_c_current is: {}'.format(b_c_all, b_c_average, b_c_current)

class b_c_OptimisationCallback(DiagnosticOptimisationCallback):
    """
    :class:`DiagnosticOptimisationCallback` that evaluates the performance of each tidal turbine farm during an optimisation.

    See the :py:mod:`optimisation` module for more info about the use of OptimisationCallbacks."""
    name = 'b_c__optimisation'
    variable_names = ['b_c_all', 'b_c_average', 'b_c_current']

    def __init__(self, solver_obj, b_c_Callback, **kwargs):
        """
        :arg solver_obj: a :class:`.FlowSolver2d` object
        :arg turbine_functional_callback: a :class:`.TurbineFunctionalCallback` used in the forward model
        :args kwargs: see :class:`.DiagnosticOptimisationCallback`"""
        self.tfc = b_c_Callback
        super().__init__(solver_obj, **kwargs)

    def compute_values(self, *args):
        b_c_currents = [b_c_current.block_variable.saved_output for b_c_current in self.tfc.b_c_current]
        b_c_all = self.tfc.b_c_all.block_variable.saved_output
        b_c_average = self.tfc.b_c_average.block_variable.saved_output
        return b_c_all, b_c_average, b_c_currents[-1]

    def message_str(self, b_c_all, b_c_average, b_c_current):
        return 'b_c_all is: {}, b_c_average is: {}, b_c_current is: {}'.format(b_c_all, b_c_average, b_c_current)