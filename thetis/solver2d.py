"""
Module for 2D depth averaged solver
"""
from __future__ import absolute_import
from .utility import *
from . import shallowwater_eq
from . import timeintegrator
from . import rungekutta
from . import implicitexplicit
from . import coupled_timeintegrator_2d
from . import tracer_eq_2d
from . import conservative_tracer_eq_2d
from . import sediment_eq_2d
from . import exner_eq
import weakref
import time as time_mod
import numpy as np
from mpi4py import MPI
from . import exporter
from .turbines import TidalTurbineFarm, DiscreteTidalTurbineFarm
from .field_defs import field_metadata
from .options import ModelOptions2d
from . import callback
from .log import *
from collections import OrderedDict
import thetis.limiter as limiter


class FlowSolver2d(FrozenClass):
    """
    Main object for 2D depth averaged solver

    **Example**

    Create mesh

    .. code-block:: python

        from thetis import *
        mesh2d = RectangleMesh(20, 20, 10e3, 10e3)

    Create bathymetry function, set a constant value

    .. code-block:: python

        fs_p1 = FunctionSpace(mesh2d, 'CG', 1)
        bathymetry_2d = Function(fs_p1, name='Bathymetry').assign(10.0)

    Create solver object and set some options

    .. code-block:: python

        solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
        options = solver_obj.options
        options.element_family = 'dg-dg'
        options.polynomial_degree = 1
        options.timestepper_type = 'CrankNicolson'
        options.simulation_export_time = 50.0
        options.simulation_end_time = 3600.
        options.timestep = 25.0

    Assign initial condition for water elevation

    .. code-block:: python

        solver_obj.create_function_spaces()
        init_elev = Function(solver_obj.function_spaces.H_2d)
        coords = SpatialCoordinate(mesh2d)
        init_elev.project(exp(-((coords[0] - 4e3)**2 + (coords[1] - 4.5e3)**2)/2.2e3**2))
        solver_obj.assign_initial_conditions(elev=init_elev)

    Run simulation

    .. code-block:: python

        solver_obj.iterate()

    See the manual for more complex examples.
    """
    def __init__(self, mesh2d, bathymetry_2d, options=None):
        """
        :arg mesh2d: :class:`Mesh` object of the 2D mesh
        :arg bathymetry_2d: Bathymetry of the domain. Bathymetry stands for
            the mean water depth (positive downwards).
        :type bathymetry_2d: :class:`Function`
        :kwarg options: Model options (optional). Model options can also be
            changed directly via the :attr:`.options` class property.
        :type options: :class:`.ModelOptions2d` instance
        """
        self._initialized = False
        self.mesh2d = mesh2d
        self.comm = mesh2d.comm

        # add boundary length info
        bnd_len = compute_boundary_length(self.mesh2d)
        self.mesh2d.boundary_len = bnd_len

        self.dt = None
        """Time step"""

        self.options = ModelOptions2d()
        """
        Dictionary of all options. A :class:`.ModelOptions2d` object.
        """
        if options is not None:
            self.options.update(options)

        # simulation time step bookkeeping
        self.simulation_time = 0
        self.iteration = 0
        self.i_export = 0
        self.next_export_t = self.simulation_time + self.options.simulation_export_time

        self.callbacks = callback.CallbackManager()
        """
        :class:`.CallbackManager` object that stores all callbacks
        """

        self.fields = FieldDict()
        """
        :class:`.FieldDict` that holds all functions needed by the solver
        object
        """

        self.function_spaces = AttrDict()
        """
        :class:`.AttrDict` that holds all function spaces needed by the
        solver object
        """

        self.fields.bathymetry_2d = bathymetry_2d

        self.export_initial_state = True
        """Do export initial state. False if continuing a simulation"""

        self.sediment_model = None
        """set up option for sediment model"""

        self.bnd_functions = {'shallow_water': {}, 'tracer': {}, 'sediment': {}}

        self._isfrozen = True

    def compute_time_step(self, u_scale=Constant(0.0)):
        r"""
        Computes maximum explicit time step from CFL condition.

        .. math :: \Delta t = \frac{\Delta x}{U}

        Assumes velocity scale :math:`U = \sqrt{g H} + U_{scale}` where
        :math:`U_{scale}` is estimated advective velocity.

        :kwarg u_scale: User provided maximum advective velocity scale
        :type u_scale: float or :class:`Constant`
        """
        csize = self.fields.h_elem_size_2d
        bath = self.fields.bathymetry_2d
        fs = bath.function_space()
        bath_pos = Function(fs, name='bathymetry')
        bath_pos.assign(bath)
        min_depth = 0.05
        bath_pos.dat.data[bath_pos.dat.data < min_depth] = min_depth
        test = TestFunction(fs)
        trial = TrialFunction(fs)
        solution = Function(fs)
        g = physical_constants['g_grav']
        u = (sqrt(g * bath_pos) + u_scale)
        a = inner(test, trial) * dx
        l = inner(test, csize / u) * dx
        solve(a == l, solution)
        return solution

    def set_time_step(self, alpha=0.05):
        """
        Sets the model the model time step

        If the time integrator supports automatic time step, and
        :attr:`ModelOptions2d.timestepper_options.use_automatic_timestep` is
        `True`, we compute the maximum time step allowed by the CFL condition.
        Otherwise uses :attr:`ModelOptions2d.timestep`.

        :kwarg float alpha: CFL number scaling factor
        """
        automatic_timestep = (hasattr(self.options.timestepper_options, 'use_automatic_timestep')
                              and self.options.timestepper_options.use_automatic_timestep)
        # TODO revisit math alpha is OBSOLETE
        if automatic_timestep:
            mesh2d_dt = self.compute_time_step(u_scale=self.options.horizontal_velocity_scale)
            dt = self.options.cfl_2d*alpha*float(mesh2d_dt.dat.data.min())
            dt = self.comm.allreduce(dt, op=MPI.MIN)
            self.dt = dt
        else:
            assert self.options.timestep is not None
            assert self.options.timestep > 0.0
            self.dt = self.options.timestep
        if self.comm.rank == 0:
            print_output('dt = {:}'.format(self.dt))
            sys.stdout.flush()

    def set_sipg_parameter(self):
        r"""
        Compute a penalty parameter which ensures stability of the Interior Penalty method
        used for viscosity and diffusivity terms, from Epshteyn et al. 2007
        (http://dx.doi.org/10.1016/j.cam.2006.08.029).

        The scheme is stable if

        ..math::
            \alpha|_K > 3*X*p*(p+1)*\cot(\theta_K),

        for all elements :math:`K`, where

        ..math::
            X = \frac{\max_{x\in K}(\nu(x))}{\min_{x\in K}(\nu(x))},

        :math:`p` the degree, and :math:`\theta_K` is the minimum angle in the element.
        """
        degree = self.function_spaces.U_2d.ufl_element().degree()
        alpha = Constant(5.0*degree*(degree+1) if degree != 0 else 1.5)
        degree_tracer = self.function_spaces.Q_2d.ufl_element().degree()
        alpha_tracer = Constant(5.0*degree_tracer*(degree_tracer+1) if degree_tracer != 0 else 1.5)

        if self.options.use_automatic_sipg_parameter:
            P0 = self.function_spaces.P0_2d
            theta = get_minimum_angles_2d(self.mesh2d)
            min_angle = theta.vector().gather().min()
            print_output("Minimum angle in mesh: {:.2f} degrees".format(np.rad2deg(min_angle)))
            cot_theta = 1.0/tan(theta)

            # Penalty parameter for shallow water
            if not self.options.tracer_only:
                nu = self.options.horizontal_viscosity
                if nu is not None:
                    alpha = alpha*get_sipg_ratio(nu)*cot_theta
                    self.options.sipg_parameter = interpolate(alpha, P0)
                    max_sipg = self.options.sipg_parameter.vector().gather().max()
                    print_output("Maximum SIPG value:        {:.2f}".format(max_sipg))
                else:
                    print_output("Using default SIPG parameter for shallow water equations")

            # Penalty parameter for tracers
            if self.options.solve_tracer or self.options.sediment_model_options.solve_suspended_sediment:
                if self.options.solve_tracer:
                    tracer_kind = 'tracer'
                elif self.options.sediment_model_options.solve_suspended_sediment:
                    tracer_kind = 'sediment'
                nu = self.options.horizontal_diffusivity
                if nu is not None:
                    alpha_tracer = alpha_tracer*get_sipg_ratio(nu)*cot_theta
                    self.options.sipg_parameter_tracer = interpolate(alpha_tracer, P0)
                    max_sipg = self.options.sipg_parameter_tracer.vector().gather().max()
                    print_output("Maximum {} SIPG value: {:.2f}".format(tracer_kind, max_sipg))
                else:
                    print_output("Using default SIPG parameter for {} equation".format(tracer_kind))
        else:
            print_output("Using default SIPG parameters")
            self.options.sipg_parameter.assign(alpha)
            self.options.sipg_parameter_tracer.assign(alpha_tracer)

    def create_function_spaces(self):
        """
        Creates function spaces

        Function spaces are accessible via :attr:`.function_spaces`
        object.
        """
        self._isfrozen = False
        # ----- function spaces: elev in H, uv in U, mixed is W
        self.function_spaces.P0_2d = get_functionspace(self.mesh2d, 'DG', 0, name='P0_2d')
        self.function_spaces.P1_2d = get_functionspace(self.mesh2d, 'CG', 1, name='P1_2d')
        self.function_spaces.P1v_2d = VectorFunctionSpace(self.mesh2d, 'CG', 1, name='P1v_2d')
        self.function_spaces.P1DG_2d = get_functionspace(self.mesh2d, 'DG', 1, name='P1DG_2d')
        self.function_spaces.P1DGv_2d = VectorFunctionSpace(self.mesh2d, 'DG', 1, name='P1DGv_2d')
        # 2D velocity space
        if self.options.element_family == 'rt-dg':
            self.function_spaces.U_2d = get_functionspace(self.mesh2d, 'RT', self.options.polynomial_degree+1, name='U_2d')
            self.function_spaces.H_2d = get_functionspace(self.mesh2d, 'DG', self.options.polynomial_degree, name='H_2d')
        elif self.options.element_family == 'dg-cg':
            self.function_spaces.U_2d = VectorFunctionSpace(self.mesh2d, 'DG', self.options.polynomial_degree, name='U_2d')
            self.function_spaces.H_2d = get_functionspace(self.mesh2d, 'CG', self.options.polynomial_degree+1, name='H_2d')
        elif self.options.element_family == 'dg-dg':
            self.function_spaces.U_2d = VectorFunctionSpace(self.mesh2d, 'DG', self.options.polynomial_degree, name='U_2d')
            self.function_spaces.H_2d = get_functionspace(self.mesh2d, 'DG', self.options.polynomial_degree, name='H_2d')
        else:
            raise Exception('Unsupported finite element family {:}'.format(self.options.element_family))
        self.function_spaces.V_2d = MixedFunctionSpace([self.function_spaces.U_2d, self.function_spaces.H_2d])

        self.function_spaces.Q_2d = get_functionspace(self.mesh2d, 'DG', 1, name='Q_2d')

        self._isfrozen = True

    def create_equations(self):
        """
        Creates shallow water equations
        """
        if not hasattr(self, 'U_2d'):
            self.create_function_spaces()
        self._isfrozen = False
        # ----- fields
        self.fields.solution_2d = Function(self.function_spaces.V_2d, name='solution_2d')
        self.fields.h_elem_size_2d = Function(self.function_spaces.P1_2d)
        get_horizontal_elem_size_2d(self.fields.h_elem_size_2d)
        self.set_sipg_parameter()
        self.depth = DepthExpression(self.fields.bathymetry_2d,
                                     use_nonlinear_equations=self.options.use_nonlinear_equations,
                                     use_wetting_and_drying=self.options.use_wetting_and_drying,
                                     wetting_and_drying_alpha=self.options.wetting_and_drying_alpha)

        # ---- tidal farms, if any
        if len(self.options.tidal_turbine_farms) + len(self.options.discrete_tidal_turbine_farms) > 0:
            self.tidal_farms = []
            p = self.function_spaces.U_2d.ufl_element().degree()
            quad_degree = 2*p + 1
            for subdomain, farm_options in self.options.tidal_turbine_farms.items():
                fdx = dx(subdomain, degree=quad_degree)
                self.tidal_farms.append(TidalTurbineFarm(farm_options.turbine_density,
                                                         fdx, farm_options))
            for subdomain, farm_options in self.options.discrete_tidal_turbine_farms.items():
                fdx = dx(subdomain, degree=farm_options.quadrature_degree)
                self.tidal_farms.append(DiscreteTidalTurbineFarm(self.mesh2d, fdx, farm_options))
        else:
            self.tidal_farms = None

        # ----- Equations
        self.eq_sw = shallowwater_eq.ShallowWaterEquations(
            self.fields.solution_2d.function_space(),
            self.depth,
            self.options,
            tidal_farms=self.tidal_farms
        )
        self.eq_sw.bnd_functions = self.bnd_functions['shallow_water']
        if self.options.solve_tracer:
            self.fields.tracer_2d = Function(self.function_spaces.Q_2d, name='tracer_2d')
            if self.options.use_tracer_conservative_form:
                self.eq_tracer = conservative_tracer_eq_2d.ConservativeTracerEquation2D(
                    self.function_spaces.Q_2d, self.depth,
                    use_lax_friedrichs=self.options.use_lax_friedrichs_tracer,
                    sipg_parameter=self.options.sipg_parameter_tracer)
            else:
                self.eq_tracer = tracer_eq_2d.TracerEquation2D(
                    self.function_spaces.Q_2d, self.depth,
                    use_lax_friedrichs=self.options.use_lax_friedrichs_tracer,
                    sipg_parameter=self.options.sipg_parameter_tracer)
            if self.options.use_limiter_for_tracers and self.options.polynomial_degree > 0:
                self.tracer_limiter = limiter.VertexBasedP1DGLimiter(self.function_spaces.Q_2d)
            else:
                self.tracer_limiter = None

        sediment_options = self.options.sediment_model_options
        if sediment_options.solve_suspended_sediment or sediment_options.solve_exner:
            uv_2d, elev_2d = self.fields.solution_2d.split()
            sediment_model_class = self.options.sediment_model_options.sediment_model_class
            self.sediment_model = sediment_model_class(
                self.options, self.mesh2d, uv_2d, elev_2d, self.depth)
        if sediment_options.solve_suspended_sediment:
            self.fields.sediment_2d = Function(self.function_spaces.Q_2d, name='sediment_2d')
            self.eq_sediment = sediment_eq_2d.SedimentEquation2D(
                self.function_spaces.Q_2d, self.depth, self.sediment_model,
                use_lax_friedrichs=self.options.use_lax_friedrichs_tracer,
                sipg_parameter=self.options.sipg_parameter_tracer,
                conservative=sediment_options.use_sediment_conservative_form)
            if self.options.use_limiter_for_tracers and self.options.polynomial_degree > 0:
                self.tracer_limiter = limiter.VertexBasedP1DGLimiter(self.function_spaces.Q_2d)
            else:
                self.tracer_limiter = None

        if sediment_options.solve_exner:
            self.eq_exner = exner_eq.ExnerEquation(
                self.fields.bathymetry_2d.function_space(), self.depth,
                depth_integrated_sediment=sediment_options.use_sediment_conservative_form, sediment_model=self.sediment_model)

        self._isfrozen = True  # disallow creating new attributes

    def get_swe_timestepper(self, integrator):
        """
        Gets shallow water timestepper object with appropriate parameters
        """
        fields = {
            'linear_drag_coefficient': self.options.linear_drag_coefficient,
            'quadratic_drag_coefficient': self.options.quadratic_drag_coefficient,
            'manning_drag_coefficient': self.options.manning_drag_coefficient,
            'nikuradse_bed_roughness': self.options.nikuradse_bed_roughness,
            'viscosity_h': self.options.horizontal_viscosity,
            'lax_friedrichs_velocity_scaling_factor': self.options.lax_friedrichs_velocity_scaling_factor,
            'coriolis': self.options.coriolis_frequency,
            'wind_stress': self.options.wind_stress,
            'atmospheric_pressure': self.options.atmospheric_pressure,
            'momentum_source': self.options.momentum_source_2d,
            'volume_source': self.options.volume_source_2d,
        }

        args = (self.eq_sw, self.fields.solution_2d, fields, self.dt, )
        kwargs = {'bnd_conditions': self.bnd_functions['shallow_water']}
        if hasattr(self.options.timestepper_options, 'use_semi_implicit_linearization'):
            kwargs['semi_implicit'] = self.options.timestepper_options.use_semi_implicit_linearization
        if hasattr(self.options.timestepper_options, 'implicitness_theta'):
            kwargs['theta'] = self.options.timestepper_options.implicitness_theta
        if self.options.timestepper_type == 'PressureProjectionPicard':
            # TODO: Probably won't work in coupled mode
            u_test = TestFunction(self.function_spaces.U_2d)
            self.eq_mom = shallowwater_eq.ShallowWaterMomentumEquation(
                u_test, self.function_spaces.U_2d, self.function_spaces.H_2d,
                self.depth,
                options=self.options,
                tidal_farms=self.tidal_farms
            )
            self.eq_mom.bnd_functions = self.bnd_functions['shallow_water']
            args = (self.eq_sw, self.eq_mom, self.fields.solution_2d, fields, self.dt, )
            kwargs['solver_parameters'] = self.options.timestepper_options.solver_parameters_pressure
            kwargs['solver_parameters_mom'] = self.options.timestepper_options.solver_parameters_momentum
            kwargs['iterations'] = self.options.timestepper_options.picard_iterations
        elif self.options.timestepper_type == 'SSPIMEX':
            # TODO meaningful solver params
            kwargs['solver_parameters'] = {
                'ksp_type': 'gmres',
                'pc_type': 'fieldsplit',
                'pc_fieldsplit_type': 'multiplicative',
            }
            kwargs['solver_parameters_dirk'] = {
                'ksp_type': 'gmres',
                'pc_type': 'fieldsplit',
                'pc_fieldsplit_type': 'multiplicative',
            }
        else:
            kwargs['solver_parameters'] = self.options.timestepper_options.solver_parameters
        return integrator(*args, **kwargs)

    def get_tracer_timestepper(self, integrator):
        """
        Gets tracer timestepper object with appropriate parameters
        """
        uv, elev = self.fields.solution_2d.split()
        fields = {
            'elev_2d': elev,
            'uv_2d': uv,
            'diffusivity_h': self.options.horizontal_diffusivity,
            'source': self.options.tracer_source_2d,
            'lax_friedrichs_tracer_scaling_factor': self.options.lax_friedrichs_tracer_scaling_factor,
            'tracer_advective_velocity_factor': self.options.tracer_advective_velocity_factor,
        }

        args = (self.eq_tracer, self.fields.tracer_2d, fields, self.dt, )
        kwargs = {
            'bnd_conditions': self.bnd_functions['tracer'],
            'solver_parameters': self.options.timestepper_options.solver_parameters_tracer,
        }
        if hasattr(self.options.timestepper_options, 'use_semi_implicit_linearization'):
            kwargs['semi_implicit'] = self.options.timestepper_options.use_semi_implicit_linearization
        if hasattr(self.options.timestepper_options, 'implicitness_theta'):
            kwargs['theta'] = self.options.timestepper_options.implicitness_theta
        return integrator(*args, **kwargs)

    def get_sediment_timestepper(self, integrator):
        """
        Gets sediment timestepper object with appropriate parameters
        """
        uv, elev = self.fields.solution_2d.split()
        fields = {
            'elev_2d': elev,
            'uv_2d': uv,
            'diffusivity_h': self.options.horizontal_diffusivity,
            'lax_friedrichs_tracer_scaling_factor': self.options.lax_friedrichs_tracer_scaling_factor,
            'tracer_advective_velocity_factor': self.sediment_model.get_advective_velocity_correction_factor(),
        }

        args = (self.eq_sediment, self.fields.sediment_2d, fields, self.dt, )
        kwargs = {
            'bnd_conditions': self.bnd_functions['sediment'],
            'solver_parameters': self.options.timestepper_options.solver_parameters_sediment,
        }
        if hasattr(self.options.timestepper_options, 'use_semi_implicit_linearization'):
            kwargs['semi_implicit'] = self.options.timestepper_options.use_semi_implicit_linearization
        if hasattr(self.options.timestepper_options, 'implicitness_theta'):
            kwargs['theta'] = self.options.timestepper_options.implicitness_theta
        return integrator(*args, **kwargs)

    def get_exner_timestepper(self, integrator):
        """
        Gets exner timestepper object with appropriate parameters
        """
        uv, elev = self.fields.solution_2d.split()
        if not self.options.sediment_model_options.solve_suspended_sediment:
            self.fields.sediment_2d = Function(elev.function_space()).interpolate(Constant(0.0))
        fields = {
            'elev_2d': elev,
            'sediment': self.fields.sediment_2d,
            'morfac': self.options.sediment_model_options.morphological_acceleration_factor,
            'porosity': self.options.sediment_model_options.porosity,
        }

        args = (self.eq_exner, self.fields.bathymetry_2d, fields, self.dt, )
        kwargs = {
            # only pass SWE bcs, used to determine closed boundaries in bedload term
            'bnd_conditions': self.bnd_functions['shallow_water'],
            'solver_parameters': self.options.timestepper_options.solver_parameters_exner,
        }
        if hasattr(self.options.timestepper_options, 'use_semi_implicit_linearization'):
            kwargs['semi_implicit'] = self.options.timestepper_options.use_semi_implicit_linearization
        if hasattr(self.options.timestepper_options, 'implicitness_theta'):
            kwargs['theta'] = self.options.timestepper_options.implicitness_theta
        return integrator(*args, **kwargs)

    def create_timestepper(self):
        """
        Creates time stepper instance
        """
        if not hasattr(self, 'eq_sw'):
            self.create_equations()

        self._isfrozen = False

        if self.options.log_output and not self.options.no_exports:
            logfile = os.path.join(create_directory(self.options.output_directory), 'log')
            filehandler = logging.logging.FileHandler(logfile, mode='w')
            filehandler.setFormatter(logging.logging.Formatter('%(message)s'))
            output_logger.addHandler(filehandler)

        self.set_time_step()

        # ----- Time integrators
        steppers = {
            'SSPRK33': rungekutta.SSPRK33,
            'ForwardEuler': timeintegrator.ForwardEuler,
            'SteadyState': timeintegrator.SteadyState,
            'BackwardEuler': rungekutta.BackwardEulerUForm,
            'DIRK22': rungekutta.DIRK22UForm,
            'DIRK33': rungekutta.DIRK33UForm,
            'CrankNicolson': timeintegrator.CrankNicolson,
            'PressureProjectionPicard': timeintegrator.PressureProjectionPicard,
            'SSPIMEX': implicitexplicit.IMEXLPUM2,
        }
        try:
            assert self.options.timestepper_type in steppers
        except AssertionError:
            raise Exception('Unknown time integrator type: {:s}'.format(self.options.timestepper_type))
        if self.options.solve_tracer:
            try:
                assert self.options.timestepper_type not in ('PressureProjectionPicard', 'SSPIMEX', 'SteadyState')
            except AssertionError:
                raise NotImplementedError("2D tracer model currently only supports SSPRK33, ForwardEuler, BackwardEuler, DIRK22, DIRK33 and CrankNicolson time integrators.")
            self.timestepper = coupled_timeintegrator_2d.CoupledMatchingTimeIntegrator2D(
                weakref.proxy(self), steppers[self.options.timestepper_type],
            )
        elif self.options.sediment_model_options.solve_suspended_sediment or self.options.sediment_model_options.solve_exner:
            try:
                assert self.options.timestepper_type not in ('PressureProjectionPicard', 'SSPIMEX', 'SteadyState')
            except AssertionError:
                raise NotImplementedError("2D sediment model currently only supports SSPRK33, ForwardEuler, BackwardEuler, DIRK22, DIRK33 and CrankNicolson time integrators.")
            self.timestepper = coupled_timeintegrator_2d.CoupledMatchingTimeIntegrator2D(
                weakref.proxy(self), steppers[self.options.timestepper_type],
            )
        else:
            self.timestepper = self.get_swe_timestepper(steppers[self.options.timestepper_type])
        print_output('Using time integrator: {:}'.format(self.timestepper.__class__.__name__))
        self._isfrozen = True  # disallow creating new attributes

    def create_exporters(self):
        """
        Creates file exporters
        """
        if not hasattr(self, 'timestepper'):
            self.create_timestepper()
        self._isfrozen = False
        # correct treatment of the split 2d functions
        uv_2d, elev_2d = self.fields.solution_2d.split()
        self.fields.uv_2d = uv_2d
        self.fields.elev_2d = elev_2d
        self.exporters = OrderedDict()
        if not self.options.no_exports:
            e = exporter.ExportManager(self.options.output_directory,
                                       self.options.fields_to_export,
                                       self.fields,
                                       field_metadata,
                                       export_type='vtk',
                                       verbose=self.options.verbose > 0)
            self.exporters['vtk'] = e
            hdf5_dir = os.path.join(self.options.output_directory, 'hdf5')
            e = exporter.ExportManager(hdf5_dir,
                                       self.options.fields_to_export_hdf5,
                                       self.fields,
                                       field_metadata,
                                       export_type='hdf5',
                                       verbose=self.options.verbose > 0)
            self.exporters['hdf5'] = e

        self._isfrozen = True  # disallow creating new attributes

    def initialize(self):
        """
        Creates function spaces, equations, time stepper and exporters
        """
        if not hasattr(self, 'U_2d'):
            self.create_function_spaces()
        if not hasattr(self, 'eq_sw'):
            self.create_equations()
        if not hasattr(self, 'timestepper'):
            self.create_timestepper()
        if not hasattr(self, 'exporters'):
            self.create_exporters()
        self._initialized = True

    def assign_initial_conditions(self, elev=None, uv=None, tracer=None, sediment=None):
        """
        Assigns initial conditions

        :kwarg elev: Initial condition for water elevation
        :type elev: scalar :class:`Function`, :class:`Constant`, or an expression
        :kwarg uv: Initial condition for depth averaged velocity
        :type uv: vector valued :class:`Function`, :class:`Constant`, or an expression
        """
        if not self._initialized:
            self.initialize()
        uv_2d, elev_2d = self.fields.solution_2d.split()
        if elev is not None:
            elev_2d.project(elev)
        if uv is not None:
            uv_2d.project(uv)
        if tracer is not None and self.options.solve_tracer:
            self.fields.tracer_2d.project(tracer)

        sediment_options = self.options.sediment_model_options
        if self.sediment_model is not None:
            # update sediment model based on initial conditions for uv and elev
            self.sediment_model.update()
        if sediment_options.solve_suspended_sediment:
            if sediment is not None:
                self.fields.sediment_2d.project(sediment)
            else:
                sediment = self.sediment_model.get_equilibrium_tracer()
                if sediment_options.use_sediment_conservative_form:
                    sediment = sediment * self.depth.get_total_depth(elev_2d)
                self.fields.sediment_2d.project(sediment)

        self.timestepper.initialize(self.fields.solution_2d)

    def add_callback(self, callback, eval_interval='export'):
        """
        Adds callback to solver object

        :arg callback: :class:`.DiagnosticCallback` instance
        :kwarg string eval_interval: Determines when callback will be evaluated,
            either 'export' or 'timestep' for evaluating after each export or
            time step.
        """
        self.callbacks.add(callback, eval_interval)

    def export(self):
        """
        Export all fields to disk

        Also evaluates all callbacks set to 'export' interval.
        """
        self.callbacks.evaluate(mode='export')
        for e in self.exporters.values():
            e.export()

    def load_state(self, i_export, outputdir=None, t=None, iteration=None):
        """
        Loads simulation state from hdf5 outputs.

        This replaces :meth:`.assign_initial_conditions` in model initilization.

        This assumes that model setup is kept the same (e.g. time step) and
        all pronostic state variables are exported in hdf5 format. The required
        state variables are: elev_2d, uv_2d

        Currently hdf5 field import only works for the same number of MPI
        processes.

        :arg int i_export: export index to load
        :kwarg string outputdir: (optional) directory where files are read from.
            By default ``options.output_directory``.
        :kwarg float t: simulation time. Overrides the time stamp stored in the
            hdf5 files.
        :kwarg int iteration: Overrides the iteration count in the hdf5 files.
        """
        if not self._initialized:
            self.initialize()
        if outputdir is None:
            outputdir = self.options.output_directory
        # create new ExportManager with desired outputdir
        state_fields = ['uv_2d', 'elev_2d']
        hdf5_dir = os.path.join(outputdir, 'hdf5')
        e = exporter.ExportManager(hdf5_dir,
                                   state_fields,
                                   self.fields,
                                   field_metadata,
                                   export_type='hdf5',
                                   verbose=self.options.verbose > 0)
        e.exporters['uv_2d'].load(i_export, self.fields.uv_2d)
        e.exporters['elev_2d'].load(i_export, self.fields.elev_2d)
        self.assign_initial_conditions()

        # time stepper bookkeeping for export time step
        self.i_export = i_export
        self.next_export_t = self.i_export*self.options.simulation_export_time
        if iteration is None:
            iteration = int(np.ceil(self.next_export_t/self.dt))
        if t is None:
            t = iteration*self.dt
        self.iteration = iteration
        self.simulation_time = t

        # for next export
        self.export_initial_state = outputdir != self.options.output_directory
        if self.export_initial_state:
            offset = 0
        else:
            offset = 1
        self.next_export_t += self.options.simulation_export_time
        for e in self.exporters.values():
            e.set_next_export_ix(self.i_export + offset)

    def print_state(self, cputime):
        """
        Print a summary of the model state on stdout

        :arg float cputime: Measured CPU time
        """
        if self.options.tracer_only:
            norm_q = norm(self.fields.tracer_2d)

            line = ('{iexp:5d} {i:5d} T={t:10.2f} '
                    'tracer norm: {q:10.4f} {cpu:5.2f}')

            print_output(line.format(iexp=self.i_export, i=self.iteration,
                                     t=self.simulation_time, q=norm_q,
                                     cpu=cputime))
        else:
            norm_h = norm(self.fields.solution_2d.split()[1])
            norm_u = norm(self.fields.solution_2d.split()[0])

            line = ('{iexp:5d} {i:5d} T={t:10.2f} '
                    'eta norm: {e:10.4f} u norm: {u:10.4f} {cpu:5.2f}')
            print_output(line.format(iexp=self.i_export, i=self.iteration,
                                     t=self.simulation_time, e=norm_h,
                                     u=norm_u, cpu=cputime))
        sys.stdout.flush()

    def iterate(self, update_forcings=None,
                export_func=None):
        """
        Runs the simulation

        Iterates over the time loop until time ``options.simulation_end_time`` is reached.
        Exports fields to disk on ``options.simulation_export_time`` intervals.

        :kwarg update_forcings: User-defined function that takes simulation
            time as an argument and updates time-dependent boundary conditions
            (if any).
        :kwarg export_func: User-defined function (with no arguments) that will
            be called on every export.
        """
        # TODO I think export function is obsolete as callbacks are in place
        if not self._initialized:
            self.initialize()

        self.options.use_limiter_for_tracers &= self.options.polynomial_degree > 0

        t_epsilon = 1.0e-5
        cputimestamp = time_mod.perf_counter()
        next_export_t = self.simulation_time + self.options.simulation_export_time

        dump_hdf5 = self.options.export_diagnostics and not self.options.no_exports
        if self.options.check_volume_conservation_2d:
            c = callback.VolumeConservation2DCallback(self,
                                                      export_to_hdf5=dump_hdf5,
                                                      append_to_log=True)
            self.add_callback(c)

        if self.options.check_tracer_conservation:
            if self.options.use_tracer_conservative_form:
                c = callback.ConservativeTracerMassConservation2DCallback('tracer_2d',
                                                                          self,
                                                                          export_to_hdf5=dump_hdf5,
                                                                          append_to_log=True)
            else:
                c = callback.TracerMassConservation2DCallback('tracer_2d',
                                                              self,
                                                              export_to_hdf5=dump_hdf5,
                                                              append_to_log=True)
            self.add_callback(c, eval_interval='export')
        if self.options.sediment_model_options.check_sediment_conservation:
            if self.options.sediment_model_options.use_sediment_conservative_form:
                c = callback.ConservativeTracerMassConservation2DCallback('sediment_2d',
                                                                          self,
                                                                          export_to_hdf5=dump_hdf5,
                                                                          append_to_log=True)
            else:
                c = callback.TracerMassConservation2DCallback('sediment_2d',
                                                              self,
                                                              export_to_hdf5=dump_hdf5,
                                                              append_to_log=True)
            self.add_callback(c, eval_interval='export')

        if self.options.check_tracer_overshoot:
            c = callback.TracerOvershootCallBack('tracer_2d',
                                                 self,
                                                 export_to_hdf5=dump_hdf5,
                                                 append_to_log=True)
            self.add_callback(c, eval_interval='export')
        if self.options.sediment_model_options.check_sediment_overshoot:
            c = callback.TracerOvershootCallBack('sediment_2d',
                                                 self,
                                                 export_to_hdf5=dump_hdf5,
                                                 append_to_log=True)
            self.add_callback(c, eval_interval='export')

        # initial export
        self.print_state(0.0)
        if self.export_initial_state:
            self.export()
            if export_func is not None:
                export_func()
            if 'vtk' in self.exporters and isinstance(self.fields.bathymetry_2d, Function):
                self.exporters['vtk'].export_bathymetry(self.fields.bathymetry_2d)

        initial_simulation_time = self.simulation_time
        internal_iteration = 0

        while self.simulation_time <= self.options.simulation_end_time - t_epsilon:
            self.timestepper.advance(self.simulation_time, update_forcings)

            # Move to next time step
            self.iteration += 1
            internal_iteration += 1
            self.simulation_time = initial_simulation_time + internal_iteration*self.dt

            self.callbacks.evaluate(mode='timestep')

            # Write the solution to file
            if self.simulation_time >= next_export_t - t_epsilon:
                self.i_export += 1
                next_export_t += self.options.simulation_export_time

                cputime = time_mod.perf_counter() - cputimestamp
                cputimestamp = time_mod.perf_counter()
                self.print_state(cputime)

                self.export()
                if export_func is not None:
                    export_func()
