"""
Thetis options for the 2D and 3D model

All options are type-checked and they are stored in traitlets Configurable
objects.
"""
from .configuration import *
from firedrake import Constant
from .sediment_model import SedimentModel


class TimeStepperOptions(FrozenHasTraits):
    """Base class for all time stepper options"""
    name = 'Time stepper'


class ExplicitTimestepperOptions(TimeStepperOptions):
    """Options for explicit time integrator"""
    use_automatic_timestep = Bool(True, help='Set time step automatically based on local CFL conditions.').tag(config=True)


class SemiImplicitTimestepperOptions2d(TimeStepperOptions):
    """Options for 2d explicit time integrator"""
    solver_parameters = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'multiplicative',
    }).tag(config=True)
    solver_parameters_tracer = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'sor',
    }).tag(config=True)
    solver_parameters_sediment = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'sor',
    }).tag(config=True)
    use_semi_implicit_linearization = Bool(
        False, help="Use linearized semi-implicit time integration").tag(config=True)
    solver_parameters_exner = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'sor',
    }).tag(config=True)


class SteadyStateTimestepperOptions2d(TimeStepperOptions):
    """Options for 2d steady state solver"""
    solver_parameters = PETScSolverParameters({
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'mat_type': 'aij'
    }).tag(config=True)


class CrankNicolsonTimestepperOptions2d(SemiImplicitTimestepperOptions2d):
    """Options for 2d Crank-Nicolson time integrator"""
    implicitness_theta = BoundedFloat(
        default_value=0.5, bounds=[0.5, 1.0],
        help='implicitness parameter theta. Value 0.5 implies Crank-Nicolson scheme, 1.0 implies fully implicit formulation.').tag(config=True)


class PressureProjectionTimestepperOptions2d(TimeStepperOptions):
    """Options for 2d pressure-projection time integrator"""
    solver_parameters_pressure = PETScSolverParameters({
        'ksp_type': 'preonly',  # we solve the full schur complement exactly, so no need for outer krylov
        'mat_type': 'matfree',
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        # velocity mass block:
        'fieldsplit_U_2d': {
            'ksp_type': 'gmres',
            'pc_type': 'python',
            'pc_python_type': 'firedrake.AssembledPC',
            'assembled_ksp_type': 'preonly',
            'assembled_pc_type': 'bjacobi',
            'assembled_sub_pc_type': 'ilu',
        },
        # schur system: explicitly assemble the schur system
        # this only works with pressureprojectionicard if the velocity block is just the mass matrix
        # and if the velocity is DG so that this mass matrix can be inverted explicitly
        'fieldsplit_H_2d': {
            'ksp_type': 'preonly',
            'pc_type': 'python',
            'pc_python_type': 'thetis.AssembledSchurPC',
            'schur_ksp_type': 'gmres',
            'schur_ksp_max_it': 100,
            'schur_pc_type': 'gamg',
        },
    }).tag(config=True)
    solver_parameters_momentum = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'sor',
    }).tag(config=True)
    implicitness_theta = BoundedFloat(
        default_value=0.5, bounds=[0.5, 1.0],
        help='implicitness parameter theta. Value 0.5 implies Crank-Nicolson scheme, 1.0 implies fully implicit formulation.').tag(config=True)
    use_semi_implicit_linearization = Bool(
        True, help="Use linearized semi-implicit time integration").tag(config=True)
    picard_iterations = PositiveInteger(
        default_value=2,
        help='number of Picard iterations to converge the nonlinearity in the equations.')


class ExplicitTimestepperOptions2d(ExplicitTimestepperOptions):
    """Options for 2d explicit time integrator"""
    solver_parameters = PETScSolverParameters({
        'snes_type': 'ksponly',
        'ksp_type': 'cg',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'ilu',
        'mat_type': 'aij',
    }).tag(config=True)
    solver_parameters_tracer = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'sor',
    }).tag(config=True)


class ExplicitTimestepperOptions3d(ExplicitTimestepperOptions):
    """Base class for all 3d time stepper options"""
    solver_parameters_2d_swe = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'multiplicative',
    }).tag(config=True)
    solver_parameters_momentum_explicit = PETScSolverParameters({
        'snes_type': 'ksponly',
        'ksp_type': 'cg',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'ilu',
    }).tag(config=True)
    solver_parameters_momentum_implicit = PETScSolverParameters({
        'snes_type': 'ksponly',
        'ksp_type': 'preonly',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'ilu',
    }).tag(config=True)
    solver_parameters_tracer_explicit = PETScSolverParameters({
        'snes_type': 'ksponly',
        'ksp_type': 'cg',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'ilu',
    }).tag(config=True)
    solver_parameters_tracer_implicit = PETScSolverParameters({
        'snes_type': 'ksponly',
        'ksp_type': 'preonly',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'ilu',
    }).tag(config=True)


class SemiImplicitTimestepperOptions3d(ExplicitTimestepperOptions3d):
    """Class for all 3d time steppers that have a configurable semi-implicit 2D solver"""
    implicitness_theta_2d = BoundedFloat(
        default_value=0.5, bounds=[0.5, 1.0],
        help='implicitness parameter theta for 2D solver. Value 1.0 implies fully implicit formulation.').tag(config=True)


class TurbulenceModelOptions(FrozenHasTraits):
    """Abstract base class for all turbulence model options"""
    name = 'Turbulence closure model'


class PacanowskiPhilanderModelOptions(TurbulenceModelOptions):
    """Options for Pacanowski-Philander turbulence model"""
    name = 'Pacanowski-Philander turbulence closure model'
    max_viscosity = PositiveFloat(5e-2, help=r"float: Constant maximum viscosity :math:`\nu_{max}`").tag(config=True)
    alpha = PositiveFloat(10.0, help="float: Richardson number multiplier").tag(config=True)
    exponent = PositiveFloat(2.0, help=r"float: Exponent of viscosity numerator :math:`n`").tag(config=True)


class GLSModelOptions(TurbulenceModelOptions):
    """Options for Generic Length Scale turbulence model"""
    name = 'Generic Length Scale turbulence closure model'
    closure_name = Enum(
        ['k-epsilon', 'k-omega', 'Generic Length Scale'],
        default_value='k-epsilon',
        help='Name of two-equation closure').tag(config=True)
    stability_function_name = Enum(
        ['Canuto A', 'Canuto B', 'Kantha-Clayson', 'Cheng'],
        default_value='Canuto A',
        help='Name of stability function family').tag(config=True)
    p = Float(3.0,
              help='float: parameter p for the definition of psi').tag(config=True)
    m = Float(1.5,
              help='float: parameter m for the definition of psi').tag(config=True)
    n = Float(-1.0,
              help='float: parameter n for the definition of psi').tag(config=True)
    schmidt_nb_tke = PositiveFloat(1.0,
                                   help='float: turbulent kinetic energy Schmidt number').tag(config=True)
    schmidt_nb_psi = PositiveFloat(1.3,
                                   help='float: psi Schmidt number').tag(config=True)
    cmu0 = PositiveFloat(0.5477,
                         help='float: cmu0 parameter').tag(config=True)
    compute_cmu0 = Bool(
        True, help="""bool: compute cmu0 from stability function parameters

        If :attr:`compute_cmu0` is True, this value will be overriden""").tag(config=True)
    c1 = Float(1.44, help='float: c1 parameter for Psi equations').tag(config=True)
    c2 = Float(1.92, help='float: c2 parameter for Psi equations').tag(config=True)
    c3_minus = Float(
        -0.52, help="""float: c3 parameter for Psi equations, stable stratification

        If :attr:`compute_c3_minus` is True this value will be overriden""").tag(config=True)
    c3_plus = Float(1.0,
                    help='float: c3 parameter for Psi equations, unstable stratification').tag(config=True)
    compute_c3_minus = Bool(True,
                            help='bool: compute :attr:`c3_minus` from :attr:`ri_st`').tag(config=True)
    f_wall = Float(1.0, help='float: wall function parameter').tag(config=True)
    ri_st = Float(0.25, help='steady state gradient Richardson number').tag(config=True)
    # FIXME should default to physical_constants['von_karman'] ??
    # FIXME remove the dualism of von Karman constant
    kappa = Float(
        0.4, help="""float: von Karman constant

        If :attr:`compute_kappa` is True this value will be overriden""").tag(config=True)
    compute_kappa = Bool(False,
                         help='bool: compute von Karman constant from :attr:`schmidt_nb_psi`').tag(config=True)
    compute_schmidt_nb_psi = Bool(True,
                                  help='bool: compute psi Schmidt number').tag(config=True)
    k_min = PositiveFloat(1.0e-6,
                          help='float: minimum value for turbulent kinetic energy').tag(config=True)
    psi_min = PositiveFloat(1.0e-14, help='float: minimum value for psi').tag(config=True)
    eps_min = PositiveFloat(1.0e-14, help='float: minimum value for epsilon').tag(config=True)
    len_min = PositiveFloat(1.0e-12,
                            help='float: minimum value for turbulent length scale').tag(config=True)
    compute_galperin_clim = Bool(True,
                                 help='bool: compute c_lim length scale limiting factor').tag(config=True)
    compute_len_min = Bool(False,
                           help='bool: compute min_len from k_min and psi_min').tag(config=True)
    compute_psi_min = Bool(False,
                           help='bool: compute psi_len from k_min and eps_min').tag(config=True)
    visc_min = PositiveFloat(1.0e-8,
                             help='float: minimum value for eddy viscosity').tag(config=True)
    diff_min = PositiveFloat(1.0e-8,
                             help='float: minimum value for eddy diffusivity').tag(config=True)
    galperin_clim = PositiveFloat(0.30,
                                  help='float: Galperin length scale limitation parameter').tag(config=True)

    limit_len = Bool(False, help='bool: apply Galperin length scale limit').tag(config=True)
    limit_psi = Bool(True,
                     help='bool: apply Galperin length scale limit on psi').tag(config=True)
    limit_eps = Bool(False,
                     help='bool: apply Galperin length scale limit on epsilon').tag(config=True)
    limit_len_min = Bool(True,
                         help='bool: limit minimum turbulent length scale to len_min').tag(config=True)

    def apply_defaults(self, closure_name):
        """
        Applies default parameters for given closure name

        :arg closure_name: name of the turbulence closure model
        :type closure_name: string

        Sets default values for parameters p, m, n, schmidt_nb_tke,
        schmidt_nb_psi, c1, c2, c3_plus, c3_minus,
        f_wall, k_min, psi_min
        """

        kepsilon = {'p': 3,
                    'm': 1.5,
                    'n': -1.0,
                    'cmu0': 0.5477,
                    'schmidt_nb_tke': 1.0,
                    'schmidt_nb_psi': 1.3,
                    'c1': 1.44,
                    'c2': 1.92,
                    'c3_plus': 1.0,
                    'c3_minus': -0.52,
                    'f_wall': 1.0,
                    'k_min': 1.0e-6,
                    'eps_min': 1.0e-14,
                    'psi_min': 1.0e-14,
                    'closure_name': 'k-epsilon',
                    }
        # k-epsilon defaults, from tables 1 and 2 in [3]
        komega = {'p': -1.0,
                  'm': 0.5,
                  'n': -1.0,
                  'cmu0': 0.5477,
                  'schmidt_nb_tke': 2.0,
                  'schmidt_nb_psi': 2.0,
                  'c1': 0.555,
                  'c2': 0.833,
                  'c3_plus': 1.0,
                  'c3_minus': -0.52,
                  'f_wall': 1.0,
                  'k_min': 7.6e-6,
                  'eps_min': 1.0e-14,
                  'psi_min': 1.0e-14,
                  'closure_name': 'k-omega',
                  }
        # k-omega defaults, from tables 1 and 2 in [3]
        gen = {'p': 2.0,
               'm': 1.0,
               'n': -0.67,
               'cmu0': 0.5477,
               'schmidt_nb_tke': 0.8,
               'schmidt_nb_psi': 1.07,
               'c1': 1.0,
               'c2': 1.22,
               'c3_plus': 1.0,
               'c3_minus': 0.05,
               'f_wall': 1.0,
               'k_min': 1.0e-6,
               'eps_min': 1.0e-14,
               'psi_min': 1.0e-14,
               'closure_name': 'Generic Length Scale',
               }
        # GLS model A defaults, from tables 1 and 2 in [3]

        if closure_name == 'k-epsilon':
            self.update(kepsilon)
        elif closure_name == 'k-omega':
            self.update(komega)
        elif closure_name == 'Generic Length Scale':
            self.update(gen)
        else:
            raise ValueError('Unknown closure name "{:}"'.format(closure_name))


class EquationOfStateOptions(FrozenHasTraits):
    """Base class of equation of state options"""
    name = 'Equation of State'


class LinearEquationOfStateOptions(EquationOfStateOptions):
    """Linear equation of state options"""
    # TODO more human readable parameter names
    # TODO document the actual equation somewhere
    name = 'Linear Equation of State'
    rho_ref = NonNegativeFloat(1000.0, help='Reference water density').tag(config=True)
    s_ref = NonNegativeFloat(35.0, help='Reference water salinity').tag(config=True)
    th_ref = Float(15.0, help='Reference water temperature').tag(config=True)
    alpha = Float(0.2, help='Thermal expansion coefficient of ocean water').tag(config=True)
    beta = Float(0.77, help='Saline contraction coefficient of ocean water').tag(config=True)


class TidalTurbineOptions(FrozenHasTraits):
    """Tidal turbine parameters"""
    name = 'Tidal turbine options'
    diameter = PositiveFloat(
        18., help='Turbine diameter').tag(config=True)
    C_support = NonNegativeFloat(
        0., help='Thrust coefficient for support structure').tag(config=True)
    A_support = NonNegativeFloat(
        0., help='Cross section of support structure').tag(config=True)


class ConstantTidalTurbineOptions(TidalTurbineOptions):
    """Options for tidal turbine with constant thrust"""
    name = 'Constant tidal turbine options'
    thrust_coefficient = PositiveFloat(
        0.8, help='Thrust coefficient C_T').tag(config=True)


class RatedTidalTurbineOptions(TidalTurbineOptions):
    """Options for tidal turbine with analytical thrust based on rated speed"""
    name = 'Rated tidal turbine options'
    rated_speed = PositiveFloat(
        3.0, help='Rated speed').tag(config=True)
    cut_in_speed = NonNegativeFloat(
        0.0, help='Cut-in speed').tag(config=True)
    # a high default, so it's not normally applied
    cut_out_speed = NonNegativeFloat(
        100.0, help='Cut-out speed').tag(config=True)


class TabulatedTidalTurbineOptions(TidalTurbineOptions):
    """Options for tidal turbine with tabulated thrust coefficient"""
    name = 'Tabulated tidal turbine options'
    thrust_coefficients = List([3.0], help='Table of thrust coefficients')
    thrust_speeds = List(
        [0.8, 0.8],
        help="""List of speeds at which thrust_coefficients are applied.
        First and last entry function as cut-in and cut-out speeds respectively""")


@attach_paired_options("turbine_type",
                       PairedEnum([('constant', ConstantTidalTurbineOptions),
                                   ('rated', RatedTidalTurbineOptions),
                                   ('table', TabulatedTidalTurbineOptions),
                                   ],
                                  "turbine_options",
                                  default_value='constant',
                                  help='Type of turbine thrust specification').tag(config=True),
                       Instance(TidalTurbineOptions, args=()).tag(config=True))
class TidalTurbineFarmOptions(FrozenHasTraits, TraitType):
    """Tidal turbine farm options"""
    name = 'Farm options'
    turbine_density = FiredrakeScalarExpression(
        Constant(0.0), help='Density of turbines within the farm')
    break_even_wattage = NonNegativeFloat(
        0.0, help='Average power production per turbine required to break even')
    considering_b_e_water_depth = Bool(False,
                             help='bool: Considering the water depth effection on the break even value').tag(config=True)
    b_e_water_depth = FiredrakeScalarExpression(
        Constant(0.0), help='Density of turbines within the farm')


class DiscreteTidalTurbineFarmOptions(TidalTurbineFarmOptions):
    """Discrete Tidal turbine farm options - defaults to 0 turbines in the field"""
    name = 'Discrete Farm options'
    turbine_coordinates = List(default_value=[], help="Coordinates for turbines").tag(config=True)
    upwind_correction = Bool(True,
                             help='bool: Apply flow correction to correct for upwind velocity').tag(config=True)
    quadrature_degree = PositiveInteger(10,
                                        help='Quadrature degree for thrust force and power output integral').tag(config=True)
    
    considering_yaw = Bool(True,
                             help='bool: consider the yaw effects for each turbine').tag(config=True)
    turbine_axis = List(default_value=[], help='The direction of turbine axis').tag(config=True)
    # farm_alpha = FiredrakeScalarExpression(Constant(0.0), help='Yaw angle function') 
    considering_individual_thrust_coefficient = Bool(False,
                             help='bool: consider thrust coefficient for each turbine').tag(config=True)
    individual_thrust_coefficient = List(default_value=[], help='Trust coefficient for each turbine').tag(config=True)


class CommonModelOptions(FrozenConfigurable):
    """Options that are common for both 2d and 3d models"""
    name = 'Model options'
    polynomial_degree = NonNegativeInteger(1, help='Polynomial degree of elements').tag(config=True)
    element_family = Enum(
        ['dg-dg', 'rt-dg', 'dg-cg'],
        default_value='dg-dg',
        help="""Finite element family

        2D solver supports 'dg-dg', 'rt-dg', or 'dg-cg' velocity-pressure pairs.
        3D solver supports 'dg-dg', or 'rt-dg' velocity-pressure pairs.""").tag(config=True)

    use_nonlinear_equations = Bool(True, help='Use nonlinear shallow water equations').tag(config=True)
    use_grad_div_viscosity_term = Bool(
        False,
        help=r"""Include :math:`\nabla (\nu_h \nabla \cdot \bar{\textbf{u}})` term in the depth-averaged viscosity

        See :class:`.shallowwater_eq.HorizontalViscosityTerm` for details.""").tag(config=True)
    use_grad_depth_viscosity_term = Bool(
        True,
        help=r"""Include :math:`\nabla H` term in the depth-averaged viscosity

        See :class:`.shallowwater_eq.HorizontalViscosityTerm` for details.""").tag(config=True)

    use_lax_friedrichs_velocity = Bool(
        True, help="use Lax Friedrichs stabilisation in horizontal momentum advection.").tag(config=True)
    lax_friedrichs_velocity_scaling_factor = FiredrakeConstantTraitlet(
        Constant(1.0), help="Scaling factor for Lax Friedrichs stabilisation term in horizontal momentum advection.").tag(config=True)
    use_lax_friedrichs_tracer = Bool(
        False, help="Use Lax Friedrichs stabilisation in tracer advection.").tag(config=True)
    lax_friedrichs_tracer_scaling_factor = FiredrakeConstantTraitlet(
        Constant(1.0), help="Scaling factor for tracer Lax Friedrichs stability term.").tag(config=True)
    use_limiter_for_tracers = Bool(
        True, help="Apply P1DG limiter for tracer fields").tag(config=True)

    check_volume_conservation_2d = Bool(
        False, help="""
        Compute volume of the 2D mode at every export

        2D volume is defined as the integral of the water elevation field.
        Prints deviation from the initial volume to stdout.
        """).tag(config=True)
    log_output = Bool(
        True, help="Redirect all output to log file in output directory").tag(config=True)
    timestep = PositiveFloat(
        10.0, help="Time step").tag(config=True)
    cfl_2d = PositiveFloat(
        1.0, help="Factor to scale the 2d time step OBSOLETE").tag(config=True)  # TODO OBSOLETE
    cfl_3d = PositiveFloat(
        1.0, help="Factor to scale the 2d time step OBSOLETE").tag(config=True)  # TODO OBSOLETE
    simulation_export_time = PositiveFloat(
        100.0, help="""
        Export interval in seconds

        All fields in fields_to_export list will be stored to disk and
        diagnostics will be computed
        """).tag(config=True)
    simulation_end_time = PositiveFloat(
        1000.0, help="Simulation duration in seconds").tag(config=True)
    horizontal_velocity_scale = FiredrakeConstantTraitlet(
        Constant(0.1), help="""
        Maximum horizontal velocity magnitude

        Used to compute max stable advection time step.
        """).tag(config=True)
    horizontal_viscosity_scale = FiredrakeConstantTraitlet(
        Constant(1.0), help="""
        Maximum horizontal viscosity

        Used to compute max stable diffusion time step.
        """).tag(config=True)
    output_directory = Unicode(
        'outputs', help="Directory where model output files are stored").tag(config=True)
    no_exports = Bool(
        False, help="""
        Do not store any outputs to disk

        Disables VTK and HDF5 field outputs. and HDF5 diagnostic outputs.
        Used in CI test suite.
        """).tag(config=True)
    export_diagnostics = Bool(
        True, help="Store diagnostic variables to disk in HDF5 format").tag(config=True)
    fields_to_export = List(
        trait=Unicode(),
        default_value=['elev_2d', 'uv_2d', 'uv_3d', 'w_3d'],
        help="Fields to export in VTK format").tag(config=True)
    fields_to_export_hdf5 = List(
        trait=Unicode(),
        default_value=[],
        help="Fields to export in HDF5 format").tag(config=True)
    verbose = Integer(0, help="Verbosity level").tag(config=True)
    linear_drag_coefficient = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        2D linear drag parameter :math:`L`

        Bottom stress is :math:`\tau_b/\rho_0 = -L \mathbf{u} H`
        """).tag(config=True)
    quadratic_drag_coefficient = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        Dimensionless 2D quadratic drag parameter :math:`C_D`

        Bottom stress is :math:`\tau_b/\rho_0 = -C_D |\mathbf{u}|\mathbf{u}`
        """).tag(config=True)
    manning_drag_coefficient = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        Manning-Strickler 2D quadratic drag parameter :math:`\mu`

        Bottom stress is :math:`\tau_b/\rho_0 = -g \mu^2 |\mathbf{u}|\mathbf{u}/H^{1/3}`
        """).tag(config=True)
    nikuradse_bed_roughness = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        Nikuradse bed roughness length used to construct the 2D quadratic drag parameter :math:`C_D`.

        In sediment transport this term is usually three times the average sediment diameter size.
        """).tag(config=True)
    norm_smoother = FiredrakeConstantTraitlet(
        Constant(0.0), help=r"""
        Coefficient used to avoid non-differentiable functions in the continuous formulation of the velocity norm in
        the quadratic bottom drag term in the momentum equation. This replaces the velocity norm in the quadratic
        bottom drag term with :math:`\|u\| \approx \sqrt{\|u\|^2 + \alpha^2}`
        """).tag(config=True)
    horizontal_viscosity = FiredrakeScalarExpression(
        None, allow_none=True, help="Horizontal viscosity").tag(config=True)
    coriolis_frequency = FiredrakeScalarExpression(
        None, allow_none=True, help="2D Coriolis parameter").tag(config=True)
    wind_stress = FiredrakeVectorExpression(
        None, allow_none=True, help="Stress at free surface (2D vector function)").tag(config=True)
    atmospheric_pressure = FiredrakeScalarExpression(
        None, allow_none=True, help="Atmospheric pressure at free surface, in pascals").tag(config=True)
    momentum_source_2d = FiredrakeVectorExpression(
        None, allow_none=True, help="Source term for 2D momentum equation").tag(config=True)
    volume_source_2d = FiredrakeScalarExpression(
        None, allow_none=True, help="Source term for 2D continuity equation").tag(config=True)
    tracer_source_2d = FiredrakeScalarExpression(
        None, allow_none=True, help="Source term for 2D tracer equation").tag(config=True)
    horizontal_diffusivity = FiredrakeCoefficient(
        None, allow_none=True, help="Horizontal diffusivity for tracers and sediment").tag(config=True)
    use_automatic_sipg_parameter = Bool(False, help=r"""
        Toggle automatic computation of the SIPG penalty parameter used in viscosity and
        diffusivity terms.

        By default, this parameter is set to

        ..math::
            \alpha = 5p(p+1),

        where :math:`p` is the polynomial degree of the velocity space.

        For anisotropic meshes, it is advisable to use the automatic SIPG parameter,
        rather than the default.
        """).tag(config=True)
    sipg_parameter = FiredrakeScalarExpression(
        Constant(10.0), help="Penalty parameter used for horizontal viscosity terms.").tag(config=True)
    sipg_parameter_tracer = FiredrakeScalarExpression(
        Constant(10.0), help="Penalty parameter used for horizontal diffusivity terms.").tag(config=True)


class SedimentModelOptions(FrozenHasTraits):
    c_bstar_constant = NonNegativeFloat(allow_none=False, help='It is 0.015 for Van Rijin').tag(config=True)
    solve_exner = Bool(False, help='Solve exner equation for bed morphology').tag(config=True)
    solve_suspended_sediment = Bool(False, help='Solve suspended sediment transport equation').tag(config=True)
    use_sediment_conservative_form = Bool(False, help='Solve 2D sediment transport in the conservative form').tag(config=True)
    use_bedload = Bool(False, help='Use bedload transport in sediment model').tag(config=True)
    use_angle_correction = Bool(True, help='Switch to use slope effect angle correction').tag(config=True)
    use_slope_mag_correction = Bool(True, help='Switch to use slope effect magnitude correction').tag(config=True)
    use_secondary_current = Bool(False, help='Switch to use secondary current for helical flow effect').tag(config=True)
    average_sediment_size = NonNegativeFloat(allow_none=False, help='Average sediment size').tag(config=True)
    bed_reference_height = NonNegativeFloat(allow_none=False, help='Bottom bed reference height').tag(config=True)
    use_advective_velocity_correction = Bool(True, help="""
        Switch to apply correction to advective velocity used in sediment equation

        Accounts for mismatch between depth-averaged product of velocity with sediment
        and product of depth-averaged velocity with depth-averaged sediment
        """).tag(config=True)
    porosity = FiredrakeCoefficient(
        Constant(0.4), help="Bed porosity for exner equation").tag(config=True)
    morphological_acceleration_factor = FiredrakeConstantTraitlet(
        Constant(1), help="""Rate at which timestep in exner equation is accelerated compared to timestep for model

        timestep in exner = morphological_acceleration_factor * timestep
        """).tag(config=True)
    morphological_viscosity = FiredrakeScalarExpression(
        None, allow_none=True, help="""Viscosity used to derive morphology terms.

        Usually equal to horizontal viscosity but can be set to have a different value""").tag(config=True)
    sediment_density = FiredrakeConstantTraitlet(
        Constant(2650), help='Density of sediment').tag(config=True)
    secondary_current_parameter = FiredrakeConstantTraitlet(
        Constant(0.75), help='Parameter controlling secondary current').tag(config=True)
    slope_effect_parameter = FiredrakeConstantTraitlet(
        Constant(1.3), help='Parameter controlling magnitude of slope effect').tag(config=True)
    slope_effect_angle_parameter = FiredrakeConstantTraitlet(
        Constant(2/3), help='Parameter controlling angle of slope effect').tag(config=True)
    check_sediment_conservation = Bool(
        False, help="""
        Compute total sediment mass at every export

        Prints deviation from the initial mass to stdout.
        """).tag(config=True)
    check_sediment_overshoot = Bool(
        False, help="""
        Compute sediment overshoots at every export

        Prints overshoot values that exceed the initial range to stdout.
        """).tag(config=True)
    sediment_model_class = Type(SedimentModel, help="""Class used to define the sediment model

    This option can be used to provide a user-defined sediment model class that should
    be a subclass of SedimentModel. For example:

    .. code-block:: python

        class UserSedimentModel(SedimentModel):
           def __init__(options, mesh2d, uv, elev, depth, extra_term):
              super().__init__(options, mesh2d, uv, elev, depth)
              self.extra_term = extra_term

           def get_bedloadterm(self, bathymetry):
              return super().get_bedloadterm(bathymetry) + self.term
    """)


# NOTE all parameters are now case sensitive
# TODO rename time stepper types? Allow capitals and spaces?
@attach_paired_options("timestepper_type",
                       PairedEnum([('SSPRK33', ExplicitTimestepperOptions2d),
                                   ('ForwardEuler', ExplicitTimestepperOptions2d),
                                   ('BackwardEuler', SemiImplicitTimestepperOptions2d),
                                   ('CrankNicolson', CrankNicolsonTimestepperOptions2d),
                                   ('DIRK22', SemiImplicitTimestepperOptions2d),
                                   ('DIRK33', SemiImplicitTimestepperOptions2d),
                                   ('SteadyState', SteadyStateTimestepperOptions2d),
                                   ('PressureProjectionPicard', PressureProjectionTimestepperOptions2d),
                                   ('SSPIMEX', SemiImplicitTimestepperOptions2d),
                                   ],
                                  "timestepper_options",
                                  default_value='CrankNicolson',
                                  help='Name of the time integrator').tag(config=True),
                       Instance(TimeStepperOptions, args=()).tag(config=True))
class ModelOptions2d(CommonModelOptions):
    """Options for 2D depth-averaged shallow water model"""
    name = 'Depth-averaged 2D model'
    sediment_model_options = Instance(SedimentModelOptions, args=()).tag(config=True)
    solve_tracer = Bool(False, help='Solve tracer transport').tag(config=True)
    use_tracer_conservative_form = Bool(False, help='Solve 2D tracer transport in the conservative form').tag(config=True)
    use_wetting_and_drying = Bool(
        False, help=r"""bool: Turn on wetting and drying

        Uses the wetting and drying scheme from Karna et al (2011).
        If ``True``, one should also set :attr:`wetting_and_drying_alpha` to control the bathymetry displacement.
        """).tag(config=True)
    wetting_and_drying_alpha = FiredrakeScalarExpression(
        Constant(0.5), help=r"""
        Coefficient: Wetting and drying parameter :math:`\alpha`.

        Used in bathymetry displacement function that ensures positive water depths. Unit is meters.
        """).tag(config=True)
    tidal_turbine_farms = Dict(trait=TidalTurbineFarmOptions(),
                               default_value={}, help='Dictionary mapping subdomain ids to the options of the corresponding farm')

    discrete_tidal_turbine_farms = Dict(trait=DiscreteTidalTurbineFarmOptions(),
                                        default_value={},
                                        help='Dictionary mapping subdomain ids to the options of the corresponding farm')

    check_tracer_conservation = Bool(
        False, help="""
        Compute total tracer mass at every export

        Prints deviation from the initial mass to stdout.
        """).tag(config=True)
    tracer_advective_velocity_factor = FiredrakeScalarExpression(
        Constant(1.0), help="""
        Custom factor multiplied to the velocity variable in tracer advection equation.

        Used to account for mismatch between depth-averaged product of velocity with tracer
        and product of depth-averaged velocity with depth-averaged tracer
        """).tag(config=True)
    check_tracer_overshoot = Bool(
        False, help="""
        Compute tracer overshoots at every export

        Prints overshoot values that exceed the initial range to stdout.
        """).tag(config=True)
    tracer_only = Bool(
        False, help="""Hold shallow water variables in initial state

        Advects tracer in the associated (constant) velocity field.
        """).tag(config=True)


@attach_paired_options("timestepper_type",
                       PairedEnum([('LeapFrog', ExplicitTimestepperOptions3d),
                                   ('SSPRK22', ExplicitTimestepperOptions3d),
                                   ],
                                  "timestepper_options",
                                  default_value='SSPRK22',
                                  help='Name of the time integrator').tag(config=True),
                       Instance(TimeStepperOptions, args=()).tag(config=True))
@attach_paired_options("turbulence_model_type",
                       PairedEnum([('gls', GLSModelOptions),
                                   ('pacanowski', PacanowskiPhilanderModelOptions)
                                   ],
                                  "turbulence_model_options",
                                  default_value='gls',
                                  help='Type of vertical turbulence model').tag(config=True),
                       Instance(TurbulenceModelOptions, args=()).tag(config=True))
@attach_paired_options("equation_of_state_type",
                       PairedEnum([('full', EquationOfStateOptions),
                                   ('linear', LinearEquationOfStateOptions)],
                                  "equation_of_state_options",
                                  default_value='full',
                                  help='Type of equation of state').tag(config=True),
                       Instance(EquationOfStateOptions, args=()).tag(config=True))
class ModelOptions3d(CommonModelOptions):
    """Options for 3D hydrostatic model"""
    name = '3D hydrostatic model'
    solve_salinity = Bool(True, help='Solve salinity transport').tag(config=True)
    solve_temperature = Bool(True, help='Solve temperature transport').tag(config=True)
    use_implicit_vertical_diffusion = Bool(True, help='Solve vertical diffusion and viscosity implicitly').tag(config=True)
    use_bottom_friction = Bool(True, help='Apply log layer bottom stress in the 3D model').tag(config=True)
    use_ale_moving_mesh = Bool(
        True, help="Use ALE formulation where 3D mesh tracks free surface").tag(config=True)
    use_baroclinic_formulation = Bool(
        False, help="Compute internal pressure gradient in momentum equation").tag(config=True)
    use_turbulence = Bool(
        False, help="Activate turbulence model in the 3D model").tag(config=True)
    use_turbulence_advection = Bool(
        False, help="Advect TKE and Psi in the GLS turbulence model").tag(config=True)
    use_smagorinsky_viscosity = Bool(
        False, help="Use Smagorinsky horisontal viscosity parametrization").tag(config=True)
    smagorinsky_coefficient = FiredrakeConstantTraitlet(
        Constant(0.1),
        help="""Smagorinsky viscosity coefficient :math:`C_S`

        See :class:`.SmagorinskyViscosity`.""").tag(config=True)

    use_limiter_for_velocity = Bool(
        True, help="Apply P1DG limiter for 3D horizontal velocity field").tag(config=True)
    check_volume_conservation_3d = Bool(
        False, help="""
        Compute volume of the 3D domain at every export

        Prints deviation from the initial volume to stdout.
        """).tag(config=True)
    check_salinity_conservation = Bool(
        False, help="""
        Compute total salinity mass at every export

        Prints deviation from the initial mass to stdout.
        """).tag(config=True)
    check_salinity_overshoot = Bool(
        False, help="""
        Compute salinity overshoots at every export

        Prints overshoot values that exceed the initial range to stdout.
        """).tag(config=True)
    check_temperature_conservation = Bool(
        False, help="""
        Compute total temperature mass at every export

        Prints deviation from the initial mass to stdout.
        """).tag(config=True)
    check_temperature_overshoot = Bool(
        False, help="""
        Compute temperature overshoots at every export

        Prints overshoot values that exceed the initial range to stdout.
        """).tag(config=True)
    timestep_2d = PositiveFloat(
        10.0, help="""
        Time step of the 2d mode

        This option is only used in the 3d solver, if 2d mode is solved
        explicitly.
        """).tag(config=True)
    vertical_velocity_scale = FiredrakeConstantTraitlet(
        Constant(1e-4), help="""
        Maximum vertical velocity magnitude

        Used to compute max stable advection time step.
        """).tag(config=True)
    use_quadratic_pressure = Bool(
        False, help="""
        Use P2DGxP2 space for baroclinic head.

        If element_family='dg-dg', P2DGxP1DG space is also used for the internal
        pressure gradient.

        This is useful to alleviate bathymetry-induced pressure gradient errors.
        If False, the baroclinic head is in the tracer space, and internal
        pressure gradient is in the velocity space.
        """).tag(config=True)
    use_quadratic_density = Bool(
        False, help="""
        Water density is projected to P2DGxP2 space.

        This reduces pressure gradient errors associated with nonlinear
        equation of state.
        If False, density is computed point-wise in the tracer space.
        """).tag(config=True)
    bottom_roughness = FiredrakeScalarExpression(
        None, allow_none=True, help="Bottom roughness length in meters.").tag(config=True)
    horizontal_diffusivity = FiredrakeScalarExpression(
        None, allow_none=True, help="Horizontal diffusivity for tracers").tag(config=True)
    vertical_diffusivity = FiredrakeScalarExpression(
        None, allow_none=True, help="Vertical diffusivity for tracers").tag(config=True)
    vertical_viscosity = FiredrakeScalarExpression(
        None, allow_none=True, help="Vertical viscosity").tag(config=True)
    momentum_source_3d = FiredrakeVectorExpression(
        None, allow_none=True, help="Source term for 3D momentum equation").tag(config=True)
    salinity_source_3d = FiredrakeScalarExpression(
        None, allow_none=True, help="Source term for salinity equation").tag(config=True)
    temperature_source_3d = FiredrakeScalarExpression(
        None, allow_none=True, help="Source term for temperature equation").tag(config=True)
    constant_temperature = FiredrakeConstantTraitlet(
        Constant(10.0), help="Constant temperature if temperature is not solved").tag(config=True)
    constant_salinity = FiredrakeConstantTraitlet(
        Constant(0.0), help="Constant salinity if salinity is not solved").tag(config=True)
    sipg_parameter_vertical = FiredrakeScalarExpression(
        Constant(10.0), help="Penalty parameter used for vertical viscosity terms.").tag(config=True)
    sipg_parameter_vertical_tracer = FiredrakeScalarExpression(
        Constant(10.0), help="Penalty parameter used for vertical diffusivity terms.").tag(config=True)
    sipg_parameter_turb = FiredrakeScalarExpression(
        Constant(1.5), help="Penalty parameter used for horizontal diffusivity terms of the turbulence model.").tag(config=True)
    sipg_parameter_vertical_turb = FiredrakeScalarExpression(
        Constant(1.0), help="Penalty parameter used for vertical diffusivity terms of the turbulence model.").tag(config=True)
