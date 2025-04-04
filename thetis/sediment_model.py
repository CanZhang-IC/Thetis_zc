from click import option
from .utility import *
from .log import warning


class CorrectiveVelocityFactor:
    def __init__(self, depth, ksp, bed_reference_height, settling_velocity, ustar):
        """
        Set up advective velocity factor `self.velocity_correction_factor`
        which accounts for mismatch between depth-averaged product of
        velocity with sediment and product of depth-averaged velocity
        with depth-averaged sediment.

        :arg depth: Depth of fluid
        :type depth: :class:`Function`
        :arg ksp: Grain roughness coefficient
        :type ksp: :class:`Constant`
        :arg bed_reference_height: Bottom bed reference height
        :type bed_reference_height: :class:`Constant`
        :arg settling_velocity: Settling velocity of the sediment particles
        :type settling_velocity: :class:`Constant`
        :arg ustar: Shear velocity
        :type ustar: :class:`Expression of functions`
        """

        a =Constant(bed_reference_height/2)

        kappa = physical_constants['von_karman']

        # correction factor to advection velocity in sediment concentration equation
        Bconv = conditional(depth > Constant(1.1)*ksp, ksp/depth, Constant(1/1.1))
        Aconv = conditional(depth > Constant(1.1)*a, a/depth, Constant(1/1.1))

        # take max of value calculated either by ksp or depth
        Amax = conditional(Aconv > Bconv, Aconv, Bconv)

        r1conv = Constant(1) - (1/kappa)*conditional(settling_velocity/ustar < Constant(1),
                                                     settling_velocity/ustar, Constant(1))

        Ione = conditional(r1conv > Constant(1e-8), (Constant(1) - Amax**r1conv)/r1conv,
                           conditional(r1conv < Constant(- 1e-8), (Constant(1) - Amax**r1conv)/r1conv, ln(Amax)))

        Itwo = conditional(r1conv > Constant(1e-8), -(Ione + (ln(Amax)*(Amax**r1conv)))/r1conv,
                           conditional(r1conv < Constant(- 1e-8), -(Ione + (ln(Amax)*(Amax**r1conv)))/r1conv,
                                       Constant(-0.5)*ln(Amax)**2))

        self.alpha = -(Itwo - (ln(Amax) - ln(30))*Ione)/(Ione * ((ln(Amax) - ln(30)) + Constant(1)))

        # final correction factor
        self.velocity_correction_factor = Function(depth.function_space(), name='velocity correction factor')
        self.velocity_correction_factor_expr = max_value(min_value(self.alpha, Constant(1)), Constant(0.))
        self.update()

    def update(self):
        """
        Update `self.velocity_correction_factor` using the updated values for velocity
        """
        # final correction factor
        self.velocity_correction_factor.interpolate(self.velocity_correction_factor_expr)


class SedimentModel(object):
    def __init__(self, options, mesh2d, uv, elev, depth):

        """
        Set up a full morphological model simulation based on provided velocity and elevation functions.

        :arg options: Model options.
        :type options: :class:`.ModelOptions2d` instance
        :arg mesh2d: :class:`Mesh` object of the 2D mesh
        :arg uv: the velocity solution during the simulation
        :type uv: :class:`Function`
        :arg elev: the elevation solution during the simulation
        :type elev: :class:`Function`
        :arg depth: a :class:`DepthExpression` instance to evaluate the current depth

        The :class:`.SedimentModel` provides various expressions to be used in a suspended sediment and/or
        the Exner equation. NOTE that the functions used in these expressions need to be updated
        with the current values of the uv and elev fields by calling :func:`.SedimentModel.update`. This is not done
        in the initialisation of the sediment model, so that the :class:`.SedimentModel` can be created before
        initial conditions have been assigned to uv and elev. After the initial conditions have been
        assigned a call to update is required to ensure that the initial values are reflected in
        the :class:`.SedimentModel` terms.
        """

        self.depth = depth
        self.options = options
        self.solve_suspended_sediment = options.sediment_model_options.solve_suspended_sediment
        self.use_bedload = options.sediment_model_options.use_bedload
        self.use_angle_correction = options.sediment_model_options.use_angle_correction
        self.use_slope_mag_correction = options.sediment_model_options.use_slope_mag_correction
        self.use_advective_velocity_correction = options.sediment_model_options.use_advective_velocity_correction
        self.use_secondary_current = options.sediment_model_options.use_secondary_current

        if not self.use_bedload:
            if self.use_angle_correction:
                warning('Slope effect angle correction only applies to bedload transport which is not used in this simulation')
            if self.use_slope_mag_correction:
                warning('Slope effect magnitude correction only applies to bedload transport which is not used in this simulation')
            if self.use_secondary_current:
                warning('Secondary current only applies to bedload transport which is not used in this simulation')

        self.average_size = options.sediment_model_options.average_sediment_size
        self.bed_reference_height = options.sediment_model_options.bed_reference_height
        self.rhos = options.sediment_model_options.sediment_density
        self.c_bstar_constant = options.sediment_model_options.c_bstar_constant

        # define function spaces
        self.P1DG_2d = get_functionspace(mesh2d, "DG", 1)
        self.P1_2d = get_functionspace(mesh2d, "CG", 1)
        self.P1v_2d = VectorFunctionSpace(mesh2d, "CG", 1)

        # define parameters
        self.g = physical_constants['g_grav']
        self.rhow = physical_constants['rho0']
        kappa = physical_constants['von_karman']

        ksp = Constant(3*self.average_size)
        self.a = Constant(self.bed_reference_height/2)
        if self.options.sediment_model_options.morphological_viscosity is None:
            self.viscosity = self.options.horizontal_viscosity
        else:
            self.viscosity = self.options.sediment_model_options.morphological_viscosity

        # magnitude slope effect parameter
        self.beta = self.options.sediment_model_options.slope_effect_parameter
        # angle correction slope effect parameters
        self.surbeta2 = self.options.sediment_model_options.slope_effect_angle_parameter
        # secondary current parameter
        self.alpha_secc = self.options.sediment_model_options.secondary_current_parameter

        # calculate critical shields parameter thetacr
        self.R = Constant(self.rhos/self.rhow - 1)

        self.dstar = Constant(self.average_size*((self.g*self.R)/(self.viscosity**2))**(1/3))
        if max(self.dstar.dat.data[:] < 1):
            raise ValueError('dstar value less than 1')
        elif max(self.dstar.dat.data[:] < 4):
            self.thetacr = Constant(0.24*(self.dstar**(-1)))
        elif max(self.dstar.dat.data[:] < 10):
            self.thetacr = Constant(0.14*(self.dstar**(-0.64)))
        elif max(self.dstar.dat.data[:] < 20):
            self.thetacr = Constant(0.04*(self.dstar**(-0.1)))
        elif max(self.dstar.dat.data[:] < 150):
            self.thetacr = Constant(0.013*(self.dstar**(0.29)))
        else:
            self.thetacr = Constant(0.055)

        # critical bed shear stress
        self.taucr = Constant((self.rhos-self.rhow)*self.g*self.average_size*self.thetacr)

        # calculate settling velocity
        if self.average_size <= 1e-04:
            self.settling_velocity = Constant(self.g*(self.average_size**2)*self.R/(18*self.viscosity))
        elif self.average_size <= 1e-03:
            self.settling_velocity = Constant((10*self.viscosity/self.average_size)
                                              * (sqrt(1 + 0.01*((self.R*self.g*(self.average_size**3))
                                                 / (self.viscosity**2)))-1))
        else:
            self.settling_velocity = Constant(1.1*sqrt(self.g*self.average_size*self.R))

        # dictionary of steps (interpolate/project) to perform in order in update()
        self.update_steps = OrderedDict()
        # fields updated in these steps
        self.fields = AttrDict()

        # first step: project velocity to CG
        self.fields.uv_cg = Function(self.P1v_2d)
        self.update_steps['uv_cg'] = Projector(uv, self.fields.uv_cg).project

        self._add_interpolation_step('old_bathymetry_2d', self.depth.bathymetry_2d)
        self._add_interpolation_step('depth', self.depth.get_total_depth(elev))

        self.u = self.fields.uv_cg[0]
        self.v = self.fields.uv_cg[1]

        # define bed friction
        hc = conditional(self.fields.depth > Constant(0.001), self.fields.depth, Constant(0.001))
        aux = conditional(11.036*hc/self.bed_reference_height > Constant(1.001),
                          11.036*hc/self.bed_reference_height, Constant(1.001))
        self.qfc = Constant(2)/(ln(aux)/kappa)**2
        # skin friction coefficient
        cfactor = conditional(self.fields.depth > ksp, Constant(2)
                              * (((1/kappa)*ln(11.036*self.fields.depth/ksp))**(-2)), Constant(0.0))
        # mu - ratio between skin friction and normal friction
        self.mu = conditional(self.qfc > Constant(0), cfactor/self.qfc, Constant(0))

        # calculate bed shear stress
        self.unorm = (self.u**2) + (self.v**2)

        self._add_interpolation_step('bed_stress', self.rhow*Constant(0.5)*self.qfc*self.unorm)

        if self.solve_suspended_sediment:
            # deposition flux - calculating coefficient to account for stronger conc at bed
            B = conditional(self.a > self.fields.depth, Constant(1.0), self.a/self.fields.depth)
            ustar = sqrt(Constant(0.5)*self.qfc*self.unorm)
            rouse_number = (self.settling_velocity/(kappa*ustar)) - Constant(1)

            intermediate_step = conditional(abs(rouse_number) > Constant(1e-04),
                                            B*(Constant(1)-B**min_value(rouse_number, Constant(3)))/min_value(rouse_number, Constant(3)), -B*ln(B))
            

            self._add_interpolation_step('integrated_rouse',
                                         min_value(max_value(conditional(intermediate_step > Constant(1e-12), Constant(1)/intermediate_step,
                                                               Constant(1e12)), Constant(1)),Constant(5)), V=self.P1DG_2d)

            # erosion flux - above critical velocity bed is eroded
            transport_stage_param = conditional(self.rhow*Constant(0.5)*self.qfc*self.unorm*self.mu > Constant(0),
                                                (self.rhow*Constant(0.5)*self.qfc*self.unorm*self.mu - self.taucr)/self.taucr,
                                                Constant(-1))

            self._add_interpolation_step('erosion_concentration', Constant(0.015)*(self.average_size/self.a)
                                         * ((max_value(transport_stage_param, Constant(0)))**1.5)
                                         / (self.dstar**0.3), V=self.P1DG_2d)


            if self.use_advective_velocity_correction:
                correction_factor_model = CorrectiveVelocityFactor(self.fields.depth, ksp,
                                                                   self.bed_reference_height, self.settling_velocity, ustar)
                self.update_steps['correction_factor'] = correction_factor_model.update
                self.fields.velocity_correction_factor = correction_factor_model.velocity_correction_factor
            self._add_interpolation_step('equilibrium_tracer', self.fields.erosion_concentration/self.fields.integrated_rouse, V=self.P1DG_2d)

            # get individual terms
            self._deposition = self.settling_velocity*self.fields.integrated_rouse
            self._erosion = self.settling_velocity*self.fields.erosion_concentration

        if self.use_bedload:
            # calculate angle of flow
            self._add_interpolation_step('calfa', self.u/sqrt(self.unorm))
            self._add_interpolation_step('salfa', self.v/sqrt(self.unorm))
            if self.use_angle_correction:
                # slope effect angle correction due to gravity
                self._add_interpolation_step('stress', self.rhow*Constant(0.5)*self.qfc*self.unorm)

    def _add_interpolation_step(self, field_name, expr, V=None):
        """Add interpolation step to update

        :arg field_name: str name of new field to project into, stored in self.fields
        :arg expr: UFL expression to interpolate
        :kwarg V: FunctionSpace for new field (default is self.P1_2d)"""
        self.fields[field_name] = Function(V or self.P1_2d, name=field_name)
        self.update_steps[field_name] = Interpolator(expr, self.fields[field_name]).interpolate

    def get_bedload_term(self, bathymetry):
        """
        Returns expression for bedload transport :math:`(qbx, qby)` to be used in the Exner equation.
        Note bathymetry is the function which is solved for in the exner equation.

        :arg bathymetry: Bathymetry of the domain. Bathymetry stands for
            the bedlevel (positive downwards).
        """

        if self.use_slope_mag_correction:
            # slope effect magnitude correction due to gravity where beta is a parameter normally set to 1.3
            # we use z_n1 and equals so that we can use an implicit method in Exner
            slopecoef = Constant(1) + self.beta*(bathymetry.dx(0)*self.fields.calfa + bathymetry.dx(1)*self.fields.salfa)
        else:
            slopecoef = Constant(1.0)

        if self.use_angle_correction:
            # slope effect angle correction due to gravity
            cparam = Constant((self.rhos-self.rhow)*self.g*self.average_size*(self.surbeta2**2))
            tt1 = conditional(self.fields.stress > Constant(1e-10), sqrt(cparam/self.fields.stress), sqrt(cparam/Constant(1e-10)))

            # define bed gradient
            dzdx = self.fields.old_bathymetry_2d.dx(0)
            dzdy = self.fields.old_bathymetry_2d.dx(1)

            # add on a factor of the bed gradient to the normal
            aa = self.fields.salfa + tt1*dzdy
            bb = self.fields.calfa + tt1*dzdx

            comb = sqrt(aa**2 + bb**2)
            angle_norm = conditional(comb > Constant(1e-10), comb, Constant(1e-10))

            # we use z_n1 and equals so that we can use an implicit method in Exner
            calfamod = (self.fields.calfa + (tt1*bathymetry.dx(0)))/angle_norm
            salfamod = (self.fields.salfa + (tt1*bathymetry.dx(1)))/angle_norm

        if self.use_secondary_current:
            # accounts for helical flow effect in a curver channel
            # use z_n1 and equals so can use an implicit method in Exner
            free_surface_dx = self.fields.depth.dx(0) - bathymetry.dx(0)
            free_surface_dy = self.fields.depth.dx(1) - bathymetry.dx(1)

            velocity_slide = (self.u*free_surface_dy)-(self.v*free_surface_dx)

            tandelta_factor = Constant(7)*self.g*self.rhow*self.fields.depth*self.qfc\
                / (Constant(2)*self.alpha_secc*((self.u**2) + (self.v**2)))

            # accounts for helical flow effect in a curver channel
            if self.use_angle_correction:
                # if angle has already been corrected we must alter the corrected angle to obtain the corrected secondary current angle
                t_1 = (self.fields.bed_stress*slopecoef*calfamod) + (self.v*tandelta_factor*velocity_slide)
                t_2 = (self.fields.bed_stress*slopecoef*salfamod) - (self.u*tandelta_factor*velocity_slide)
            else:
                t_1 = (self.fields.bed_stress*slopecoef*self.fields.calfa) + (self.v*tandelta_factor*velocity_slide)
                t_2 = ((self.fields.bed_stress*slopecoef*self.fields.salfa) - (self.u*tandelta_factor*velocity_slide))

            # calculated to normalise the new angles
            t4 = sqrt((t_1**2) + (t_2**2))

            # updated magnitude correction and angle corrections
            slopecoef_secc = t4/self.fields.bed_stress

            calfanew = t_1/t4
            salfanew = t_2/t4

        # implement meyer-peter-muller bedload transport formula
        thetaprime = self.mu*(self.rhow*Constant(0.5)*self.qfc*self.unorm)/((self.rhos-self.rhow)*self.g*self.average_size)

        # if velocity above a certain critical value then transport occurs
        phi = conditional(thetaprime < self.thetacr, 0, Constant(8)*(thetaprime-self.thetacr)**1.5)

        # bedload transport flux with magnitude correction
        if self.use_secondary_current:
            qb_total = slopecoef_secc*phi*sqrt(self.g*self.R*self.average_size**3)
        else:
            qb_total = slopecoef*phi*sqrt(self.g*self.R*self.average_size**3)

        # formulate bedload transport flux with correct angle depending on corrections implemented
        if self.use_angle_correction and self.use_secondary_current is False:
            qbx = qb_total*calfamod
            qby = qb_total*salfamod
        elif self.use_secondary_current:
            qbx = qb_total*calfanew
            qby = qb_total*salfanew
        else:
            qbx = qb_total*self.fields.calfa
            qby = qb_total*self.fields.salfa

        return qbx, qby

    def get_deposition_coefficient(self):
        """Returns coefficient :math:`C` such that :math:`C/H*sediment` is deposition term in sediment equation

        If sediment field is depth-averaged, :math:`C*sediment` is (total) deposition (over the column)
        as it appears in the Exner equation, but deposition term in sediment equation needs
        averaging: :math:`C*sediment/H`
        If sediment field is depth-integrated, :math:`C*sediment/H` is (total) deposition (over the column)
        as it appears in the Exner equation, and is the same in the sediment equation."""
        return self._deposition

    def get_erosion_term(self):
        """Returns expression for (depth-integrated) erosion."""
        return self._erosion

    def get_equilibrium_tracer(self):
        """Returns expression for (depth-averaged) equilibrium tracer."""
        return self.fields.equilibrium_tracer

    def get_advective_velocity_correction_factor(self):
        """Returns correction factor for the advective velocity in the sediment equations

        With :attr:`.SedimentModelOptions.use_advective_velocity_correction`, this applies a correction to
        the supplied velocity solution `uv` to take into account the mismatch between
        depth-averaged product of velocity with sediment and product of depth-averaged
        velocity with depth-averaged sediment.
        """
        if self.use_advective_velocity_correction:
            return self.fields.velocity_correction_factor
        else:
            return 1

    def update(self):
        """Update all functions used by :class:`.SedimentModel`

        This repeats all projection and interpolations steps based on the current values
        of the `uv` and `elev` functions, provided in __init__."""
        for step in self.update_steps.values():
            step()
            # print_output('D:B1'+'\t'+str(self.fields.integrated_rouse.at((442416.06750038936, 3323689.309654328))))
            # print_output('E:B1'+'\t'+str(self.fields.erosion_concentration.at((442416.06750038936, 3323689.309654328))))
            
