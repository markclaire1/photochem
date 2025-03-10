submodule(photochem_atmosphere) photochem_atmosphere_utils
  implicit none
  
  ! Contains routines utility routines for returning or saving 
  ! model output
  
contains
  
  module subroutine out2atmosphere_txt(self, filename, overwrite, clip, err)
    use photochem_common, only: out2atmosphere_txt_base
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: rhs(self%var%neqs)  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ! update wrk variables
    call self%right_hand_side_chem(wrk%usol, rhs, err)
    if (allocated(err)) return

    call out2atmosphere_txt_base(dat, var, &
                                 wrk%pressure, wrk%density, wrk%densities, wrk%molecules_per_particle, &
                                 filename, overwrite, clip, err)
    if (allocated(err)) return

  end subroutine
  
  module subroutine out2in(self, err)
    class(Atmosphere), intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    
    if (self%var%at_photo_equilibrium) then
      self%var%usol_init = self%var%usol_out
    else
      err = "Can not set output to input without first converging to photochemical equilibrium."
      return
    endif
  end subroutine
  
  module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(out) :: surf_fluxes(:)
    real(dp), intent(out) :: top_fluxes(:)
    character(:), allocatable, intent(out) :: err
  
    real(dp) :: rhs(self%var%neqs)  
    real(dp) :: diffusive_production
    real(dp) :: chemical_production
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
  
    integer :: i
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(surf_fluxes) /= dat%nq .or. size(top_fluxes) /= dat%nq) then
      err = "Input fluxes to gas_fluxes has the wrong dimensions"
      return
    endif
  
    call self%right_hand_side_chem(wrk%usol, rhs, err)
    if (allocated(err)) return
    
    ! surface flux is molecules required to sustain the lower boundary
    ! chemical production + diffusion production = total change in lower cell    
    do i = 1,dat%nq
      diffusive_production = (wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                            + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1)) &
                              *wrk%density(1)*var%dz(1)
      chemical_production = rhs(i)*wrk%density(1)*var%dz(1)
      surf_fluxes(i) = -(diffusive_production + chemical_production)
    enddo
    
    ! fluxes going into or out of the top of the atmosphere.
    do i = 1,dat%nq
      diffusive_production = &
         (wrk%DD(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADD(i,var%nz)*wrk%usol(i,var%nz) &
        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1)) &
          *wrk%density(var%nz)*var%dz(var%nz)
    
      chemical_production = rhs(i + (var%nz-1)*dat%nq)*wrk%density(var%nz)*var%dz(var%nz)
      top_fluxes(i) = diffusive_production + chemical_production
    enddo
    
  end subroutine
  
  module function atom_conservation(self, atom, err) result(con)
    use photochem_enum, only: VelocityDistributedFluxBC
    use photochem_eqns, only: damp_condensation_rate
    use photochem_types, only: AtomConservation
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: atom
    character(:), allocatable, intent(out) :: err
    type(AtomConservation) :: con
    
    real(dp) :: surf_fluxes(self%dat%nq)
    real(dp) :: top_fluxes(self%dat%nq)
    real(dp) :: integrated_rainout(self%dat%nq)
    real(dp) :: cond_rate, con_evap_rate
    
    integer :: ind(1), i, j, kk
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ind = findloc(dat%atoms_names,atom)
    kk = ind(1)
    if (ind(1) == 0) then
      err = "Atom "//trim(atom)//" is not in the list of atoms."
      return
    endif
    if (dat%species_composition(ind(1),dat%nsp) /= 0) then
      err = "Atom "//trim(atom)//" makes up the background gas"// &
            " which is not conserved."
      return 
    endif
    
    call self%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (allocated(err)) return
    
    con%in_surf = 0
    con%in_top = 0
    con%in_dist = 0
    con%in_other = 0
    con%out_surf = 0
    con%out_top = 0
    con%out_rain = 0
    con%out_other = 0
    
    ! Upper and lower boundary
    do i = 1,dat%nq
      if (surf_fluxes(i) > 0) then
        con%in_surf = con%in_surf + surf_fluxes(i)*dat%species_composition(kk,i)
      else
        con%out_surf = con%out_surf + (-1.0_dp)*surf_fluxes(i)*dat%species_composition(kk,i)
      endif
      if (top_fluxes(i) > 0) then
        con%out_top = con%out_top + top_fluxes(i)*dat%species_composition(kk,i)
      else
        con%in_top = con%in_top + (-1.0_dp)*top_fluxes(i)*dat%species_composition(kk,i)
      endif
    enddo
    
    ! distributed fluxes
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        con%in_dist = con%in_dist + var%lower_flux(i)*dat%species_composition(kk,i)
      endif
    enddo
    
    ! rainout
    if (dat%gas_rainout) then
      integrated_rainout = 0.0_dp
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          integrated_rainout(i) = integrated_rainout(i) + &
                wrk%rainout_rates(i,j)*wrk%usol(i,j)*wrk%density(j)*var%dz(j)
        enddo
      enddo
      
      do i = 1,dat%nq
        con%out_rain = con%out_rain + integrated_rainout(i)*dat%species_composition(kk,i)
      enddo
    endif

    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        con_evap_rate = var%fast_arbitrary_rate*(wrk%H2O_sat_mix(j)*wrk%H2O_rh(j) - wrk%usol(dat%LH2O,j)) &
                        *wrk%density(j)*var%dz(j)*dat%species_composition(kk,dat%LH2O)
        if (con_evap_rate > 0.0_dp) then
          con%in_other = con%in_other + con_evap_rate
        else
          con%out_other = con%out_other + (-1.0_dp)*con_evap_rate
        endif
      enddo
    endif
    
    if (dat%water_cond) then
      if (dat%fix_water_in_trop) then
        i = var%trop_ind+1
      else
        i = 1
      endif
      do j = i,var%nz
        if (wrk%usol(dat%LH2O,j) >= var%H2O_condensation_rate(2)*wrk%H2O_sat_mix(j)) then
          cond_rate = damp_condensation_rate(var%H2O_condensation_rate(1), &
                                             var%H2O_condensation_rate(2), &
                                             var%H2O_condensation_rate(3), &
                                             wrk%usol(dat%LH2O,j)/wrk%H2O_sat_mix(j))
          con%out_other = con%out_other + cond_rate*(wrk%usol(dat%LH2O,j)  &
                        - var%H2O_condensation_rate(2)*wrk%H2O_sat_mix(j)) &
                        *wrk%density(j)*var%dz(j)*dat%species_composition(kk,dat%LH2O)
           
        endif
      enddo
    endif

    ! custom rate functions
    ! NOTE, This might be wrong. The model does not seem to conserve well
    ! for these custom functions.
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then; block
        real(dp) :: tmp_rate
        ! note that we use the time set in wrk%tn
        call var%rate_fcns(i)%fcn(wrk%tn, var%nz, wrk%xp) ! using wrk%xp space.
        tmp_rate = sum(wrk%xp*var%dz)*dat%species_composition(kk,i) ! atoms/cm^2/s
        if (tmp_rate > 0.0_dp) then
          con%in_other = con%in_other + tmp_rate
        else
          con%out_other = con%out_other + (-1.0_dp)*tmp_rate
        endif
      endblock; endif
    enddo

    con%net = con%in_surf + con%in_top + con%in_dist + con%in_other &
            - con%out_surf - con%out_top - con%out_rain - con%out_other
    
    con%factor = abs(con%net/maxval([con%in_surf, con%in_top, con%in_dist, con%in_other, &
                                     con%out_surf, con%out_top, con%out_rain, con%out_other]))
    
  end function
  
  module function redox_conservation(self, err) result(redox_factor)
    use photochem_enum, only: VelocityDistributedFluxBC
    class(Atmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    real(dp) :: redox_factor
    
    real(dp) :: surf_fluxes(self%dat%nq)
    real(dp) :: top_fluxes(self%dat%nq)
    real(dp) :: integrated_rainout(self%dat%nq)
    
    real(dp) :: oxi_in_surf, oxi_out_surf
    real(dp) :: red_in_surf, red_out_surf
    real(dp) :: oxi_in_top, oxi_out_top
    real(dp) :: red_in_top, red_out_top
    real(dp) :: oxi_in_dist
    real(dp) :: red_in_dist
    real(dp) :: oxi_out_rain
    real(dp) :: red_out_rain
    real(dp) :: oxi_in_other, oxi_out_other
    real(dp) :: red_in_other, red_out_other
    
    real(dp) :: oxi_in, oxi_out, red_in, red_out
    real(dp) :: net_redox
    
    integer :: i, j
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    call self%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (allocated(err)) return
    
    oxi_in_surf = 0.0_dp
    oxi_out_surf = 0.0_dp
    red_in_surf = 0.0_dp
    red_out_surf = 0.0_dp
    oxi_in_top = 0.0_dp
    oxi_out_top = 0.0_dp
    red_in_top = 0.0_dp
    red_out_top = 0.0_dp
    oxi_in_dist = 0.0_dp
    red_in_dist = 0.0_dp
    oxi_out_rain = 0.0_dp
    red_out_rain = 0.0_dp
    oxi_in_other = 0.0_dp
    oxi_out_other = 0.0_dp
    red_in_other = 0.0_dp
    red_out_other = 0.0_dp
    
    ! All numbers will be treated as positive.
    ! Later on, we implement signs
    
    ! boundary fluxes
    do i = 1,dat%nq
      if (dat%species_redox(i) > 0) then
        
        if (surf_fluxes(i) > 0) then
          oxi_in_surf = oxi_in_surf + surf_fluxes(i)*dat%species_redox(i)
        else
          oxi_out_surf = oxi_out_surf + (-1.0_dp)*surf_fluxes(i)*dat%species_redox(i)
        endif
        
        if (top_fluxes(i) > 0) then
          oxi_out_top = oxi_out_top + top_fluxes(i)*dat%species_redox(i)
        else
          oxi_in_top = oxi_in_top + (-1.0_dp)*top_fluxes(i)*dat%species_redox(i)
        endif
        
      elseif (dat%species_redox(i) < 0) then
        
        if (surf_fluxes(i) > 0) then
          red_in_surf = red_in_surf + surf_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        else
          red_out_surf = red_out_surf + (-1.0_dp)*surf_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        endif
        
        if (top_fluxes(i) > 0) then
          red_out_top = red_out_top + top_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        else
          red_in_top = red_in_top + (-1.0_dp)*top_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        endif
        
      endif
    enddo
    
    ! distributed fluxes
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        if (dat%species_redox(i) > 0) then
          oxi_in_dist = oxi_in_dist + var%lower_flux(i)*dat%species_redox(i)
        elseif (dat%species_redox(i) < 0) then
          red_in_dist = red_in_dist + var%lower_flux(i)*dat%species_redox(i)*(-1.0_dp)
        endif
      endif
    enddo
    
    ! rainout
    if (dat%gas_rainout) then
      integrated_rainout = 0.0_dp
      ! rhs_chem already got layer 1. So
      ! we start at layer 2
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          integrated_rainout(i) = integrated_rainout(i) + &
                wrk%rainout_rates(i,j)*wrk%usol(i,j)*wrk%density(j)*var%dz(j)
        enddo
      enddo
      
      do i = 1,dat%nq
        if (dat%species_redox(i) > 0) then
          oxi_out_rain = oxi_out_rain + integrated_rainout(i)*dat%species_redox(i)
        elseif (dat%species_redox(i) < 0) then
          red_out_rain = red_out_rain + integrated_rainout(i)*dat%species_redox(i)*(-1.0_dp)
        endif
      enddo
    endif

    ! custom rate functions
    ! NOTE, This might be wrong. The model does not seem to conserve well
    ! for these custom functions.
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then; block
        real(dp) :: tmp_rate
        ! note that we use the time set in wrk%tn
        call var%rate_fcns(i)%fcn(wrk%tn, var%nz, wrk%xp) ! using wrk%xp space.
        tmp_rate = sum(wrk%xp*var%dz)

        if (dat%species_redox(i) > 0) then
          if (tmp_rate > 0) then
            oxi_in_other = oxi_in_other + tmp_rate*dat%species_redox(i)
          else
            oxi_out_other = oxi_out_other + (-1.0_dp)*tmp_rate*dat%species_redox(i)
          endif
        elseif (dat%species_redox(i) < 0) then
          if (tmp_rate > 0) then
            red_in_other = red_in_other + tmp_rate*dat%species_redox(i)*(-1.0_dp)
          else
            red_out_other = red_out_other + (-1.0_dp)*tmp_rate*dat%species_redox(i)*(-1.0_dp)
          endif
        endif
      endblock; endif
    enddo
    
    ! total fluxes going in and out
    oxi_in = oxi_in_surf + oxi_in_top + oxi_in_dist + oxi_in_other
    red_in = red_in_surf + red_in_top + red_in_dist + red_in_other
    oxi_out = oxi_out_surf + oxi_out_top + oxi_out_rain + oxi_out_other
    red_out = red_out_surf + red_out_top + red_out_rain + red_out_other
    ! Net redox. Should be close to zero
    net_redox = oxi_in - red_in - oxi_out + red_out
    ! compute how close net_redox is to zero, relative to redox fluxes going in and out
    redox_factor = abs(net_redox/maxval([oxi_in, red_in, oxi_out, red_out]))

  end function
  
  module subroutine set_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC, VelocityDistributedFluxBC
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(dp), optional, intent(in) :: vdep
    real(dp), optional, intent(in) :: mix
    real(dp), optional, intent(in) :: flux
    real(dp), optional, intent(in) :: height
    character(:), allocatable, intent(out) :: err
    
    integer :: ind(1)
    
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species))
    if (ind(1) == 0) then
      err = "Can not change boundary condntion of '"//trim(species)// &
            "' because it is not in the list of species"
      return
    endif
    
    if (self%dat%fix_water_in_trop) then
      if (species == "H2O") then
        err = "You can not change the boundary condition for H2O because"// &
              " you have water fixed in the troposphere."
        return
      endif
    endif
    
    if (bc_type == 'vdep') then
      if (.not. present(vdep)) then
        err = "To change boundary condition to deposition"// &
              " velocity must supply the 'vdep' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = VelocityBC
      self%var%lower_vdep(ind(1)) = vdep
      
    elseif (bc_type == 'mix') then
      if (.not. present(mix)) then
        err = "To change boundary condition to fixed mixing"// &
              " ratio must supply the 'mix' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = MixingRatioBC
      self%var%lower_fix_mr(ind(1)) = mix
      
    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a surface flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = FluxBC
      self%var%lower_flux(ind(1)) = flux

    elseif (bc_type == 'vdep + dist flux') then
      if (.not.present(vdep) .or. .not.present(flux) .or. .not.present(height)) then
        err = "To change boundary condition to deposition velocity with"// &
              " a distributed flux, must supply the 'vdep', 'flux', and 'height' arguments"
        return
      endif
      self%var%lowerboundcond(ind(1)) = VelocityDistributedFluxBC
      self%var%lower_vdep(ind(1)) = vdep
      self%var%lower_flux(ind(1)) = flux
      self%var%lower_dist_height(ind(1)) = height
      
    elseif (bc_type == 'Moses') then
      self%var%lowerboundcond(ind(1)) = MosesBC
    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " boundary condition type"
      return
    endif
    
  end subroutine
    
  module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
    use photochem_enum, only: VelocityBC, FluxBC
    use photochem_enum, only: DiffusionLimHydrogenEscape
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(dp), optional, intent(in) :: veff
    real(dp), optional, intent(in) :: flux
    character(:), allocatable, intent(out) :: err
    
    integer :: ind(1)
    
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species))
    if (ind(1) == 0) then
      err = "Can not change boundary condition of '"//trim(species)// &
            "' because it is not in the list of species"
      return
    endif
    
    if (self%dat%H_escape_type == DiffusionLimHydrogenEscape) then
      if (species == "H2") then
        err = "You can not change the boundary condition for H2 because"// &
              " diffusion limited H2 escape is on."
        return
      endif
      if (species == "H") then
        err = "You can not change the boundary condition for H because"// &
              " diffusion limited H escape is on."
        return
      endif
    endif
    
    if (bc_type == 'veff') then
      if (.not. present(veff)) then
        err = "To change boundary condition to effusion"// &
              " velocity must supply the 'veff' argument"
        return
      endif
      self%var%upperboundcond(ind(1)) = VelocityBC
      self%var%upper_veff(ind(1)) = veff

    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%upperboundcond(ind(1)) = FluxBC
      self%var%upper_flux(ind(1)) = flux

    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " upper boundary condition type"
      return
    endif
  
  end subroutine
  
  module subroutine set_temperature(self, temperature, trop_alt, err)
    use photochem_input, only: interp2xsdata, compute_gibbs_energy
    
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: temperature(:)
    real(dp), optional, intent(in) :: trop_alt
    character(:), allocatable, intent(out) :: err
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemVars) :: var_save
    
    dat => self%dat
    var => self%var
    
    if (size(temperature) /= var%nz) then
      err = "temperature has the wrong input dimension"
      return
    endif
    
    ! save in case there is an issue
    var_save = var
    
    var%temperature = temperature
    
    ! xsections and gibbs energy needs updating
    call interp2xsdata(dat, var, err)
    if (allocated(err)) then
      var = var_save
      return
    endif
    if (dat%reverse) then
      call compute_gibbs_energy(dat, var, err)
      if (allocated(err)) then
        var = var_save
        return
      endif
    endif
    
    ! if water is fixed in troposhere or gas rainout, and trop_alt present
    ! then we need to change trop_ind, reallocate some stuff
    ! in wrk, then we will re-prep the atmosphere
    if ((dat%fix_water_in_trop .or. dat%gas_rainout) .and. present(trop_alt)) then
      if (trop_alt < var%bottom_atmos .or. trop_alt > var%top_atmos) then
        var = var_save
        err = "trop_alt is above or bellow the atmosphere!"
        return
      endif
      
      var%trop_alt = trop_alt
      var%trop_ind = max(minloc(abs(var%z - var%trop_alt), 1) - 1, 1)

      if (var%trop_ind < 3) then
        var = var_save
        err = 'Tropopause is too low.'
        return
      elseif (var%trop_ind > var%nz-2) then
        var = var_save
        err = 'Tropopause is too high.'
        return
      endif
                      
      call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                         self%var%nz, self%dat%nrT, self%dat%kj, &
                         self%dat%nw)

    endif
    
  end subroutine

  module subroutine set_photon_flux_fcn(self, photon_flux_fcn)
    use photochem_types, only: time_dependent_flux_fcn
    class(Atmosphere), target, intent(inout) :: self
    procedure(time_dependent_flux_fcn), pointer :: photon_flux_fcn
    self%var%photon_flux_fcn => photon_flux_fcn
  end subroutine

  module subroutine set_rate_fcn(self, species, fcn, err)
    use photochem_types, only: time_dependent_rate_fcn
    class(Atmosphere), target, intent(inout) :: self
    character(*), intent(in) :: species
    procedure(time_dependent_rate_fcn), pointer :: fcn
    character(:), allocatable, intent(inout) :: err
    
    integer :: ind

    ind = findloc(self%dat%species_names(1:self%dat%nq), species, 1)
    if (ind == 0) then
      err = 'Species "'//species//'" is not in the list of species, '// &
            'or is a background or short-lived species.'
      return
    endif

    self%var%rate_fcns(ind)%fcn => fcn

  end subroutine
  
end submodule