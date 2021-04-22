
module photochem_types ! make a giant IO object
  implicit none
  
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
  public PhotoMechanism, PhotoSettings, PhotoPlanet
  
  type :: PhotoSettings
    real(real_kind) :: bottom_atmosphere
    real(real_kind) :: top_atmosphere 
    integer :: nz = 0
    
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw !number of bins
  end type
  
  type :: PhotoPlanet
    real(real_kind) :: gravity
    real(real_kind) :: surface_pressure
    real(real_kind) :: planet_radius
    real(real_kind) :: surface_albedo
    logical :: water_sat_trop
    real(real_kind) :: trop_alt
    logical :: lightning
    real(real_kind) :: lightning_NO_production
    logical :: rainout
    real(real_kind) :: rainout_multiplier
  end type
  
  ! type :: PhotoMolecules
  !   integer :: nsp, natoms
  !   character(len=8), allocatable :: atoms_names(:)
  !   character(len=8), allocatable :: species_names(:)
  !   integer, allocatable :: species_composition(:,:)
  !   integer, allocatable :: lowerboundcond(:)
  !   real(real_kind), allocatable :: lower_vdep(:)
  !   real(real_kind), allocatable :: lower_flux(:)
  !   real(real_kind), allocatable :: lower_distributed_height(:)
  !   real(real_kind), allocatable :: lower_fixed_mr(:)
  !   integer, allocatable :: upperboundcond(:)
  !   real(real_kind), allocatable :: upper_veff(:)
  !   real(real_kind), allocatable :: upper_flux(:)
  !   real(real_kind), allocatable :: thermo_data(:,:,:)
  !   real(real_kind), allocatable :: thermo_temps(:,:)  
  ! end type
  ! 
  ! type :: PhotoReactions
  !   integer :: nr
  !   character(len=8), allocatable :: reactions_names(:,:)
  !   integer, allocatable :: reactions_indices(:,:)
  !   character(len=15), allocatable :: rxtypes(:)
  !   real(real_kind), allocatable :: rateparams(:,:)
  ! 
  !   integer, allocatable :: nump(:)
  !   integer, allocatable :: numl(:)
  !   integer, allocatable :: iprod(:,:)
  !   integer, allocatable :: iloss(:,:,:)
  ! end type
  
  type :: PhotoMechanism
    type(PhotoSettings) :: settings
    type(PhotoPlanet) :: planet
    ! type(PhotoMolecules) :: molecules
    ! type(PhotoReactions) : reactions
    
    integer :: nsp
    integer :: nrF
    integer :: nrR
    integer :: nrT
    integer :: natoms
    character(len=8), allocatable :: atoms_names(:)
    character(len=8), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    integer, allocatable :: lowerboundcond(:)
    real(real_kind), allocatable :: lower_vdep(:)
    real(real_kind), allocatable :: lower_flux(:)
    real(real_kind), allocatable :: lower_distributed_height(:)
    real(real_kind), allocatable :: lower_fixed_mr(:)
    integer, allocatable :: upperboundcond(:)
    real(real_kind), allocatable :: upper_veff(:)
    real(real_kind), allocatable :: upper_flux(:)
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    
    integer :: max_num_reactants
    integer :: max_num_products
    character(len=8), allocatable :: reactants_names(:,:) ! not really needed.
    character(len=8), allocatable :: products_names(:,:)
    integer, allocatable :: reactants_sp_inds(:,:)
    integer, allocatable :: products_sp_inds(:,:)
    integer, allocatable :: nreactants(:)
    integer, allocatable :: nproducts(:)
    
    integer, allocatable :: reverse_info(:,:) ! all for calculating rates
    character(len=15), allocatable :: rxtypes(:)
    real(real_kind), allocatable :: rateparams(:,:)
    
    integer, allocatable :: nump(:) ! length nsp. number of 
    integer, allocatable :: numl(:)
    integer, allocatable :: iprod(:,:)
    integer, allocatable :: iloss(:,:)
    
  end type
  
end module


! module photochem_vars ! unpack the object to plain variables 
!   use photochem_types, only : PhotoMechanism
!   implicit none
!   integer, private, parameter :: real_kind = kind(1.0d0)
!   public
! end module

! module settings
! end module
! 
! module planet
! end module
! 
! module molecules
! end module
! 
! module reactions
! end module

module photochem_io
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use stringifor, only : string
  use photochem_types, only : PhotoMechanism, PhotoSettings, PhotoPlanet
  implicit none

  private 
  
  public get_photodata
  
  integer,parameter :: real_kind = kind(1.0d0)
  
contains
  
  subroutine get_photodata(infile,photomech) ! read yaml and make giant IO object
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(out) :: photomech
    
    character(error_length) :: error
    class (type_node), pointer :: root
    
    root => parse(infile,unit=100,error=error)
    if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
      stop
    end if
    
    select type (root)
      class is (type_dictionary)
        call parserootdict(root,infile,photomech)
        call root%finalize()
        deallocate(root)
      class default
        print*,"yaml file must have dictionaries at root level"
        stop
    end select
  end subroutine
  
  subroutine parserootdict(mapping, infile,photomech)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(out) :: photomech
    
    class (type_dictionary), pointer :: settings, planet
    class (type_list), pointer :: atoms, species, reactions
    type (type_error), pointer :: config_error
    class (type_list_item), pointer :: item
    class (type_dictionary), pointer :: dict
    class (type_key_value_pair), pointer :: key_value_pair

    ! temporary work variables
    type(string) :: tmp
    type(string), allocatable :: tmps(:), eqr(:), eqp(:)
    character(len=8) :: outstr(5)
    integer :: outarr(5)
    integer :: i, j, ind(1)
    logical :: reverse
    
    settings => mapping%get_dictionary('settings',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    planet => mapping%get_dictionary('planet',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    atoms => mapping%get_list('atoms',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    species => mapping%get_list('species',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    reactions => mapping%get_list('reactions',.true.,error = config_error) 
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! settings
    call get_settings(settings, infile, photomech%settings)
    
    ! planet
    call get_planet(planet, infile, photomech%planet)
    
    ! atoms
    item => atoms%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_scalar)
        tmp = tmp//" "//trim(element%string)
      class default
        print*,"IOError: Problem reading in atoms."
        stop
      end select
      item => item%next
    enddo
    call tmp%split(tokens=tmps) ! list of atoms
    photomech%natoms = size(tmps)
    allocate(photomech%atoms_names(photomech%natoms))
    do i = 1,photomech%natoms
      photomech%atoms_names(i) = tmps(i)%chars()
    enddo
    ! done with atoms
    
    ! now do species
    photomech%nsp = 0 ! count number of species
    item => species%first
    do while (associated(item))
      item => item%next
      photomech%nsp = photomech%nsp + 1
    enddo
        
    allocate(photomech%species_composition(photomech%natoms,photomech%nsp+2))
    photomech%species_composition = 0
    allocate(photomech%species_names(photomech%nsp+2))
    photomech%species_names(photomech%nsp+1) = "M" ! always add these guys
    photomech%species_names(photomech%nsp+2) = "hv"
    allocate(photomech%lowerboundcond(photomech%nsp))
    allocate(photomech%lower_vdep(photomech%nsp))
    allocate(photomech%lower_flux(photomech%nsp))
    allocate(photomech%lower_distributed_height(photomech%nsp))
    allocate(photomech%lower_fixed_mr(photomech%nsp))
    allocate(photomech%upperboundcond(photomech%nsp))
    allocate(photomech%upper_veff(photomech%nsp))
    allocate(photomech%upper_flux(photomech%nsp))
    allocate(photomech%thermo_data(7,2,photomech%nsp))
    allocate(photomech%thermo_temps(3,photomech%nsp))
    j = 1
    item => species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        dict => element%get_dictionary("composition",.true.,error = config_error)  ! get composition
        if (associated(config_error)) call handleerror(config_error%message,infile)
        key_value_pair => dict%first ! to see if
        do while (associated(key_value_pair))
          ind = findloc(photomech%atoms_names,trim(key_value_pair%key))
          if (ind(1) == 0) then
            print*,'IOError: The atom "', trim(key_value_pair%key), '" is not in the list of atoms.'
            stop
          endif
          key_value_pair =>key_value_pair%next
        enddo
        
        do i=1,photomech%natoms
          photomech%species_composition(i,j) = dict%get_integer(photomech%atoms_names(i),0,error = config_error) ! no error possible.
        enddo
        photomech%species_names(j) = trim(element%get_string("name",error = config_error)) ! get name
        if (associated(config_error)) call handleerror(config_error%message,infile)
        call get_boundaryconds(element,photomech%species_names(j), infile, &
                               photomech%lowerboundcond(j), photomech%lower_vdep(j), &
                               photomech%lower_flux(j), photomech%lower_distributed_height(j), &
                               photomech%lower_fixed_mr(j), &
                               photomech%upperboundcond(j), photomech%upper_veff(j), photomech%upper_flux(j))! get boundary conditions
        
        call get_thermodata(element,photomech%species_names(j), infile,photomech%thermo_temps(:,j),photomech%thermo_data(:,:,j)) ! get thermodynamic data
        
      class default
        print*,"IOError: Problem with species number ", j,"  in the input file"
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    
    
    ! reactions
    photomech%nrF = 0 ! count forward reactions
    item => reactions%first
    do while (associated(item))
      item => item%next
      photomech%nrF = photomech%nrF + 1
    enddo
    
    allocate(photomech%reverse_info(2,photomech%nrF))
    photomech%reverse_info = 0
    ! determine which reactions to reverse. Determine maximum number of reactants, and productants
    photomech%max_num_reactants = 1
    photomech%max_num_products = 1
    photomech%nrR = 0
    j = 1
    item => reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = config_error))
        call parse_reaction(tmp, reverse, eqr, eqp)
        if (reverse) then
          photomech%nrR = photomech%nrR + 1
          photomech%reverse_info(1,j) = 1  ! whether the reaction is reversed
          photomech%reverse_info(2,j) = photomech%nrR + photomech%nrF ! the reaction number of reversed reaction
          if (size(eqr) > photomech%max_num_products) photomech%max_num_products = size(eqr)
          if (size(eqp) > photomech%max_num_reactants) photomech%max_num_reactants = size(eqp)
        endif
        if (size(eqr) > photomech%max_num_reactants) photomech%max_num_reactants = size(eqr)
        if (size(eqp) > photomech%max_num_products) photomech%max_num_products = size(eqp)
      class default
        print*,"IOError: Problem with reaction number ",j," in the input file."
        stop
      end select
      item => item%next
      j = j+1
    enddo
    photomech%nrT = photomech%nrR + photomech%nrF

    ! allocate stuff and loop through reactions again
    allocate(photomech%nreactants(photomech%nrT))
    allocate(photomech%nproducts(photomech%nrT))
    allocate(photomech%reactants_sp_inds(photomech%max_num_reactants,photomech%nrT))
    allocate(photomech%products_sp_inds(photomech%max_num_products,photomech%nrT))
    allocate(photomech%reactants_names(photomech%max_num_reactants,photomech%nrF))
    allocate(photomech%products_names(photomech%max_num_products,photomech%nrF))
    allocate(photomech%rateparams(6,photomech%nrF)) ! This will need knowledge of which rection it is reversing
    allocate(photomech%rxtypes(photomech%nrF))
    j = 1
    item => reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = config_error))
        if (associated(config_error)) call handleerror(config_error%message,infile)
        
        call parse_equation(tmp, photomech%max_num_reactants, photomech%max_num_products, &
                            photomech%nreactants(j), photomech%nproducts(j), &
                            photomech%reactants_names(:,j), photomech%products_names(:,j))
        if (photomech%reverse_info(1,j) == 1) then
          ! reaction has a reverse
          i = photomech%reverse_info(2,j)
          photomech%nreactants(i) = photomech%nproducts(j)
          photomech%nproducts(i) = photomech%nreactants(j)
        endif

        call species_name2number(tmp, photomech%max_num_reactants, photomech%max_num_products, &
                                 photomech%reactants_names(:,j), photomech%products_names(:,j), &
                                 photomech%species_names, photomech%species_composition, &
                                 photomech%natoms, photomech%nsp, &
                                 photomech%reactants_sp_inds(:,j), photomech%products_sp_inds(:,j))
        if (photomech%reverse_info(1,j) == 1) then
          ! reaction has a reverse
          i = photomech%reverse_info(2,j)
          photomech%reactants_sp_inds(:,i) = photomech%products_sp_inds(:,j)
          photomech%products_sp_inds(:,i) = photomech%reactants_sp_inds(:,j)
        endif

        call get_rateparams(element, infile, photomech%rxtypes(j), photomech%rateparams(:,j))
      class default
        print*,"IOError: Problem with reaction number ",j," in the input file."
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! now get nump and numl
    
    ! check for inconsistencies
    ! if rainout is on, water must be a reactant
    ind = findloc(photomech%species_names,'H2O')
    if ((photomech%planet%rainout) .and. (ind(1)==0)) then
      print*,'IOError: H2O must be a species when rainout is turned on'
      stop
    endif
    if ((photomech%planet%rainout) .and. (ind(1)==0)) then
      print*,'IOError: H2O must be a species when rainout is turned on'
      stop
    endif
    

  end subroutine
  
  
  subroutine get_settings(filesettings, infile, outsettings)
    class(type_dictionary), intent(in) :: filesettings
    character(len=*), intent(in) :: infile
    type(PhotoSettings), intent(out) :: outsettings
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    
    tmpdict => filesettings%get_dictionary("atmosphere-grid",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outsettings%bottom_atmosphere = tmpdict%get_real("bottom",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%top_atmosphere = tmpdict%get_real("top",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%nz = tmpdict%get_integer("number-of-layers",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    tmpdict => filesettings%get_dictionary("photo-grid",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outsettings%lower_wavelength = tmpdict%get_real("lower-wavelength",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%upper_wavelength = tmpdict%get_real("upper-wavelength",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%nw = tmpdict%get_integer("number-of-bins",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
  
  end subroutine
  
  subroutine get_planet(fileplanet, infile, outplanet)
    class(type_dictionary), intent(in) :: fileplanet
    character(len=*), intent(in) :: infile
    type(PhotoPlanet), intent(out) :: outplanet
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    
    outplanet%gravity = fileplanet%get_real("gravity",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outplanet%surface_pressure = fileplanet%get_real("surface-pressure",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outplanet%planet_radius = fileplanet%get_real("planet-radius",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)  
    outplanet%surface_albedo = fileplanet%get_real("surface-albedo",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)  
    outplanet%water_sat_trop = fileplanet%get_logical("water-saturated-troposhere",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    if (outplanet%water_sat_trop) then
      outplanet%trop_alt = fileplanet%get_real("tropopause-altitude",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    else
      outplanet%trop_alt = 0.d0 ! no tropopause.
    endif
    
    tmpdict => fileplanet%get_dictionary("lightning",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outplanet%lightning = tmpdict%get_logical("on-off",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    if (.not. outplanet%lightning) then
      outplanet%lightning_NO_production = 0.d0
    else
      outplanet%lightning_NO_production = tmpdict%get_real("NO-production",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    endif
    tmpdict => fileplanet%get_dictionary("rainout",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outplanet%rainout= tmpdict%get_logical("on-off",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    if (.not.outplanet%rainout) then
      outplanet%rainout_multiplier = 1.d0
    else
      outplanet%rainout_multiplier = tmpdict%get_real("multiplier",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    endif
  end subroutine
  
  subroutine get_boundaryconds(molecule, molecule_name, infile, &
                               lowercond, Lvdep, Lflux, LdistH, Lmr, &
                               uppercond, Uveff, Uflux)
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    integer, intent(out) :: lowercond
    real(real_kind), intent(out) :: Lvdep, Lflux, LdistH, Lmr
    integer, intent(out) :: uppercond
    real(real_kind), intent(out) :: Uveff, Uflux
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    character(len=:), allocatable :: bctype
    
    tmpdict => molecule%get_dictionary("lower-boundary",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    bctype = tmpdict%get_string("type",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! constant deposition velocity = vdep
    ! constant mixing-ratio = mixing-ratio
    ! constant flux = flux
    ! deposition velocity + distributed flux = vdep + flux + distributed-height
    if (bctype == "constant deposition velocity") then
      lowercond = 0
      Lvdep = tmpdict%get_real("vdep",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "constant mixing-ratio") then
      lowercond = 1
      Lvdep = 0.d0
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = tmpdict%get_real("mixing-ratio",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    elseif (bctype == "constant flux") then
      lowercond = 2
      Lvdep = 0.d0
      Lflux = tmpdict%get_real("flux",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "deposition velocity + distributed flux") then
      lowercond = 3
      Lvdep = tmpdict%get_real("vdep",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Lflux = tmpdict%get_real("flux",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      LdistH = tmpdict%get_real("distributed-height",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Lmr = 0.d0 
    else
      print*,'IOError: "',trim(bctype),'" is not a valid lower boundary condition for ',trim(molecule_name)
      stop
    endif
    
    tmpdict => molecule%get_dictionary("upper-boundary",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    bctype = tmpdict%get_string("type",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! constant effusion velocity = veff
    ! constant flux = flux
    if (bctype == "constant effusion velocity") then
      uppercond = 0
      Uveff = tmpdict%get_real("veff",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Uflux = 0.d0
    elseif (bctype == "constant flux") then
      uppercond = 2
      Uveff = 0.d0
      Uflux = tmpdict%get_real("flux",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    else
      print*,'IOError: "',trim(bctype),'" is not a valid upper boundary condition for ',trim(molecule_name)
      stop
    endif
    
  end subroutine
    
  
  subroutine get_thermodata(molecule, molecule_name, infile, thermo_temps_entry, thermo_data_entry)
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    real(real_kind), intent(out) :: thermo_temps_entry(3)
    real(real_kind), intent(out) :: thermo_data_entry(7,2)
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict, tmpdict1
    class(type_list), pointer :: tmplist
    class(type_list_item), pointer :: item
    logical :: success
    
    integer :: j, i
    
    thermo_temps_entry = -1.d0
    thermo_data_entry = -1.d0
    
    tmpdict => molecule%get_dictionary("thermo",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! check thermodynamic model
    if (tmpdict%get_string("model",error = config_error) /= "Shomate") then
      print*,"IOError: Thermodynamic data must be in Shomate format for ",trim(molecule_name)
      stop
    endif
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    j = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_scalar)
        if (j > 3) then
          print*,"IOError: Too many temperature ranges for ",trim(molecule_name)
          stop
        endif
        thermo_temps_entry(j) = listitem%to_real(-1.d0,success)
        if (.not. success) then
          print*,"IOError: Problem reading thermodynamic data for  ",trim(molecule_name)
          stop
        endif
      class default
        print*,"IOError: Problem reading thermodynamic data for ",trim(molecule_name)
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! check amount of thermodynamic data
    if (thermo_temps_entry(3) == -1.d0) then
      i = 1
    else
      i = 2
    endif
    
    ! get data
    tmpdict1 => tmpdict%get_dictionary("data",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    tmplist =>tmpdict1%get_list("poly1",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    j = 1
    item => tmplist%first
    do while (associated(item)) 
      select type (listitem => item%node)
      class is (type_scalar)
        if (j > 7) then
          print*,"IOError: Too much thermodynamic data for ",trim(molecule_name)
          stop
        endif
        thermo_data_entry(j, 1) = listitem%to_real(-1.d0,success)
        if (.not. success) then
          print*,"IOError: Problem reading thermodynamic data for  ",trim(molecule_name)
          stop
        endif
      class default
        print*,"IOError: Problem reading thermodynamic data for ",trim(molecule_name)
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    if (j-1 /= 7) then
      print*,"IOError: Missing thermodynamic data for ",trim(molecule_name)
      stop
    endif

    tmplist =>tmpdict1%get_list("poly2",.false.,error = config_error)
    if ((.not.associated(tmplist)) .and. (i == 1)) then
      ! nothing happens
    elseif ((.not.associated(tmplist)) .and. (i == 2)) then
      print*,'IOError: More temperature ranges than thermodynamic data for ',trim(molecule_name)
      stop
    else
      j = 1
      item => tmplist%first
      do while (associated(item)) 
        select type (listitem => item%node)
        class is (type_scalar)
          if (j > 7) then
            print*,"IOError: Too much thermodynamic data for ",trim(molecule_name)
            stop
          endif
          thermo_data_entry(j, 2) = listitem%to_real(-1.d0,success)
          if (.not. success) then
            print*,"IOError: Problem reading thermodynamic data for  ",trim(molecule_name)
            stop
          endif
        class default
          print*,"IOError: Problem reading thermodynamic data for ",trim(molecule_name)
          stop
        end select
        item => item%next
        j = j + 1
      enddo
      if (j-1 /= 7) then
        print*,"IOError: Missing thermodynamic data for ",trim(molecule_name)
        stop
      endif
    endif
    
  end subroutine
  
  subroutine get_rateparams(reaction, infile, rxtype, rateparam)
    class(type_dictionary), intent(in) :: reaction
    character(len=*), intent(in) :: infile
    character(len=15), intent(out) :: rxtype
    real(real_kind), intent(out) :: rateparam(6)
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    class(type_scalar), pointer :: tmpscalar
    
    rateparam = 0.d0
    
    rxtype= reaction%get_string("type","",error = config_error) ! no error possible
    if (trim(rxtype) == '') rxtype = "elementary"
    
    ! get params
    if ((trim(rxtype) == 'elementary') .or. (trim(rxtype) == 'three-body')) then
      tmpdict => reaction%get_dictionary('rate-constant',.true.,error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(1) = tmpdict%get_real('A',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(2) = tmpdict%get_real('b',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(3) = tmpdict%get_real('Ea',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    elseif (trim(rxtype) == 'falloff') then
      tmpdict => reaction%get_dictionary('low-P-rate-constant',.true.,error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(1) = tmpdict%get_real('A',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(2) = tmpdict%get_real('b',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(3) = tmpdict%get_real('Ea',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      tmpdict => reaction%get_dictionary('high-P-rate-constant',.true.,error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(4) = tmpdict%get_real('A',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(5) = tmpdict%get_real('b',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      rateparam(6) = tmpdict%get_real('Ea',error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    elseif (trim(rxtype) == 'photolysis') then
      ! nothing
    else
      print*,'IOError: reaction type ',trim(rxtype),' is not a valid reaction type.'
      stop
    endif
    
  end subroutine
  
  subroutine parse_reaction(instring, reverse, eqr, eqp)
    type(string), intent(in) :: instring
    logical, intent(out) :: reverse
    type(string), allocatable, intent(out) :: eqr(:), eqp(:)
    
    type(string) :: string1, string2, string3, string4
    type(string), allocatable :: eq1(:)
    integer i
    string1 = instring%replace(old='+', new=' ')
    string2 = string1%replace(old='(', new=' ')
    string3 = string2%replace(old=')', new=' ')
    if (index(instring%chars(), "<=>") /= 0) then
      call string3%split(eq1, sep="<=>")
      reverse = .true.
    elseif (index(instring%chars(), " =>") /= 0) then
      call string3%split(eq1, sep=" =>")
      reverse = .false.
    else
      print*,"IOError: Invalid reaction arrow in reaction ",instring
      stop
    endif
    call eq1(1)%split(eqr, sep=" ")
    call eq1(2)%split(eqp, sep=" ")
  end subroutine
  
  subroutine parse_equation(instring, max_num_react, max_num_prod, numr, nump, outreact, outprod)
    type(string), intent(in) :: instring
    integer, intent(in) :: max_num_react, max_num_prod
    integer, intent(out) :: numr, nump
    character(len=8), intent(out) :: outreact(max_num_react), outprod(max_num_prod)
    logical :: reverse
    type(string), allocatable :: eqr(:), eqp(:)
    integer :: i
    
    call parse_reaction(instring, reverse, eqr, eqp)
    
    numr = size(eqr)
    nump = size(eqp)
    
    outreact = ''
    outprod = ''
    do i=1,numr
      outreact(i) = eqr(i)%chars()
    enddo
    do i=1,nump
      outprod(i) = eqp(i)%chars()
    enddo
  end subroutine
  
  subroutine species_name2number(reaction, max_num_react, max_num_prod, reacts, prods, &
                                 species_names, species_composition, natoms, nsp, &
                                 react_sp_nums, prod_sp_nums)
    type(string) :: reaction
    integer, intent(in) :: max_num_react, max_num_prod
    
    character(len=8), intent(in) :: reacts(max_num_react)
    character(len=8), intent(in) :: prods(max_num_prod)
    
    character(len=8), intent(in) :: species_names(nsp+2)
    integer, intent(in) :: species_composition(natoms,nsp+2)
    integer, intent(in) :: nsp, natoms
    
    integer, intent(out) :: react_sp_nums(max_num_react)
    integer, intent(out) :: prod_sp_nums(max_num_prod)
    
    integer :: i, ind(1)
    integer :: reactant_atoms(natoms), product_atoms(natoms)
    reactant_atoms = 0
    product_atoms = 0
    
    do i = 1,max_num_react
      ind = findloc(species_names,reacts(i))
      react_sp_nums(i) = ind(1)
      if ((reacts(i) /= '') .and. (ind(1) == 0)) then
        print*,"IOError: ", & 
               "Species ",trim(reacts(i))," in reaction ",reaction, &
               "is not in the list of species."
        stop
      endif
    enddo
    
    do i = 1,max_num_prod
      ind = findloc(species_names,prods(i))
      prod_sp_nums(i) = ind(1)
      if ((prods(i) /= '') .and. (ind(1) == 0)) then
        print*,"IOError: ", & 
               "Species ",trim(reacts(i))," in reaction ",reaction, &
               "is not in the list of species."
        stop
      endif
    enddo
  
    do i=1,max_num_react
      if (react_sp_nums(i) /= 0) then
        reactant_atoms = reactant_atoms + species_composition(:,react_sp_nums(i))
      endif
    enddo
    do i=1,max_num_prod
      if (prod_sp_nums(i) /= 0) then
        product_atoms = product_atoms + species_composition(:,prod_sp_nums(i))
      endif
    enddo
    if (.not. all(reactant_atoms == product_atoms)) then
      print*,"IOError: ", & 
             "Bad mass balance in reaction",reaction
      print*,"You could have messed up how many atoms one of the species has."
      stop
    endif
  end subroutine
  
  subroutine handleerror(message,infile)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: infile
    print*,"IOError: ",trim(infile),trim(message)
    stop
  end subroutine

end module

program main
  use photochem_io, only: get_photodata
  use photochem_types, only: PhotoMechanism
  implicit none
  type(PhotoMechanism) :: photomech
  integer i

  call get_photodata("../zahnle.yaml", photomech)


  do i = 1,photomech%nrT
    print*,photomech%reactants_sp_inds(:,i), "=>",photomech%products_sp_inds(:,i), "|", &
    photomech%nreactants(i),photomech%nproducts(i)
  enddo
end program






