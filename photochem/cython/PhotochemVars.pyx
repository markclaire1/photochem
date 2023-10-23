cimport PhotochemVars_pxd as var_pxd

cdef class PhotochemVars:
  """This class contains data that can change between independent
  model integrations.
  """

  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      var_pxd.allocate_photochemvars(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      var_pxd.deallocate_photochemvars(&self._ptr)
      self._ptr = NULL
  
  property nz:
    "The number of vertical atmospheric layers"
    def __get__(self):
      cdef int nz
      var_pxd.photochemvars_nz_get(&self._ptr, &nz)
      return nz

  property top_atmos:
    "The top of the model domain (cm)"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_top_atmos_get(&self._ptr, &val)
      return val

  property bottom_atmos:
    "The bottom of the model domain (cm)"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_bottom_atmos_get(&self._ptr, &val)
      return val
      
  property usol_init:
    """ndarray[double,dim=2], shape (nq,nz). Contains the initial concentration
    of atmospheric species. For the model `Atmosphere` then units are mixing ratios,
    and if the model is `EvoAtmosphere` then the units are molecules/cm^3.
    """
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_usol_init_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_usol_init_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property z:
    "ndarray[double,dim=1], shape (nz). The altitude of the center of each atmopsheric layer (cm)"
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_z_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_z_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property temperature:
    "ndarray[double,dim=1], shape (nz). The temperature of each atmospheric layer (K)"
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_temperature_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_temperature_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property edd:
    "ndarray[double,dim=1], shape (nz). The eddy diffusion of each atmospheric layer (cm^2/s)"
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_edd_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_edd_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property photon_flux:
    "ndarray[double,dim=1], shape (nw). photon/cm^2/s in each wavelength bin hitting planet."
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_photon_flux_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_photon_flux_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property grav:
    "ndarray[double,dim=1], shape (nz). The gravitational acceleration at the center of each grid cell (cm/s^2)."
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_grav_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_grav_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property at_photo_equilibrium:
    "bool. If True, then the model is at photochemical equilibrium."
    def __get__(self):
      cdef bool at_photo_equilibrium
      var_pxd.photochemvars_at_photo_equilibrium_get(&self._ptr, &at_photo_equilibrium)
      return at_photo_equilibrium
      
  property surface_pressure:
    "double. The surface pressure in bars."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_surface_pressure_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_surface_pressure_set(&self._ptr, &val)

  property max_error_reinit_attempts:
    """int. number of times to reinitialize CVODE when it returns
    a potentially recoverable error. Only used in `EvoAtmosphere` (not `Atmosphere`)
    """
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_max_error_reinit_attempts_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_max_error_reinit_attempts_set(&self._ptr, &val)
  
  property rtol:
    "double. Integration relative tolerance."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_rtol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_rtol_set(&self._ptr, &val)
      
  property atol:
    "double. Integration absolute tolerance."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_atol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_atol_set(&self._ptr, &val)

  property mxsteps:
    "int. Max number of steps before integrator will give up."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_mxsteps_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_mxsteps_set(&self._ptr, &val)
      
  property equilibrium_time:
    "double. Atomsphere considered in equilibrium if integrations reaches this time (seconds)"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_equilibrium_time_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_equilibrium_time_set(&self._ptr, &val)
  
  property verbose:
    "int. 0 == no printing. 1 == some printing. 2 == bunch of printing."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_verbose_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_verbose_set(&self._ptr, &val)

  property fast_arbitrary_rate:
    "double. arbitrary rate that is fast (1/s). Used for keeping H2O at saturation in troposphere"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_fast_arbitrary_rate_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_fast_arbitrary_rate_set(&self._ptr, &val)

  property clipval:
    "clipping value used in out2atmosphere_txt"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_clipval_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_clipval_set(&self._ptr, &val)
    
  