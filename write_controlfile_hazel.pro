

pro write_controlfile_hazel,file_name,input_model,file_obs,file_inv_profiles,file_inv_parameters,file_error_inv,action,verbose,linear_solver,stopping_volume,synthesis_mode,stimulated_emission,magnetic_field,Paschen_Back,atomic_polarization,magneto_optical,stimulated_emission_RT,multiplet,LOS_angles,Wavelength_axis,number_slabs,boundary_condition,aa,height,beta,ff,BB1,thetaB1,chiB1,vdopp1,tau1,vmac1,BB2,thetaB2,chiB2,vdopp2,tau2,vmac2,aa_ranges,beta_ranges,ff_ranges,BB1_ranges,thetaB1_ranges,chiB1_ranges,vdopp1_ranges,tau1_ranges,vmac1_ranges,BB2_ranges,thetaB2_ranges,chiB2_ranges,vdopp2_ranges,tau2_ranges,vmac2_ranges,iterationsLM,ncycles,inversion_modes,aa_cycles,beta_cycles,ff_cycles,BB1_cycles,thetaB1_cycles,chiB1_cycles,vdopp1_cycles,tau1_cycles,vmac1_cycles,BB2_cycles,thetaB2_cycles,chiB2_cycles,vdopp2_cycles,tau2_cycles,vmac2_cycles,Weights_I,Weights_Q,Weights_U,Weights_V




OPENW,1,file_name, WIDTH=150


PRINTF,1,'# Hazel configuration File'
PRINTF,1,' '
PRINTF,1,'#####################'
PRINTF,1,'# General information'
PRINTF,1,'#####################'
PRINTF,1,' '
PRINTF,1,'[Files]'
PRINTF,1,'Input model file = ',input_model
PRINTF,1,'File with observations = ',file_obs
PRINTF,1,'File with inverted profiles = ',file_inv_profiles
PRINTF,1,'File with inverted parameters = ',file_inv_parameters
PRINTF,1,'File with errors in inverted parameters = ',file_error_inv

PRINTF,1,' '
PRINTF,1,'[Working mode]'
PRINTF,1,'Action = ',action
PRINTF,1,'Verbose = ',verbose
PRINTF,1,'Linear system solver = ',linear_solver
PRINTF,1,'Stopping volume for DIRECT = ',stopping_volume

PRINTF,1,' '
PRINTF,1,'[General parameters]'
PRINTF,1,'Synthesis mode = ',synthesis_mode
PRINTF,1,'Include stimulated emission = ',stimulated_emission
PRINTF,1,'Include magnetic field = ',magnetic_field
PRINTF,1,'Include Paschen-Back effect = ',Paschen_Back
PRINTF,1,'Include atomic level polarization = ',atomic_polarization
PRINTF,1,'Include magneto-optical effects in the RT = ',magneto_optical
PRINTF,1,'Include stimulated emission in the RT = ',stimulated_emission_RT
PRINTF,1,'Multiplet = ',multiplet
PRINTF,1,'Line-of-sight angles = ',LOS_angles
PRINTF,1,'Wavelength axis = ',Wavelength_axis

PRINTF,1,' '
PRINTF,1,'#####################'
PRINTF,1,'# General information'
PRINTF,1,'#####################'
PRINTF,1,'[Synthesis]'
PRINTF,1,'Number of slabs = ',number_slabs
PRINTF,1,'Boundary condition = ',boundary_condition
PRINTF,1,'a = ',aa
PRINTF,1,'height = ',height
PRINTF,1,'beta = ',beta
PRINTF,1,'ff = ',ff


PRINTF,1,'	[[Slab 1]]'
PRINTF,1,'	B = ',BB1
PRINTF,1,'	thetaB = ',thetaB1
PRINTF,1,'	chiB = ',chiB1
PRINTF,1,'	vdopp = ',vdopp1
PRINTF,1,'	tau = ',tau1
PRINTF,1,'	vmac = ',vmac1

PRINTF,1,'	[[Slab 2]]'
PRINTF,1,'	B = ',BB2
PRINTF,1,'	thetaB = ',thetaB2
PRINTF,1,'	chiB = ',chiB2
PRINTF,1,'	vdopp = ',vdopp2
PRINTF,1,'	tau = ',tau2
PRINTF,1,'	vmac = ',vmac2




PRINTF,1,' '
PRINTF,1,'#####################'
PRINTF,1,'# Ranges for the DIRECT method [min, max]'
PRINTF,1,'#####################'


PRINTF,1,'[Ranges]'
PRINTF,1,'a = ',aa_ranges
PRINTF,1,'beta = ',beta_ranges
PRINTF,1,'ff = ',ff_ranges


PRINTF,1,'	[[Slab 1]]'
PRINTF,1,'	B = ',BB1_ranges
PRINTF,1,'	thetaB = ',thetaB1_ranges
PRINTF,1,'	chiB = ',chiB1_ranges
PRINTF,1,'	vdopp = ',vdopp1_ranges
PRINTF,1,'	tau = ',tau1_ranges
PRINTF,1,'	vmac = ',vmac1_ranges

PRINTF,1,'	[[Slab 2]]'
PRINTF,1,'	B = ',BB2_ranges
PRINTF,1,'	thetaB = ',thetaB2_ranges
PRINTF,1,'	chiB = ',chiB2_ranges
PRINTF,1,'	vdopp = ',vdopp2_ranges
PRINTF,1,'	tau = ',tau2_ranges
PRINTF,1,'	vmac = ',vmac2_ranges



PRINTF,1,' '
PRINTF,1,'#####################'
PRINTF,1,'# Parameters to invert'
PRINTF,1,'#####################'
PRINTF,1,'[Inversion]'
PRINTF,1,'Iterations in LM = ',iterationsLM
PRINTF,1,'Number of cycles = ',ncycles
PRINTF,1,'Inversion modes = ',inversion_modes



PRINTF,1,'	[[Cycles]]'
PRINTF,1,'	a = ',aa_cycles
PRINTF,1,'	beta = ',beta_cycles
PRINTF,1,'	ff = ',ff_cycles

PRINTF,1,'		[[[Slab 1]]]'
PRINTF,1,'		B = ',BB1_cycles
PRINTF,1,'		thetaB = ',thetaB1_cycles
PRINTF,1,'		chiB = ',chiB1_cycles
PRINTF,1,'		vdopp = ',vdopp1_cycles
PRINTF,1,'		tau = ',tau1_cycles
PRINTF,1,'		vmac = ',vmac1_cycles

PRINTF,1,'		[[[Slab 2]]]'
PRINTF,1,'		B = ',BB2_cycles
PRINTF,1,'		thetaB = ',thetaB2_cycles
PRINTF,1,'		chiB = ',chiB2_cycles
PRINTF,1,'		vdopp = ',vdopp2_cycles
PRINTF,1,'		tau = ',tau2_cycles
PRINTF,1,'		vmac = ',vmac2_cycles
  



PRINTF,1,'	[[Weights]]'
PRINTF,1,'		Stokes I = ',Weights_I
PRINTF,1,'		Stokes Q = ',Weights_Q
PRINTF,1,'		Stokes U = ',Weights_U
PRINTF,1,'		Stokes V = ',Weights_V








  CLOSE,1  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

end