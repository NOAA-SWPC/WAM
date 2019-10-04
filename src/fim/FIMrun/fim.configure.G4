
#TBH:  Note that settings in this file work for a 
#TBH:  G5 24-hour run.  We need a way to generate this file 
#TBH:  automatically from FIMnamelist.  And a way to set 
#TBH:  "dt_int" and its duplicate "deltim", during the model 
#TBH:  run since FIM computes these on-the-fly.  .  

#TBH:  I have commented out any field not needed by FIM.  
#TBH:  Eventually remove them after debugging.  

###############################
######## Debug Prints #########
###############################
# turns on "iprint" in module_ERR_MSG.F90 for loquacious prints
iprint: false

###############################
######## Model Core ###########
###############################

core: fim                   # The dynamic core options:
                                # nmm
                                # gfs
                                # fim
                                # arw
                                                                                                                                              
                                                                                                                                              
###############################
##### Grid Specifications #####
###############################

#im:  230                   # I gridpoints
#jm:  141                   # J gridpoints
#lm:  60                    # Number of atmospheric layers

#tph0d:  40.                # Central geographic latitude of grid (degrees)
#tlm0d: -100.               # Central geographic longitude of grid (degrees, positive east)

#wbd:  -22.9                # Grid's western boundary (rotated degrees)
#sbd:  -14.0                # Grid's southern boundary (rotated degrees)


##############################
##### Grid Decomposition #####
##############################

#inpes:  04                 # Number of compute tasks in the I direction
#jnpes:  03                 # Number of compute tasks in the J direction


###############################
#### Specify the I/O tasks ####
###############################
                                                                                                                                             
#quilting:             false    # Do you want asynchronous quilting/history writes?
#read_groups:           0
#read_tasks_per_group:  0
#write_groups:          3
#write_tasks_per_group: 1


##########################################
##### Fundamental Timestep (seconds) #####
##########################################

dt_int: 720                 # Integer seconds
dt_num: 0                   # Numerator of fractional second
dt_den: 1                   # Denominator of fractional second


###############################
##### Forecast Start Time #####
###############################

start_year:   2007   
###start_month:  04
###start_day:    15
###start_hour:   12
start_month:  07
start_day:    17
start_hour:   00
start_minute: 0
start_second: 0
                                                                                                                                              
tstart:       0.            # The starting forecast hour
                                                                                                                                              
                                                                                                                                              
###########################
##### Forecast Length #####
###########################

nhours_fcst:     24          # Length of the forecast (hours)

#TBH:  this replaces nhours_history in nems r11258
#nhours_history: 01          # Frequency of history output (hours)
nfhout:          1
#TBH:  specify output frequency in time steps if nsout>0
nsout:           0
#TBH:  seems redundant with dt_int, see job/exglobal_fcst.sh.sms_nems
deltim: 720.0                # real seconds, deltim="${dt_int}.0"

#restart:        false       # True--> A restarted run


#########################
##### General modes #####
#########################

#global:  true               # True--> Global ; False--> Regional
#hydro:   false              # True--> Hydrostatic ;  False--> Nonhydrostatic
#JR adiabatic true not yet implemented
adiabatic: false            # True--> Adiabatic ;  False--> Diabatic


##############################
##### Lateral Boundaries #####
##############################

#specified:  true
#nested:     false


####################
##### Dynamics #####
####################
 
#secadv:  true               # True--> 2nd order advection ; False--> 4th order advection
                                                  
#secdif:  true               # True--> 2nd order diffusion ; False--> 4th order diffusion
#smag2:   0.2                # Smagorinsky constant for 2nd order diffusion
#smag4:   1.3                # Smagorinsky constant for 4th order diffusion
                                        
#codamp:  9.0                # Divergence damping constant
#wcor:    0.0075             # Divergence correction factor
 
#pwrc:    1.0                # Fourier filter factor
                
#idtad:   2                  # Number of adjustment timesteps between passive advection calls
                                                                                                                                              
#advect_tracers:  false      # Will tracers be advected?
#idtadt:  4                  # Number of adjustment timesteps between tracer advection calls (normally 1)

#num_tracers_met:   4        # Number of specified meterological "tracer" scalars (e.g., water)
#num_tracers_chem:  0        # Number of specified chem/aerosol "tracer" scalars
                 
########################################
##### Read/Write Global Summations #####
########################################

#read_global_sums:  false    # Read in global summations or not
#write_global_sums: false    # Write out global summations or not


##################################
### Precipitation Assimilation####
##################################

#pcpflg:  false              # True--> assimilation on ; False--> assimilation off


##################################
##### Physics Scheme Options #####
##################################
#---------------------------------
#shortwave: gfdl             # Radiation schemes
#longwave:  gfdl        
                                # gfdl--> GFDL (shortwave and longwave)
                                # gsfc--> Goddard shortwave
                                # dudh--> Dudhia shortwave
                                # rrtm--> RRTM longwave

#co2tf:  1                      # co2tf--> Read (0) or generate internally (1)
                                #          the GFDL CO2 transmission functions

#---------------------------------
#convection:  bmj           # Convection scheme
                                # bmj--> Betts-Miller-Janjic
                                #  kf--> Kain-Fritsch
                                # sas--> Simplified Arakawa-Schubert
                                #  gd--> Grell-Devenyi

#---------------------------------
#microphysics:  fer         # Microphysics scheme
                                # fer--> Ferrier
                                # kes--> Kessler
                                # lin--> Lin
                                # tho--> Thompson
                                                                                                                                              
#---------------------------------
#turbulence:  myj           # Turbulence schemes
                                # fer--> Ferrier
                                # kes--> Kessler
                                # lin--> Lin
                                # tho--> Thompson
                                                                                                                                              
#---------------------------------
                                                                                                                                              
#sfc_layer:  myj             # Surface layer schemes
                                # myj--> Mellor-Yamada-Jamnjic
                                # mm5--> MM5 sfc layer
                                                                                                                                              
#---------------------------------
                                                                                                                                              
#land_surface:  noah        # Land surface schemes
                                # nmm-->  LSM NMM
                                # noah--> NOAH scheme

#ucmcall: 0                  # Use (1) or do not use (0) the Urban Canopy Model

#---------------------------------
                                                                                                                                              
#gfs:  false                 # Select entire GFS physics suite


#############################
##### Physics Timesteps #####
#############################

#nrads:    100               # Number of dynamics timesteps between calls to shortwave
#nradl:    100               # Number of dynamics timesteps  between calls to longwave

#nphs:       4               # Number of dynamics timesteps between calls to landsurface and turbulence

#nprecip:    4               # Number of dynamics timesteps between calls to convection and microphysics

#nhrs_prec:  6               # Frequency in hours between times when precip arrays are emptied
#nhrs_heat:  3               # Frequency in hours between times when heating arrays are emptied
#nhrs_clod:  3               # Frequency in hours between times when cloud arrays are emptied
#nhrs_rdlw:  3               # Frequency in hours between times when LW radiation arrays are emptied
#nhrs_rdsw:  3               # Frequency in hours between times when SW radiation arrays are emptied
#nhrs_srfc:  3               # Frequency in hours between times when sfc evap/flux arrays are emptied


#######################################
####  Write History Specifications ####
#######################################
 
#write_flag:              .TRUE.       # True--> Yes write history files
 
#io_unit:                 -999
#io_file:                 'DEFERRED'
 
#-------------------------------------------------
# If IO_FILE == 'DEFERRED', use FILENAME_BASE
# as the first part of the history file names.
#-------------------------------------------------
#filename_base:           'nmm_b_history.'
#io_status:               'REPLACE'
#io_access:               'SEQUENTIAL'
#io_form:                 'UNFORMATTED'
#io_recl:                 100
#io_position:             ' '
#io_action:               'WRITE'
#io_delim:                ' '
#io_pad:                  ' '

#######################################
####  timing info                  ####
#######################################
nhours_clocktime:   1        # Frequency in hours between clocktime diagnostic prints.
npe_print:          1        # The MPI task that will provide the clocktimes.


#
# New config flags from nems r11258
# TODO:  eliminate obsolete flags above...
#

##########################
####  Ensemble flags  ####
##########################

ENS_SPS:            .FALSE.
total_member:       1
RUN_CONTINUE:       .FALSE.


