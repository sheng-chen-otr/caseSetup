#This is a config file where you can turn on and off certain settings for the post processing
#do not make comments on the same line as any of the variables!


BODY_GEOM = True
BODY_CP = True
BODY_CPX = True
BODY_CPZ = True
BODY_WSS = True
BODY_UMEANNEAR = True
BODY_CTPISO = True
BODY_QISO = True
SLICES = True
MOVIES = True
SLICE_LIC = True 
#Will only run if SLICES = True and CTP_SLICES = True
CP_SLICES = True
CTP_SLICES = True
PERCENT_U_SLICES = True


X_SLICES = True
Y_SLICES = True
Z_SLICES = True

HIGH_RES = True
LIC_HIGH_RES = False

#[PLOT_RANGES]
#### COEFFICIENT OF PRESSURE ####
CP_RANGE = [-1,1] 
#DEFAULT = [-1,1]
CPX_RANGE = [-0.7,0.7] 
#DEFAULT = [-1,1]
CPZ_RANGE = [-0.7,0.7] 
#DEFAULT = [-1,1]

#### COEFFICIENT OF TOTAL PRESSURE ####
CTP_RANGE = [-0.5,1] 
#DEFAULT = [-1.5, 1]

#### PERCENT U RANGE ####
PER_U_RANGE = [0,200]


#### PERCENT U RANGE ####
UMEANNEAR_RANGE = [0,200]

#### WSS RANGE ####
WSS_RANGE = [0,2]

#### UMEAN RANGE ###
UMEAN_RANGE = [0,200]

#### VORTICITY RANGE ###
VORTICITY_RANGE = [0,1000]


#FOR VIEW SCALING, LARGER NUMBER, CAMERA FARTHER BACK
#### X View Scaling ####
XSCALE = 0.940 
# default 0.940

#### Y View Scaling ####
YSCALE = 1.13775

#### Z View Scaling ####
ZSCALE = 1.13775

#### X CENTER ####
XCENTER =  0.48676 
# x coordinate of center of wheel base

#### Z COORD ####
ZCOORD = 0.5605
