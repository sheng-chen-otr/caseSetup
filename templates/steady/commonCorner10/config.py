#This is a config file where you can turn on and off certain settings for the post processing
#do not make comments on the same line as any of the variables!

#### POST PRO FUNCTIONS ####
BODY_GEOM = True
BODY_CP = True
BODY_WSS = True
BODY_UMEANNEAR = True
BODY_CTPISO = True
BODY_QISO = True
SLICES = False
SLICE_LIC = True #Will only run if SLICES = True and CTP_SLICES = True
CP_SLICES = True
CTP_SLICES = True
PERCENT_U_SLICES = True


X_SLICES = True
Y_SLICES = True
Z_SLICES = True

#### RESOLUTIONS ####
HIGH_RES = True
LIC_HIGH_RES = False

#### COEFFICIENT OF PRESSURE ####
CP_RANGE = [-1,1] #DEFAULT = [-1,1]

#### COEFFICIENT OF TOTAL PRESSURE ####
CTP_RANGE = [-1.5,1] #DEFAULT = [-1.5, 1]

#### PERCENT U RANGE ####
PER_U_RANGE = [0,200]


#### PERCENT U RANGE ####
UMEANNEAR_RANGE = [0,30]

#### WSS RANGE ####
WSS_RANGE = [0,2]

#### UMEAN RANGE ###
UMEAN_RANGE = [0,30]

#### VORTICITY RANGE ###
VORTICITY_RANGE = [0,1000]


