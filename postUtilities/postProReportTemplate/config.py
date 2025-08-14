#This is a config file where you can turn on and off certain settings for the post processing
#do not make comments on the same line as any of the variables!

#### POST PRO FUNCTIONS ####
BODY_GEOM = False
BODY_CP = False
BODY_WSS = False
BODY_UMEANNEAR = False
BODY_CTPISO = True
BODY_QISO = True
SLICES = False
SLICE_LIC = False #Will only run if SLICES = True and CTP_SLICES = True
CP_SLICES = False
CTP_SLICES = False
PERCENT_U_SLICES = False


X_SLICES = False
Y_SLICES = False
Z_SLICES = False

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


