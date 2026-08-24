import os
import numpy as np
import fileinput
import os.path
import sys
import fileinput
import glob
from bcParser_v1_0 import bcParser
from matplotlib import pyplot as plt
ver=1.0

print('#### forceBinPlot version %s ####' % (ver))


path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]

if len(sys.argv) > 1:
	comp = sys.argv[1]
	print("Comparing Trial%s with Trial %s" % (case,comp))
	if "half" in comp:
		scale2 = 2
	else:
		scale2 = 1
else:
	comp = ""
	pass
COMPARE = False

if "half" in case:
	scale = 2
else:
	scale = 1





[inletMag,lastTime,yaw,wheelRotation,simType,turbModel] = bcParser(path,case)

surfaceLoc = glob.glob("%s/%s/postProcessing/binForceCoeffs/*" % (path,case))
lastTime = os.path.basename(surfaceLoc[0])
#print(lastTime)

#checking for bin coeffs

if not os.path.isdir("%s/%s/postProcessing/binForceCoeffs/%s" % (path,case,lastTime)):
	sys.exit("Force bin coefficients NOT found! Exiting!")


#reads coefficient.dat file
coeffs=np.loadtxt('%s/%s/postProcessing/binForceCoeffs/%s/forceCoeffBin.dat' % (path,case,lastTime), dtype='float', comments='#', delimiter=None, converters=None, skiprows=10, unpack=False, ndmin=0)




forceMultiplier = (inletMag ** 2) * 0.5
coeffs = coeffs[1:]
totalCoeffs = np.zeros(300)
xCoeffs = np.zeros(100)
yCoeffs = np.zeros(100)
zCoeffs = np.zeros(100)
xForce = np.zeros(100)
yForce = np.zeros(100)
zForce = np.zeros(100)

n = list(range(100))
#print(coeffs.shape)

for i in n:
	start = (9*i)+0
	end = (9*i)+3
	totalCoeffs[(3*i):(3*i)+3] = coeffs[start:end] * scale
	xCoeffs[i] = totalCoeffs[3*i] 
	yCoeffs[i] = totalCoeffs[(3*i)+1] 
	zCoeffs[i] = totalCoeffs[(3*i)+2] * -1
	xForce[i] = xCoeffs[i] * forceMultiplier * 1.225
	yForce[i] = yCoeffs[i] * forceMultiplier * 1.225
	zForce[i] = zCoeffs[i] * forceMultiplier * 1.225
	

#### Comparison Trial ####
if comp == "":
	pass

else:
	
	if os.path.isdir("%s/%s" % (path,comp)):
		[inletMag2,lastTime2,yaw2,wheelRotation2,simType2,turbModel2] = bcParser(path,comp)
		#print(lastTime2)
		#print("%s/%s/postProcessing/binForcesCoeffs/%s/forceCoeffBin.dat" % (path,comp,lastTime2))
		surfaceLoc2 = glob.glob("%s/%s/postProcessing/binForceCoeffs/*" % (path,comp))
		lastTime2 = os.path.basename(surfaceLoc2[0])

		if os.path.isfile("%s/%s/postProcessing/binForceCoeffs/%s/forceCoeffBin.dat" % (path,comp,lastTime2)):
			
			COMPARE = True
			coeffs2=np.loadtxt('%s/%s/postProcessing/binForceCoeffs/%s/forceCoeffBin.dat' % (path,comp,lastTime2), dtype='float', comments='#', delimiter=None, converters=None, skiprows=10, unpack=False, ndmin=0)
			forceMultiplier2 = (inletMag2 ** 2) * 0.5
			coeffs2 = coeffs2[1:]
			totalCoeffs2 = np.zeros(300)
			xCoeffs2 = np.zeros(100)
			yCoeffs2 = np.zeros(100)
			zCoeffs2 = np.zeros(100)
			xForce2 = np.zeros(100)
			yForce2 = np.zeros(100)
			zForce2 = np.zeros(100)

			n = list(range(100))
			

			for i in n:
				start = (9*i)+0
				end = (9*i)+3
				totalCoeffs2[(3*i):(3*i)+3] = coeffs2[start:end] * scale2
				xCoeffs2[i] = totalCoeffs2[3*i] 
				yCoeffs2[i] = totalCoeffs2[(3*i)+1] 
				zCoeffs2[i] = totalCoeffs2[(3*i)+2] * -1
				xForce2[i] = xCoeffs2[i] * forceMultiplier2 * 1.225
				yForce2[i] = yCoeffs2[i] * forceMultiplier2 * 1.225
				zForce2[i] = zCoeffs2[i] * forceMultiplier2 * 1.225
		else:
			sys.exit("Compared trial does not have binned forces! Exiting!")
	else:
		sys.exit("There is no trial by the name %s!" % (comp))

allForces = np.vstack((n,xCoeffs,yCoeffs,zCoeffs,xForce,yForce,zForce)).T


np.savetxt("trial%s_binForces.csv" % (case), allForces, delimiter=",",header = "xCoeffs,yCoeffs,zCoeffs,xForce,yForce,zForce")


fig, ax = plt.subplots()

if os.path.isfile("%s/%s/postProcessing/images/%s_geom_left.png" % (path,case,case)):
	img = plt.imread("%s/%s/postProcessing/images/%s_geom_left.png" % (path,case,case))
else:
	img = plt.imread("/home/openfoam/openFoam/templates/postPro/binForceImage/092_half_geom_left.png")

ax.plot(zForce,'-', linewidth=1, color='firebrick',label = 'Downforce - Trial%s' % (case))

if COMPARE == True:
	ax.plot(zForce2,'--', linewidth=1, color='firebrick',label = 'Downforce - Trial%s' % (comp))

ax.plot(xForce,'-', linewidth=1, color='blue',label = 'Drag - Trial%s' % (case))
if COMPARE == True:
	ax.plot(xForce2,'--', linewidth=1, color='blue',label = 'Drag - Trial%s' % (comp))
ax.set_xlim(0,99)
ax.imshow(img,extent = [-28.5,129.5,0,int(max(zForce))*1.2],aspect = 'auto',alpha = 0.5)


#ax.set_xlim(0,99)
ax.set_ylim(-20,int(max(zForce))*1.1)
#size = fig.get_size_inches()*fig.dpi # size in pixels
ax.set_xlabel('Percent Length of Car')
ax.set_ylabel('Force (N)')
ax.set_title('Binned Forces')
ax.legend()

if COMPARE == True:
	print("Printing binned forces plot to %s_binnedForces_compare_trial%s.png" % (case,comp))
	plt.savefig('%s_binnedForces_compare_trial%s.png' % (case,comp), dpi=600)
else:
	print("Printing binned forces plot to %s_binnedForces.png" % (case))
	plt.savefig('%s_binnedForces.png' % (case), dpi=600)



	




