# Standard modules
import os, time, sys, platform

# Numeca modules
from fidelity import project
from fidelity import gui
from fidelity import geometry
from fidelity import domain
from fidelity import meshing
from fidelity import material
from fidelity import simulation
from fidelity import analysis
from fidelity.plugins import hexstream_parameters
from fidelity.base import Point

path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]


# ================================================================================
#  Functions
# ================================================================================

def sortSurfaces(geomAssy):

	#defining the strings that mark the levels
	LVL9 = ['lev9']
	LVL10 = ['fine','lev10']
	LVL11 = ['extra-fine','lev11','super-fine']
	
	#initializing the assembly level lists
	assyLevel8 = []
	assyLevel9 = []
	assyLevel10 = []
	assyLevel11 = []
	#getting surfaces in the assembly
	
	#checking if geom assy is a list, if list it will loop through, if not it will just continue as regular
	#list may occur when multiple geometry files are grouped together
	
	if type(geomAssy) is list:
		for geom in geomAssy:
			surfaces = geom.get_surfaces()	
			assyList.append(geom)
			geomList.append(geom)
			for surf in surfaces:
				surfName = surf.get_name()
				if any(x in surfName for x in LVL11):
					print("	Adding %s to level 11 group..." % (surfName))
					globalLev11.append(surf)
					
				elif any(x in surfName for x in LVL10):
					print("	Adding %s to level 10 group..." % (surfName))
					globalLev10.append(surf)
					
				elif any(x in surfName for x in LVL9):
					print("	Adding %s to level 9 group..." % (surfName))
					globalLev9.append(surf)
				else:
					print("	Adding %s to level 8 group..." % (surfName))
					globalLev8.append(surf)
	else:
		assyList.append(geomAssy)
		geomList.append(geomAssy)
		surfaces = geomAssy.get_surfaces()
		for surf in surfaces:
			surfName = surf.get_name()
			if any(x in surfName for x in LVL11):
				print("	Adding %s to level 11 group..." % (surfName))
				globalLev11.append(surf)
				
			elif any(x in surfName for x in LVL10):
				print("	Adding %s to level 10 group..." % (surfName))
				globalLev10.append(surf)
				
			elif any(x in surfName for x in LVL9):
				print("	Adding %s to level 9 group..." % (surfName))
				globalLev9.append(surf)
			else:
				print("	Adding %s to level 8 group..." % (surfName))
				globalLev8.append(surf)
	
	
		





print("\n\n\n\n\n")
print("##################################################################################")
print("##################################################################################")
print("##################################################################################")
print("OTR FIDELITY MESHING - EXTERNAL AERO - FULL CAR					 ")
print("##################################################################################")
print("##################################################################################")
print("##################################################################################")
print("\n\n")

print("Opening DBS template...")
project.open_template('DBS')
print("Saving file as fidelityMesh_%s..." % (case))
project_file="%s/%s/fidelityMesh_%s" % (path,case,case)
project.save_project_as(project_file)

# ================================================================================
#  Geometry
# ================================================================================

print("Getting geometry...")

gList = os.listdir("%s/%s/constant/triSurface/" % (path,case))

fw = []
rw = []
body = []
fr = []
rr = []
#grouping all the similar parts into one group
for part in gList:
	print("	Found geometry: %s" %(part))
	partPath = "%s/%s/constant/triSurface/%s" % (path,case,part)
	if "fw" in part:
		print("		Adding to fw assembly")
		fw.append(partPath)
	elif "rw" in part:
		print("		Adding to rw assembly")
		rw.append(partPath)
	elif "trial" in part:
		print("		Adding to body assembly")
		body.append(partPath)
	elif "fr" in part:
		print("		Adding to fr assembly")
		fr.append(partPath)
	elif "rr" in part:
		print("		Adding to rr assembly")
		rr.append(partPath)
	else:
		print("		WARNING! %s is not identified, skipping..." % (part))

r_options = geometry.ImportResourceOptions()
r_options._threads = [4,4]
r_options._parallel = True
assyList = []
geomList = []

#initializing global refinement arrays
globalLev8 = []
globalLev9 = []
globalLev10 = []
globalLev11 = []

print("\n")
#putting all similar geometries into one assembly
if len(fw) > 0:
	print("Creating front wing assembly...")
	frontWing = geometry.import_geometry(fw,resource_options = r_options)
	sortSurfaces(frontWing)
	
if len(rw) > 0:
	print("Creating rear wing assembly...")
	rearWing = geometry.import_geometry(rw,resource_options = r_options)
	sortSurfaces(rearWing)
	
if len(body) > 0:
	print("Creating body assembly...")
	bodyAssy = geometry.import_geometry(body,resource_options = r_options)
	sortSurfaces(bodyAssy)
	
if len(rr) > 0:
	print("Creating rear wheel assembly...")
	rearWheels = geometry.import_geometry(rr,resource_options = r_options)
	sortSurfaces(rearWheels)
	
if len(fr) > 0:
	print("Creating front wheel assembly...")
	frontWheels = geometry.import_geometry(fr,resource_options = r_options)
	sortSurfaces(frontWheels)

print("\n\n")
print("Grouped surfacees:")
print("	Level 8 Surfaces:")
for surf in globalLev8:	
	print("		" + surf.get_name())
	
print("	Level 9 Surfaces:")	
for surf in globalLev9:
	print("		" + surf.get_name())
	
print("	Level 10 Surfaces:")
for surf in globalLev10:
	print("		" + surf.get_name())
	
print("	Level 11 Surfaces:")	
for surf in globalLev11:
	print("		" + surf.get_name())

#Creating external bounding box
print("Creating external bounding box...")
extBox = geometry.create_box([-16,-16,0], [32,16,16], 'boundingBox')
extBox.get_boundary('Back').rename('x-min')
extBox.get_boundary('Front').rename('x-max')
extBox.get_boundary('Left').rename('y-min')
extBox.get_boundary('Right').rename('y-max')
extBox.get_boundary('Bottom').rename('z-min')
extBox.get_boundary('Top').rename('z-max')
assyList.append(extBox)

#creating wake boxes
print("Creating wake boxes...")
wb1 = geometry.create_box([-1.5,-1,-0.29], [3.5,1,1.5], 'wake1')
wb2 = geometry.create_box([-2.5,-1.4,-0.29], [5,1.4,2], 'wake2')
wb3 = geometry.create_box([-3,-2,-0.29], [7,2,2.4], 'wake3')
wb4 = geometry.create_box([-5,-3,-0.29], [100,3,3], 'wake4')
# ================================================================================
#  Domain Creation
# ================================================================================

print("Creating domain: fluidDomain")
fD = domain.create_domain(assyList,'fluidDomain')
fD.set_seed_point(-1,-1,4)
print("	Setting domain boundary conditions...")
fD.set_boundary_type('boundingBox|x-min','Inlet')
fD.set_boundary_type('boundingBox|x-max','Outlet')
fD.set_boundary_type('boundingBox|y-min','Inlet')
fD.set_boundary_type('boundingBox|y-max','Outlet')
fD.set_boundary_type('boundingBox|z-min','Wall')
fD.set_boundary_type('boundingBox|z-max','Wall')

# ================================================================================
#  Mesh Setups
# ================================================================================

print("Creating mesh setup instance...")
mesh_setup = meshing.create_mesh_setup(fD, 'meshSetup', mesher_name='Hexpress')
mesh_setup.set_cells_type('Full hexa')
mesh_setup.set_initial_cell_size(1.6)
mesh_setup.set_far_field_treatment('Auto')
mesh_setup.set_buffer_insertion('Cross through', 130)
mesh_setup.set_initial_mesh_type('Cartesian')
#mesh_setup.set_close_viscous_layers_based_on_openfoam_skewness(True, value=4.0)
#mesh_setup.enable_compute_OpenFOAM_quality_criteria()
mesh_setup.enable_compute_Fidelity_quality_criteria()
mesh_setup.enable_capturing_of_corners_of_degree_3_4()
#mesh_setup.set_improve_quality_based_on(True,value='OpenFOAM')
#mesh_setup.set_improve_quality_based_on(True,value='PBS', allow_deviation_from_geometry=False)
mesh_setup.set_auto_edge_proximity_mode('Fine')
mesh_setup.set_auto_edge_proximity(True,min_cell_size=0.001,diffusion=2)

#mesh_setup.enable_improve_sharp_angle_narrow_gap() #Doesn't work for Full Hexa


#creating fluid domain
mesh_domain = mesh_setup.get_domain('fluidDomain')

#setting up inflation layers
print("Setting up inflation layers... \n\n")
mesh_domain.set_viscous_layers_enlarge('boundingBox|z-min', 
                                        True, 
                                        insertion_method = 'Inflation',
                                        first_layer_thickness=0.0007, 
                                        stretching_ratio=1.3,
                                        num_layers=5)
if not geomList == []:
	mesh_domain.set_viscous_layers_enlarge(geomList, 
                                        True, 
                                        insertion_method = 'Inflation',
                                        first_layer_thickness=0.0007, 
                                        stretching_ratio=1.3,
                                        num_layers=5)


#surface refinements
print("	Setting up refinements...")
print("		Setting refinement for level 8...")
if len(globalLev8) > 0:
	mesh_domain.set_uniform_refinement(globalLev8, True, 'Refinement level', 8)

print("		Setting refinement for level 9...")
if len(globalLev9) > 0:
	mesh_domain.set_uniform_refinement(globalLev9, True, 'Refinement level', 9)

print("		Setting refinement for level 10...")
if len(globalLev10) > 0:
	mesh_domain.set_uniform_refinement(globalLev10, True, 'Refinement level', 9)

print("		Setting refinement for level 11...")
if len(globalLev11) > 0:
	mesh_domain.set_uniform_refinement(globalLev11, True, 'Refinement level', 10)

mesh_domain.set_curvature_refinement(geomList,True,max_division=5,min_cell_size=0.001,diffusion=2)

#edge capture settings
mesh_domain.enable_surface_edges_capture(geomList)
mesh_domain.enable_viscous_layers_stop_propagation(assyList)
mesh_domain.enable_improve_capturing(geomList)
mesh_domain.enable_boundary_edges_capture(geomList)
mesh_domain.set_overwrite_buffer_angle('boundingBox|y-max',True,angle=175)

# refinement box
print("	Setting refinement for wake 1 BOI...")
wakeVol_1 = mesh_setup.add_refinement_volume(wb1)
wakeVol_1.set_refinement_region('Volume')
wakeVol_1.set_isotropic_refinement('Refinement level', 7, diffusion=2)

print("	Setting refinement for wake 2 BOI...")
wakeVol_2 = mesh_setup.add_refinement_volume(wb2)
wakeVol_2.set_refinement_region('Volume')
wakeVol_2.set_isotropic_refinement('Refinement level', 6, diffusion=2)

print("	Setting refinement for wake 3 BOI...")
wakeVol_3 = mesh_setup.add_refinement_volume(wb3)
wakeVol_3.set_refinement_region('Volume')
wakeVol_3.set_isotropic_refinement('Refinement level', 5, diffusion=2)

print("	Setting refinement for wake 4 BOI...")
wakeVol_4 = mesh_setup.add_refinement_volume(wb4)
wakeVol_4.set_refinement_region('Volume')
wakeVol_4.set_isotropic_refinement('Refinement level', 4, diffusion=2)
			
# # # ================================================================================
# # #  Running the Mesher
# # # ================================================================================
domainsInProject = domain.get_number_of_domains()
print("The number of domains in this project is: %s" % (domainsInProject))
numMeshSetups = meshing.get_number_of_mesh_setup()
print("The number of mesh setups is: %s" % (numMeshSetups))

print("Creating mesh!")
mesh_setup.start_mesh_generation(num_proc=32)
if mesh_domain.mesh_exists_for_domain():
	mesh_setup.print_mesh_statistics()
else:
	print("Error in mesh generation!")


# # ================================================================================
# #  Saving Project
# # ================================================================================
print(	"Saving project...")
project.save_project()


# # ================================================================================
# #  Exporting Mesh
# # ================================================================================

if mesh_domain.mesh_exists_for_domain():
	mesh_setup.print_mesh_statistics()
	for key in mesh_setup.get_mesh_quality_field_keys():
		print(key)
		print(mesh_setup.get_mesh_quality_statistics(key))
else:
	print("No mesh to export!")
	
print(	"Saving project...")
project.save_project()
				
		






