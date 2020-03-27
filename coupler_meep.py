
# Import meep and mpb (from meep)
import meep as mp
from meep import mpb

# arrays
import numpy as np

# plotting
import matplotlib.pyplot as plt

ring_radius = 5 # um
ring_width = 0.5 # um
coupler_gap = 0.2 # um

res = 128        # pixels/μm

gdsII_file = 'coupler_zeropdk.gds'
CELL_LAYER = 0
SOURCE_LAYER = 10
Si_LAYER = 1
PORT1_LAYER = 20
PORT2_LAYER = 21
PORT3_LAYER = 22
PORT4_LAYER = 23

t_oxide = 1.0
t_Si = 0.22
t_SiO2 = 0.78

dpml = 1
cell_thickness = dpml+t_oxide+t_Si+t_SiO2+dpml
si_zmin = 0

oxide = mp.Medium(epsilon=2.25)
silicon=mp.Medium(epsilon=12)

lcen = 1.55
fcen = 1/lcen
df = 0.2*fcen

cell_zmax =  0
cell_zmin =  0
si_zmax = 0

# read cell size, volumes for source region and flux monitors,
# and coupler geometry from GDSII file

# Clear MEEP memory (otherwise will just add stuff over old stuff)
try: del si_layer
except NameError: x = None
    
# Manual fix
# geometry = []
# geometry.append(mp.Cylinder(material=silicon, center=mp.Vector3(0,-1*ring_radius,0), radius=ring_radius + ring_width/2, height=0))
# geometry.append(mp.Cylinder(material=oxide, center=mp.Vector3(0,-1*ring_radius,0), radius=ring_radius - ring_width/2, height=0))
# geometry.append(mp.Block(material=oxide, center=mp.Vector3(-0.5*ring_radius,0,0), size=mp.Vector3(ring_radius,2*ring_radius,0)))
# geometry.append(mp.Block(material=oxide, center=mp.Vector3(ring_radius,-1.5*ring_radius,0), size=mp.Vector3(2*ring_radius,ring_radius,0)))
    
si_layer = mp.get_GDSII_prisms(silicon, gdsII_file, Si_LAYER, si_zmin, si_zmax)

# # Later objects get priority : fix
final_geometry = []
# for fix in geometry:
#     final_geometry.append(fix)
for fix in si_layer:
    final_geometry.append(fix)

cell = mp.GDSII_vol(gdsII_file, CELL_LAYER, cell_zmin, cell_zmax)
src_vol = mp.GDSII_vol(gdsII_file, SOURCE_LAYER, si_zmin, si_zmax)
p1 = mp.GDSII_vol(gdsII_file, 20, si_zmin, si_zmax)
p2 = mp.GDSII_vol(gdsII_file, 21, si_zmin, si_zmax)
p3 = mp.GDSII_vol(gdsII_file, 22, si_zmin, si_zmax)
p4 = mp.GDSII_vol(gdsII_file, 23, si_zmin, si_zmax)


sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                              size=src_vol.size,
                              center=src_vol.center,
                              eig_band=1,
                              eig_parity=mp.EVEN_Y+mp.ODD_Z,
                              eig_match_freq=True)]

# Display simulation object
sim = mp.Simulation(resolution=res,
                    default_material=oxide,
                    eps_averaging=False,
                    subpixel_maxeval=1,
                    subpixel_tol=1,
                    cell_size=cell.size,
                    boundary_layers=[mp.PML(dpml)],
                    sources=sources,
                    geometry=final_geometry,
                    geometry_center=mp.Vector3(ring_radius/2, -ring_radius/2))

# Could add monitors at many frequencies by looping over fcen
# Means one FDTD for many results!
mode1 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p1))
mode2 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p2))
mode3 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p3))
mode4 = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(volume=p4))

# Setup and run the simulation
f = plt.figure(dpi=100)
animate = mp.Animate2D(sim,mp.Ez,f=f,normalize=True)
sim.run(mp.at_every(1,animate),until_after_sources=100)
#sim.run(until_after_sources=100)
plt.close()

# Do the analysis we want
# S parameters
print(sim.get_eigenmode_coefficients(mode1, [1], eig_parity=mp.NO_PARITY))
print(sim.get_eigenmode_coefficients(mode2, [1], eig_parity=mp.NO_PARITY))

# Save a video
filename = 'media/coupler.mp4'
animate.to_mp4(10,filename)
