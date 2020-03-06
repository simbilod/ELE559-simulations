from __future__ import division

import meep as mp
from meep import mpb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-res', type=int, default=64, help='resolution (default: 64 pixels/um)')
parser.add_argument('-w', type=float, default=0.5, help='branch separation (default: 0.5 um)')
parser.add_argument('-wavelength', type=float, default=0.5, help='wavelength for fields (default: 1.55 um)')
args = parser.parse_args()

w = args.w

h = 0.22  # Si height (um)

Si = mp.Medium(index=3.45)
SiO2 = mp.Medium(index=1.45)

# Define the computational cell.  We'll make x the propagation direction.
# the other cell sizes should be big enough so that the boundaries are
# far away from the mode field.
sc_y = 2  # supercell width (um)
sc_z = 2  # supercell height (um)
geometry_lattice = mp.Lattice(size=mp.Vector3(0, sc_y, sc_z))

# define the 2d blocks for the strip and substrate
geometry = [mp.Block(size=mp.Vector3(mp.inf, mp.inf, 0.5 * (sc_z - h)),
                     center=mp.Vector3(z=0.25 * (sc_z + h)), material=SiO2),
            mp.Block(size=mp.Vector3(mp.inf, w, h), material=Si)]

# The k (i.e. beta, i.e. propagation constant) points to look at, in
# units of 2*pi/um.  We'll look at num_k points from k_min to k_max.
num_k = 9
k_min = 0.1
k_max = 3.0
k_points = mp.interpolate(num_k, [mp.Vector3(k_min), mp.Vector3(k_max)])

resolution = args.res  # pixels/um

# Increase this to see more modes.  (The guided ones are the ones below the
# light line, i.e. those with frequencies < kmag / 1.45, where kmag
# is the corresponding column in the output if you grep for "freqs:".)
num_bands = 4

filename_prefix = 'strip-script-'  # use this prefix for output files

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands,
    filename_prefix=filename_prefix
)


ms.run(mpb.display_yparities, mpb.display_zparities)
omega = 1 / args.wavelength # frequency corresponding to 1.55um

# Output the x component of the Poynting vector for num_bands bands at omega
ms.find_k(mp.NO_PARITY, omega, 1, num_bands, mp.Vector3(1), 1e-3, omega * 3.45,
          omega * 0.1, omega * 4, mpb.output_poynting_x, mpb.display_yparities,
          mpb.display_group_velocities)