import numpy as np
import sys
from analytic import *

print('Example usage: python run_analytic.py theta')

# setting up model parameters
a1 = -0.5 # bottom depth of the anisotropic layer
a2 = -0.1 # top depth
mid_layer_depth = (a1 + a2)/2
w  = 1
ux0 = 1
theta = float(sys.argv[1]) # 
n1 = np.cos((theta+90)/180*np.pi)
n2 = np.sin((theta+90)/180*np.pi)
es = 0.1 # weak viscosity
e  = 1 # strong viscosity
dd = 0.1 # spatial grids along the depth profile.
mid_layer_depth = (a1 + a2)/2
id = round((mid_layer_depth-(-w))/dd+1) # counting from -w to 0, idx 8 is inside a1 ~ a2 depths of the anisotropic zone.

# compute
d, sig11, sig12, sig22, str11, str12, str22, u1, p = analytic(a1, a2, w, ux0, n1, n2, es, e, dd)

smax, smin, smax_eigvec, smin_eigvec, s_J2 = calc_principal(sig11[id-1], sig22[id-1], sig12[id-1])
srmax, srmin, srmax_eigvec, srmin_eigvec, srJ2 = calc_principal(str11[id-1], str22[id-1], str12[id-1])
srmax_i, srmin_i, _, _, srJ2_iso = calc_principal(str11[2], str22[2], str12[2])

theta_sigma_max = np.degrees(np.arctan(smax_eigvec[1]/smax_eigvec[0]))
theta_eps_min = np.degrees(np.arctan(srmin_eigvec[1]/srmin_eigvec[0]))

print(' ')
print('********* results *********')
print(f"1. Normal director to weak anisotropy direction rotates counterclockwisely from y+ at {theta:.2f} degs")
print('2. Stress tensor [sxx, syy, sxy] = ', [sig11[id-1], sig22[id-1], sig12[id-1]])
print(f"3. Principal stresses smax and smin are ", smax, smin)
print(f"4. Rotate from the x axis clockwisely {-theta_sigma_max:.2f} to get the maximum compressive s_max")
print(f"5. Principal strain rate eps_max and eps_min are ", srmax, srmin)
print(f"6. Rotate from the x axis counterclockwisely {theta_eps_min:.2f} to get the minimum/extensional eps_min")

print(f"Strain localization is {srJ2/srJ2_iso:.2f}")
