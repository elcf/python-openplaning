# OpenPlaning

[![Documentation Status](https://readthedocs.org/projects/openplaning/badge/?version=latest)](https://openplaning.readthedocs.io/en/latest/?badge=latest)

OpenPlaning is a Python library for the hydrodynamic evaluation of planing hulls based on the Savitsky empirical methods.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install openplaning.

```bash
pip install openplaning
```

## Examples

You can run the example below, plus an optimization case study, online with Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/elcf/binder-openplaning/main?filepath=OpenPlaningExamples.ipynb)

```python
from openplaning import PlaningBoat

#Vessel particulars (from the Savitsky '76 example)
speed = 13.07 #m/s
weight = 827400 #N
beam = 7.315 #m
length = 24.38 #m, vessel LOA
lcg = 10.67 #m, long. center of gravity
vcg = beam/7 #m, vert. center of gravity
r_g = 0.25*length #m, radius of gyration
beta = 15 #deg, deadrise

#Propulsion
epsilon = 0 #deg, thrust angle w.r.t. keel
vT = vcg #m, thrust vertical distance
lT = lcg #m, thrust horizontal distance

#Trim tab particulars
sigma = 1.0 #flap span-hull beam ratio
delta = 5 #deg, flap deflection
Lf = 0.3048 #m, flap chord

#Seaway
H_sig = 1.402 #m, significant wave height

#Create boat object
boat = PlaningBoat(speed, weight, beam, lcg, vcg, r_g, beta, epsilon, vT, lT, length, H_sig, Lf=Lf, sigma=sigma, delta=delta, wetted_lengths_type=3)

#Calculates the equilibrium trim and heave, and updates boat.tau and boat.z_wl
boat.get_steady_trim()

boat.print_description()
# RETURNS:
# ---VESSEL---
# Speed            13.07 m/s
# V_k              25.40808 knot
# Fn (beam)        1.543154 
# Fn (volume)      2.001392 
# 
# Weight           827400 N
# Mass             84371.75 kg
# Volume           82.24409 m³
# Beam             7.315 m
# LCG              10.67 m from stern
# VCG              1.045 m from keel
# R_g              6.095 m
# Deadrise         15 deg
# 
# LOA              24.38 m
# AHR              0.00015 m, average hull roughness
#
# ---ATTITUDE---
# z_wl             0.1384811 m, vertical distance of center of gravity to the calm water line
# tau              2.878945 deg, trim angle
# η₃               0 deg, additional heave
#η₅               0 deg, additional trim
#Transom draft    1.441111 m, draft of keel at transom
#
# ---PROPULSION---
# Thrust angle     0 deg w.r.t. keel (CCW with body-fixed origin at 9 o'clock)
# LCT              10.67 m from stern, positive forward
# VCT              1.045 m from keel, positive up
# 
# ---FLAP---
# Chord            0.3048 m
# Span/Beam        1 
# Angle            5 deg w.r.t. keel (CCW with body-fixed origin at 9 o'clock)
# 
# ---AIR DRAG---
# l_air            0 m, distance from stern to center of air pressure
# h_air            0 m, height from keel to top of square which bounds the air-drag-inducing shape
# b_air            0 m, transverse width of square which bounds the air-drag-inducing shape
# C_shape          0 area coefficient for air-drag-inducing shape. C_shape = 1 means the air drag reference area is h_air*b_air
# C_D              0.7 air drag coefficient
# 
# ---ENVIRONMENT---
# ρ                1025.87 kg/m³, water density
# ν                1.19e-06 m²/s, water kinematic viscosity
# ρ_air            1.225 kg/m³, air density
# g                9.8066 m/s², gravitational acceleration
# 
# ---WETTED LENGTH OPTIONS---
# LC_type          3 (1 = Use Faltinsen 2010 wave rise approximation, 2 = Use Savitsky's '64 approach, 3 = Use Savitsky's '76 approach)
# zmax_type        1 (1 = Uses 3rd order polynomial fit (faster, recommended), 2 = Use cubic interpolation)
# 
# ---WETTED LENGTHS---
# L_K              28.69256 m, keel wetted length
# L_C              17.67617 m, chine wetted length
# λ                3.199428 mean wetted-length to beam ratio (L_K+L_C)/(2*beam)
# x_s              11.0164 m, distance from keel/water-line intersection to start of wetted chine
# z_max            0.7704615 maximum presssure coordinate coefficient (z_max/Ut)
# 
# ---FORCES [F_x (N, +aft), F_z (N, +up), M_cg (N*m, +pitch up)]---
# Hydrodynamic Force =
# [39245.86 780400.3 301189.8]
#
# Skin Friction =
# [31893.98 -1603.929 -18962.39]
#
# Air Resistance =
# [0 0 0]
#
# Flap Force =
# [1840.949 44933.51 -282227.4]
# 
# Net Force =
# [72980.79 2.608704e-08 2.55066e-07]
# 
# Resultant Thrust =
# [-72980.79 3670.16 0]
# 
# 
# ---THURST & POWER---
# Thrust Magnitude 73073.02 N
# Effective Thrust 72980.79 N
# Eff. Power       953.859 kW
# Eff. Horsepower  1279.146 hp
# 
# ---EOM MATRICES---
# Mass matrix, [kg, kg*m/rad; kg*m, kg*m²/rad] =
# [[501800.8 67648.82]
#  [67648.82 2.046024e+07]]
# 
# Damping matrix, [kg/s, kg*m/(s*rad); kg*m/s, kg*m²/(s*rad)] =
# [[447299.8 -8182636]
#  [3078703 2.909537e+07]]
#
# Restoring matrix, [N/m, N/rad; N, N*m/rad] =
# [[1325673 -2390482]
#  [4940227 5.375431e+07]]
#
#
# ---PORPOISING---
# [[Eigenvalue check result, Est. pitch settling time (s)],
#  [Savitsky chart result, Critical trim angle (deg)]] =
# [[0 7.097941]
#  [0 9.955598]]
#
#
# ---BEHAVIOR IN WAVES---
# H_sig            1.402 m, significant wave heigth
# R_AW             38406.03 N, added resistance in waves
# Average impact acceleration [n_cg, n_bow] (g's) =
# [0.3082269 0.754686]
```

## Dependencies

* [NumPy](https://numpy.org/)
* [SciPy](https://www.scipy.org/)
* [ndmath](https://github.com/elcf/python-ndmath)

## Contributing
Contributions and feedback are welcome and greatly appreciated. Feel free to open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

## Citing
This package is scheduled to be presented as a conference paper at the SNAME FAST Conference 2021:
* Castro-Feliciano, E. L., 2021, "OpenPlaning: Open-Source Framework for the Hydrodynamic Design of Planing Hulls," SNAME FAST '21 Conference Proceedings

## References
* Castro-Feliciano, E. L., Sun, J., and Troesch, A. W., 2017, "First Step Toward the Codesign of Planing Craft and Active Control Systems," J. Offshore Mech. Arct. Eng., 139(1)
* Faltinsen, O. M., 2005, "Planing Vessels," Hydrodynamics of High-Speed Marine Vehicles, Cambridge University Press, New York, p. 342.
* Fridsma, G., 1971, "A Systematic Study of the Rough-Water Performance of Planing Boats (Irregular Waves - Part II)," Tech. Rep. 1495, Stevens Institute of Technology
* Hadler, J. B., 1966, "The Prediction of Power Performance on Planing Craft," SNAME Trans., 74, pp. 563–610.
* Savitsky, D., 1964, "Hydrodynamic Design of Planing Hulls," Mar. Technol., 1(1), pp. 71–94.
* Savitsky, D., and Brown, P. W., 1976, "Procedures for Hydrodynamic Evaluation of Planing Hulls in Smooth and Rough Water," Mar. Technol., 13(4), pp. 381-400
