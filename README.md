# Astrodynamics Projects

## Table of Contents
### Kepler's Equation solution
**Description**  
Kepler’s equation
**M = E − e sin E**
is transcendental in E, i.e. given M, a solution in E cannot be expressed in closed form. However, during the years
a variety of numerical and analytical methods have been developed to compute the solutions.  
* **Newton Rhapson Method**  
  Kepler's Equation Solution E(M) is computed using Newton Rhapson method for different values of eccentricity e. 
* **Lagrange's Theorem**  
  In 1770, Lagrange’s work on Kepler’s equation led to a generally useful expansion theorem. Given an equation
of the form y = x + αφ(y), its solution is approximated by a series expansion. Langrage’s theorem is applied to Kepler’s equation and the series solution is computed. 
The series expression E(M) is computed by keeping terms up to order n = 3 and n = 10.  
* **Fourier Bessel Expansion**  
Another representation of the solution in E is given by re-arranging the terms as a Fourier-Bessel series expansion. The Bessel functions are computed
manually. The series is computed by keeping terms up to order j = 3 and j = 10.


### Astrodynamics Project
* **Time in Earth’s shadow**  
The time a satellite spends in Earth's shadow (which is considered cylindrical) is computed from this program. The satellite moves in an 
elliptical orbit around Earth with eccentricity e and semi-major axis a. Moreover, the Sun lies on the (1, 0, 0) direction in the satellite’s perifocal reference frame.  
* **Lunar Orbiter Hill's radius**  
The semi major axis of a lunar-sychronous orbiter is computed. That distance is compared to Moon's Hill's radius. 
The Hill’s radius is the distance from the Moon within which smaller bodies would tend to orbit around it. 
* **Lunar Orbiter elliptical orbit**  
The Keplerian ellipse that the satellite describes about the Moon is drawn in a three dimensional plot. The distance of the satellite from
Moon’s surface is computed for an orbital period assuming the Moon is a spherical surface.  The orbital stability of the orbiter is verified by orbit propagation.  
* **Lunar Orbiter Groundtracks**  
The Groundtracks of the lunar orbiter are comptuted for 3 and 30 days time spans. 
* **Earth satellite maneuvering**  
From an initial orbital state, the satellite is transferred to another orbital location using
change of plane, change of perigee and change of shape maneuvers.  
The time of flight needed to perform the transfer as well as the ΔV cost are computed. 
### ISS Orbit Integration without pertubations  
In this script our goal is to numerically integrate the differential equations of a satellite in cartesian coordinates.
There are no pertubations in the system we are simulating.

The only force affecting the satellite is considered to be Earth's gravitational pull.

Thus, we expect the orbital elements to remain unchanged during the orbit propagation, since no pertubative forces affect the movement of the satellite.
### ISS Orbit Integration considering Stoke's Drag  
In this script our goal is to numerically integrate the differential equations of a satellite in cartesian coordinates
The system is considered pertrubed from Stoke's drag.
The integration time span is 6000 orbital periods.
The orbit of the satellite is plotted. As we see, the orbital radius progressively decreases.
We compute the time needed for the satellite to impact the earth's surface due to the pertubative force of Stoke's drag.
Next, we calculate the theoretical time according to the prediction from Gauss Equations.
We find that the simulation time coincides with the theoretical time.
### ISS Orbit Integration considering Solar Radiation Pressure
### Satellite orbit in oblate Earth potential
* Part A  
* Part B
