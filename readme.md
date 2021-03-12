# Launch Trajectory
Using Python to simulate a rocket launch with constant acceleration (3DOF).

## Objective
Simulate a rocket's flight path through atmosphere. This is accomplished by creating an acceleration model using SciPy's ODEINT and integrating to calculate velocity and position.

Features include an indexed atmospheric model, coarse mass model, and drag vector.

##  Units: SI
- pos:   altitude           m
- vel:   velocity           m/s
- a:     speed of sound     m/s
- acc:   acceleration       m/s^2
- g0:    gravity            m/s^2
- u:     kinematic viscosity   m^2/s
- T:     temperature        °K
- P:     pressure           Pa
- m:     mass               kg
- s:     frontal area       m^2
- rho:   density            kg/m^3
- R:     gas constant       N-m/kg-K
- d:     drag               N

# Constants
- gamma = 1.4
- R = 287.05287 N-m/kg-K
- g0 = 9.80665 m/s^2
- Re = 6356766 m
- azimuth = 45 degrees
- angleOfAttack = 45 degrees

## Performance Parameters
- coefficient of drag = random float in range .1-.7
- frontal area = random float in range .1-.7 m^2
- thrust = 50,000 N
- burnTime = 1.8 s

## Initial Values
- position x = 0 m
- position y = 0 m
- position z = 210 m
