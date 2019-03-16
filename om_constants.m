function om_constants

% astrodynamic and utility constants

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global emu mmu smu aunit omega req

% earth gravitational constant (DE421 value; km**3/sec**2)

emu = 398600.436233;

% moon gravitational constant (DE421 value; km**3/sec**2)

mmu = 4902.800076;

% sun gravitational constant (DE421 value; km**3/sec**2)

smu = 132712440040.944;

% lunar inertial rotation rate (radians/second)

omega = 2.6621916d-6;

% astronomical unit (kilometers)

aunit = 149597870.691;

% radius of the moon (kilometers)

req = 1738.0;

