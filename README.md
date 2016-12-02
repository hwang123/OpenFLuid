# OpenFLuid
6.837 final project particle-based fluid simulation, William Caruso &amp; Harrison Wang

## Smooth Particle Hydrodynamics (SPH) 

SPH is a computational method used for simulating fluid flows. It was developed by Gingold and Monaghan (1977) and Lucy (1977) initially for astrophysical problems. It is a mesh-free Lagrangian method (where the coordinates move with the fluid), and the resolution of the method can easily be adjusted with respect to variables such as the density.

The smoothed-particle hydrodynamics (SPH) method works by dividing the fluid into a set of discrete elements, referred to as particles. These particles have a spatial distance (known as the "smoothing length", typically represented in equations by `h`, over which their properties are "smoothed" by a kernel function. This means that the physical quantity of any particle can be obtained by summing the relevant properties of all the particles which lie within the range of the kernel. 

## Plan of Action

### Phase 1

The goal of Phase 1 is to show fluid-like particle interaction with a small amount of particles (100). We will need to calculate particle-particle interactions over all of the remaining particles (expensive but small amount of particles). 

This involves the following:

* 10 x 10 grid of particles (Particle, Particle System classes)
* Container (Wall class))
* compute acceleration
* compute density pressure

> It should be clear that we are not taking Surface Tension into account in Phase 1

### Phase 2

The goal of Phase 2 is to include more factors in the particle simulation. 

* Surface Tension

