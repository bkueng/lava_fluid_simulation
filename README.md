### SPH Lava Simulation on CPU ###

This is a 3D lava particle simulator based on Smoothed-Particle Hydrodynamics.
It is a project developed in the ETH course Physically-based Simulation.


#### Compilation ####
There are no external dependencies. The code works under Linux & Mac.

Release build:
`$ make`

Debug build:
`$ make debug`

Then to render you need to first compile the shaders in data/shaders. Eg for
Pixie:
`$ cd data/shaders`
`$ sdrc *.sl`


#### Usage ####
Input is an XML scene configuration file under config/ together with a height
field (see include/simulation.h for documentation of the config parameters).
The output is a RenderMan Scene File (RIB), which can be rendered with a
RenderMan compliant renderer. We used Pixie (http://www.renderpixie.com).

##### Example #####
`$ ./simulator -c config/grid_simple.xml -f 100`
`$ ./scripts/render.sh config/grid_simple.xml -r <pixie_bin>/rndr`
`$ ./scripts/video.sh output/rendering/grid_simple 20`

##### Configuration #####
How to setup & tune configuraton parameters:
* Set `lookup_dist` equal to `smoothing_kernel_size`, then make sure that the
  simulation has around 30-40 average neighbors (avg\_nei).
* Then slightly increase the `lookup_dist` to have around 20% neighbor list
  updates (nei\_upd) for optimal performance.
* Make sure the surface particles are correctly calculated by changing
  surface\_air\_threshold. To visualize, set color="surface" in output.


#### Literature ####
Parts of the implementation are based on these papers:
* Particle-Based Fluid Simulation for Interactive Applications, Matthias Mueller
  et al., 2003
* Animating Lava Flows, Dan Stora et al., 1999
* Numerical simulation of lava flow using a GPU SPH model, Alexis Herault et
  al., 2011


Copyright 2014 Hans Hardmeier <hanshardmeier@gmail.com>
Copyright 2014 Beat Küng <beat-kueng@gmx.net>

