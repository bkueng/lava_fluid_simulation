### PCISPH Lava Simulation on CPU ###

This is a 3D lava particle simulator based on Smoothed-Particle Hydrodynamics
(PCISPH).
It is a project developed in the ETH course Physically-based Simulation.


#### Features ####
* Predictive-Corrective Incompressible SPH simulation based on 3D grid and
  neighbor lists (for pure SPH version, see git tag pure\_SPH)
* Temperature diffusion inside volume and to the ground and air
* Spring and elastic impact based ground collisions
* Viscosity calculation based on temperature
* Multithreading using OpenMP
* Configurable neighbor lookup distance so that neighbor lists can be used for
  multiple timesteps
* Multiple rendering modes using a RIB (Renderman) scene output file:
  * Point rendering (with different color modes: density, temperature, surface
    or shader)
  * Spheres or oriented disks
  * Surface rendering using disk splatting


#### Compilation ####
There are no external dependencies for the simulator. The code works under Linux
and Mac.

Release build:  
`$ make`

Debug build:  
`$ make debug`

Then to render you need to first compile the shaders in data/shaders. Eg for
Pixie:  
`$ cd data/shaders`  
`$ sdrc *.sl`

The simulation output can get considerably large, so there is an option to
compile with compressed output (gzip files, requires libz). Pixie can directly
read it, however simulation performance suffers quite a bit.  
`$ make clean && make COMPRESSION=1`


##### Script Dependencies #####
* render script: Renderman compliant renderer (eg Pixie), ImageMagick for
  surface rendering (spatial and temporal blurring steps)
* video script: ffmpeg


#### Usage ####
Input is an XML scene configuration file under config/ together with a height
field (see include/simulation.h for documentation of the config parameters).
The output is a RenderMan Scene File (RIB), which can be rendered with a
RenderMan compliant renderer. We used Pixie (http://www.renderpixie.com).

Simulator output:  
`#P: 13656, T:  1.463/10.000, steps:  1.8/s (549ms), avg_nei: 41, nei_upd:50% temp:[617 1000] avg_it:12`  
With current number of particles, current and total simulation time, simulation
steps per second and average time of one timestep, average number of neighbors,
percentage of steps where neighbor lists needed to be updated, [min max]
temperature and average number of iterations for PCI loop.

##### Example #####
`$ ./simulator -c config/grid_simple.xml -f 100`  
`$ ./scripts/render.sh config/grid_simple.xml -r <pixie_bin>/rndr`  
`$ ./scripts/video.sh output/rendering/grid_simple --fps 20`  


##### Configuration #####
Configuration files are under ./config using xml files.
How to setup & tune configuraton parameters:
* To figure out the particle mass, first use a particles\_grid with the desired
  particle density and calc\_mass=true. Running the simulation will output the
  particle mass.
* Set `lookup_dist` equal to `smoothing_kernel_size`, then make sure that the
  simulation has around 30-40 average neighbors (avg\_nei).
* Then slightly increase the `lookup_dist` to have around 20% neighbor list
  updates (nei\_upd) for optimal performance.
* To get a stable simulation, the time\_step needs to be low enough and the
  viscosity high enough
* Make sure the surface particles are correctly calculated by changing
  surface\_air\_threshold. To visualize, set color="surface" in output.
* The rate how fast the temperature diffuses within the fluid and to the air,
  can be changed with the parameters temperature\_diffusion\_coeff and
  temperature\_diffusion\_coeff\_air, temperature\_diffusion\_coeff\_ground.
  Increasing these will result in faster diffusion.
* Multipass Rendering:
  * Image blur: use pass{i}="blur 0xX file" where file is one of the outputs of
    the previous rendering pass (eg normal). blurring will be done in-place, so
    the same file can be used in following rendering passes.
  * Temporal blur: use pass{i}="temporalblur size [name]" where size is the
    neighborhood size in the past and future (1 means to filter over previous,
    current and next image). If name is not given, the final output is blurred,
    otherwise it's the same convention as for the blur.
    Note: this is not in-place, the new filename is
    tblur\_[name]\_out\_%06i.tif. So the following passes need to use the
    prefix tblur\_ to use the output in a .rib file.
    To render the video with final temporal blur use --prefix tblur
* High quality rendering:
  * Increase image resolution (Format x y)
  * Decrease ShadingRate
  * Enable PixelSamples and PixelFilter
  * Change the amount of blur in the config file according to the changed image
    resolution


#### Literature ####
Parts of the implementation are based on these papers:
* Particle-Based Fluid Simulation for Interactive Applications, Matthias Mueller
  et al., 2003
* Animating Lava Flows, Dan Stora et al., 1999
* Numerical simulation of lava flow using a GPU SPH model, Alexis Herault et
  al., 2011
* Predictive-Corrective Incompressible SPH, B. Solenthaler et al., 2009


Copyright 2014 Hans Hardmeier <hanshardmeier@gmail.com>  
Copyright 2014 Beat Küng <beat-kueng@gmx.net>

