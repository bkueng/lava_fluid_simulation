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
field.
The output is a RenderMan Scene File (RIB), which can be rendered with a
RenderMan compliant renderer. We used Pixie (http://www.renderpixie.com).

##### Example #####
`$ ./simulator -c config/sample.xml -f 100`
`$ ./scripts/render.sh config/sample.xml -r <pixie_bin>/rndr`
`$ ./scripts/video.sh output/rendering/sample 20`


Copyright 2014 Hans Hardmeier <hanshardmeier@gmail.com>
Copyright 2014 Beat Küng <beat-kueng@gmx.net>

