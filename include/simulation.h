/*
 * Copyright (C) 2014 Beat Küng <beat-kueng@gmx.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "height_field.h"
#include "particle.h"
#include "grid.h"
#include "global.h"
#include "kernel.h"

#include <vector>

struct ErruptionConfig {
	dfloat start_time;
	//...
};

struct SimulationConfig {
	dfloat simulation_time = 0.1; /** total time to simulate */

	dfloat neighbor_lookup_dist = 0.08;

	dfloat cell_size = 0.08; /** cell size in [0,1] space */
	int num_y_cells = 25; /** number of cells in y direction: increase this for smaller cell_size! */

	//pressure constants
	dfloat k = 1000; /** stiffness parameter: higher stiffness needs smaller timesteps */
	dfloat rho0 = 1000; /** rest density of a particle */

	Math::Vec3f g = Math::Vec3f(0, -0.981, 0); /** gravity acceleration */

	dfloat viscosity = 1; //TODO: make temp-dependent

	dfloat particle_mass = 0.0072;

	std::vector<ErruptionConfig> erruptions;
};

/**
 ** class Simulation
 *  main class for a simulation
 */
class Simulation {
public:
	Simulation(HeightField& height_field, const SimulationConfig& config);
	~Simulation();

	void run();

private:

	/** add a new particle. Note that this (can) invalidate all existing particle
	 * pointers!
	 */
	inline void addParticle(const Math::Vec3f& position, const Math::Vec3f& velocity,
			dfloat temperature);

	inline dfloat pressure(const Particle& particle) {
		return std::max((dfloat)0., m_config.k * (particle.density - m_config.rho0));
	}

	SimulationConfig m_config;

	std::vector<Particle> m_particles; /** all particles: can dynamically change */

	Grid<Particle>* m_grid = NULL;
	HeightField& m_height_field;

	KernelPoly6 m_kernel;
	KernelSpiky m_kernel_pressure;
	KernelViscosity m_kernel_viscosity;
};

#endif /* _SIMULATION_H_ */
