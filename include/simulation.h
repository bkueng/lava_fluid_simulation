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
#include <string>

struct ErruptionConfig {
	dfloat start_time;
	//...
};

struct SimulationConfig {
	dfloat simulation_time = 60; /** total time to simulate */
	int num_frames = -1; /** limit number of frames (-1 = no limit) */

	dfloat neighbor_lookup_dist = 0.03;

	dfloat cell_size = 0.03; /** cell size in [0,1] space */
	int num_y_cells = 65; /** number of cells in y direction: increase this for smaller cell_size! */

	//pressure constants
	dfloat k = 1000; /** stiffness parameter: higher stiffness needs smaller timesteps */
	dfloat rho0 = 1000; /** rest density of a particle */

	Math::Vec3f g = Math::Vec3f(0, -0.981, 0); /** gravity acceleration */

	dfloat viscosity = 1; //TODO: make temp-dependent

	dfloat particle_mass = 0.0072;

	std::vector<ErruptionConfig> erruptions;

	/* output */
	std::string output_dir; /** output directory without ending '/' */
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

	/**
	 * add particles on a grid defined by min_pos and max_pos. y values will
	 * be offset by the height field
	 * @param min_pos
	 * @param max_pos
	 * @param counts              how many particles in each direction
	 * @param initial_velocity
	 * @param temperature
	 * @param calc_mass           whether particle mass should be calculated
	 *                            (from the volume). if false, use the config value
	 */
	void addParticlesOnGrid(const Math::Vec3f& min_pos, const Math::Vec3f& max_pos,
        const Math::Vec3i& counts, const Math::Vec3f& initial_velocity,
		dfloat temperature, bool calc_mass);
private:

	/** add a new particle. Note that this (can) invalidate all existing particle
	 * pointers!
	 */
	inline void addParticle(const Math::Vec3f& position, const Math::Vec3f& velocity,
			dfloat temperature);

	inline dfloat pressure(const Particle& particle) {
		return std::max((dfloat)0., m_config.k * (particle.density - m_config.rho0));
	}

	void initOutput();
	void writeOutput(int frame);

	SimulationConfig m_config;

	std::vector<Particle> m_particles; /** all particles: can dynamically change */

	Grid<Particle>* m_grid = NULL;
	HeightField& m_height_field;

	KernelPoly6 m_kernel;
	KernelSpiky m_kernel_pressure;
	KernelViscosity m_kernel_viscosity;
};

#endif /* _SIMULATION_H_ */
