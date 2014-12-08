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

#include "Math/Vec3.h"
#include "Math/Vec2.h"
#include "Math/Rand.h"

#include <vector>
#include <string>
#include <memory>

class ErruptionSourceBase {
public:
	/**
	 * get a random position
	 * @param height_field
	 * @param u, v             [0,1] random values
	 * @return                 point above the heightfield
	 */
	virtual Math::Vec3f getPosition(const HeightField& height_field, dfloat u, dfloat v) = 0;
private:
};

class ErruptionSourceLineSegment : public ErruptionSourceBase {
public:
	ErruptionSourceLineSegment(const Math::Vec2f& start, const Math::Vec2f& end, dfloat y_offset)
		: m_start(start), m_end(end), m_y_offset(y_offset) {}

	virtual Math::Vec3f getPosition(const HeightField& height_field, dfloat u, dfloat v);
private:
	Math::Vec2f m_start;
	Math::Vec2f m_end;
	dfloat m_y_offset;
};

//TODO: more source types: rectangle source...


struct ErruptionConfig {
	dfloat start_time;
	dfloat duration;
	dfloat particles_per_sec;
	dfloat init_temperature = 1000;
	Math::Vec3f init_velocity;
	std::shared_ptr<ErruptionSourceBase> source;

	//additional non-config data
	dfloat particles_left = 0.;

	bool operator<(const ErruptionConfig& other) const {
		return start_time > other.start_time; //smallest is last entry
	}
};

struct SimulationConfig {
	dfloat simulation_time = 60; /** total time to simulate */
	dfloat time_step = 0.001;
	int num_frames = -1; /** limit number of frames (-1 = no limit) */

	dfloat smoothing_kernel_size = 0.03; /** the kernel size for smoothing neighbor
				particles. should be such that avg neighbors is within ~30-40 */
	dfloat neighbor_lookup_dist = 0.03; /** must be >= smoothing_kernel_size.
				this defines the max lookup distance for neighbor search.
				a greater distance means the same set of neighbors can be used
				for multiple timesteps, which improves performance */
	dfloat cell_size = 0.03; /** 3D grid cell size in [0,1] space (normally set
				equal to neighbor_lookup_dist) */
	int num_y_cells = 65; /** number of cells in y direction: increase this for smaller cell_size! */

	//pressure constants
	dfloat rho0 = 1000; /** rest density of a particle */

	dfloat viscosity_coeff_a = 0.002;
	dfloat viscosity_coeff_b = 10; /** temperature T to viscosity vi formula:
				vi = b*exp(-a*T). a,b > 0. a smaller = flatter curve */

	dfloat temperature_diffusion_coeff = 90000; /** diffusion coefficient inside volume */
	dfloat temperature_diffusion_coeff_air = 0.2; /** diffusion coefficient volume <--> air */
	dfloat temperature_diffusion_coeff_ground = 9000000; /** diffusion coefficient volume <--> ground */
	dfloat temperature_air = 20; /** constant air temperature in degrees */
	dfloat temperature_ground = 10;

	dfloat surface_air_threshold = 0.4; /** threshold to determine whether a particle
				is at the surface (to the air). [0, 1], lower value means less surface
				particles (best values in [0.1, 0.5]. use 'color="surface"' to validate */
	dfloat surface_ground_radius_factor = 0.1; /** defines how far away a
				particle can be from the ground, and still be considered to touch it.
				max distance from ground is = surface_ground_radius_factor * particle_radius */

	Math::Vec3f g = Math::Vec3f(0, -9.81, 0); /** gravity acceleration */

	dfloat particle_mass = 0.0072;

	dfloat init_velocity_perturb_angle = 0.; /** randomly perturb the initial
				velocity of an errupted particle by maximally this angle (in degrees) */

	enum GroundMethod {
		GroundForceSpring = 0, /** "spring" */
		GroundElastic          /** "elastic" */
	};
	GroundMethod ground_method = GroundElastic;
	dfloat ground_spring = 1000.; /** ground spring constant */

	std::vector<ErruptionConfig> erruptions;

	/* output */
	std::string output_dir; /** output directory without ending '/' */
	int output_rate = 1; /** write output data every x timestep */
	enum OutputColor {
		ColorPressure = 0,  /** "pressure" */
		ColorTemperature,  /** "temperature" (lowest=blue, middle=green, highest=red) */
		ColorSurface,      /** "surface" indicate whether particles are on
							   ground and/or on surface */
		ColorShader        /** "shader" additional normalized temperature array
		                       for shader. only with FormatSurface */
	};
	OutputColor output_color = ColorPressure;
	enum OutputFormat {
		FormatPoint = 0, /** "point". render flat points */
		FormatSphere,    /** "sphere". render spheres */
		FormatDisk,      /** "disk". render (flat) disks */
		FormatSurface    /** "surface". render a surface */
	};
	OutputFormat output_format = FormatPoint;
	float output_constantwidth = 0.003; /** output sphere width/point width */
	int min_neighborhood_size = 1; /** minimum number of neighbors for a particle
	                           to be outputed (including the particle itself */

	/* PCISPH */
	dfloat density_error_threshold = 0.01; /** acceptable error threshold for
								density update iterations */
	int min_density_iterations = 3; /** minimum number of density update iterations */
	int max_density_iterations = 200; /** maximum number of density update iterations */
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
	 * @param print_avg_density   if true calculate average density & print (may be slow)
	 */
	void addParticlesOnGrid(const Math::Vec3f& min_pos, const Math::Vec3f& max_pos,
        const Math::Vec3i& counts, const Math::Vec3f& initial_velocity,
		dfloat temperature, bool calc_mass, bool print_avg_density = false);
private:

	/** add a new particle. Note that this (can) invalidate all existing particle
	 * pointers!
	 */
	inline void addParticle(const Math::Vec3f& position, const Math::Vec3f& velocity,
			dfloat temperature);

	inline dfloat viscosity(const Particle& particle) const {
		return m_config.viscosity_coeff_b *
				exp(-m_config.viscosity_coeff_a * particle.temperature);
	}

	/** check whether a position is within the grid and adjust position &
	 * velocity if it's outside.
	 */
	inline void checkGridBoundary(Math::Vec3f& position, Math::Vec3f& velocity) const;

	/**
	 * find a position above ground with conditions:
	 * @param position     current position below ground, will be changed
	 * @param dir          opposite velocity direction where (position+dir) is above ground
	 */
	inline void findPosAboveGround(Math::Vec3f& position, const Math::Vec3f& dir) const;

	/**
	 * @param try_to_defer     try not to add any particles in this step, but
	 *                         accumulate and add in the future. this is for
	 *                         performance.
	 * @return true if particle array was changed (ie particles added/removed)
	 */
	bool handleErruptions(dfloat simulation_time, dfloat dt, bool try_to_defer);

	Math::Vec3f getInitVelocity(const ErruptionConfig& config);


	void initOutput();
	void writeOutput(int frame);

	SimulationConfig m_config;

	std::vector<Particle> m_particles; /** all particles: can dynamically change */

	Grid<Particle>* m_grid = NULL;
	HeightField& m_height_field;

	KernelPoly6 m_kernel;
	KernelSpiky m_kernel_pressure;
	KernelViscosity m_kernel_viscosity;

	dfloat m_particle_radius = 0.;
	dfloat m_max_temperature = 0.;

	std::list<ErruptionConfig> m_active_erruptions;
	int m_deferred_erruption_steps = 0;
	static constexpr dfloat m_max_deferred_erruption_steps = 3;


	Math::RandMT m_rand;
};

#endif /* _SIMULATION_H_ */
