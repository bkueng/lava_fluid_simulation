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

#include "simulation.h"
#include "timer.h"
#include "memory_pool.h"

#include <cstdio>

using namespace Math;

Simulation::Simulation(HeightField& height_field, const SimulationConfig& config)
	: m_config(config), m_height_field(height_field),
	m_kernel(config.neighbor_lookup_dist),
	m_kernel_pressure(config.neighbor_lookup_dist),
	m_kernel_viscosity(config.neighbor_lookup_dist) {

	m_grid = new Grid<Particle>(height_field, config.cell_size, config.num_y_cells);
}
Simulation::~Simulation() {
	if (m_grid) delete (m_grid);
}

void Simulation::run() {

	initOutput();

	dfloat simulation_time = 0;
	int num_timesteps = 0;
	int num_timesteps_total = 0;
	int64_t start_time = getTickCount();
	dfloat particle_mass = m_config.particle_mass;

	const int min_neighbors_array_size = 500;
	MemoryPool<Particle*> memory_pool(min_neighbors_array_size);
	bool simulation_running = true;

	while(simulation_running) {

		dfloat dt = 0.003; //TODO: dynamic??

		/* erruptions: add particles */
		//TODO


		/* neighbor search */
		memory_pool.reset();
		m_grid->updateEntries(m_particles);

		long total_neighbors=0;
		for(auto& particle : m_particles) {
			int num_neighbors = 0;
			Particle** neighbor_array = memory_pool.next();
			particle.neighbors = neighbor_array;

			auto neighbor_cb = [&](Particle* neighbor, dfloat dist2) {
				DEBUG_ASSERT1(num_neighbors < min_neighbors_array_size);
				neighbor_array[num_neighbors++] = neighbor;
			};
			m_grid->iterateNeighbors(particle.position,
					m_config.neighbor_lookup_dist, neighbor_cb);

			memory_pool.setNumUsedElements(num_neighbors);
			particle.num_neighbors = num_neighbors;
			total_neighbors += num_neighbors;
		}
		int avg_neighbors = (int)(total_neighbors / m_particles.size());


		/* density update */
		for (auto& particle : m_particles) {
			dfloat kernel_sum = 0;
			for (int i = 0; i < particle.num_neighbors; ++i) {
				kernel_sum += m_kernel.eval(particle.position, particle.neighbors[i]->position);
			}
			particle.density = m_kernel.kernelWeight() * particle_mass * kernel_sum;
		}


		/* force computation */
		for (auto& particle : m_particles) {
			dfloat particle_pressure = pressure(particle);

			Math::Vec3f force_pressure(0.);
			Math::Vec3f force_viscosity(0.);
			for (int i = 0; i < particle.num_neighbors; ++i) {
				Particle* neighbor_part = particle.neighbors[i];

				force_pressure += ((particle_pressure + pressure(*neighbor_part)) /
					neighbor_part->density) *
					m_kernel_pressure.evalGrad(particle.position, neighbor_part->position);

				force_viscosity += (neighbor_part->velocity - particle.velocity) *
					(m_kernel_viscosity.evalLaplace(particle.position, neighbor_part->position) /
					neighbor_part->density);
			}
			force_pressure *= -0.5 * particle_mass * m_kernel_pressure.kernelWeightGrad();
			force_viscosity *= m_config.viscosity * particle_mass *
					m_kernel_viscosity.kernelWeightLaplace();

			Math::Vec3f force_gravity = particle.density*m_config.g;

			particle.forces = force_pressure + force_viscosity + force_gravity;
		}


		/* update positions, velocities & handle collisions */
		for (auto& particle : m_particles) {
			//TODO: better integration scheme? -> leapfrog?

			particle.velocity += dt * particle.forces / particle.density;
			particle.position += dt * particle.velocity;

			//TODO: collisions: use grid

			Vec3f& pos = particle.position;
			if(pos.x < Math::FEQ_EPS)
				pos.x = Math::FEQ_EPS;
			else if(pos.x > m_height_field.fieldWidth()-Math::FEQ_EPS)
				pos.x = m_height_field.fieldWidth()-Math::FEQ_EPS;
			if(pos.z < Math::FEQ_EPS)
				pos.z = Math::FEQ_EPS;
			else if(pos.z > m_height_field.fieldDepth()-Math::FEQ_EPS)
				pos.z = m_height_field.fieldDepth()-Math::FEQ_EPS;

			dfloat height_field_val = m_height_field.lookup(pos.x, pos.z);
			if(pos.y < height_field_val) pos.y = height_field_val;

			//test//////////
			dfloat y_max = m_config.cell_size*(m_config.num_y_cells-5) + height_field_val;
			if(pos.y > y_max) pos.y = y_max;
		}


		/* write output */
		writeOutput(num_timesteps_total);


		simulation_time += dt;
		++num_timesteps;
		++num_timesteps_total;

		simulation_running = simulation_time < m_config.simulation_time &&
				(m_config.num_frames == -1 || num_timesteps_total < m_config.num_frames);

		/* statistics */
		int64_t cur_time = getTickCount();
		double elapsed_time = getTickSeconds(cur_time-start_time);
		if(elapsed_time >= 1. || !simulation_running) {
			printf("particles:%6i, time:%7.3f / %.3f, steps:%6.2f/s (%4.0fms/step), avg_neighbors:%3i\n",
					(int)m_particles.size(), (float)simulation_time,
					(float)m_config.simulation_time,
					(float)num_timesteps / elapsed_time,
					(float)elapsed_time/num_timesteps*1000.f,
					avg_neighbors);
			start_time = cur_time;
			num_timesteps = 0;
		}
	}
}

inline void Simulation::addParticle(const Math::Vec3f& position,
		const Math::Vec3f& velocity, dfloat temperature) {
	Particle particle;
	particle.position = position;
	particle.velocity = velocity;
	particle.temperature = temperature;
	m_particles.push_back(particle);
}

void Simulation::initOutput() {
	printf("Saving output to %s\n", m_config.output_dir.c_str());
}

void Simulation::writeOutput(int frame) {
	char buffer[32];
	/* RIB output */
	sprintf(buffer, "/frame_%06i.rib", frame);
	string file_name = m_config.output_dir + buffer;

	FILE* file = fopen(file_name.c_str(), "w");
	if (!file) throw EXCEPTION(EFILE_ERROR);
	//positions
	fprintf(file, "Points \"P\" [ ");
	for(const auto& particle : m_particles) {
		fprintf(file, "%.7lf %.7lf %.7lf ",
			(double)particle.position.x, (double)particle.position.y, (double)particle.position.z);
	}

	fprintf(file, "]\n \"Cs\" [ ");
	//colors
	dfloat max_density = -1e12;
	dfloat min_density =  1e12;
	for(const auto& particle : m_particles) {
		if(particle.density > max_density) max_density = particle.density;
		if(particle.density < min_density) min_density = particle.density;
	}
	for(const auto& particle : m_particles) {
		float r = 1;
		float g = (particle.density-min_density)/(max_density-min_density);
		float b = 0;
		fprintf(file, "%.3f %.3f %.3f ", r, g, b);
	}
	fprintf(file, "]\n");

	fclose(file);
}

void Simulation::addParticlesOnGrid(const Math::Vec3f& min_pos,
		const Math::Vec3f& max_pos, const Math::Vec3i& counts,
		const Math::Vec3f& initial_velocity, dfloat temperature, bool calc_mass) {
	dfloat x_min = min_pos.x, x_max = max_pos.x;
	dfloat y_min = min_pos.y, y_max = max_pos.y;
	dfloat z_min = min_pos.z, z_max = max_pos.z;
	dfloat dx, dy, dz;
	Vec3f p;
	for (dfloat z = z_min; z <= z_max; z+=dz=(z_max - z_min) / counts.z) {
		for (dfloat x = x_min; x <= x_max; x+=dx=(x_max - x_min) / counts.x) {
			for (dfloat y = y_min; y <= y_max; y+=dy=(y_max - y_min) / counts.y) {
				p.x = x;
				p.z = z;
				p.y = m_height_field.lookup(p.x, p.z) + y;
				addParticle(p, initial_velocity, temperature);
			}
		}
	}
	if(calc_mass) {
		m_config.particle_mass = dx * dy * dz * m_config.rho0;
		printf("Particle Mass=%lf\n", m_config.particle_mass);
	}
}
