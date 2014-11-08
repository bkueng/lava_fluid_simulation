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
	dfloat simulation_time = 0;
	/*
	//test: add particles & test grid neighbor iter
	Vec3f p, v(0);
	p.x = 0.3; p.z = 0.23;
	p.y = m_height_field.lookup(p.x, p.z)+0.01;
	addParticle(p, v, 10);
	p += Vec3f(0.001, 0., 0);
	p.y = m_height_field.lookup(p.x, p.z)+0.01;
	addParticle(p, v, 10);

	m_grid->updateEntries(m_particles);
	printf("lookup pos=%.3f %.3f %.3f\n", p.x, p.y, p.z);
	auto cb = [&](Particle* p, dfloat dist2) {
		Vec3f& pos = p->position;
		printf("neighbor: %.3f, %.3f, %.3f, d=%.4f\n", pos.x, pos.y, pos.z, sqrt(dist2));
	};
	m_grid->iterateNeighbors(p, 0.0013, cb);
	//*/

	//add particles over a grid
	dfloat x_min = 0.3, x_max = 0.6;
	dfloat z_min = 0.3, z_max = 0.6;
	dfloat y_min = 0.01, y_max = 0.05;
	dfloat dx, dy, dz;
	for(dfloat z = z_min; z<=z_max; z+=dz=(z_max-z_min)/10) {
		for(dfloat x = x_min; x<=x_max; x+=dx=(x_max-x_min)/10) {
			for(dfloat y = y_min; y<=y_max; y+=dy=(y_max-y_min)/5) {
				Vec3f p, v(0);
				p.x = x;
				p.z = z;
				p.y = m_height_field.lookup(p.x, p.z)+y;
				addParticle(p, v, 100);
			}
		}
	}
	m_config.particle_mass = dx*dy*dz*m_config.rho0;
	printf("particle mass=%lf\n", m_config.particle_mass);


	int num_timesteps = 0;
	int num_timesteps_total = 0;
	int64_t start_time = getTickCount();
	dfloat particle_mass = m_config.particle_mass;

	MemoryPool<Particle*> memory_pool(500);

	while(simulation_time < m_config.simulation_time) {
		dfloat dt = 0.00001; //TODO: dynamic??

		/* erruptions: add particles */


		/* neighbor search */
		memory_pool.reset();
		m_grid->updateEntries(m_particles);

		long total_neighbors=0;
		for(auto& particle : m_particles) {
			int num_neighbors = 0;
			Particle** neighbor_array = memory_pool.next();
			particle.neighbors = neighbor_array;

			auto neighbor_cb = [&](Particle* neighbor, dfloat dist2) {
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
			force_viscosity *= m_config.viscosity * particle_mass * m_kernel_viscosity.kernelWeightLaplace();

			Math::Vec3f force_gravity = particle.density*m_config.g;

			particle.forces = force_pressure + force_viscosity + force_gravity;
		}


		/* update positions, velocities & handle collisions */


		/* write output */
		//writeOutput(num_timesteps_total)


		simulation_time += dt;
		++num_timesteps;
		++num_timesteps_total;


		/* statistics */
		int64_t cur_time = getTickCount();
		double elapsed_time = getTickSeconds(cur_time-start_time);
		if(elapsed_time >= 1.) {
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
