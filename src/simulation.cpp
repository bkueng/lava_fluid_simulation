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
	m_kernel(config.smoothing_kernel_size),
	m_kernel_pressure(config.smoothing_kernel_size),
	m_kernel_viscosity(config.smoothing_kernel_size) {

	m_grid = new Grid<Particle>(height_field, config.cell_size, config.num_y_cells);
}
Simulation::~Simulation() {
	if (m_grid) delete (m_grid);
}

void Simulation::run() {

	initOutput();

	dfloat simulation_time = 0;
	int num_timesteps_statistics = 0;
	int num_neighbor_updates_statistics = 1;
	int num_timesteps_total = 0;
	int64_t start_time = getTickCount();
	dfloat particle_mass = m_config.particle_mass;
	dfloat max_velocity_dt = 0; /* max position change of a particle */
	dfloat smoothing_kernel_size2 = m_config.smoothing_kernel_size
			* m_config.smoothing_kernel_size;
	dfloat cur_neighbors_dist = m_config.neighbor_lookup_dist;
	bool need_neighbors_update = true;

	const int min_neighbors_array_size = 500;
	MemoryPool<Particle*> memory_pool(min_neighbors_array_size);
	bool simulation_running = true;
	long global_num_neighbors = 0; /** total number of neighbors including all threads */
	const int chunk_size = 1000; /** amount of particles to iterate per chunk for each thread */

	dfloat dt = m_config.time_step; //TODO: dynamic??

	m_grid->updateEntries(m_particles);

#pragma omp parallel firstprivate(memory_pool)
	while(simulation_running) {

		memory_pool.reset();

		/* neighbor search */
		long local_total_neighbors=0;
		if (need_neighbors_update) {
#pragma omp for schedule(dynamic, chunk_size)
			for(int particle_idx = 0; particle_idx < (int)m_particles.size(); ++particle_idx) {
				Particle& particle = m_particles[particle_idx];
				int num_neighbors = 0, num_neighbors_statistics = 0;
				Particle** neighbor_array = memory_pool.next();
				particle.neighbors = neighbor_array;

				auto neighbor_cb = [&](Particle* neighbor, dfloat dist2) {
					DEBUG_ASSERT1(num_neighbors < min_neighbors_array_size);
					neighbor_array[num_neighbors++] = neighbor;
					if (dist2 < smoothing_kernel_size2) ++num_neighbors_statistics;
				};
				m_grid->iterateNeighbors(particle.position,
						m_config.neighbor_lookup_dist, neighbor_cb);

				memory_pool.setNumUsedElements(num_neighbors);
				particle.num_neighbors = num_neighbors;
				local_total_neighbors += num_neighbors_statistics;
			}
#pragma omp atomic
			global_num_neighbors += local_total_neighbors;
		}



		/* density update */
#pragma omp for schedule(dynamic, chunk_size)
		for(int particle_idx = 0; particle_idx < (int)m_particles.size(); ++particle_idx) {
			Particle& particle = m_particles[particle_idx];
			dfloat kernel_sum = 0;
			for (int i = 0; i < particle.num_neighbors; ++i) {
				kernel_sum += m_kernel.eval(particle.position, particle.neighbors[i]->position);
			}
			particle.density = m_kernel.kernelWeight() * particle_mass * kernel_sum;
		}



		/* force computation */
#pragma omp for schedule(dynamic, chunk_size)
		for(int particle_idx = 0; particle_idx < (int)m_particles.size(); ++particle_idx) {
			Particle& particle = m_particles[particle_idx];
			dfloat particle_pressure = pressure(particle);

			Vec3f force_pressure(0.);
			Vec3f force_viscosity(0.);
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

			Vec3f force_gravity = particle.density*m_config.g;

			particle.forces = force_pressure + force_viscosity + force_gravity;
		}


		/* update positions, velocities & handle collisions */
		dfloat local_max_velocity2 = 0;
#pragma omp for schedule(dynamic, chunk_size) nowait
		for(int particle_idx = 0; particle_idx < (int)m_particles.size(); ++particle_idx) {
			Particle& particle = m_particles[particle_idx];
			Vec3f& pos = particle.position;

			//ground contact: spring force
			if(m_config.ground_method == SimulationConfig::GroundForceSpring) {
				dfloat height_field_val = m_height_field.lookup(pos.x, pos.z);
				if(pos.y < height_field_val) {
					Vec3f n;
					m_height_field.normal(pos.x, pos.z, n);
					particle.forces += n * ((height_field_val - pos.y)*m_config.ground_spring);
				}
			}

			//symplectic euler
			particle.velocity += dt * particle.forces / particle.density;
			particle.position += dt * particle.velocity;

			dfloat velocity2 = particle.velocity.length2();
			if(velocity2 > local_max_velocity2) local_max_velocity2 = velocity2;


			//grid boundaries: assume perfect elastic (should not occur anyway)
			if(pos.x < Math::FEQ_EPS) {
				pos.x = Math::FEQ_EPS;
				particle.velocity.x = -particle.velocity.x;
			} else if(pos.x > m_height_field.fieldWidth()-Math::FEQ_EPS) {
				pos.x = m_height_field.fieldWidth()-Math::FEQ_EPS;
				particle.velocity.x = -particle.velocity.x;
			}
			if(pos.z < Math::FEQ_EPS) {
				pos.z = Math::FEQ_EPS;
				particle.velocity.z = -particle.velocity.z;
			} else if(pos.z > m_height_field.fieldDepth()-Math::FEQ_EPS) {
				pos.z = m_height_field.fieldDepth()-Math::FEQ_EPS;
				particle.velocity.z = -particle.velocity.z;
			}

			//ground contact: elastic
			if(m_config.ground_method == SimulationConfig::GroundElastic) {
				dfloat height_field_val = m_height_field.lookup(pos.x, pos.z);
				if(pos.y < height_field_val) {
					Vec3f n;
					m_height_field.normal(pos.x, pos.z, n);
					pos.y = height_field_val + Math::FEQ_EPS;
					//reflect velocity
					particle.velocity = particle.velocity - (2.*dot(particle.velocity, n))*n;
				}
			}

			//make sure no particle leaves the field in y direction
			//This should never happen with correct scene config!
			//TODO: turn this off once we have a stable version??
			m_grid->moveInsideGrid(pos);
		}

		//reduce max displacement from each thread
		dfloat local_max_velocity_dt = sqrt(local_max_velocity2) * dt;
		if (local_max_velocity_dt > max_velocity_dt) {
#pragma omp critical
			{
				if (local_max_velocity_dt > max_velocity_dt)
					max_velocity_dt = local_max_velocity_dt;
			}
		}
#pragma omp flush(max_velocity_dt)
#pragma omp barrier /* wait until max_velocity_dt is updated by all threads */


		//serial part: only one thread should update this. the others will wait
#pragma omp single
		{
			simulation_time += dt;
			++num_timesteps_statistics;
			++num_timesteps_total;

			simulation_running = simulation_time < m_config.simulation_time &&
					(m_config.num_frames == -1 || num_timesteps_total
						< m_config.num_frames * m_config.output_rate);

			/* statistics */
			int64_t cur_time = getTickCount();
			double elapsed_time = getTickSeconds(cur_time-start_time);
			if(elapsed_time >= 1. || !simulation_running) {
				int avg_neighbors = (int)(global_num_neighbors / m_particles.size());
				printf("#ele:%6i, T:%7.3f / %.3f, steps:%6.2f/s (%4.0fms/step), avg_nei:%3i, nei_upd:%2i%%\n",
						(int)m_particles.size(), (float)simulation_time,
						(float)m_config.simulation_time,
						(float)num_timesteps_statistics / elapsed_time,
						(float)elapsed_time/num_timesteps_statistics*1000.f,
						avg_neighbors,
						100*num_neighbor_updates_statistics/num_timesteps_statistics);
				start_time = cur_time;
				num_timesteps_statistics = 0;
				num_neighbor_updates_statistics = 0;
			}

			//is neighbor update needed?
			cur_neighbors_dist -= 2. * max_velocity_dt;
			max_velocity_dt = 0;
			need_neighbors_update = cur_neighbors_dist < m_config.smoothing_kernel_size;


			/* erruptions: add particles */
			//TODO
			// if added/removed -> need_neighbors_update = true;


			if (need_neighbors_update) {
				global_num_neighbors = 0;
				cur_neighbors_dist = m_config.neighbor_lookup_dist;
				++num_neighbor_updates_statistics;
			}
		}
#pragma omp flush(need_neighbors_update)


		//the following sections can be done in parallel
#pragma omp sections
		{
#pragma omp section
			{
				/* write output */
				if ((num_timesteps_total - 1) % m_config.output_rate == 0) {
					writeOutput((num_timesteps_total - 1) / m_config.output_rate);
				}
			}
#pragma omp section
			{
				//TODO: try to parallelize this more (this makes about half of the serial
				//part, the other is the file output) ...
				if (need_neighbors_update)
					m_grid->updateEntries(m_particles);
			}
		}
	}

	printf("Done: calculated %i timesteps, wrote %i files\n", num_timesteps_total,
			num_timesteps_total/m_config.output_rate);
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
	dfloat dx=0, dy=0, dz=0;
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
