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
#include "output_file.h"
#include "Math/Rand.h"
#include "Math/Warp.h"

#include <cstdio>
#include <algorithm>

using namespace Math;
using namespace std;

Simulation::Simulation(HeightField& height_field, const SimulationConfig& config)
	: m_config(config), m_height_field(height_field),
	m_kernel(config.smoothing_kernel_size),
	m_kernel_pressure(config.smoothing_kernel_size),
	m_kernel_viscosity(config.smoothing_kernel_size) {

	m_grid = new Grid<Particle>(height_field, config.cell_size, config.num_y_cells);
	sort(m_config.erruptions.begin(), m_config.erruptions.end());
	m_rand.init(0);
	for (auto& erruption : m_config.erruptions) {
		if (erruption.init_temperature > m_max_temperature)
			m_max_temperature = erruption.init_temperature;
	}
}
Simulation::~Simulation() {
	if (m_grid) delete (m_grid);
}
void Simulation::checkGridBoundary(Math::Vec3f& position, Math::Vec3f& velocity)
 const {
	//grid boundaries: assume perfect elastic (should not occur anyway)
	if(position.x < Math::FEQ_EPS) {
		position.x = Math::FEQ_EPS;
		velocity.x = -velocity.x;
	} else if(position.x > m_height_field.fieldWidth()-Math::FEQ_EPS) {
		position.x = m_height_field.fieldWidth()-Math::FEQ_EPS;
		velocity.x = -velocity.x;
	}
	if(position.z < Math::FEQ_EPS) {
		position.z = Math::FEQ_EPS;
		velocity.z = -velocity.z;
	} else if(position.z > m_height_field.fieldDepth()-Math::FEQ_EPS) {
		position.z = m_height_field.fieldDepth()-Math::FEQ_EPS;
		velocity.z = -velocity.z;
	}
}

void Simulation::findPosAboveGround(Math::Vec3f& position, const Math::Vec3f& dir) const {
	dfloat a=0, b=1, dheight=1;
	int iters = 10;
	while (--iters >= 0 && dheight > Math::FEQ_EPS) {
		dfloat middle = (a+b)/2.;
		Vec3f middle_pos = position + middle * dir;
		dfloat height = m_height_field.lookup(middle_pos.x, middle_pos.z);
		if (middle_pos.y <= height) {
			a = middle;
		} else {
			dheight = middle_pos.y - height;
			b = middle;
		}
	}
	position = position + b*dir;
}

void Simulation::run() {

	initOutput();

	/* global particle attributes */
	dfloat particle_mass = m_config.particle_mass;
	dfloat rho0 = m_config.rho0;
	m_particle_radius = pow(particle_mass / rho0 * 3./4./M_PI, 1./3.);

	printf("Particle mass: %.10f, radius: %.10f, max temperature: %.1f\n",
			particle_mass, m_particle_radius, m_max_temperature);

	dfloat max_ground_distance = m_particle_radius * m_config.surface_ground_radius_factor;
	dfloat simulation_time = 0;
	int num_timesteps_statistics = 0;
	int num_neighbor_updates_statistics = 1;
	int num_timesteps_total = 0;
	int num_density_iters_statistics = 0;
	int64_t start_time = getTickCount();
	dfloat max_velocity_dt = 0; /* max position change of a particle */
	dfloat smoothing_kernel_size2 = m_config.smoothing_kernel_size
			* m_config.smoothing_kernel_size;
	dfloat cur_neighbors_dist = m_config.neighbor_lookup_dist;
	bool need_neighbors_update = true;

	const int min_neighbors_array_size = 500;
	MemoryPool<Particle*> memory_pool(min_neighbors_array_size);
	bool simulation_running = true;
	long global_num_neighbors = 0; /** total number of neighbors including all threads */
	int global_max_num_neighbors = 0, global_max_num_neighbors_idx = 0;
	const int chunk_size = 1000; /** amount of particles to iterate per chunk for each thread */

	dfloat dt = m_config.time_step;
	dfloat beta = dt*dt * particle_mass*particle_mass * 2. / (rho0*rho0); /* PCISPH update factor */
	dfloat delta = 1.;
	int num_prediction_correction_iters = 0;
	dfloat max_density_error = 0;

	m_grid->updateEntries(m_particles);

#pragma omp parallel firstprivate(memory_pool)
	while(simulation_running) {

		memory_pool.reset();

		/* neighbor search */
		long local_total_neighbors = 0;
		int local_max_num_neighbors = 0, local_max_num_neighbors_idx = 0;
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
				if (num_neighbors > local_max_num_neighbors) {
					local_max_num_neighbors = num_neighbors;
					local_max_num_neighbors_idx = particle_idx;
				}
				local_total_neighbors += num_neighbors_statistics;
			}
#pragma omp critical
			{
				global_num_neighbors += local_total_neighbors;
				if (local_max_num_neighbors > global_max_num_neighbors) {
					global_max_num_neighbors = local_max_num_neighbors;
					global_max_num_neighbors_idx = local_max_num_neighbors_idx;
				}
			}
		}



		/* force computation */
#pragma omp for schedule(dynamic, chunk_size)
		for(int particle_idx = 0; particle_idx < (int)m_particles.size(); ++particle_idx) {
			Particle& particle = m_particles[particle_idx];

			Vec3f force_viscosity(0.);
			Vec3f mass_density_gradient(0.); /* direction of highest neighbor mass */
			dfloat dtemperature = 0.;
			for (int i = 0; i < particle.num_neighbors; ++i) {
				Particle* neighbor_part = particle.neighbors[i];

				Vec3f gradient = m_kernel_pressure.evalGrad(
						particle.position, neighbor_part->position);

				dfloat laplacian_weight = m_kernel_viscosity.evalLaplace(
						particle.position, neighbor_part->position);
				force_viscosity += (neighbor_part->velocity - particle.velocity) *
					(laplacian_weight / rho0);

				dtemperature += (neighbor_part->temperature - particle.temperature) *
						(laplacian_weight / rho0);
				mass_density_gradient += gradient;
			}
			mass_density_gradient *= particle_mass * m_kernel_pressure.kernelWeightGrad();
			particle.density_gradient = mass_density_gradient;
			//surface particles: count number of particles on the half-space
			//pointed to by the negative mass_density_gradient
			int num_neighbors = 0, num_neighbors_half_space = 0;
			for (int i = 0; i < particle.num_neighbors; ++i) {
				Particle* neighbor_part = particle.neighbors[i];
				Vec3f dir = neighbor_part->position - particle.position;
				if (dir.length2() <= smoothing_kernel_size2) {
					++num_neighbors;
					if (dot(dir, mass_density_gradient) < 0.)
						++num_neighbors_half_space;
				}
			}
			particle.changeFlag(Particle::FLAG_IS_AT_AIR,
				(dfloat)num_neighbors_half_space / num_neighbors < m_config.surface_air_threshold);

			force_viscosity *= viscosity(particle) * particle_mass *
					m_kernel_viscosity.kernelWeightLaplace();

			Vec3f force_gravity = rho0*m_config.g;

			particle.forces = force_viscosity + force_gravity;
			particle.dtemperature = dtemperature;

			particle.pressure = 0;
			particle.force_pressure = 0;
		}


		/* prediction-correction iteration */
		while ((max_density_error >= rho0*m_config.density_error_threshold
				|| num_prediction_correction_iters < m_config.min_density_iterations)
				&& num_prediction_correction_iters < m_config.max_density_iterations) {

			dfloat local_max_density_error = 0;

			/* predict velocity & position */
#pragma omp for schedule(dynamic, chunk_size)
			for(int particle_idx = 0; particle_idx < (int)m_particles.size();
					++particle_idx) {
				Particle& particle = m_particles[particle_idx];
				Vec3f& pos = particle.predicted_position;
				//FIXME: spring ground force is not handled here

				particle.predicted_velocity = particle.velocity +
						dt * (particle.forces+particle.force_pressure)/rho0;
				particle.predicted_position = particle.position +
						dt * particle.predicted_velocity;

				checkGridBoundary(particle.predicted_position, particle.predicted_velocity);

				//ground contact: elastic
				if(m_config.ground_method == SimulationConfig::GroundElastic) {
					dfloat height_field_val = m_height_field.lookup(pos.x, pos.z);
					if(pos.y < height_field_val) {
						Vec3f dir = -dt * particle.predicted_velocity;
						findPosAboveGround(pos, dir);
						//no need to update velocity
					}
				}
			}

#pragma omp atomic
			//reset, in between 2 barriers to make sure all threads see this consistently
			max_density_error *= 0.;


			/* predict density & update pressure */
#pragma omp for schedule(dynamic, chunk_size)
			for(int particle_idx = 0; particle_idx < (int)m_particles.size();
					++particle_idx) {
				Particle& particle = m_particles[particle_idx];
				//Note: here we use the 'old' neighborhood, which could be inaccurate.
				dfloat predicted_density = 0;
				for (int i = 0; i < particle.num_neighbors; ++i) {
					Particle* neighbor_part = particle.neighbors[i];
					predicted_density += m_kernel.eval(particle.predicted_position,
							neighbor_part->predicted_position);
				}
				predicted_density *= m_kernel.kernelWeight() * particle_mass;
				dfloat density_err = predicted_density - rho0;
				particle.pressure += delta * density_err;
				if (particle.pressure < 0.) particle.pressure = 0.;

				if (density_err > local_max_density_error)
					local_max_density_error = density_err;
			}

			/* compute pressure force at current time */
#pragma omp for schedule(dynamic, chunk_size) nowait
			for(int particle_idx = 0; particle_idx < (int)m_particles.size();
					++particle_idx) {
				Particle& particle = m_particles[particle_idx];

				Vec3f force_pressure(0.);
				for (int i = 0; i < particle.num_neighbors; ++i) {
					Particle* neighbor_part = particle.neighbors[i];
					Vec3f gradient = m_kernel_pressure.evalGrad(
							particle.position, neighbor_part->position);
					force_pressure += (particle.pressure + neighbor_part->pressure) * gradient;
				}
				particle.force_pressure = force_pressure *
						(-0.5*particle_mass / rho0 * m_kernel_pressure.kernelWeightGrad());
			}

			//reduce max_density_error
#pragma omp critical
			{
				if (local_max_density_error > max_density_error)
					max_density_error = local_max_density_error;
			}
#pragma omp single
			++num_prediction_correction_iters;
#pragma omp flush(max_density_error, num_prediction_correction_iters)
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
				//ground particle flag
				particle.changeFlag(Particle::FLAG_IS_ON_GROUND,
						pos.y < height_field_val + max_ground_distance);
			}

			particle.temperature += particle.dtemperature * dt
					* m_config.temperature_diffusion_coeff;

			//symplectic euler
			particle.velocity += dt * (particle.forces+particle.force_pressure)/rho0;
			particle.position += dt * particle.velocity;

			dfloat velocity2 = particle.velocity.length2();
			if(velocity2 > local_max_velocity2) local_max_velocity2 = velocity2;


			checkGridBoundary(pos, particle.velocity);


			//ground contact: elastic
			if(m_config.ground_method == SimulationConfig::GroundElastic) {
				dfloat height_field_val = m_height_field.lookup(pos.x, pos.z);
				if(pos.y < height_field_val) {
					Vec3f dir = -dt * particle.velocity;
					findPosAboveGround(pos, dir);
					Vec3f n;
					m_height_field.normal(pos.x, pos.z, n);
					//reflect velocity
					particle.velocity = particle.velocity - (2.*dot(particle.velocity, n))*n;
					height_field_val = pos.y; //make sure to set the ground flag
				}
				//ground particle flag
				particle.changeFlag(Particle::FLAG_IS_ON_GROUND,
						pos.y < height_field_val + max_ground_distance);
			}

			//FIXME: here we assume a particle is either at air or touches ground.
			// if the volume is only 1 particle think, it can be both, ground & air.
			// but simply testing for ground & air flag is not enough, because
			// usually ground particles also have the air flag set.
			if (particle.isOnGround()) {
				dfloat r2_over_rho = m_particle_radius*m_particle_radius / rho0;
				particle.temperature += (m_config.temperature_ground - particle.temperature)
						* r2_over_rho * dt * m_config.temperature_diffusion_coeff_ground;
			} else if (particle.isAtAir()) {
				dfloat r2_over_rho = m_particle_radius*m_particle_radius / rho0;
				particle.temperature += (pow(m_config.temperature_air, 4) - pow(particle.temperature, 4))
						* r2_over_rho * dt * m_config.temperature_diffusion_coeff_air;
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
			num_density_iters_statistics += num_prediction_correction_iters;
			num_prediction_correction_iters = 0; //reset for next iteration...

			/* calculate delta for PCISPH */
			global_max_num_neighbors = 0;
			dfloat grad_dot_sum = 0;
			Vec3f grad_sum(0);
			if (m_particles.size() > 0) {
				Particle& max_particle = m_particles[global_max_num_neighbors_idx];
				for (int i = 0; i < max_particle.num_neighbors; ++i) {
					Particle* neighbor_part = max_particle.neighbors[i];
					Vec3f gradient = m_kernel_pressure.evalGrad(
							max_particle.position, neighbor_part->position);
					grad_sum += gradient;
					grad_dot_sum += dot(gradient, gradient);
				}
			}
			grad_sum *= m_kernel_pressure.kernelWeight();
			grad_dot_sum *= m_kernel_pressure.kernelWeight() * m_kernel_pressure.kernelWeight();
			delta = 1./(beta * (dot(grad_sum, grad_sum) + grad_dot_sum));
			if (grad_dot_sum < Math::FEQ_EPS2) delta = 1;
			//this seems to be needed to make it converge at all. not sure if
			//something else is wrong...
			delta /= 5.;


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
				dfloat min_temp = 1e20, max_temp = 0;
				for(const auto& particle : m_particles) {
					if(particle.temperature > max_temp) max_temp = particle.temperature;
					if(particle.temperature < min_temp) min_temp = particle.temperature;
				}
				int avg_neighbors = 0;
				if (m_particles.size() > 0)
					avg_neighbors = (int) (global_num_neighbors / m_particles.size());
				printf("#P:%6i, T:%7.3f/%.3f, steps:%5.1f/s (%3.0fms), avg_nei:%3i,"
						" nei_upd:%2i%% temp:[%.0f %.0f] avg_it:%2i\n",
						(int)m_particles.size(), (float)simulation_time,
						(float)m_config.simulation_time,
						(float)num_timesteps_statistics / elapsed_time,
						(float)elapsed_time/num_timesteps_statistics*1000.f,
						avg_neighbors,
						100*num_neighbor_updates_statistics/num_timesteps_statistics,
						min_temp, max_temp, num_density_iters_statistics / num_timesteps_statistics);
				start_time = cur_time;
				num_timesteps_statistics = 0;
				num_neighbor_updates_statistics = 0;
				num_density_iters_statistics = 0;
			}

			//is neighbor update needed?
			cur_neighbors_dist -= 2. * max_velocity_dt;
			max_velocity_dt = 0;
			need_neighbors_update = cur_neighbors_dist < m_config.smoothing_kernel_size;


			/* erruptions: add particles */
			if (handleErruptions(simulation_time, dt, !need_neighbors_update))
				need_neighbors_update = true;


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
	dfloat temperature_scaling = 1./m_max_temperature;
	dfloat max_pressure = -1e12, min_pressure = 1e12, pressure_span = 1;
	float r = 0, g = 0, b = 0;
	/* init color */
	switch (m_config.output_color) {
	case SimulationConfig::ColorPressure:
		for (const auto& particle : m_particles) {
			if (particle.pressure > max_pressure) max_pressure = particle.pressure;
			if (particle.pressure < min_pressure) min_pressure = particle.pressure;
		}
		pressure_span = max_pressure - min_pressure;
		if (pressure_span < Math::FEQ_EPS)
			pressure_span = 1;
		break;
	default:
		break;
	}

	auto f_color_pressure = [&](const Particle& particle) {
		r = 1;
		g = (particle.pressure - min_pressure) / pressure_span;
		b = 0;
	};
	auto f_color_temperature = [&](const Particle& particle) {
		float color_temp = (float)particle.temperature * temperature_scaling;
		if (color_temp > 0.5) {
			r = (color_temp - 0.5)*2.;
			g = 1. - (color_temp - 0.5)*2.;
			b = 0.;
		} else {
			r = 0.;
			g = color_temp*2.;
			b = 1. - color_temp*2.;
		}
	};
	auto f_color_surface = [&](const Particle& particle) {
		r = particle.isOnGround() ? 1. : 0.;
		g = particle.isAtAir() ? 1. : 0.;
		b = 0.8;
	};


	char buffer[32];
	/* RIB output */
	sprintf(buffer, "/frame_%06i.rib", frame);
	string file_name = m_config.output_dir + buffer;
	OutputFile file(file_name);


	switch (m_config.output_format) {
	case SimulationConfig::FormatPoint:
		//positions
		file.printf("Points \"P\" [ ");
		for(const auto& particle : m_particles) {
			if (particle.num_neighbors >= m_config.min_neighborhood_size) {
				file.printf("%.7lf %.7lf %.7lf ",
						(double)particle.position.x, (double)particle.position.y,
						(double)particle.position.z);
			}
		}
		file.printf("]\n \"Cs\" [ ");
		switch (m_config.output_color) {
		case SimulationConfig::ColorPressure:
			for (const auto& particle : m_particles) {
				if (particle.num_neighbors >= m_config.min_neighborhood_size) {
					f_color_pressure(particle);
					file.printf("%.3f %.3f %.3f ", r, g, b);
				}
			}
			break;
		case SimulationConfig::ColorTemperature:
			for (const auto& particle : m_particles) {
				if (particle.num_neighbors >= m_config.min_neighborhood_size) {
					f_color_temperature(particle);
					file.printf("%.3f %.3f %.3f ", r, g, b);
				}
			}
			break;
		case SimulationConfig::ColorSurface:
			for (const auto& particle : m_particles) {
				if (particle.num_neighbors >= m_config.min_neighborhood_size) {
					f_color_surface(particle);
					file.printf("%.3f %.3f %.3f ", r, g, b);
				}
			}
			break;
		default:
			break;
		}
		file.printf("]\n");
		break;

	case SimulationConfig::FormatSphere:
	{
		float radius = m_config.output_constantwidth / 2.;

		for(const auto& particle : m_particles) {
			if (particle.num_neighbors >= m_config.min_neighborhood_size) {
				file.printf("AttributeBegin\n");
				if (m_config.output_format != SimulationConfig::FormatSurface) {
					switch (m_config.output_color) {
					case SimulationConfig::ColorPressure:
						f_color_pressure(particle);
						break;
					case SimulationConfig::ColorTemperature:
						f_color_temperature(particle);
						break;
					case SimulationConfig::ColorSurface:
						f_color_surface(particle);
						break;
					default:
						break;
					}
					file.printf("Color[%.3f %.3f %.3f]\n", r, g, b);
				}
				if (m_config.output_color == SimulationConfig::ColorShader) {
					file.printf("Translate %.7lf %.7lf %.7lf\n"
							"Sphere %.4f -%.4f %.4f 360\n\"Temp\"[%.4f]\nAttributeEnd\n",
							(double)particle.position.x, (double)particle.position.y,
							(double)particle.position.z,
							radius, radius, radius, (float)particle.temperature * temperature_scaling);
				} else {
					file.printf("Translate %.7lf %.7lf %.7lf\n"
							"Sphere %.4f -%.4f %.4f 360\nAttributeEnd\n",
							(double)particle.position.x, (double)particle.position.y,
							(double)particle.position.z,
							radius, radius, radius);
				}
			}
		}
	}
		break;
	case SimulationConfig::FormatDisk:
	case SimulationConfig::FormatSurface: /* surface also uses disks */
	{
		float radius = m_config.output_constantwidth / 2.;

		for(const auto& particle : m_particles) {
			if (particle.num_neighbors >= m_config.min_neighborhood_size) {
				file.printf("AttributeBegin\n");
				if (m_config.output_format != SimulationConfig::FormatSurface) {
					switch (m_config.output_color) {
					case SimulationConfig::ColorPressure:
						f_color_pressure(particle);
						break;
					case SimulationConfig::ColorTemperature:
						f_color_temperature(particle);
						break;
					case SimulationConfig::ColorSurface:
						f_color_surface(particle);
						break;
					default:
						break;
					}
					file.printf("Color[%.3f %.3f %.3f]\n", r, g, b);
				}
				//rotate z axis to -particle.gradient
				Vec3f zaxis(0, 0, 1), grad(0, 1, 0);
				if (particle.density_gradient.length2() > Math::FEQ_EPS)
					grad = -particle.density_gradient.normalized();
				float angle = acos(dot(zaxis, grad))*(180./M_PI);
				Vec3f rot_dir = cross(zaxis, grad);
				if (m_config.output_color == SimulationConfig::ColorShader) {
					file.printf("Translate %.7lf %.7lf %.7lf\nRotate %f %f %f %f\n"
							"Disk %.4f %.4f 360\n\"Temp\"[%.4f]\nAttributeEnd\n",
							(double)particle.position.x, (double)particle.position.y,
							(double)particle.position.z,
							angle, rot_dir.x, rot_dir.y, rot_dir.z, radius/3.f, radius,
							(float)particle.temperature * temperature_scaling);
				} else {
					file.printf("Translate %.7lf %.7lf %.7lf\nRotate %f %f %f %f\n"
							"Disk %.4f %.4f 360\nAttributeEnd\n",
							(double)particle.position.x, (double)particle.position.y,
							(double)particle.position.z,
							angle, rot_dir.x, rot_dir.y, rot_dir.z, radius/3.f, radius);
				}
			}
		}
	}
		break;
	}
}

void Simulation::addParticlesOnGrid(const Math::Vec3f& min_pos,
		const Math::Vec3f& max_pos, const Math::Vec3i& counts,
		const Math::Vec3f& initial_velocity, dfloat temperature, bool calc_mass,
		bool print_avg_density) {
	dfloat x_min = min_pos.x, x_max = max_pos.x;
	dfloat y_min = min_pos.y, y_max = max_pos.y;
	dfloat z_min = min_pos.z, z_max = max_pos.z;
	dfloat dx=0, dy=0, dz=0;
	Vec3f p;
	int num_particles_before = (int)m_particles.size();
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
	int num_particles_added = (int)m_particles.size()-num_particles_before;
	if(calc_mass) {
		m_config.particle_mass = dx * dy * dz * m_config.rho0;
	}

	if (print_avg_density) {
		dfloat avg_density = 0;
		for(int particle_idx = num_particles_before; particle_idx < (int)m_particles.size();
				++particle_idx) {
			Particle& particle = m_particles[particle_idx];
			dfloat kernel_sum = 0;
			for (int i = 0; i < (int)m_particles.size(); ++i) {
				kernel_sum += m_kernel.eval(particle.position, m_particles[i].position);
			}
			avg_density += m_kernel.kernelWeight() * m_config.particle_mass * kernel_sum;
		}
		avg_density /= num_particles_added;
		printf("Added %i grid particles with average density %.3f\n",
				num_particles_added, avg_density);
	}

	if (temperature > m_max_temperature) m_max_temperature = temperature;
}

bool Simulation::handleErruptions(dfloat simulation_time, dfloat dt, bool try_to_defer) {
	bool changed_particles = false;

	//new active
	auto& erruptions = m_config.erruptions;
	while (!erruptions.empty() && erruptions.back().start_time <= simulation_time) {
		m_active_erruptions.push_back(erruptions.back());
		erruptions.pop_back();
		LOG(DEBUG, "Starting Erruption");
	}
	int deferred = 0;
	if (try_to_defer && m_deferred_erruption_steps < m_max_deferred_erruption_steps)
		++deferred;

	//iterate active & check finished
	for (auto active = m_active_erruptions.begin(); active != m_active_erruptions.end(); ) {
		if (active->start_time + active->duration < simulation_time) {
			active = m_active_erruptions.erase(active);
			LOG(DEBUG, "Ending Erruption");
		} else {
			//add particles if necessary
			if (deferred > 0) {
				active->particles_left += dt * active->particles_per_sec;
				++deferred;
			} else {
				dfloat num_particles = dt * active->particles_per_sec
						+ active->particles_left;
				int inum_particles = (int)num_particles;
				active->particles_left = num_particles - (dfloat)inum_particles;
				for (int i = 0; i < inum_particles; ++i) {
					dfloat u = m_rand.nextf(), v = m_rand.nextf();
					Vec3f position = active->source->getPosition(m_height_field, u, v);
					Vec3f velocity = getInitVelocity(*active);
					addParticle(position, velocity, active->init_temperature);
				}
				changed_particles = true;
			}
			++active;
		}
	}
	if (deferred > 1) ++m_deferred_erruption_steps;
	else m_deferred_erruption_steps = 0;

	return changed_particles;
}
Math::Vec3f Simulation::getInitVelocity(const ErruptionConfig& config) {
	Vec3f velocity = config.init_velocity;
	if (m_config.init_velocity_perturb_angle > 0.) {
		dfloat angle = m_rand.nextf() * m_config.init_velocity_perturb_angle
				* (M_PI / 180.);
		Vec3f dir;
		Warp::uniformSphere(&dir, m_rand.nextf(), m_rand.nextf());
		velocity.rotate(angle, dir);
	}
	return velocity;
}

Math::Vec3f ErruptionSourceBase::interpolateBetween(const Math::Vec2f& start,
		const Math::Vec2f& end, const HeightField& height_field,
		dfloat u, dfloat v) {
	Vec3f pos;
	pos.x = Math::lerp(start.x, end.x, u);
	pos.z = Math::lerp(start.y, end.y, v);
	pos.y =  Math::FEQ_EPS + m_y_offset;
	if (!m_absolute_offset) pos.y += height_field.lookup(pos.x, pos.z);
	return pos;
}

Math::Vec3f ErruptionSourceLineSegment::getPosition(
		const HeightField& height_field, dfloat u, dfloat v) {
	return interpolateBetween(m_start, m_end, height_field, u, u);
}
Math::Vec3f ErruptionSourceGrid::getPosition(
		const HeightField& height_field, dfloat u, dfloat v) {
	return interpolateBetween(m_start, m_end, height_field, u, v);
}
