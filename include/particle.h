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

#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "Math/Vec3.h"
#include "global.h"

struct Particle {

	Math::Vec3f position; /** in [0,1] coordinate space */
	Math::Vec3f velocity;

	dfloat temperature;

	dfloat density;
	Math::Vec3f forces; /** sum of all forces */

	Particle* next_in_grid = NULL; /** linked list of particles in a grid cell */

	Particle** neighbors = NULL; /** pointer to an array with its neighbors */
	unsigned short num_neighbors = 0;
};


#endif /* _PARTICLE_H_ */
