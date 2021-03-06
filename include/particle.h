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
	dfloat dtemperature; /** temperature update */

	Math::Vec3f forces; /** sum of all forces */

	/* PCISPH pressure loop */
	dfloat predicted_density;
	Math::Vec3f predicted_velocity;
	Math::Vec3f predicted_position;
	Math::Vec3f force_pressure;
	dfloat pressure;

	Math::Vec3f density_gradient;

	Particle* next_in_grid = NULL; /** linked list of particles in a grid cell */

	Particle** neighbors = NULL; /** pointer to an array with its neighbors */
	unsigned short num_neighbors = 0;

	unsigned char flags = 0;

	static constexpr unsigned char FLAG_IS_ON_GROUND = 1 << 0;
	static constexpr unsigned char FLAG_IS_AT_AIR = 1 << 1;

	bool isOnGround() const { return flags & FLAG_IS_ON_GROUND; }
	bool isAtAir() const { return flags & FLAG_IS_AT_AIR; }
	void setFlag(unsigned char flag) { flags |= flag; }
	void clearFlag(unsigned char flag) { flags &= ~flag; }
	void changeFlag(unsigned char flag, bool set)
		{ if (set) setFlag(flag); else clearFlag(flag); }
};


#endif /* _PARTICLE_H_ */
