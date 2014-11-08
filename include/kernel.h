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

#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "global.h"
#include "Math/Vec3.h"
#include <cmath>

/*! @file different weighting kernels for the SPH simulation */

/**
 ** class KernelPoly6
 *  (no gradient or laplacian)
 */
class KernelPoly6 {
public:
	KernelPoly6(dfloat max_dist) {
		m_max_dist2 = max_dist*max_dist;
		m_kernel_weight = 315./(64.*M_PI*pow(max_dist, 9));
	}

	/** evaluate the kernel weight at a given position */
	inline dfloat eval(const Math::Vec3f& center, const Math::Vec3f& pos) {
		dfloat dist2 = (center - pos).length2();
		if (dist2 > m_max_dist2) return 0.;
		dfloat diff = m_max_dist2 - dist2;
		return diff * diff * diff;
	}

	/** kernel weight applied after summation (indep of particles) */
	inline dfloat kernelWeight() const { return m_kernel_weight; }
private:
	dfloat m_kernel_weight;
	dfloat m_max_dist2;
};


/**
 ** class KernelSpiky
 *  (no laplacian)
 */
class KernelSpiky {
public:
	KernelSpiky(dfloat max_dist) {
		m_max_dist = max_dist;
		m_kernel_weight = 15./(M_PI*pow(max_dist, 6));
	}

	/** evaluate the kernel weight at a given position */
	inline dfloat eval(const Math::Vec3f& center, const Math::Vec3f& pos) {
		dfloat dist = (center - pos).length();
		if (dist > m_max_dist) return 0.;
		dfloat diff = m_max_dist - dist;
		return diff*diff*diff;
	}
	/** kernel weight applied after summation (indep of particles) */
	inline dfloat kernelWeight() const { return m_kernel_weight; }


	/** evaluate the gradient of the kernel weight at a given position */
	inline Math::Vec3f evalGrad(const Math::Vec3f& center, const Math::Vec3f& pos) {
		Math::Vec3f r = center - pos;
		//gradient is: -3(max_dist - ||r||)^2 * r / ||r||
		dfloat dist = r.length();
		if (dist > m_max_dist || dist < Math::FEQ_EPS) return Math::Vec3f(0.);
		dfloat diff = m_max_dist - dist;
		return diff*diff/dist * r;
	}

	/** kernel weight applied after summation of gradient (indep of particles) */
	inline dfloat kernelWeightGrad() const { return -3.*m_kernel_weight; }

private:
	dfloat m_kernel_weight;
	dfloat m_max_dist;
};


/**
 ** class KernelViscosity
 *  (no gradient)
 */
class KernelViscosity {
public:
	KernelViscosity(dfloat max_dist) {
		m_max_dist = max_dist;
		m_kernel_weight = 15./(2.*M_PI*pow(max_dist, 3));
		m_kernel_weight_laplace = 45./(M_PI*pow(max_dist, 6));
	}

	/** evaluate the kernel weight at a given position */
	inline dfloat eval(const Math::Vec3f& center, const Math::Vec3f& pos) {
		dfloat dist = (center - pos).length();
		if (dist > m_max_dist) return 0.;
		dfloat doverh = dist / m_max_dist;
		return -0.5*doverh*doverh*doverh + doverh*doverh + 0.5/doverh - 1.;
	}
	/** kernel weight applied after summation (indep of particles) */
	inline dfloat kernelWeight() const { return m_kernel_weight; }


	/** evaluate the laplacian of the kernel weight at a given position */
	inline dfloat evalLaplace(const Math::Vec3f& center, const Math::Vec3f& pos) {
		dfloat dist = (center - pos).length();
		return std::max((dfloat)0., m_max_dist - dist);
	}

	/** kernel weight applied after summation of the laplacian (indep of particles) */
	inline dfloat kernelWeightLaplace() const { return m_kernel_weight_laplace; }

private:
	dfloat m_kernel_weight;
	dfloat m_kernel_weight_laplace;
	dfloat m_max_dist;
};

#endif /* _KERNEL_H_ */
