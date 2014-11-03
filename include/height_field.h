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

#ifndef _HEIGHT_FIELD_H_
#define _HEIGHT_FIELD_H_

#include <tiny_obj_loader.h>
#include "global.h"
#include <iostream>
#include <string>


/**
 * Base class to store a height field. it is in the range [0,1] x [0,z]
 * (z is such that it matches the aspect ratio of the input field)
 * y axis is height & points upwards.
 */
class HeightField {
public:
	virtual ~HeightField();

	/**
	 * write the height field to a RIB (RenderMan Scene Description file).
	 * The height field is consistent to how it is used in the simulation.
	 * Form:
        PointsPolygons
            [ 3 ]            # 1 triangle
            [ 0 1 2 ]        # vertex indexes
            "P" [ 1.0 0.9999999403953552 -1.0  1.0 -1.0 -1.0  -1.0000001192092896 -0.9999998211860657 -1.0 ]
			"st" [ 0.5 0  0.5 1  0 1 ]
	 *
	 * @param file opened file
	 */
	void writeRIBFile(FILE* file);


	int width() const { return m_width; }
	int depth() const { return m_depth; }

	dfloat fieldWidth() const { return 1.; }
	dfloat fieldDepth() const { return m_field_depth; }


	//TODO: intersection test:
	//  - height check (mirror speed)
	//  - also check for particles leaving the field


	/**
	 * @param x [0, 1]
	 * @param z [0, fieldDepth()]
	 * @return the y value at position [x,z]
	 * it uses the same calculation as for rendering (RIB output)
	 */
	inline dfloat lookup(dfloat x, dfloat z);

	const dfloat& field(int x, int z) const { return m_field[x + m_width*z]; }
	dfloat& getField(int x, int z) { return m_field[x + m_width*z]; }
protected:
	HeightField() {}


	virtual void writeRIBTexCoords(FILE* file) = 0;

	/**
	 * create the field: called by the subclass
	 */
	void init(int width, int depth, dfloat max_height);

	dfloat* m_field = NULL;
	int m_width = 0;
	int m_depth = 0;

	dfloat m_field_depth; //depth in field coordinates
	dfloat m_max_height; //max y value: field will be scaled to be in this range (min is 0)
};


dfloat HeightField::lookup(dfloat x, dfloat z) {
	DEBUG_ASSERT(x >= 0. && x < 1., "x out of range: %f", (float)x);
	DEBUG_ASSERT(z >= 0. && z < m_field_depth, "z out of range: %f", (float)z);

	dfloat fx = x*(dfloat)(m_width-1);
	dfloat fz = z/m_field_depth*(dfloat)(m_depth-1);
	dfloat fidx_x = floor(fx);
	dfloat fidx_z = floor(fz);
	int idx_x = (int)fidx_x;
	int idx_z = (int)fidx_z;
	dfloat u = fx - fidx_x;
	dfloat v = fz - fidx_z;

	//interpolate triangle: test whether to use lower or upper
	//TODO: use this if triangle meshes used for rendering
	/*
	dfloat a, b, c;
	if (u > v) { //lower
		a = field(idx_x, idx_z);
		b = field(idx_x+1, idx_z);
		c = field(idx_x+1, idx_z+1);
		return a + u*(b-a) + v*(c-b);
	}
	//upper
	a = field(idx_x, idx_z);
	b = field(idx_x+1, idx_z+1);
	c = field(idx_x, idx_z+1);
	return a + u*(b-c) + v*(c-a);
	*/

	//bilinear interpolation
	return (1.-v)*((1.-u)*field(idx_x, idx_z) + u* field(idx_x+1, idx_z))
			+ v*((1-u)* field(idx_x, idx_z+1) + u* field(idx_x+1, idx_z+1));
}



/**
 * Height field from a wavefront object file
 */
class HeightFieldObj : public HeightField {
public:
	/**
	 * @param shapes       wavefront objects. assumes y axis points upwards
	 * @param field_size   field resolution. shapes will be resamples with this.
	 * @param max_height   rescale the field to y values between [0, max_height]
	 */
	HeightFieldObj(const std::vector<tinyobj::shape_t>& shapes, int field_size=2000,
			dfloat max_height=0.6);

	virtual ~HeightFieldObj() {}
protected:
};


/**
 * Height field from a TIFF image
 */
class HeightFieldTiff : public HeightField {
public:
	/**
	 * @param tiff_file    a tiff file with one channel, 8 or 16 bits
	 *                     image size defines the field resolution
	 * @param max_height   rescale the field to y values between [0, max_height]
	 * @param step_x       use decreased resolution of tiff file (1 for full resolution)
	 * @param step_y
	 */
	HeightFieldTiff(const std::string& tiff_file, dfloat max_height,
			int step_x=1, int step_y=1);

	virtual ~HeightFieldTiff() {}
protected:
	virtual void writeRIBTexCoords(FILE* file);
};


#endif /* _HEIGHT_FIELD_H_ */
