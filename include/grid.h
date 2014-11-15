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

#ifndef _GRID_H_
#define _GRID_H_

#include <vector>
#include <functional>
#include <array>
#include <stddef.h>
#include <string.h>

#include "height_field.h"
#include "Math/Vec3.h"

/**
 ** class Grid
 *  3D grid where each grid cell contains a list of Entry types.
 *
 *  details: x, z dimensions are given by the height_field & together with
 *  cell_size define how many grid points each dimension has (actual cell_size
 *  will be rounded s.t. ie nr x cells  * cell_size = x dimension)
 *  y dimension is special: at each (x,z) grid location, height_field determines
 *  at which y location the grid starts in an upwards direction (since the
 *  height_field is a constant surface, no particles can go below)
 *
 *  updateEntries must be called before any other manipulations!
 */
template<class Entry>
class Grid {
public:
	/**
	 * @param height_field
	 * @param cell_size         desired cell_size in all 3 dimensions
	 * @param num_y_cells       number of stacked cells in y direction
	 */
	Grid(HeightField& height_field, dfloat cell_size, int num_y_cells);
	~Grid();

	/**
	 * update all entries: this can include changed positions but also
	 * additions or removals of entries.
	 * Not thread-safe!
	 * @param entries    iterable (array) of updated entries
	 */
	template<class Iterable=std::vector<Entry> >
	void updateEntries(Iterable& entries);


	/**
	 * iterate through all neighbors and call the callback for the neighbor.
	 * thread-safe
	 * @param position    where to search
	 * @param max_dist    maximal distance
	 * @param callback    called with the entry & squared distance to position
	 */
	inline void iterateNeighbors(const Math::Vec3f& position, dfloat max_dist,
			const std::function<void (Entry*, dfloat dist2)>& callback) const;

	/**
	 * make sure a position is within the grid. this only checks the y coordinate!
	 * @param position
	 * @return          true if position changed & was outside the grid
	 */
	inline bool moveInsideGrid(Math::Vec3f& position) const;

private:
	/** convert 3d position in scale [0,1] -> grid (index) */
	inline Entry*& grid(const Math::Vec3f& pos) const;
	inline Entry*& grid(int x, int y, int z) const;
	inline void gridIndex(const Math::Vec3f& pos, int& x, int& y, int& z) const;
	inline int gridIndex(int x, int y, int z) const;

	inline dfloat& elevation(const Math::Vec3f& pos) const;
	inline dfloat& elevation(int x, int z) const;
	inline void elevationIndex(const Math::Vec3f& pos, int& x, int& z) const;
	inline int elevationIndex(int x, int z) const;

	Entry** m_grid = NULL; /** 3d grid: linked list of entries */
	int m_width;
	int m_height;
	int m_depth;
	dfloat* m_elevations = NULL; /** 2d grid with elevations at each (x,z) location */

	Math::Vec3f m_cell_size;
};

template<class Entry>
Grid<Entry>::Grid(HeightField& height_field, dfloat cell_size, int num_y_cells) {

	//cell sizes & dimension
	m_width = (int)(height_field.fieldWidth() / cell_size);
	m_cell_size.x = height_field.fieldWidth() / (dfloat)m_width;
	m_height = num_y_cells;
	m_cell_size.y = cell_size;
	m_depth = (int)(height_field.fieldDepth() / cell_size);
	m_cell_size.z = height_field.fieldDepth() / (dfloat)m_depth;

	LOG(DEBUG, "3D Grid Dimension: w=%i h=%i d=%i cell_size=(%.5f, %.5f, %.5f)",
		m_width, m_height, m_depth,
		(float)m_cell_size.x, (float)m_cell_size.y, (float)m_cell_size.z);

	//init elevations
	m_elevations = new dfloat[m_width * m_depth];
	dfloat cell_y_offset = -2.*m_cell_size.y; //safety gap, to make sure no
								//particle position will ever be below the grid
	//query the heightfield in a (coarse) grid to find the lowest value inside
	//each grid cell
	std::array<dfloat, 3> query_locations{{ 0.1, 0.5, 0.9 }};
	for (int z = 0; z < m_depth; ++z) {
		for (int x = 0; x < m_width; ++x) {
			dfloat min_y_val = 1.e12;
			for (std::size_t zz = 0; zz < query_locations.size(); ++zz) {
				dfloat z_pos = ((dfloat) z + query_locations[zz]) * m_cell_size.z;
				for (std::size_t xx = 0; xx < query_locations.size(); ++xx) {
					dfloat x_pos = ((dfloat) x + query_locations[xx]) * m_cell_size.x;
					dfloat val = height_field.lookup(x_pos, z_pos);
					if (val < min_y_val) min_y_val = val;
				}
			}
			elevation(x, z) = min_y_val + cell_y_offset;
		}
	}

	m_grid = new Entry*[m_width * m_height * m_depth];
}

template<class Entry>
Grid<Entry>::~Grid() {
	if (m_grid) delete[] (m_grid);
	if (m_elevations) delete[] (m_elevations);
}

template<class Entry> template<class Iterable>
void Grid<Entry>::updateEntries(Iterable& entries) {
	memset(m_grid, 0, sizeof(Entry*) * m_width * m_height * m_depth);
	for(auto& entry : entries) {
		Entry*& grid_entry = grid(entry.position);
		entry.next_in_grid = grid_entry;
		grid_entry = &entry;
	}
}

template<class Entry>
inline void Grid<Entry>::iterateNeighbors(const Math::Vec3f& position,
		dfloat max_dist, const std::function<void(Entry*, dfloat)>& callback) const {

	int x_min, y_min, z_min, x_max, y_max, z_max;
	dfloat max_dist2 = max_dist*max_dist;

	x_min = std::max(0, (int)((position.x-max_dist) / m_cell_size.x));
	z_min = std::max(0, (int)((position.z-max_dist) / m_cell_size.z));
	x_max = std::min(m_width-1, (int)((position.x+max_dist) / m_cell_size.x));
	z_max = std::min(m_depth-1, (int)((position.z+max_dist) / m_cell_size.z));

	for (int z = z_min; z <= z_max; ++z) {
		for (int x = x_min; x <= x_max; ++x) {
			dfloat y_val = position.y - elevation(x, z);
			y_min = std::max(0, (int)((y_val-max_dist) / m_cell_size.y));
			y_max = std::min(m_height-1, (int)((y_val+max_dist) / m_cell_size.y));
			for (int y = y_min; y <= y_max; ++y) {
				Entry* e = grid(x, y, z);
				while (e) {
					dfloat dist2 = (e->position-position).length2();
					if (dist2 <= max_dist2)
						callback(e, dist2);
					e = e->next_in_grid;
				}
			}
		}
	}
}

template<class Entry>
inline void Grid<Entry>::gridIndex(const Math::Vec3f& pos, int& x, int& y, int& z) const {

	x = (int)(pos.x / m_cell_size.x);
	z = (int)(pos.z / m_cell_size.z);
	y = (int)((pos.y - elevation(x, z)) / m_cell_size.y);
	//TODO: try w/o elevation & larger y grid -> faster??

	DEBUG_ASSERT1(x>=0 && x<m_width && y>=0 && y<m_height && z>=0 && z<m_depth);
}

template<class Entry>
inline bool Grid<Entry>::moveInsideGrid(Math::Vec3f& position) const {
	int x = (int) (position.x / m_cell_size.x);
	int z = (int) (position.z / m_cell_size.z);
	dfloat elev = elevation(x, z);
	int y = (int) ((position.y - elev) / m_cell_size.y);
	if (y < 0) {
		position.y = Math::FEQ_EPS + elev;
		return true;
	}
	if (y >= m_height) {
		position.y = m_cell_size.y * (dfloat)(m_height - 1) - Math::FEQ_EPS + elev;
		return true;
	}
	return false;
}


template<class Entry>
inline int Grid<Entry>::gridIndex(int x, int y, int z) const {
	return x + (y + z * m_height) * m_width;
}
template<class Entry>
inline Entry*& Grid<Entry>::grid(int x, int y, int z) const {
	return m_grid[gridIndex(x, y, z)];
}

template<class Entry>
inline Entry*& Grid<Entry>::grid(const Math::Vec3f& pos) const {
	int x, y, z;
	gridIndex(pos, x, y, z);
	return m_grid[gridIndex(x, y, z)];
}

template<class Entry>
inline void Grid<Entry>::elevationIndex(const Math::Vec3f& pos, int& x, int& z) const {

	x = (int)(pos.x / m_cell_size.x);
	z = (int)(pos.z / m_cell_size.z);
	DEBUG_ASSERT1(x>=0 && x<m_width && z>=0 && z<m_depth);
}

template<class Entry>
inline int Grid<Entry>::elevationIndex(int x, int z) const {
	return x + z * m_width;
}
template<class Entry>
inline dfloat& Grid<Entry>::elevation(int x, int z) const {
	return m_elevations[elevationIndex(x, z)];
}

template<class Entry>
inline dfloat& Grid<Entry>::elevation(const Math::Vec3f& pos) const {
	int x, z;
	elevationIndex(pos, x, z);
	return m_elevations[elevationIndex(x, z)];
}

#endif /* _GRID_H_ */
