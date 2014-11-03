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

#include "height_field.h"
#include <tiffio.h>
#include <cmath>

void HeightField::writeRIBFile(FILE* file) {

	fprintf(file, "PointsPolygons\n");
	int num_triangles = (m_width-1)*(m_depth-1)*2;
	fprintf(file, " [ ");
	// triangles
	for (int i = 0; i < num_triangles; ++i) {
		fprintf(file, "3 ");
	}
	fprintf(file, "]\n [ ");

	//vertex indices
	for (int z = 0; z < m_depth - 1; ++z) {
		for (int x = 0; x < m_width - 1; ++x) {
			int base_idx = z*m_width + x;
			//2 triangles, CCW
			fprintf(file, "%i %i %i %i %i %i ",
					base_idx, base_idx+1, base_idx+m_width+1,
					base_idx, base_idx+m_width+1, base_idx+m_width);
		}
	}
	fprintf(file, "]\n \"P\" [ ");

	//vertex positions
	for (int z = 0; z < m_depth; ++z) {
		for (int x = 0; x < m_width; ++x) {
			fprintf(file, "%.7lf %.7lf %.7lf ",
					(double)x / (double)(m_width - 1), field(x, z),
					(double)z / (double)(m_depth - 1) * m_field_depth);
		}
	}

	fprintf(file, "]\n");

	writeRIBTexCoords(file);
}

void HeightField::init(int width, int depth, dfloat max_height) {
	m_width = width;
	m_depth = depth;
	m_max_height = max_height;
	m_field_depth = (dfloat)depth / width;
	if (m_field) delete[] m_field;
	m_field = new dfloat[m_width * m_depth];
}

HeightField::~HeightField() {
	if (m_field) delete[] m_field;
}


HeightFieldObj::HeightFieldObj(const std::vector<tinyobj::shape_t>& shapes,
		int field_size, dfloat max_height) {
	//TODO

	/*
	//test code
std::string inputfile = "test.obj";
std::vector<tinyobj::material_t> materials;

std::string err = tinyobj::LoadObj(shapes, materials, inputfile.c_str());

if (!err.empty()) {
  std::cerr << err << std::endl;
  exit(1);
}

std::cout << "# of shapes    : " << shapes.size() << std::endl;
std::cout << "# of materials : " << materials.size() << std::endl;

for (size_t i = 0; i < shapes.size(); i++) {
  printf("shape[%ld].name = %s\n", i, shapes[i].name.c_str());
  printf("Size of shape[%ld].indices: %ld\n", i, shapes[i].mesh.indices.size());
  printf("Size of shape[%ld].material_ids: %ld\n", i, shapes[i].mesh.material_ids.size());
//assert((shapes[i].mesh.indices.size() % 3) == 0);
  for (size_t f = 0; f < shapes[i].mesh.indices.size() / 3; f++) {
    printf("  idx[%ld] = %d, %d, %d. mat_id = %d\n", f, shapes[i].mesh.indices[3*f+0], shapes[i].mesh.indices[3*f+1], shapes[i].mesh.indices[3*f+2], shapes[i].mesh.material_ids[f]);
  }

  printf("shape[%ld].vertices: %ld\n", i, shapes[i].mesh.positions.size());
  //assert((shapes[i].mesh.positions.size() % 3) == 0);
  for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
    printf("  v[%ld] = (%f, %f, %f)\n", v,
      shapes[i].mesh.positions[3*v+0],
      shapes[i].mesh.positions[3*v+1],
      shapes[i].mesh.positions[3*v+2]);
  }
}

for (size_t i = 0; i < materials.size(); i++) {
  printf("material[%ld].name = %s\n", i, materials[i].name.c_str());
  printf("  material.Ka = (%f, %f ,%f)\n", materials[i].ambient[0], materials[i].ambient[1], materials[i].ambient[2]);
  printf("  material.Kd = (%f, %f ,%f)\n", materials[i].diffuse[0], materials[i].diffuse[1], materials[i].diffuse[2]);
  printf("  material.Ks = (%f, %f ,%f)\n", materials[i].specular[0], materials[i].specular[1], materials[i].specular[2]);
  printf("  material.Tr = (%f, %f ,%f)\n", materials[i].transmittance[0], materials[i].transmittance[1], materials[i].transmittance[2]);
  printf("  material.Ke = (%f, %f ,%f)\n", materials[i].emission[0], materials[i].emission[1], materials[i].emission[2]);
  printf("  material.Ns = %f\n", materials[i].shininess);
  printf("  material.Ni = %f\n", materials[i].ior);
  printf("  material.dissolve = %f\n", materials[i].dissolve);
  printf("  material.illum = %d\n", materials[i].illum);
  printf("  material.map_Ka = %s\n", materials[i].ambient_texname.c_str());
  printf("  material.map_Kd = %s\n", materials[i].diffuse_texname.c_str());
  printf("  material.map_Ks = %s\n", materials[i].specular_texname.c_str());
  printf("  material.map_Ns = %s\n", materials[i].normal_texname.c_str());
  std::map<std::string, std::string>::const_iterator it(materials[i].unknown_parameter.begin());
  std::map<std::string, std::string>::const_iterator itEnd(materials[i].unknown_parameter.end());
  for (; it != itEnd; it++) {
    printf("  material.%s = %s\n", it->first.c_str(), it->second.c_str());
  }
  printf("\n");
}
	 */
}


HeightFieldTiff::HeightFieldTiff(const std::string& tiff_file,
		dfloat max_height, int step_x, int step_y) {

	TIFF* tif = TIFFOpen(tiff_file.c_str(), "r");
	if(!tif) throw EXCEPTION_s(EFILE_ERROR, "failed to open TIFF file %s\n", tiff_file.c_str());

	uint32 image_height, image_width;
	tdata_t buf;

	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &image_height);
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &image_width);

	uint16 samples_per_pix;
	TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pix); //1 meaning grayscale/monochrome
	if(samples_per_pix != 1) {
		TIFFClose(tif);
		throw EXCEPTION_s(EFILE_PARSING_ERROR, "TIFF file contains more than one channel!\n");
	}

	init(image_width/step_x, image_height/step_y, max_height);

	uint16 bits_per_sample;
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);

	tsize_t scanline = TIFFScanlineSize(tif);
	buf = _TIFFmalloc(scanline);


	dfloat max_value = 0;
	if (bits_per_sample == 16) {
		for (int row = 0; row < (int) image_height; row+=step_y) {
			int z_idx = (image_height - row - 1)/step_y;
			if (z_idx >= m_depth) continue;
			TIFFReadScanline(tif, buf, row);
			uint16* b = (uint16*) buf;
			for (int col = 0; col < scanline / ((int) bits_per_sample / 8); col+=step_x) {
				if(col/step_x >= m_width) break;
				dfloat val = (dfloat) b[col];
				getField(col/step_x, z_idx) = val;
				if (val > max_value) max_value = val;
			}
		}
	} else if (bits_per_sample == 8) {
		for (int row = 0; row < (int) image_height; row+=step_y) {
			int z_idx = (image_height - row - 1)/step_y;
			if (z_idx >= m_depth) continue;
			TIFFReadScanline(tif, buf, row);
			uint8* b = (uint8*) buf;
			for (int col = 0; col < scanline / ((int) bits_per_sample / 8); col+=step_x) {
				if(col/step_x >= m_width) break;
				dfloat val = (dfloat) b[col];
				getField(col/step_x, z_idx) = val;
				if (val > max_value) max_value = val;
			}
		}
	} else {
		_TIFFfree(buf);
		TIFFClose(tif);
		throw EXCEPTION_s(EFILE_PARSING_ERROR, "TIFF file contains unsupported depth (%i)\n", (int)bits_per_sample);
	}

	_TIFFfree(buf);
	TIFFClose(tif);

	//normalize height
	if (max_value > 0.) {
		dfloat normalize_factor = m_max_height / max_value;
		for (int z = 0; z < m_depth; ++z) {
			for (int x = 0; x < m_width; ++x) {
				getField(x, z) *= normalize_factor;
			}
		}
	}
}

void HeightFieldTiff::writeRIBTexCoords(FILE* file) {
	fprintf(file, " \"st\" [ ");
	for (int z = 0; z < m_depth; ++z) {
		for (int x = 0; x < m_width; ++x) {
			//u, v coords
			fprintf(file, "%.5lf %.5lf ", (double)x/(m_width-1), 1.-(double)z/(m_depth-1));
		}
	}
	fprintf(file, "]\n");
}
