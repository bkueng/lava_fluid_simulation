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

#include "output_file.h"

OutputFile::OutputFile(const std::string& file_name) {

#ifdef USE_COMPRESSION
	m_file_handle = gzopen(file_name.c_str(), "wb");
#else
	m_file_handle = fopen(file_name.c_str(), "w");
#endif
	if (!m_file_handle) throw EXCEPTION(EFILE_ERROR);
}

OutputFile::~OutputFile() {
	if (m_file_handle) {
#ifdef USE_COMPRESSION
		gzclose(m_file_handle);
#else
		fclose(m_file_handle);
#endif
	}
}
