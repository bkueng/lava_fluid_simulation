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

#ifndef _OUTPUT_FILE_H_
#define _OUTPUT_FILE_H_

#include <string>
#include <stdarg.h>
#include <stdio.h>

#ifdef USE_COMPRESSION
#include <zlib.h>
typedef gzFile FileHandle;
#else
typedef FILE* FileHandle;
#endif

#include "exception.h"

/**
 ** class OutputFile
 * write output files, that can be compressed using gzip (defined at compile time)
 */
class OutputFile {
public:
	OutputFile(const std::string& file_name);
	~OutputFile();

	template <int buffer_len=32>
	inline void printf(const char* format, ...);
private:
	FileHandle m_file_handle = NULL;
};


template<int buffer_len>
void OutputFile::printf(const char* format, ...) {
   va_list arg;
   char buffer[buffer_len];

   va_start(arg, format);
   int len = vsprintf(buffer, format, arg);
   va_end(arg);
   if (len < 0) THROW(EFILE_ERROR);

#ifdef USE_COMPRESSION
	gzwrite(m_file_handle, buffer, len);
#else
	fwrite(buffer, sizeof(char), len, m_file_handle);
#endif
}

#endif /* _OUTPUT_FILE_H_ */
