/*
 * Copyright (C) 2010-2011 Beat Küng <beat-kueng@gmx.net>
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

#ifndef MAIN_CLASS_H_
#define MAIN_CLASS_H_

#include "global.h"
#include "command_line.h"
#include "height_field.h"

#include <string>

/*********************************************************************//*
 * class CMain
 *
 * main class with the main task and command line parser
 *//*********************************************************************/


class CMain
{
public:
	CMain();
	~CMain();
	
	/* parse the command line parameters */
	void init(int argc, char* argv[]);
	
	void exec();
	
private:
	void parseCommandLine(int argc, char* argv[]);
	void printHelp();
	void wrongUsage(const char* fmt, ...);
	void printVersion();
	
	void processArgs();
	
	/** create a directory if it does not exist */
	void createDirectory(const std::string& directory);

	/** get base name of a file name (strip preceding directory & file extension */
	std::string fileBaseName(const std::string& directory);
	std::string removeExtension(std::string const& filename);

	HeightField* m_height_field = NULL;
	
	CCommandLineParser* m_parameters;
	ECLParsingResult m_cl_parse_result;
};



#endif /* MAIN_CLASS_H_ */
