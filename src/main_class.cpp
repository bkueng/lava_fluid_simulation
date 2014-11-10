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

#include "main_class.h"
#include "pugixml.hpp"
#include "height_field.h"
#include "version.h"

using namespace pugi;

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>

/*********************************************************************//*
 * class CMain
 *//*********************************************************************/


CMain::CMain() : m_parameters(NULL), m_cl_parse_result(Parse_none_found)
{

}

CMain::~CMain()
{
	SAVE_DEL(m_parameters);
	SAVE_DEL(m_height_field);
}



void CMain::init(int argc, char* argv[])
{
	parseCommandLine(argc, argv);
	
}

void CMain::parseCommandLine(int argc, char* argv[])
{

	SAVE_DEL(m_parameters);
	m_parameters = new CCommandLineParser(argc, argv);
	
	//init known arguments
	m_parameters->addSwitch("help", 'h');
	m_parameters->addSwitch("version");
	m_parameters->addSwitch("verbose", 'v');
	//control the logging
	m_parameters->addParam("log");
	m_parameters->addSwitch("no-log");
	m_parameters->addParam("file-log");
	m_parameters->addSwitch("no-file-log");
	
	
	m_parameters->addParam("config", 'c');

	m_parameters->addSwitch("simulate", 's');
	m_parameters->addParam("generate-rib-heightfield", 'g');
	
	
	m_cl_parse_result = m_parameters->parse();
	
}

void CMain::printHelp()
{
	printf("Usage:\n"
		   " " APP_NAME " [-v] -c <config> [--simulate]\n"
		   " " APP_NAME " [-v] -c <config> [--generate-rib-heightfield <file>]\n"
		   " " APP_NAME " --version\n"
		   "\n"
		   "  -c, --config <config>           specify configuration file\n"
		   "                                  under ./config/\n"
		   " tasks\n"
		   "  -s, --simulate                  start simulation.\n"
		   "                                  this is the default if no task is given\n"
		   "  -g, --generate-rib-heightfield  <file>\n"
		   "                                  create a RIB file for rendering\n"
		   "                                  from the configured height field\n"
		   "                                  and write the result to <file> (eg height_field.rib)\n"
		   "\n"
		   "  -v, --verbose                   print debug messages\n"
		   "                                  (same as --log debug)\n"
		   "  -h, --help                      print this message\n"
		   "  --version                       print the version\n"
		   "\n"
		   " logging\n"
		   "  --log <level>                   set console log level\n"
		   "  --file-log <level>              set file log level\n"
		   "   <level>                        none, error, warn, info, debug\n"
		   "  --no-log                        no console logging (--log none)\n"
		   "  --no-file-log                   no file logging (--file-log none)\n"
		  );
}


void CMain::exec()
{

	ASSERT_THROW(m_parameters, ENOT_INITIALIZED);
	
	switch (m_cl_parse_result) {
	case Parse_none_found:
		processArgs();
		break;
	case Parse_unknown_command:
		wrongUsage("Unknown command: %s",
				   m_parameters->getUnknownCommand().c_str());
		break;
	case Parse_success:
		if (m_parameters->getSwitch("help")) {
			printHelp();
		} else if (m_parameters->getSwitch("version")) {
			printVersion();
		} else {
			processArgs();
		}
		break;
	}
}

void CMain::printVersion()
{
	printf("%s\n", getAppVersion().toStr().c_str());
}

void CMain::wrongUsage(const char* fmt, ...)
{

	printHelp();
	
	printf("\n ");
	
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	
	printf("\n");
	
}


void CMain::processArgs()
{

	//set console log level
	string level;
	if (m_parameters->getSwitch("verbose"))
		CLog::getInstance().setConsoleLevel(DEBUG);
	else if (m_parameters->getSwitch("no-log"))
		CLog::getInstance().setConsoleLevel(NONE);
	else if (m_parameters->getParam("log", level)) {
		ELOG log_level;
		if (CLog::parseLevel(level, log_level))
			CLog::getInstance().setConsoleLevel(log_level);
	}
	//set file log level
	if (m_parameters->getSwitch("no-file-log"))
		CLog::getInstance().setFileLevel(NONE);
	else if (m_parameters->getParam("file-log", level)) {
		ELOG log_level;
		if (CLog::parseLevel(level, log_level))
			CLog::getInstance().setFileLevel(log_level);
	}



	bool do_simulate = true;
	string config_file = "";
	if (!m_parameters->getParam("config", config_file)) {
		printf("Error: no config file given. use --help for the usage\n");
		return;
	}

	xml_document doc;
	xml_parse_result result = doc.load_file(config_file.c_str());
	if(!result)
		throw EXCEPTION_s(EINVALID_PARAMETER, "Error: Failed to open config file %s (%s)",
				config_file.c_str(), result.description());


	xpath_node height_field_node = doc.select_single_node("/config/simulation/heightfield");
	if(!height_field_node) throw EXCEPTION_s(EFILE_PARSING_ERROR, "Error: config parsing error");
	string height_field_file = height_field_node.node().attribute("file").as_string();

	float max_height;
	if(hasSuffix(toLower(height_field_file), ".tif")
			|| hasSuffix(toLower(height_field_file), ".tiff")) {
		max_height = height_field_node.node().attribute("max_height").as_float(1.);
		xml_node tiff_child = height_field_node.node().child("tiff");
		int step_x = tiff_child.attribute("step_x").as_int(1);
		int step_y = tiff_child.attribute("step_y").as_int(1);
		m_height_field = new HeightFieldTiff(height_field_file, max_height, step_x, step_y);
	} else {
		//TODO: obj file??
		throw EXCEPTION_s(EINVALID_PARAMETER, "Error: Unsupported height field %s",
			height_field_file.c_str());
	}

	LOG(DEBUG, "height field: w=%i, d=%i, field depth=%f, max_height=%f",
		m_height_field->width(), m_height_field->depth(),
		(float)m_height_field->fieldDepth(), max_height);


	string rib_file;
	if(m_parameters->getParam("generate-rib-heightfield", rib_file)) {
		do_simulate = false;
		FILE* file = fopen(rib_file.c_str(), "w");
		if (!file) throw EXCEPTION(EFILE_ERROR);

		LOG(DEBUG, "writing RIB heightfield");
		m_height_field->writeRIBFile(file);
		fclose(file);
	}

	if(do_simulate || m_parameters->getSwitch("simulate")) {
		LOG(DEBUG, "do simulation");
		//TODO
	}
}










