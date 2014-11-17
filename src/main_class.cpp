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
#include "simulation.h"
#include "version.h"

using namespace Math;
using namespace pugi;

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>

/* creating directories */
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

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

	m_parameters->addParam("frames", 'f');
	
	
	m_cl_parse_result = m_parameters->parse();
	
}

void CMain::printHelp()
{
	printf("Usage:\n"
		   " " APP_NAME " [-v] -c <config> [--simulate] [-f <num>]\n"
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
		   " options\n"
		   "  -f, --frames <num_frames>       limit number of simulated frames\n"
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

	float height_scaling;
	if(hasSuffix(toLower(height_field_file), ".tif")
			|| hasSuffix(toLower(height_field_file), ".tiff")) {
		height_scaling = height_field_node.node().attribute("scaling").as_float(1.);
		xml_node tiff_child = height_field_node.node().child("tiff");
		int step_x = tiff_child.attribute("step_x").as_int(1);
		int step_y = tiff_child.attribute("step_y").as_int(1);
		m_height_field = new HeightFieldTiff(height_field_file, height_scaling, step_x, step_y);
	} else {
		//TODO: obj file??
		throw EXCEPTION_s(EINVALID_PARAMETER, "Error: Unsupported height field %s",
			height_field_file.c_str());
	}

	LOG(DEBUG, "height field: w=%i, d=%i, field depth=%.3f, height_scaling=%.3f",
		m_height_field->width(), m_height_field->depth(),
		(float)m_height_field->fieldDepth(), height_scaling);


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
		SimulationConfig config;
		xml_node config_node = doc.child("config");
		xml_node sim_node = config_node.child("simulation");
		//xml_node output_node = config_node.child("output");

		/* read the config */
		config.simulation_time = (dfloat)sim_node.attribute("simulation_time").as_double(60);
		config.smoothing_kernel_size = (dfloat)sim_node.attribute("smoothing_kernel_size").as_double(0.03);
		config.neighbor_lookup_dist =
			(dfloat)sim_node.attribute("lookup_dist").as_double(config.smoothing_kernel_size);
		if(config.neighbor_lookup_dist < config.smoothing_kernel_size)
			THROW_s(EINVALID_PARAMETER, "Config Error: lookup_dist < smoothing_kernel_size");
		config.cell_size = (dfloat)sim_node.attribute("cell_size").as_double(config.neighbor_lookup_dist);
		config.num_y_cells = sim_node.attribute("num_y_cells").as_int(20);
		config.k = (dfloat)sim_node.attribute("pressure_k").as_double(1000);
		config.rho0 = (dfloat)sim_node.attribute("rho0").as_double(1000);
		config.particle_mass = (dfloat)sim_node.attribute("particle_mass").as_double(0.01);
		istringstream(sim_node.attribute("gravity").as_string("0 -9.81 0")) >> config.g;
		config.time_step = (dfloat)sim_node.attribute("time_step").as_double(0.001);
		config.init_velocity_perturb_angle = (dfloat)sim_node.attribute(
				"init_velocity_perturb_angle").as_double(0.);
		config.output_rate = (int)sim_node.attribute("output_rate").as_int(1);
		string ground_method = sim_node.attribute("ground_method").as_string("elastic");
		if (ground_method == "elastic") {
			config.ground_method = SimulationConfig::GroundElastic;
		} else if (ground_method == "spring") {
			config.ground_method = SimulationConfig::GroundForceSpring;
			config.ground_spring = (dfloat) sim_node.attribute("ground_spring").as_double(1000.);
		} else {
			LOG(ERROR, "unknown ground_method (valid are: 'elastic', 'spring')");
			return;
		}

		//erruptions
		xml_node erruptions_node = sim_node.child("erruptions");
		for (xml_node erruption = erruptions_node.child("erruption"); erruption;
				erruption = erruption.next_sibling("erruption")) {
			ErruptionConfig erruption_config;
			erruption_config.start_time =
					(dfloat)erruption.attribute("start_time").as_double(0.);
			erruption_config.duration =
					(dfloat)erruption.attribute("duration").as_double(0.1);
			erruption_config.particles_per_sec =
					(dfloat)erruption.attribute("particles_per_sec").as_double(10000);
			erruption_config.init_temperature =
					(dfloat)erruption.attribute("init_temperature").as_double(1000);
			istringstream(erruption.attribute("init_velocity").as_string("0 0.1 0"))
					>> erruption_config.init_velocity;
			xml_node source_node;
			if ((source_node = erruption.child("source-line"))) {
				Vec2f start, end;
				istringstream(source_node.attribute("start").as_string("0 0")) >> start;
				istringstream(source_node.attribute("end").as_string("0 0")) >> end;
				dfloat y_offset = (dfloat)source_node.attribute("y_offset").as_double(0.);
				erruption_config.source = std::make_shared<ErruptionSourceLineSegment>(start, end, y_offset);
			} //else if: TODO: other types...

			config.erruptions.push_back(erruption_config);
		}

		//output directory
		string output_dir = "output/simulation";
		createDirectory(output_dir);
		string config_file_base = fileBaseName(config_file);
		config.output_dir = output_dir + "/" + config_file_base;
		createDirectory(config.output_dir);

		string sframes;
		if(m_parameters->getParam("frames", sframes)) {
			config.num_frames = atoi(sframes.c_str());
		}

		Simulation* simulation = new Simulation(*m_height_field, config);

		xml_node particles_grid = sim_node.child("particles_grid");
		if(!particles_grid.empty()) {
			Vec3f min_pos, max_pos, velocity;
			Vec3i counts;
			istringstream(particles_grid.attribute("min_pos").as_string("0 0 0")) >> min_pos;
			istringstream(particles_grid.attribute("max_pos").as_string("1 0.05 1")) >> max_pos;
			istringstream(particles_grid.attribute("count").as_string("100 5 100")) >> counts;
			istringstream(particles_grid.attribute("velocity").as_string("0 0 0")) >> velocity;
			dfloat temperature = (dfloat)particles_grid.attribute("temperature").as_double(100);
			bool calc_mass = particles_grid.attribute("calc_mass").as_bool(false);
			simulation->addParticlesOnGrid(min_pos, max_pos, counts, velocity,
					temperature, calc_mass);
		}

		simulation->run();
		LOG(DEBUG, "simulation ended");
		delete simulation;
	}
}

void CMain::createDirectory(const std::string& directory) {
	struct stat st;
	if (stat(directory.c_str(), &st) == 0) return;

	int status = mkdir(directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status) {
		THROW_s(EFILE_ERROR, "Failed to create directory %s (%s)",
				directory.c_str(), strerror(errno));
	}
}

std::string CMain::fileBaseName(const std::string& directory) {
	string tmp = directory;
	auto pivot = std::find(directory.rbegin(), directory.rend(), '/');
	if(pivot != directory.rend()) tmp = string(pivot.base(), directory.end());
	return removeExtension(tmp);
}

string CMain::removeExtension(std::string const& filename) {
	auto pivot = std::find(filename.rbegin(), filename.rend(), '.');
	return pivot == filename.rend() ?
			filename : std::string(filename.begin(), pivot.base() - 1);
}
