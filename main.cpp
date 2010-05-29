/*
 * Copyright (C) 2010 Beat Küng <beat-kueng@gmx.net>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#include "main_class.h"


int main(int argc, char *argv[]) {
	try {
		
#ifdef _DEBUG
		CLog::getInstance().setConsoleLevel(DEBUG);
#else
		CLog::getInstance().setConsoleLevel(WARN);
#endif
		CLog::getInstance().setFileLevel(WARN);
		CLog::getInstance().setLogDateTime(false);
		CLog::getInstance().setLogSourceFileAll(false);
		CLog::getInstance().setLogSourceFile(ERROR, true);
		
		CMain main;
		main.init(argc, argv);
		main.exec();
	} catch(Exception& e) {
		return(-1);
	}
	return(0);
}

