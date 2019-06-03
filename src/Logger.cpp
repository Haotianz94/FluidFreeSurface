//#####################################################################
// Copyright 2017. Haotian Zhang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Logger.h"
#include "Configer.h"
#include <ctime>
#include <cstdio>

using namespace std;
using namespace cv;

Logger* Logger::_logger = NULL;
bool Logger::out_to_file = true;
bool Logger::out_to_console = true;

Logger::Logger()
{
	out_to_file = false;
	out_to_console = true;

	if(out_to_file)
	{
		time_t now = time(0);
		char logname[50];
		sprintf(logname, "../log/%ld.log", now);
		_logStream.open(logname);
	}
}

Logger::~Logger()
{
	_logStream << flush;
	if(out_to_file)
		_logStream.close();
}

Logger* Logger::getLogger()
{
	if (!_logger)
	{
	  	_logger = new Logger();
	}
	return _logger;
}

Logger& operator << (Logger& logger, const int num)
{
	if(Logger::out_to_console)
		cout << num;
	if(Logger::out_to_file)
	{
		logger._logStream << num;
		logger._logStream << flush;
	}
	return logger;
}

Logger& operator << (Logger& logger, const char* str)
{
	if(Logger::out_to_console)
		cout << str;
	if(Logger::out_to_file)
	{
		logger._logStream << str;
		logger._logStream << flush;
	}
	return logger;
}

Logger& operator << (Logger& logger, const char c)
{
	if(Logger::out_to_console)
		cout << c;
	if(Logger::out_to_file)
	{
		logger._logStream << c;
		logger._logStream << flush;
	}
	return logger;
}

Logger& operator << (Logger& logger, const double num)
{
	if(Logger::out_to_console)
		cout << num;
	if(Logger::out_to_file)
	{
		logger._logStream << num;
		logger._logStream << flush;
	}
	return logger;
}

Logger& operator << (Logger& logger, const float num)
{
	if(Logger::out_to_console)
		cout << num;
	if(Logger::out_to_file)
	{
		logger._logStream << num;
		logger._logStream << flush;
	}
	return logger;
}

Logger& operator << (Logger& logger, const CvPoint pos)
{
	if(Logger::out_to_console)
		cout << '(' << pos.x << ',' << pos.y << ')';
	if(Logger::out_to_file)
	{
		logger._logStream << '(' << pos.x << ',' << pos.y << ')';
		logger._logStream << flush;
	}
	return logger;
}

