#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <opencv2/opencv.hpp>

#define LOG *Logger::getLogger()
#define REPORT(X) *Logger::getLogger() << "[Logger]  : " << #X << " = " << X << "\n"
#define PRINT(X) *Logger::getLogger() << X << "\n"
#define PRINT_ONELINE(X) *Logger::getLogger() << X
#define LOGSEG *Logger::getLogger() << "================================================================\n"

class Logger
{
public:
	Logger();
	~Logger();
	static Logger* getLogger();

	friend Logger& operator << (Logger& logger, const int num);
	friend Logger& operator << (Logger& logger, const char* str);
	friend Logger& operator << (Logger& logger, const char c);
	friend Logger& operator << (Logger& logger, const double num);
	friend Logger& operator << (Logger& logger, const float num);
	friend Logger& operator << (Logger& logger, const CvPoint pos);
	template<typename T>
	friend Logger& operator << (Logger& logger, const cv::Mat_<T> m);

private:
	std::ofstream _logStream;
	static Logger* _logger;
	static bool out_to_file;
	static bool out_to_console;
};





template<typename T>
Logger& operator << (Logger& logger, const cv::Mat_<T> m)
{
	int w = m.cols;
	int h = m.rows;

	if(Logger::out_to_console)
	{
		std::cout << "\n";
		for(int y = 0; y < h; y++)
			for(int x = 0; x < w; x++)
			{
				std::cout << m(x, y) << ' ';
				if(x == w-1)
					std::cout << "\n";
			}
	}
	if(Logger::out_to_file)
	{
		logger._logStream << "\n";
		for(int y = 0; y < h; y++)
			for(int x = 0; x < w; x++)
			{
				logger._logStream << m(x, y) << ' ';
				if(x == w-1)
					logger._logStream << "\n";
			}
		logger._logStream << std::flush;
	}
	return logger;
}


#endif