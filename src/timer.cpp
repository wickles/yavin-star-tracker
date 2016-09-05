/* timer.cpp
Alexander Wickes
September 20, 2010
*/

#include "timer.h"

Timer::Timer()
{
	StartTimer();
	StopTimer();
}

void Timer::StartTimer()
{
	start_time = clock();
}

void Timer::StopTimer()
{
	end_time = clock();
}

clock_t Timer::GetTime()
{
	if (end_time <= start_time)
		StopTimer();
	return (end_time - start_time);
}

double Timer::GetTimeSeconds()
{
	return GetTime()/double(CLOCKS_PER_SEC);
}