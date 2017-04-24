/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
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