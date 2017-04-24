/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#include <ctime>

#pragma once

class Timer {
	clock_t start_time;
	clock_t end_time;
public:
	Timer();
	void StartTimer();
	void StopTimer();
	clock_t GetTime();
	double GetTimeSeconds();
};