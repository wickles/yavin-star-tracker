#include <ctime>

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