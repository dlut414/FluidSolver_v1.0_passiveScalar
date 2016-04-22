/*
*/
#pragma once
#include <sys/timeb.h>

class Timer {
public:
	Timer() {
	}

	~Timer() {
	}

	float now() { return getSystemTime() / 1000.f; }
	float cost() { return stop - start; }

public:
	float start, stop;

private:
	long getSystemTime() {
		timeb t;
		ftime(&t);
		return long(1000 * t.time + t.millitm);
	}
};

