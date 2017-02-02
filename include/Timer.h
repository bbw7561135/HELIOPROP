#ifndef HELIOPROP_TIMER_H_
#define HELIOPROP_TIMER_H_

#include <ctime>

class Timer
{
public:
	Timer() {
		reset();
	}
	void reset() {
		beg_ = clock();
	}
	double elapsed() {
		end_ = clock();
		return double(end_ - beg_) / CLOCKS_PER_SEC;
	}
private:
	clock_t beg_, end_;
};

#endif /* HELIOPROP_TIMER_H_ */
