#ifndef INCLUDE_ENERGYAXIS_H_
#define INCLUDE_ENERGYAXIS_H_

#include <cmath>
#include <vector>

class EnergyAxis {
public:
	EnergyAxis(const Input& input) {
		size = input.kenergyn_size;
		double logKE_min = log10(input.kenergyn_min);
		double logKE_max = log10(input.kenergyn_max);
		for (size_t i = 0; i < size; ++i) {
			double value = pow(10, logKE_min + double(i) / double(size - 1) * (logKE_max - logKE_min));
			axis.push_back(value);
		}
	}
	~EnergyAxis() {
		axis.clear();
	}
	size_t get_size() {
		return size;
	}
	double at(const size_t& i) {
		return axis.at(i);
	}
private:
	std::vector<double> axis;
	size_t size;
};

#endif /* INCLUDE_ENERGYAXIS_H_ */
