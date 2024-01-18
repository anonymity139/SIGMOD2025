#ifndef RUN_PIECEWISE_EXHAUSTION_H
#define RUN_PIECEWISE_EXHAUSTION_H
#include "hyperplane_set.h"
#include "u_vector.h"

void Piecewise_exhaustion(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound);

void Piecewise_exhaustion_robust(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound);

#endif //RUN_PIECEWISE_EXHAUSTION_H
