#ifndef RUN_ACTIVERANKING_H
#define RUN_ACTIVERANKING_H
#include "point_set.h"
#include "hyperplane_set.h"
#include "Partition.h"

int ActiveRanking(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound);

int ActiveRanking_robust(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound);
#endif //RUN_ACTIVERANKING_H
