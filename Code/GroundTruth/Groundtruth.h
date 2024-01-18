#ifndef REGRET_GROUNDTRUTH_H
#define REGRET_GROUNDTRUTH_H
#include "point_set.h"
#include "structure/define.h"

point_set* ground_truth(point_set *pSet, point_t *u, int k, std::vector<std::string> &category, std::vector<std::vector<double>> &bound);

#endif //REGRET_GROUNDTRUTH_H
