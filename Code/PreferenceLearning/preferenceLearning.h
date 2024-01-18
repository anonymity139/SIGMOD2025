#ifndef RUN_PREFERENCELEARNING_H
#define RUN_PREFERENCELEARNING_H
#include "point_set.h"
#include "Partition.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include "s_node.h"
#include "cluster_t.h"
#include "../Others/QuadProg++.hh"
using namespace quadprogpp;

int PreferenceLearning(point_set *original_set, point_t *u, int k, point_set *realResult, std::vector<std::string> &category,
                       std::vector<std::vector<double>> &bound, bool error);


#endif //RUN_PREFERENCELEARNING_H
