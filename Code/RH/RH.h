#ifndef RUN_RH_H
#define RUN_RH_H
#include "hyperplane_set.h"
#include "Partition.h"
#include "choose_item.h"


int RH(point_set* pset, point_t* u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound, bool error);


#endif //RUN_RH_H
