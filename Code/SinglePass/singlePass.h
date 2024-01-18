#ifndef RUN_SINGLEPASS_H
#define RUN_SINGLEPASS_H
#include "hyperplane_set.h"
#include "Partition.h"
void singlePass(point_set *pset, point_t *u, int k, point_set *realResult,std::vector<std::string> &category,
                    std::vector<std::vector<double>> &bound, bool error);


#endif //RUN_SINGLEPASS_H
