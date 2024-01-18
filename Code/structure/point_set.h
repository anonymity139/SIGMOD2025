#ifndef POINT_SET_H
#define POINT_SET_H

#include "point_t.h"
#include <string>


class point_set
{
public:
    std::vector<point_t*> points;

    point_set();
    explicit point_set(point_set *p_set);
    explicit point_set(const std::string input);
    ~point_set();

    void print();
    void print_with_score(point_t *u);
    void random(double RandRate);
    point_set* sort(point_t *u);
    void findTopk(point_t *u, int k, point_set *topSet);
    void findRanking(point_t *u, point_set *rankingSet);
    void sort_point(std::vector<point_t*> &return_point, point_t *u);
    void sort_point(std::vector<point_t*> &return_point);
    void fairTopk(const int k, const std::vector<std::string> &category, std::vector<std::vector<double>> bound, point_set *topSet);
    void markFairTopk(const int k, const std::vector<std::string> &category, std::vector<std::vector<double>> bound, int *numOfEachCluster);
    void write(std::string fileName);
    void prunePt(point_t* p);
    bool is_prune(point_t *p);
    bool isSame(point_set *pset);
    bool isSame_exact(point_set *pset);
    point_set* findsame(point_set *pset);
    bool checkExist(point_t *p);
    void printResult(char *name, int Qcount, int s,timeval t1, double preTime, long mem_baseline);
    void printResult(char *name, int Qcount, point_set* realResultSet, point_t *u, timeval t1, double preTime, long mem_baseline);
    double PercentageOfSamePoint(point_set* pset);
    double PercentageOfUtility(point_set* pset, point_t *u);
    void skyband(point_set *returnSet, int k);
    void skyband_with_clusters(std::vector<std::string> &category, int k, std::vector<std::vector<double>> &bound);
    int findBest(point_t *u);
    void findTopk_sampling(std::vector<point_set*> &topSet, point_t *u, int k, int level, int used_seg);
    void findRanking_sampling(std::vector<point_set*> &rankSet, point_t *u, int level, int used_seg);
    void findTopk_sampling(std::vector<point_set*> &topSet, double *max, double *min, point_t *u, int k, int level, int used_seg);
    void subtract(point_set *pset);
};


#endif //U_2_POINT_SET_H
