#include "Groundtruth.h"

/**
 * @brief Find the nearest point of e
 * @param p_skyline The point set
 * @param e         The expected point
 */
point_set* ground_truth(point_set *pSet, point_t *u, int k, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    point_set *rankingSet = new point_set, *resultSet = new point_set();
    pSet->findRanking(u, rankingSet);
    //rankingSet->print_with_score(u);
    rankingSet->fairTopk(k, category, bound,resultSet);

    //rankingSet->print();
    std::cout << "-----------------------------------------------------------------------------------\n";
    printf("|%15s |%15s |%15s |%15s |%10d |\n", "Ground Truth", "-", "-", "-", resultSet->points[0]->id);
    for(int i = 1; i < k; ++i)
        printf("|%15s |%15s |%15s |%15s |%10d |\n", "-", "-", "-", "-", resultSet->points[i]->id);
    std::cout << "-----------------------------------------------------------------------------------\n";

    return resultSet;
}
