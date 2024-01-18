#include "hyperplane_set.h"

/**
 * @brief Constructor
 */
hyperplane_set::hyperplane_set(){}


/**
 * @brief Constructor: Initialize all the possible hyperplanes as questions
 * @param pset
 */
hyperplane_set::hyperplane_set(point_set *pset, point_t *exp, double RandRate)
{
    int M = pset->points.size();
    for(int i = 0; i < M; ++i)
    {
        for(int j = i + 1; j < M; ++j)
        {
            hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
            hyperplanes.push_back(h);
        }
    }

    //reinsert
    for (int i = 0; i < M * RandRate; i++)
    {
        int n = ((int) rand()) % M;
        hyperplane *h = hyperplanes[n];
        hyperplanes.erase(hyperplanes.begin() + n);
        hyperplanes.push_back(h);
    }
}


/**
 * @brief Check whether the hyperplane is in the set
 * @param h     The hyperplane
 * @return      True if the hyperplane is in the set
 *              False if the hyperplane is not in the set
 */
bool hyperplane_set::is_hyperplane_exist(hyperplane *h)
{
    for(int i = 0; i < hyperplanes.size(); ++i)
    {
        if (hyperplanes[i]->p_1->id == h->p_1->id && hyperplanes[i]->p_2->id == h->p_2->id)
            return true;
    }
    return false;
}




