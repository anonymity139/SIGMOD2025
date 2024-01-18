#include "PI.h"
#include <iostream>
using namespace std;

struct MinHeapCmp
{
    inline bool operator()(const u_vector* y, const u_vector* z)const
    {
        double v1=y->x;
        double v2=z->x;
        //printf("v1%lf v2%lf\n", v1, v2);
        if(v1 - v2<0.000000001&&v1 - v2>-0.000000001)
        {
            if(y->place == z->place)
            {
                return y->time < z->time;
            }
            else
            {
                return y->place > z->place;
            }

        }
        else if(v1 > v2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

/*
 * @brief Used for 2D dimensional datasets.
 * Start from x-axis to y-axis, record the range of utility and the corresponding top-k record.
 * Then find one of the top-k point by asking questions
 * @param original_set       The original dataset
 * @param u                  User's utility function
 * @param k                  The threshold
 */
void twoPI(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    timeval t1; gettimeofday(&t1, 0);

    point_set *rankingSet = new point_set();
    int *numOfEachCluster = new int[category.size()];
    for(int i = 0; i < category.size(); ++i)
        numOfEachCluster[i] = 0;
    point_t* u0 = new point_t(2);
    u0->attr[0] = 0; u0->attr[1] = 1;
    pset->findRanking(u0, rankingSet); //rankingSet is sorted by y-axis
    rankingSet->markFairTopk(k, category, bound, numOfEachCluster); //mark the fair top-k points
    std::vector<point_t *> Q = rankingSet->points;
    int M = Q.size();
    //rankingSet->print_with_score(u0);

    //initial the utility vectors for points order change
    std::vector<u_vector *> H;
    for (int i = 0; i < M - 1; i++)
    {
        double x0 = (Q[i]->attr[1] - Q[i + 1]->attr[1]) /
                    (Q[i + 1]->attr[0] - Q[i + 1]->attr[1] - Q[i]->attr[0] + Q[i]->attr[1]);
        //note that x0 in [0,1]
        if (0 <= x0 && x0 <= 1)
        {
            u_vector *ut = new u_vector(x0, Q[i], Q[i + 1]);
            ut->inserted(H);
        }
    }

    std::vector<point_set *> representSet;
    std::vector<u_vector *> border;
    representSet.push_back(new point_set());
    int j = 0, count = 0;
    while (count < k)
    {
        if (Q[j]->topk == 1)
        {
            representSet[representSet.size() - 1]->points.push_back(Q[j]);
            ++count;
        }
        ++j;
    }

    //representSet[0]->print();

    while (!H.empty()) //&& H[0]->x < 1)
    {
        //select one vector
        u_vector *current_u = H[H.size() - 1];
        H.pop_back();

        if(H.size() <= 0 || current_u->x != H[H.size() - 1]->x)
        {
            //change the location of point_up, point_down
            //One is in top-k, the other is not
            if (current_u->point_up->place <= k && current_u->point_down->place > k)
            {
                int upCluID = current_u->point_up->clusterID;
                int downCluID = current_u->point_down->clusterID;
                if (upCluID == downCluID ||
                    numOfEachCluster[upCluID] > bound[upCluID][0] && numOfEachCluster[downCluID] < bound[downCluID][1])
                {
                    current_u->point_down->place = current_u->point_up->place;
                    current_u->point_down->topk = 1;
                    numOfEachCluster[downCluID] += 1;
                    current_u->point_up->place = 2 * k;
                    current_u->point_up->topk = 0;
                    numOfEachCluster[upCluID] -= 1;

                    representSet.push_back(new point_set());
                    int j = 0, count = 0;
                    while (count < k)
                    {
                        if (Q[j]->topk == 1)
                        {
                            representSet[representSet.size() - 1]->points.push_back(Q[j]);
                            ++count;
                        }
                        ++j;
                    }
                    border.push_back(current_u);
                }
            }
            else if (current_u->point_up->place <= k && current_u->point_down->place <= k) //both are in fair top-k
            {
                int place_up = current_u->point_up->place;
                current_u->point_up->place = current_u->point_down->place;
                current_u->point_down->place = place_up;
            }

            //change the location of point_up, point_down
            point_t *p = current_u->point_down;
            std::vector<point_t *>::iterator iter = std::find(Q.begin(), Q.end(), p);
            int downPlace = std::distance(Q.begin(), iter);

            Q.erase(Q.begin() + downPlace);
            Q.insert(Q.begin() + downPlace - 1, p);
            //insert new utility vector of point: 1. place_down - 2 and 2. place_down - 1
            if (downPlace > 1)
            {
                double x0 = (Q[downPlace - 2]->attr[1] - Q[downPlace - 1]->attr[1]) /
                            (Q[downPlace - 1]->attr[0] - Q[downPlace - 1]->attr[1] -
                             Q[downPlace - 2]->attr[0] + Q[downPlace - 2]->attr[1]);

                if (current_u->x < x0 && x0 <= 1)
                {
                    u_vector *ut = new u_vector(x0, Q[downPlace - 2], Q[downPlace - 1]);
                    ut->inserted(H);
                }

            }
            //insert new utility vector of point: 1. place_down and 2. place_down + 1
            if (downPlace + 1 < Q.size())
            {
                double x0 = (Q[downPlace]->attr[1] - Q[downPlace + 1]->attr[1]) /
                            (Q[downPlace + 1]->attr[0] - Q[downPlace + 1]->attr[1] -
                             Q[downPlace]->attr[0] + Q[downPlace]->attr[1]);

                if (current_u->x <= x0 && x0 <= 1)
                {
                    u_vector *ut = new u_vector(x0, Q[downPlace], Q[downPlace + 1]);
                    ut->inserted(H);
                }

            }
        }
        else
        {
            H.pop_back();
            while (H.size() > 0 && current_u->x == H[H.size() - 1]->x)
                H.pop_back();

            u0->attr[0] = current_u->x + EQN2;
            u0->attr[1] = 1 - u0->attr[0];
            rankingSet->points.clear();
            pset->findRanking(u0, rankingSet); //rankingSet is sorted by y-axis
            Q = rankingSet->points;
            H.clear();
            for (int i = 0; i < M - 1; i++)
            {
                double x0 = (Q[i]->attr[1] - Q[i + 1]->attr[1]) /
                            (Q[i + 1]->attr[0] - Q[i + 1]->attr[1] - Q[i]->attr[0] + Q[i]->attr[1]);
                //note that x0 in [0,1]
                if (current_u->x + EQN2 < x0 && x0 <= 1)
                {
                    u_vector *ut = new u_vector(x0, Q[i], Q[i + 1]);
                    ut->inserted(H);
                }
            }

            point_set *resultSet = new point_set();
            rankingSet->fairTopk(k, category, bound, resultSet);
            rankingSet->markFairTopk(k, category, bound, numOfEachCluster); //mark the fair top-k points
            if (!resultSet->isSame(representSet[representSet.size() - 1]))
            {
                border.push_back(current_u);
                representSet.push_back(resultSet);
            }

        }

        /*
        for (int i = 0; i < Q.size(); i++)
        {
            std::cout<< Q[i]->id <<"  ";
            for (int j = 0; j < Q[i]->dim; j++)
                std::cout<< Q[i]->attr[j] << "  ";
            std::cout << Q[i]->cateAttr << "  ";
            std::cout << "score: " << Q[i]->dot_product(u) << "  top-k: " <<
                      Q[i]->topk << "  place: " << Q[i]->place << "\n";
        }
        std::cout << "\n";
        */
    }

    /*
    for(int i=0; i<border.size();i++)
    {
        printf("border %d u %lf\n", i, border[i]->x);
        cout << border[i]->point_up->place << "  " << border[i]->point_down->place << "\n";
    }
    for(int i=0; i < representSet.size();i++)
    {
        cout << "PresentSet" << i << "  " << endl;
        representSet[i]->print();
    }
    */

    //Interaction
    int numOfQuestion = 0;
    int b_left = 0, b_right = border.size() - 1;
    while (b_left <= b_right)
    {
        numOfQuestion++;
        int index = (b_left + b_right) / 2;
        double v1 = border[index]->point_up->dot_product(u);
        double v2 = border[index]->point_down->dot_product(u);
        if (v1 > v2)
        {
            b_right = index - 1;
        } else
        {
            b_left = index + 1;
        }
    }

    representSet[b_left]->printResult("TDIA", numOfQuestion, realResult, u, t1, 0, 0);


}





