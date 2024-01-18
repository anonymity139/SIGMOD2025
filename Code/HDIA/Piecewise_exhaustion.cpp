#include "Piecewise_exhaustion.h"
#include "Partition.h"
#include "halfspace_tree.h"
#include "interaction_tree.h"
#define eps 0.2

/*

void find_shortest_tree(hyperplane_set* crt, hyperplane_set *rem, std::vector<point_t*> partitionPt, int &curMinHeight)
{
    int M = rem->hyperplanes.size();
    if(M == 0)
        return;
    else
    {
        for (int i = 0; i < M; ++i)
        {
            hyperplane *h = rem->hyperplanes[i];
            crt->hyperplanes.push_back(h);
            rem->hyperplanes.erase(rem->hyperplanes.begin() + i);

            find_shortest_tree(crt, rem, partitionPt, curMinHeight);

            crt->hyperplanes.pop_back();
            rem->hyperplanes.insert(rem->hyperplanes.begin() + i, h);
        }
    }
}
*/

void Piecewise_exhaustion(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int dim = pset->points[0]->dim;
    int piece_size = 40;
    Partition *R = new Partition(dim);
    point_set *resultSet = new point_set();
    int numOfQuestion = 0;

    while(1)
    {
        //randomly build some hyper-planes based on the points without overlapping
        hyperplane_set *hset = new hyperplane_set();
        int count = 0;
        while (hset->hyperplanes.size() < piece_size && count < pset->points.size() * 10)
        {
            count++;
            int pindex1 = rand() % pset->points.size();
            int pindex2 = rand() % pset->points.size();
            while (pindex1 == pindex2)
                pindex2 = rand() % pset->points.size();
            hyperplane *h = new hyperplane(pset->points[pindex1], pset->points[pindex2]);

            if (R->check_relation(h) == 0 && !hset->is_hyperplane_exist(h))
                hset->hyperplanes.push_back(h);
        }

        //obtain the partitions divided by the hyper-planes
        halfspace_tree *hTree = new halfspace_tree(R);
        for(int i = 0; i < hset->hyperplanes.size(); ++i)
        {
            std::cout << i << "\n";
            hTree->insert(hset->hyperplanes[i]);
            //hTree->print();
        }
        //hTree->print_leaves();

        //find the shortest tree based on the partitions
        std::vector<point_t*> partitionPt;
        int currentMinRound = hTree->find_leafPt_height(partitionPt);

        if(partitionPt.size() <= 1)
        {
            point_set* rankingSet1 = new point_set(), *resultSet1 = new point_set();
            pset->findRanking(R->ext_pts[0], rankingSet1);
            rankingSet1->fairTopk(k, category, bound,resultSet1);
            resultSet1->printResult("HDIA", numOfQuestion, realResult, u, t1, 0, 0);
            return;
        }


        while(partitionPt.size() > 1)
        {
            int partitionLeft = partitionPt.size() / 10;
            if(partitionLeft < 1)
                partitionLeft = 1;
            it_tree *itTree = new it_tree(hset, partitionPt, partitionLeft);
            //itTree->print();

            //interaction based on the tree
            it_node *curNode = itTree->root;
            while (curNode->bestpt != NULL && curNode->bestng != NULL) //if it is an internal node
            {
                numOfQuestion++;
                double v1 = curNode->bestdivideHyper->p_1->dot_product(u);
                double v2 = curNode->bestdivideHyper->p_2->dot_product(u);
                if (v1 > v2)
                {
                    hyperplane *h = new hyperplane(curNode->bestdivideHyper->p_2, curNode->bestdivideHyper->p_1);
                    R->hyperplanes.push_back(h);
                    curNode = curNode->bestpt;
                }
                else
                {
                    hyperplane *h = new hyperplane(curNode->bestdivideHyper->p_1, curNode->bestdivideHyper->p_2);
                    R->hyperplanes.push_back(h);
                    curNode = curNode->bestng;
                }
                R->set_ext_pts();
                R->center = R->average_point();

                //stopping condition
                point_set *rankingSet1 = new point_set(), *resultSet1 = new point_set();
                pset->findRanking(R->ext_pts[0], rankingSet1);
                rankingSet1->fairTopk(k, category, bound, resultSet1);
                bool flag = true;
                for (int i = 1; i < R->ext_pts.size(); ++i)
                {
                    point_set *rankingSet2 = new point_set(), *resultSet2 = new point_set();
                    pset->findRanking(R->ext_pts[i], rankingSet2);
                    rankingSet2->fairTopk(k, category, bound, resultSet2);
                    if (!resultSet1->isSame(resultSet2))
                    {
                        flag = false;
                        break;
                    }
                }
                if (flag)
                {
                    resultSet = resultSet1;
                    resultSet->printResult("HDIA", numOfQuestion, realResult, u, t1, 0, 0);
                    return;
                }
            }

            partitionPt = curNode->ptSet;
        }
    }

}


void Piecewise_exhaustion_robust(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int dim = pset->points[0]->dim;
    int piece_size = 7;
    Partition *R = new Partition(dim);
    point_set *resultSet = new point_set();
    int numOfQuestion = 0;
    int count = 0;

    while(1)
    {
        //randomly build some hyper-planes based on the points without overlapping
        hyperplane_set *hset = new hyperplane_set();
        while (hset->hyperplanes.size() < piece_size)
        {
            int pindex1 = rand() % pset->points.size();
            int pindex2 = rand() % pset->points.size();
            while (pindex1 == pindex2)
                pindex2 = rand() % pset->points.size();
            hyperplane *h = new hyperplane(pset->points[pindex1], pset->points[pindex2]);

            if (R->check_relation(h) == 0 && !hset->is_hyperplane_exist(h))
                hset->hyperplanes.push_back(h);
        }

        //obtain the partitions divided by the hyper-planes
        halfspace_tree *hTree = new halfspace_tree(R);
        for(int i = 0; i < hset->hyperplanes.size(); ++i)
        {
            //std::cout << i << "\n";
            hTree->insert(hset->hyperplanes[i]);
            //hTree->print();
        }
        R->set_ext_pts();
        //hTree->print_leaves();

        //find the shortest tree based on the partitions
        std::vector<point_t*> partitionPt;
        int currentMinRound = hTree->find_leafPt_height(partitionPt);

        /*
        if(partitionPt.size() <= 1)
        {
            point_set* rankingSet1 = new point_set(), *resultSet1 = new point_set();
            pset->findRanking(R->ext_pts[0], rankingSet1);
            rankingSet1->fairTopk(k, category, bound,resultSet1);
            bool flag = true;
            resultSet1->printResult("R-HDIA", numOfQuestion, t1, 0, 0);
            return;
        }
        */

        it_tree *itTree; //= new it_tree(hset, partitionPt);
        //itTree->print();

        //interaction based on the tree
        it_node *curNode = itTree->root;
        while(curNode->bestpt != NULL && curNode->bestng != NULL) //if it is an internal node
        {
            hyperplane *htest = new hyperplane(curNode->bestdivideHyper->p_2, curNode->bestdivideHyper->p_1);
            if(R->check_relation(htest) != 0)
                break;
            double v1 = curNode->bestdivideHyper->p_1->dot_product(u);
            double v2 = curNode->bestdivideHyper->p_2->dot_product(u);
            count++;

            //probability of succeed
            double simProb1 = R->largest_utility(new point_t(curNode->bestdivideHyper->p_1, curNode->bestdivideHyper->p_2));
            double simProb2 = R->largest_utility(new point_t(curNode->bestdivideHyper->p_2, curNode->bestdivideHyper->p_1));
            double simProb;
            if(simProb1 > simProb2)
                simProb = simProb1;
            else
                simProb = simProb2;
            simProb = 1.0 / (1.0 + exp(-1.0 * simProb));
            simProb = (1-2*simProb)*(1-2*simProb);
            simProb = (-1 / simProb) * log(eps/2);
            if(simProb > 200)
                simProb = 200;
            double prob = 0;
            if(v1 > v2)
                prob = v1 - v2;
            else
                prob = v2 - v1;
            prob = 1.0 / (1.0 + exp(-1.0 * prob));

            double sumr = 0;
            for(int t = 0; t < simProb; ++t)
            {
                //a random number from 0 to 1
                double r = (double)rand() / RAND_MAX;
                if(r < prob)
                    sumr++;
                numOfQuestion++;
            }
            sumr = sumr / simProb;
            if (v1 > v2 && sumr > 0.5 || v1 < v2 && sumr < 0.5)
            {
                hyperplane *h = new hyperplane(curNode->bestdivideHyper->p_2, curNode->bestdivideHyper->p_1);
                R->hyperplanes.push_back(h);
                curNode = curNode->bestpt;
            }
            else
            {
                hyperplane *h = new hyperplane(curNode->bestdivideHyper->p_1, curNode->bestdivideHyper->p_2);
                R->hyperplanes.push_back(h);
                curNode = curNode->bestng;
            }
            R->set_ext_pts();
            R->center = R->average_point();

            //stopping condition
            point_set* rankingSet1 = new point_set(), *resultSet1 = new point_set();
            pset->findRanking(R->ext_pts[0], rankingSet1);
            rankingSet1->fairTopk(k, category, bound,resultSet1);
            bool flag = true;
            for(int i = 1; i < R->ext_pts.size(); ++i)
            {
                point_set* rankingSet2 = new point_set(), *resultSet2 = new point_set();
                pset->findRanking(R->ext_pts[i], rankingSet2);
                rankingSet2->fairTopk(k, category, bound,resultSet2);
                if(!resultSet1->isSame(resultSet2))
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
            {
                resultSet = resultSet1;
                resultSet->printResult("R-HDIA", count, realResult, u, t1, 0, 0);
                return;
            }
        }
    }

}