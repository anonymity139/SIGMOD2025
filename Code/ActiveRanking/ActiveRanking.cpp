#include "ActiveRanking.h"
#define Repetition 15

/**
 * @brief Ask user questions and give a ranking
 * @param original_set 		The original dataset
 * @param u 				The linear function
 * @param k 				The threshold top-k
 */
int ActiveRanking(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int dim = pset->points[0]->dim, numOfQuestion = 0, M = pset->points.size();;
    pset->random(0.5);

    //initialization
    Partition *R = new Partition(dim);
    point_set *current = new point_set();
    current->points.push_back(pset->points[0]);


    //store all the points in order
    for (int i = 1; i < M; i++) //compare: p_set contains all the points
    {
        int num_point = current->points.size();
        int place = 0; //the place of the point inserted into the current_use
        //find the question asked user
        for (int j = 0; j < num_point; j++)
        {
            hyperplane *h = new hyperplane(pset->points[i], current->points[j]);
            int relation = R->check_relationlose(h);
            delete h;
            //if intersect, calculate the distance
            if (relation == 0)
            {
                numOfQuestion++;
                double v1 = pset->points[i]->dot_product(u);
                double v2 = current->points[j]->dot_product(u);
                if (v1 > v2)
                {
                    hyperplane *h = new hyperplane(current->points[j], pset->points[i]);
                    R->hyperplanes.push_back(h);
                    if(!R->set_ext_pts())
                        R->hyperplanes.pop_back();
                    break;

                }
                else
                {
                    hyperplane *h = new hyperplane(pset->points[i], current->points[j]);
                    R->hyperplanes.push_back(h);
                    if(!R->set_ext_pts())
                        R->hyperplanes.pop_back();
                    place = j + 1;
                }
                //R->print();
            }
            else if (relation == -1)
            {
                place = j + 1;
            }
            else
            {
                break;
            }
        }
        current->points.insert(current->points.begin() + place, pset->points[i]);
    }

    point_set *resultSet = new point_set();
    current->fairTopk(k, category, bound,resultSet);
    resultSet->printResult("ActiveRanking", numOfQuestion, realResult, u, t1, 0, 0);



    return numOfQuestion;

}





/**
 * @brief Ask user questions and give a ranking
 * @param original_set 		The original dataset
 * @param u 				The linear function
 * @param k 				The threshold top-k
 */
int ActiveRanking_robust(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int dim = pset->points[0]->dim, numOfQuestion = 0, M = pset->points.size();;
    pset->random(0.5);

    //initialization
    Partition *R = new Partition(dim);
    point_set *current = new point_set();
    current->points.push_back(pset->points[0]);


    //store all the points in order
    for (int i = 1; i < M; i++) //compare: p_set contains all the points
    {
        int num_point = current->points.size();
        int place = 0; //the place of the point inserted into the current_use
        //find the question asked user
        for (int j = 0; j < num_point; j++)
        {
            hyperplane *h = new hyperplane(pset->points[i], current->points[j]);
            int relation = R->check_relationlose(h);
            delete h;
            //if it intersects, calculate the distance
            if (relation == 0)
            {
                numOfQuestion++;
                double v1 = pset->points[i]->dot_product(u);
                double v2 = current->points[j]->dot_product(u);
                double prob = 0;
                if(v1 > v2)
                    prob = v1 - v2;
                else
                    prob = v2 - v1;
                prob = 1.0 / (1.0 + exp(-1.0 * prob));

                double sumr = 0;
                for(int t = 0; t < Repetition; ++t)
                {
                    //a random number from 0 to 1
                    double r = (double)rand() / RAND_MAX;
                    if(r < prob)
                        sumr++;
                }
                sumr = sumr / Repetition;

                if (v1 > v2 && sumr > 0.5 || v1 < v2 && sumr < 0.5)
                {
                    hyperplane *h = new hyperplane(current->points[j], pset->points[i]);
                    R->hyperplanes.push_back(h);
                    if(!R->set_ext_pts())
                        R->hyperplanes.pop_back();
                    break;

                }
                else
                {
                    hyperplane *h = new hyperplane(pset->points[i], current->points[j]);
                    R->hyperplanes.push_back(h);
                    if(!R->set_ext_pts())
                        R->hyperplanes.pop_back();
                    place = j + 1;
                }
                //R->print();
            }
            else if (relation == -1)
            {
                place = j + 1;
            }
            else
            {
                break;
            }
        }
        current->points.insert(current->points.begin() + place, pset->points[i]);
    }

    point_set *resultSet = new point_set();
    current->fairTopk(k, category, bound,resultSet);
    resultSet->printResult("ActiveRanking", numOfQuestion, realResult, u, t1, 0, 0);



    return numOfQuestion;

}










