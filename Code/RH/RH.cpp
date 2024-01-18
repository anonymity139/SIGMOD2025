#include "RH.h"
#include <sys/time.h>
#define Beta 0.1
/**
 * @brief Build all the partition(intersection of the halfspace), each partition corresponds to a top-1 point
 * @param p_set 				The set containing all the points which are possible to be the top-1 point
 *                              for some utility vectors
 * @param choose_item_points		The set containing points used to build choose_item
 * @param half_set_set			The returned partitions
 * @param considered_half_set	The set recorded all the partitions in the form of index
 */
void construct_partitions(Partition *R, std::vector<point_set*> &rankSet, std::vector<hyperplane*> &hyperSet,
                          std::vector<Partition*> &partitionSet, std::vector<int> &considered_halfset)
{
    int M = rankSet[0]->points.size(), dim = rankSet[0]->points[0]->dim;
    Partition *halfset;
    //record we use i-th point as the pivot(p[i]>p)
    for(int i = 0; i < rankSet.size(); ++i)
    {
        std::cout << i << "\n";
        halfset = new Partition(dim);
        for (int j = 0; j < M; j++)
        {
            for(int l = j + 1; l < M; ++l)
            {
                hyperplane *h = new hyperplane(rankSet[i]->points[l], rankSet[i]->points[j]);
                halfset->hyperplanes.push_back(h);
            }
        }
        for(int j = 0; j < R->hyperplanes.size(); ++j)
            halfset->hyperplanes.push_back(R->hyperplanes[j]);

        if (halfset->set_ext_pts())
        {
            halfset->topSet = rankSet[i];
            halfset->insertHyperplane(hyperSet);
            partitionSet.push_back(halfset);
            halfset->ID = partitionSet.size() - 1;
            considered_halfset.push_back(partitionSet.size() - 1);
        }
    }
}

/**
 * @brief Build the choose_item table used for selecting the hyperplane(question) asked user
 * @param half_set_set 		All the partitions
 * @param p_set 			All the points which used to build choose_item
 * @param choose_item_set 	The returned choose_item
 * @return The index of the choose_item(hyperplane) which divide the half_set most evenly
 */
int build_choose_item_table1(std::vector<Partition*> PartitionSet, std::vector<hyperplane*> hset,
                            std::vector<choose_item*> &choose_item_set)
{
    int HM = hset.size(); //The number of points used to construct hyperplane(questions)
    int PM = PartitionSet.size(); //The number of halfspace sets
    double ES_h = -1000; //Even score for finding hyperplane(question)
    int num_hyperplane = 0; //The index of the chosen hyperplane(question)
    for (int i = 0; i < HM; ++i)
    {
        //hset[i]->print();
        if(hset[i]->p_1 != NULL && hset[i]->p_2 != NULL)
        {
            choose_item *c_item = new choose_item();
            choose_item_set.push_back(c_item);
            for (int j = 0; j < PM; ++j)
            {
                c_item->h = hset[i];
                int which_side = PartitionSet[j]->check_relation(hset[i]);
                if (which_side == 1)
                {
                    c_item->positiveSide.push_back(PartitionSet[j]->ID);
                } else if (which_side == -1)
                {
                    c_item->negativeSide.push_back(PartitionSet[j]->ID);
                } else if (which_side == 0)
                {
                    c_item->intersectCase.push_back(PartitionSet[j]->ID);
                } else
                {
                    printf("Error: check side failed.\n");
                    return -1;
                }
            }

            //calculate the even score
            double ES_i;
            int p_s = c_item->positiveSide.size();
            int n_s = c_item->negativeSide.size();
            int i_s = c_item->intersectCase.size();
            if (p_s > n_s)
            {
                ES_i = n_s - (Beta * i_s);
            } else
            {
                ES_i = p_s - (Beta * i_s);
            }
            //renew even score if necessary
            if (ES_h < ES_i)
            {
                ES_h = ES_i;
                num_hyperplane = choose_item_set.size() - 1;
            }
            //std::cout << ES_i <<"\n";
        }
    }
    //std::cout << ES_h <<"\n";
    return num_hyperplane;
}





/**
 * @brief Based on the answer from the user, modify the choose_item table
 * @param choose_item_set 		The choose_item table which will be modified
 * @param half_set_set 			All the partitions
 * @param considered_half_set 	Based on the user's answer, the partitions still considered
 * @param H_num					The choose_item used for asking question
 * @param direction 			The user's answer(which side)
 *                              -true on the positive side
 *                              -false on the negative side
 * @return The choose_item(hyperplane) which divide the half_set most evenly
 */
int modify_choose_item_table1(std::vector<choose_item*> &choose_item_set, std::vector<Partition*> PartitionSet,
                             std::vector<int> &considered_halfset, int H_num, bool direction)
{
    int M = choose_item_set.size();
    //refine the considered half_set
    if (direction == true)//delete all the negative side partitions
    {
        int consider_count = 0; //index for scanning the considered_half_set
        for (int i = 0; i < choose_item_set[H_num]->negativeSide.size(); i++)
        {
            while (consider_count < considered_halfset.size() &&
                   considered_halfset[consider_count] <= choose_item_set[H_num]->negativeSide[i])
            {

                if (considered_halfset[consider_count] == choose_item_set[H_num]->negativeSide[i])
                    considered_halfset.erase(considered_halfset.begin() + consider_count);
                else
                    consider_count++;
            }
        }
    }
    else
    {
        int consider_count = 0;
        for (int i = 0; i < choose_item_set[H_num]->positiveSide.size(); i++)
        {
            while (consider_count < considered_halfset.size() &&
                   considered_halfset[consider_count] <= choose_item_set[H_num]->positiveSide[i])
            {
                if (considered_halfset[consider_count] == choose_item_set[H_num]->positiveSide[i])
                {
                    considered_halfset.erase(considered_halfset.begin() + consider_count);
                }
                else
                {
                    consider_count++;
                }
            }
        }
    }


    //build a halfspace based on choose_item[H_num]
    hyperplane *half;
    if (direction == true)
        half = new hyperplane(choose_item_set[H_num]->h->p_2, choose_item_set[H_num]->h->p_1);
    else
        half = new hyperplane(choose_item_set[H_num]->h->p_1, choose_item_set[H_num]->h->p_2);

    std::vector<int> intersectDelete;
    //refine all the half_set in the intersect.case
    for (int i = 0; i < choose_item_set[H_num]->intersectCase.size(); i++)
    {
        //num_set: index of the halfspace_set in half_set_set
        int num_set = choose_item_set[H_num]->intersectCase[i];
        PartitionSet[num_set]->hyperplanes.push_back(half);
        if(!PartitionSet[num_set]->set_ext_pts())
        {
            for(int j = 0; j < considered_halfset.size(); ++j)
            {
                if(considered_halfset[j] == num_set)
                    considered_halfset.erase(considered_halfset.begin() + j);
            }
            intersectDelete.push_back(num_set);
        }
    }


    //delete all the half_sets which are on the non-direction side
    if (direction)
    {
        //delete negative side
        for (int j = 0; j < M; j++)//scan all the choose_item, j-th choose_item
        {
            if (H_num != j)
            {
                //deal with one item, delete all the half_set on the negative side of the hyperplane
                int positive_count = 0, negative_count = 0, intersect_count = 0;
                for (int i = 0; i < choose_item_set[H_num]->negativeSide.size(); i++)
                {
                    //deal with vector positive_side
                    while (positive_count < choose_item_set[j]->positiveSide.size() &&
                           choose_item_set[j]->positiveSide[positive_count] <=
                           choose_item_set[H_num]->negativeSide[i])
                    {
                        if (choose_item_set[j]->positiveSide[positive_count] ==
                            choose_item_set[H_num]->negativeSide[i])
                        {
                            choose_item_set[j]->positiveSide.erase(
                                    choose_item_set[j]->positiveSide.begin() + positive_count);
                        }
                        else
                        {
                            positive_count++;
                        }
                    }

                    //deal with vector negative_side
                    while (negative_count < choose_item_set[j]->negativeSide.size() &&
                           choose_item_set[j]->negativeSide[negative_count] <=
                           choose_item_set[H_num]->negativeSide[i])
                    {
                        if (choose_item_set[j]->negativeSide[negative_count] ==
                            choose_item_set[H_num]->negativeSide[i])
                        {
                            choose_item_set[j]->negativeSide.erase(
                                    choose_item_set[j]->negativeSide.begin() + negative_count);
                        }
                        else
                        {
                            negative_count++;
                        }
                    }

                    //deal with vector intersect_side
                    while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                           choose_item_set[j]->intersectCase[intersect_count] <=
                           choose_item_set[H_num]->negativeSide[i])
                    {
                        if (choose_item_set[j]->intersectCase[intersect_count] ==
                            choose_item_set[H_num]->negativeSide[i])
                        {
                            choose_item_set[j]->intersectCase.erase(
                                    choose_item_set[j]->intersectCase.begin() + intersect_count);
                        }
                        else
                        {
                            intersect_count++;
                        }
                    }
                }
            }
        }
    }
    else
    {
        //delete positive side
        for (int j = 0; j < M; j++)//scan all the choose_item, j-th choose_item
        {
            if (H_num != j)
            {
                //deal with one item, delete all the half_set on the positive side of the hyperplane
                int positive_count = 0, negative_count = 0, intersect_count = 0;
                for (int i = 0; i < choose_item_set[H_num]->positiveSide.size(); i++)
                {
                    //deal with vector positive_side
                    while (positive_count < choose_item_set[j]->positiveSide.size() &&
                           choose_item_set[j]->positiveSide[positive_count] <=
                           choose_item_set[H_num]->positiveSide[i])
                    {
                        if (choose_item_set[j]->positiveSide[positive_count] ==
                            choose_item_set[H_num]->positiveSide[i])
                        {
                            choose_item_set[j]->positiveSide.erase(
                                    choose_item_set[j]->positiveSide.begin() + positive_count);
                        }
                        else
                        {
                            positive_count++;
                        }
                    }

                    //deal with vector negative_side
                    while (negative_count < choose_item_set[j]->negativeSide.size() &&
                           choose_item_set[j]->negativeSide[negative_count] <=
                           choose_item_set[H_num]->positiveSide[i])
                    {
                        if (choose_item_set[j]->negativeSide[negative_count] ==
                            choose_item_set[H_num]->positiveSide[i])
                        {
                            choose_item_set[j]->negativeSide.erase(
                                    choose_item_set[j]->negativeSide.begin() + negative_count);
                        }
                        else
                        {
                            negative_count++;
                        }
                    }

                    //deal with vector intersect_side
                    while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                           choose_item_set[j]->intersectCase[intersect_count] <=
                           choose_item_set[H_num]->positiveSide[i])
                    {
                        if (choose_item_set[j]->intersectCase[intersect_count] ==
                            choose_item_set[H_num]->positiveSide[i])
                        {
                            choose_item_set[j]->intersectCase.erase(
                                    choose_item_set[j]->intersectCase.begin() + intersect_count);
                        }
                        else
                        {
                            intersect_count++;
                        }
                    }
                }
            }
        }
    }

    for (int j = 0; j < M; j++)//scan all the choose_item, j-th choose_item
    {
        if (H_num != j)
        {
            //deal with one item, delete all the half_set on the negative side of the hyperplane
            int positive_count = 0, negative_count = 0, intersect_count = 0;
            for (int i = 0; i < intersectDelete.size(); i++)
            {
                //deal with vector positive_side
                while (positive_count < choose_item_set[j]->positiveSide.size() &&
                       choose_item_set[j]->positiveSide[positive_count] <= intersectDelete[i])
                {
                    if (choose_item_set[j]->positiveSide[positive_count] == intersectDelete[i])
                    {
                        choose_item_set[j]->positiveSide.erase(
                                choose_item_set[j]->positiveSide.begin() + positive_count);
                    }
                    else
                    {
                        positive_count++;
                    }
                }

                //deal with vector negative_side
                while (negative_count < choose_item_set[j]->negativeSide.size() &&
                       choose_item_set[j]->negativeSide[negative_count] <= intersectDelete[i])
                {
                    if (choose_item_set[j]->negativeSide[negative_count] == intersectDelete[i])
                    {
                        choose_item_set[j]->negativeSide.erase(
                                choose_item_set[j]->negativeSide.begin() + negative_count);
                    }
                    else
                    {
                        negative_count++;
                    }
                }

                //deal with vector intersect_side
                while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                       choose_item_set[j]->intersectCase[intersect_count] <= intersectDelete[i])
                {
                    if (choose_item_set[j]->intersectCase[intersect_count] == intersectDelete[i])
                    {
                        choose_item_set[j]->intersectCase.erase(
                                choose_item_set[j]->intersectCase.begin() + intersect_count);
                    }
                    else
                    {
                        intersect_count++;
                    }
                }
            }
        }
    }

    //deal with the refined half_set
    for (int i = 0; i < choose_item_set[H_num]->intersectCase.size(); i++)
    {
        int num_set = choose_item_set[H_num]->intersectCase[i];
        //scan all the choose_item, j-th choose_item
        for (int j = 0; j < M; j++)
        {
            if (H_num != j)
            {
                //deal with vector intersect_case
                int intersect_count = 0;
                while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                       choose_item_set[j]->intersectCase[intersect_count] <= num_set)
                {
                    if (choose_item_set[j]->intersectCase[intersect_count] == num_set)
                    {
                        int which_side = PartitionSet[num_set]->check_relation(choose_item_set[j]->h);
                        if (which_side == -2)
                        {
                            printf("%s\n", "Error: check side failed.");
                            return -1;
                        }
                        if (which_side != 0)
                        {
                            choose_item_set[j]->intersectCase.erase(
                                    choose_item_set[j]->intersectCase.begin() + intersect_count);
                            if (which_side == 1)//the half_set is on the positive side
                            {
                                bool is_insert = false;
                                for (int q = 0; q < choose_item_set[j]->positiveSide.size(); q++)
                                {
                                    if (num_set < choose_item_set[j]->positiveSide[q])
                                    {
                                        choose_item_set[j]->positiveSide.insert(
                                                choose_item_set[j]->positiveSide.begin() + q, num_set);
                                        is_insert = true;
                                        break;
                                    }
                                }
                                if (is_insert == false)
                                {
                                    choose_item_set[j]->positiveSide.push_back(num_set);
                                }

                            }
                            else//the half_set is on the negative side
                            {
                                bool is_insert = false;
                                for (int q = 0; q < choose_item_set[j]->negativeSide.size(); q++)
                                {
                                    if (num_set < choose_item_set[j]->negativeSide[q])
                                    {
                                        choose_item_set[j]->negativeSide.insert(
                                                choose_item_set[j]->negativeSide.begin() + q, num_set);
                                        is_insert = true;
                                        break;
                                    }
                                }
                                if (is_insert == false)
                                {
                                    choose_item_set[j]->negativeSide.push_back(num_set);
                                }
                            }

                        }
                        break;
                    }
                    intersect_count++;
                }
            }
        }
    }

    //refine the choose_item_set and select the hyperplane asked user
    choose_item_set.erase(choose_item_set.begin() + H_num);

    int item_count = 0;                         //the index for scanning the choose_item_set
    int ES_h = -1000;                           //Even score
    int num_hyperplane = 0;                     //index of the chosen hyperplane(question)
    int choose_size = choose_item_set.size();   //the original size of the choose_item_set
    for (int i = 0; i < choose_size; i++)
    {
        choose_item *c_item = choose_item_set[item_count];
        int M_positive = choose_item_set[item_count]->positiveSide.size();
        int M_negative = choose_item_set[item_count]->negativeSide.size();
        int M_intersect = choose_item_set[item_count]->intersectCase.size();
        if (M_negative + M_intersect == 0 || M_positive + M_intersect == 0)// || M_negative + M_positive == 0)
        {
            choose_item_set.erase(choose_item_set.begin() + item_count);
        }
        else
        {
            if (M_positive < M_negative)
            {
                double ES_i = M_positive - (Beta * M_intersect);
                if (ES_h < ES_i)
                {
                    ES_h = ES_i;
                    num_hyperplane = item_count;
                }
            }
            else
            {
                double ES_i = M_negative - (Beta * M_intersect);
                if (ES_h < ES_i)
                {
                    ES_h = ES_i;
                    num_hyperplane = item_count;
                }
            }
            item_count++;
        }
    }
    return num_hyperplane;
}


/**
 * @brief Ask user question and find one of the top-k point
 *        The input set is refined by skyband
 *        For find the point which is the top-k point w.r.t any u in R, we search the whole set
 * @param original_set 		The refined dataset
 * @param u 				The linear function
 * @param k 				The threshold top-k
 */
int RH(point_set* pset, point_t* u, int k, point_set *realResult, std::vector<std::string> &category, std::vector<std::vector<double>> &bound, bool error)
{
    timeval t1;
    gettimeofday(&t1, 0);

    pset->random(0.5);
    int dim = pset->points[0]->dim, M = pset->points.size(), numOfQuestion = 0;

    Partition *R = new Partition(dim);
    std::vector<point_t*> current, considered;
    considered.push_back(pset->points[0]);

    //compare: p_set contains all the points
    for(int i = 1; i < M; ++i)
    {
        for(int j = 0; j < considered.size(); j++)
            current.push_back(considered[j]);
        considered.push_back(pset->points[i]);

        while(current.size() > 0)
        {
            bool needAsk = false;
            double distance = 999999;
            int p_index = -1;// the index where the point we select to ask question
            //find the question to be asked the user
            for(int j = 0; j < current.size(); j++)
            {
                hyperplane *h = new hyperplane(pset->points[i], current[j]);
                int relation = R->check_relationlose(h); //check_situation_accelerate(h, R_half_set, 0);
                //if intersect, calculate the distance
                if(relation == 0)
                {
                    needAsk = true;
                    double d_h = h->check_distance(R->center);
                    if(d_h < distance)
                    {
                        distance = d_h;
                        p_index = j;
                    }
                }
                else
                {
                    current.erase(current.begin() + j);
                    --j;
                }
            }
            if(needAsk)
            {
                numOfQuestion++;
                double v1 = pset->points[i]->dot_product(u);
                double v2 = current[p_index]->dot_product(u);
                double prob = 2;
                if(error)
                {
                    if (v1 > v2)
                        prob = v1 - v2;
                    else
                        prob = v2 - v1;
                    prob = 1.0 / (1.0 + exp(-1.0 * prob));
                }
                double sumr = (double)rand() / RAND_MAX;
                if(v1 > v2 && sumr < prob || v1 < v2 && sumr > prob)
                {
                    hyperplane *half = new hyperplane(current[p_index], pset->points[i]);
                    R->hyperplanes.push_back(half);
                }
                else
                {
                    hyperplane *half = new hyperplane(pset->points[i], current[p_index]);
                    R->hyperplanes.push_back(half);
                }
                current.erase(current.begin() + p_index);
                R->set_ext_pts();
            }
        }
    }


    point_set *rankingSet = new point_set, *resultSet = new point_set();
    pset->findRanking(R->average_point(), rankingSet);
    rankingSet->fairTopk(k, category, bound,resultSet);
    resultSet->printResult("RH", numOfQuestion, realResult, u, t1, 0, 0);

    return numOfQuestion;
}


