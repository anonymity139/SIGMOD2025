#include "UH.h"
#include <sys/time.h>



/**
 * @brief Get the index of the "current best" point.
 *        Find the average utility vector u of the extreme points of range R. Find the point with the maximum utility w.r.t u
 * @param P         the input car set
 * @param C_idx     the indexes of the current candidate favorite car in P
 * @param R         the utility range
 * @return          the index of the point in C_idx
 */
int get_current_best_pt(point_set *P, vector<int> &C_idx, Partition *R)
{
    int dim = P->points[0]->dim;

    // use the "mean" utility vector in R (other strategies could also be used)
    point_t *mean = R->average_point();

    // look for the maximum utility point w.r.t. the "mean" utility vector
    int best_pt_idx;
    double max = 0;
    for (int i = 0; i < C_idx.size(); i++)
    {
        point_t *pt = P->points[C_idx[i]];

        double v = mean->dot_product(pt);
        if (v > max)
        {
            max = v;
            best_pt_idx = C_idx[i];
        }
    }
    return best_pt_idx;
}


/**
 * @brief Generate s cars for selection in a round
 * @param P                     the input car set
 * @param C_idx                 the indexes of the current candidate favorite car in P
 * @param s                     the number of cars for user selection
 * @param current_best_idx      the current best car
 * @param last_best             the best car in previous interaction
 * @param frame                 the frame for obtaining the set of neigbouring vertices of the current best vertex
 *                              (used only if cmp_option = SIMPLEX)
 * @param cmp_option            the car selection mode, which must be either SIMPLEX or RANDOM
 * @return The set of indexes of the points chosen for asking questions
 */
vector<int> generate_S(point_set *P, vector<int> &C_idx, int s, int current_best_idx, int &last_best,
                       vector<int> &frame, int cmp_option)
{
    // the set of s cars for selection
    vector<int> S;
    if (cmp_option == RANDOM) // RANDOM car selection mode
    {
        // randoming select at most s non-overlaping cars in the candidate set
        while (S.size() < s && S.size() < C_idx.size())
        {
            int idx = rand() % C_idx.size();

            bool isNew = true;
            for (int i = 0; i < S.size(); i++)
            {
                if (S[i] == idx)
                {
                    isNew = false;
                    break;
                }
            }
            if (isNew)
            {
                S.push_back(idx);
            }
        }
    }
    else if (cmp_option == SIMPLEX) // SIMPLEX car selection mode
    {
        if (last_best != current_best_idx || frame.size() == 0) // the new frame is not computed before (avoid duplicate frame computation)
        {
            // create one ray for each car in P for computing the frame
            vector<point_t *> rays;
            int best_i = -1;
            for (int i = 0; i < P->points.size(); i++)
            {
                if (i == current_best_idx)
                {
                    best_i = i;
                    continue;
                }

                point_t *best = P->points[current_best_idx];
                point_t *newRay = P->points[i]->sub(best);
                rays.push_back(newRay);
            }

            // frame compuatation
            frameConeFastLP(rays, frame);
            // update the indexes lying after current_best_idx
            for (int i = 0; i < frame.size(); i++)
            {
                if (frame[i] >= current_best_idx)
                {
                    frame[i]++;
                }
                //S[i] = C_idx[S[i]];
            }
        }
        for (int j = 0; j < C_idx.size(); j++)
        {
            if (C_idx[j] ==
                current_best_idx)// it is possible that current_best_idx is no longer in the candidate set, then no need to compare again
            {
                S.push_back(j);
                break;
            }
        }
        // select at most s non-overlaping cars in the candidate set based on "neighboring vertices" obtained via frame compuation
        for (int i = 0; i < frame.size() && S.size() < s; i++)
        {
            for (int j = 0; j < C_idx.size() && S.size() < s; j++)
            {
                if (C_idx[j] == current_best_idx)
                {
                    continue;
                }

                if (C_idx[j] == frame[i])
                {
                    S.push_back(j);
                }
            }
        }

        // if less than s car are selected, fill in the remaining one
        if (S.size() < s && C_idx.size() > s)
        {
            for (int j = 0; j < C_idx.size(); j++)
            {
                bool noIn = true;
                for (int i = 0; i < S.size(); i++)
                {
                    if (j == S[i])
                    {
                        noIn = false;
                    }
                }
                if (noIn)
                {
                    S.push_back(j);
                }

                if (S.size() == s)
                {
                    break;
                }
            }
        }
    }
    else // for testing only. Do not use this!
    {
        vector<point_t *> rays;

        int best_i = -1;
        for (int i = 0; i < C_idx.size(); i++)
        {
            if (C_idx[i] == current_best_idx)
            {
                best_i = i;
                continue;
            }

            point_t *best = P->points[current_best_idx];

            point_t *newRay = P->points[C_idx[i]]->sub(best);

            rays.push_back(newRay);
        }

        partialConeFastLP(rays, S, s - 1);
        if (S.size() > s - 1)
        {
            S.resize(s - 1);
        }
        for (int i = 0; i < S.size(); i++)
        {
            if (S[i] >= best_i)
            {
                S[i]++;
            }

            //S[i] = C_idx[S[i]];
        }
        S.push_back(best_i);
    }
    return S;
}

/**
 * @brief Generate the options for user selection. Update the extreme vectors based on the user feedback
 *        Prune points which are impossible to be top-1
 * @param P                     The skyline dataset
 * @param C_idx                 the indexes of the current candidate favorite car in P
 * @param u                     the utility vector
 * @param s                     the number of cars for user selection
 * @param ext_vec               the set of extreme vecotr
 * @param current_best_idx      the current best car
 * @param last_best             the best car in previous interaction
 * @param frame                 the frame for obtaining the set of neigbouring vertices of the current best vertex
 *                              (used only if cmp_option = SIMPLEX)
 * @param cmp_option            the car selection mode, which must be either SIMPLEX or RANDOM
 */
void update_ext_vec(point_set *P, vector<int> &C_idx, point_t *u, int s, vector<point_t *> &ext_vec, int &current_best_idx,
               int &last_best, vector<int> &frame, int cmp_option, Partition *R, bool error)
{
    // generate s cars for selection in a round
    vector<int> S = generate_S(P, C_idx, s, current_best_idx, last_best, frame, cmp_option);


    int max_i = -1;
    double max = -1;

    /*
    for (int i = 0; i < S.size(); i++)
    {
        point_t *p = P->points[C_idx[S[i]]];
        double v = p->dot_product(u);
        if (v > max)
        {
            max = v;
            max_i = i;
        }
    }
    */

    double v1 = P->points[C_idx[S[0]]]->dot_product(u);
    double v2 = P->points[C_idx[S[1]]]->dot_product(u);
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
    if (v1 >= v2 && sumr < prob || v1 < v2 && sumr > prob)
    {
        max_i = 0;
    }
    else
    {
        max_i = 1;
    }



    // get the better car among those from the user
    last_best = current_best_idx;
    current_best_idx = C_idx[S[max_i]];

    // for each non-favorite car, create a new extreme vector
    for (int i = 0; i < S.size(); i++)
    {
        if (max_i == i)
        {
            continue;
        }
        point_t *tmp = P->points[C_idx[S[i]]]->sub(P->points[C_idx[S[max_i]]]);
        point_t *new_ext_vec = tmp->scale(1 / tmp->cal_len());
        ext_vec.push_back(new_ext_vec);
        R->hyperplanes.push_back(new hyperplane(P->points[C_idx[S[i]]], P->points[C_idx[S[max_i]]]));
        C_idx[S[i]] = -1;
    }
    R->set_ext_pts();

    // directly remove the non-favorite car from the candidate set
    vector<int> newC_idx;
    for (int i = 0; i < C_idx.size(); i++)
    {
        if (C_idx[i] >= 0)
        {
            newC_idx.push_back(C_idx[i]);
        }
    }
    C_idx = newC_idx;
}


/**
 * @brief The interactive algorithm UH-Random-Adapt.
 *        Find a points/return a point which satisfy the regret ratio
 *        Modified stopping condition, prune points which are not able to be top-k
 * @param original_P     the original dataset
 * @param u              the unkonwn utility vector
 * @param s              the question size
 * @param k              The threshold top-k
 */
int UHRandom(point_set *pset, point_t *u, int k, std::vector<std::string> &category, std::vector<std::vector<double>> &bound)
{
    //initialize time
    timeval t1;
    gettimeofday(&t1, 0);

    // M: the number of cars, dim: the dimension of the utility vector
    int M = pset->points.size(), dim = pset->points[0]->dim, numOfQuestion = 0;
    Partition *R = new Partition(dim);
    for(int i = 0; i < M; ++i)
        pset->points[i]->place = i;
    std::vector<std::vector<int>> is_considered(M, std::vector<int>(M, 0));
    point_t *pidx1, *pidx2;

    // interactively reduce the candidate set and shrink the candidate utility range
    while (pset->points.size() > k)
    {
        M = pset->points.size();
        bool choose = false;
        int considered_count = 0;

        //point selection
        while (!choose)
        {
            pidx1 = pset->points[rand() % pset->points.size()];
            pidx2 = pset->points[rand() % pset->points.size()];

            if (is_considered[pidx1->place][pidx2->place] == 0)
            {
                is_considered[pidx1->place][pidx2->place] = numOfQuestion;
                is_considered[pidx2->place][pidx1->place] = numOfQuestion;
                considered_count++;
                if(!pidx1->is_same(pidx2))
                {
                    if (R->check_relation(new hyperplane(pidx1, pidx2)) == 0)
                        choose = true;
                }
            }
            else if(is_considered[pidx1->place][pidx2->place] < numOfQuestion)
            {
                is_considered[pidx1->place][pidx2->place] = numOfQuestion;
                is_considered[pidx2->place][pidx1->place] = numOfQuestion;
                considered_count++;
            }
            if (considered_count >= (M * (M + 1)) / 2 && (!choose))
                goto stop;
        }

        //interaction
        double v1 = pidx1->dot_product(u);
        double v2 = pidx2->dot_product(u);
        hyperplane *h;
        if (v1 > v2)
            h = new hyperplane(pidx2, pidx1);
        else
            h = new hyperplane(pidx1, pidx2);
        R->hyperplanes.push_back(h);
        R->set_ext_pts();
        numOfQuestion++;
        R->prune_points(pset, k, category, bound);
    }

    stop:
    if(pset->points.size() > k)
    {
        point_set *rankingSet = new point_set, *resultSet = new point_set();
        pset->findRanking(R->average_point(), rankingSet);
        rankingSet->fairTopk(k, category, bound,resultSet);
        pset = resultSet;
    }
    //pset->printResult("UH-Random", numOfQuestion, t1, 0, 0);

    return numOfQuestion;
}






/**
 * @brief The interactive algorithm Random/Simplex. Find best points/return a point which satisfy the regret ratio
 * @param original_P     the original dataset
 * @param u              the unkonwn utility vector
 * @param s              the question size
 * @param epsilon        the required regret ratio
 * @param maxRound       the maximum number of rounds of interactions
 * @param cmp_option     the car selection mode
 *                       -SIMPLEX UH-Simplex method
 *                       -RANDOM UH-Simplex method
 * @param stop_option    the stopping condition
 *                       -NO_BOUND
 *                       -EXACT_BOUND
 *                       -APRROX_BOUND
 * @param prune_option   the skyline algorithm
 *                       -SQL
 *                       -RTREE
 * @param dom_option     the domination checking mode
 *                       -HYPER_PLANE
 *                       -CONICAL_HULL
 */
int max_utility(point_set *pSet, point_t *u, int k, point_set *realResult, std::vector<std::string> &category,
                    std::vector<std::vector<double>> &bound, int s, int cmp_option, bool error)
{
    timeval t1, t2;
    gettimeofday(&t1, 0);
    int prune_option = RTREE;
    int dom_option = HYPER_PLANE;
    int stop_option = EXACT_BOUND;


    int dim = pSet->points[0]->dim;
    int numOfQuestion = 0;
    double rr = 1;
    Partition *R = new Partition(dim);
    vector<point_t *> ext_vec;
    for (int i = 0; i < dim; i++)
    {
        point_t *e = new point_t(dim);
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
                e->attr[j] = -1;
            else
                e->attr[j] = 0;
        }
        ext_vec.push_back(e);
    }
    point_set *rankingSet = new point_set(), *resultSet = new point_set();
    rankingSet->fairTopk(k, category, bound,resultSet);
    vector<int> frame;

    while (resultSet->points.size() < k)
    {
        resultSet->points.clear();
        vector<int> C_idx;
        for (int i = 0; i < pSet->points.size(); i++)
            C_idx.push_back(i);
        rtree_pruning(pSet, C_idx, R);

        int current_best_idx = -1;
        int last_best = -1;
        current_best_idx = get_current_best_pt(pSet, C_idx, R);

        // interactively reduce the candidate set and shrink the candidate utility range
        while (C_idx.size() > 1) // while none of the stopping condition is true
        {
            ++numOfQuestion;

            sort(C_idx.begin(), C_idx.end()); // prevent select two different points after different skyline algorithms

            // generate the options for user selection and update the extreme vectors based on the user feedback
            update_ext_vec(pSet, C_idx, u, s, ext_vec, current_best_idx, last_best, frame, cmp_option, R, error);
            //R->print();
            //std::cout << C_idx.size() << "\n";
            if (C_idx.size() == 1)
                break;

            //update candidate set
            rtree_pruning(pSet, C_idx, R);
        }

        int pindex = C_idx[0];
        rankingSet->points.push_back(pSet->points[pindex]);
        pSet->points.erase(pSet->points.begin() + pindex);
        rankingSet->fairTopk(k, category, bound,resultSet);
    }

    if(cmp_option == RANDOM)
        resultSet->printResult("UHRandom", numOfQuestion, realResult, u, t1, 0, 0);
    else
        resultSet->printResult("UHSimplex", numOfQuestion, realResult, u, t1, 0, 0);

}










/**
 * @brief The interactive algorithm UH-Simplex-Adapt.
 *        Find a points/return a point which satisfy the regret ratio
 *        Modified stopping condition, prune points which are not able to be top-k
 * @param original_P     the original dataset
 * @param u              the unkonwn utility vector
 * @param s              the question size
 * @param k              The threshold top-k
 */
int UHSimplex0(point_set *originalSet, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    point_set *pset = new point_set(originalSet);
    int dim = pset->points[0]->dim, numOfQuestion = 0, ccount; //the number of questions asked

    bool orderd = false, last_exist = true;
    std::vector<point_t *> near_set;
    std::vector<int> considered_id;
    Partition *R = new Partition(dim);
    int idx1 = pset->findBest(R->average_point());
    int idx2;
    // interactively reduce the candidate set and shrink the candidate utility range
    while (pset->points.size() > k && (!orderd))
    {
        R->print();
        numOfQuestion++;
        bool choose = false;
        int M = pset->points.size(), considered_count = 0;
        considered_id.clear();
        //considered_id.shrink_to_fit();
        bool **is_considered = new bool *[M];
        for (int i = 0; i < M; i++)
            is_considered[i] = new bool[M];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < M; j++)
                is_considered[i][j] = false;
        }
        if (!last_exist)
        {
            while (!choose)
            {
                idx1 = rand() % pset->points.size();
                idx2 = rand() % pset->points.size();
                if (!is_considered[idx1][idx2] && !is_considered[idx2][idx1])
                {
                    is_considered[idx1][idx2] = true;
                    is_considered[idx2][idx1] = true;
                    considered_count++;
                    if (!pset->points[idx1]->is_same(pset->points[idx2]))
                    {
                        hyperplane *h = new hyperplane(pset->points[idx1], pset->points[idx2]);
                        if (R->check_relation(h) == 0)
                            choose = true;
                    }
                }
                if (considered_count >= (M * (M + 1)) / 2 && (!choose))
                {
                    orderd = true;
                    break;
                }
            }
        }
        else
        {
            considered_id.push_back(pset->points[idx1]->id); //record the point id, we have used for pixel
            ccount = 0;
            near_set.clear();
            //near_set.shrink_to_fit(); //record the nearest point of p_set[idx1]
            for (int q = 0; q < pset->points.size(); q++)
            {
                if (q != idx1)
                {
                    pset->points[q]->value = 0;
                    for (int j = 0; j < dim; j++)
                    {
                        pset->points[q]->value += (pset->points[q]->attr[j] - pset->points[idx1]->attr[j]) *
                                           (pset->points[q]->attr[j] - pset->points[idx1]->attr[j]);
                    }
                    pset->points[q]->value = sqrt(pset->points[q]->value);
                    pset->points[q]->place = q;
                    near_set.push_back(pset->points[q]);
                }
            }
            sort(near_set.begin(), near_set.end(), valueCmp());
            int idx = -1;
            while (!choose)
            {
                ++idx;
                idx2 = (near_set[idx]->place) % pset->points.size();
                ccount++;
                if (!is_considered[idx1][idx2] && !is_considered[idx2][idx1])
                {
                    is_considered[idx1][idx2] = true;
                    is_considered[idx2][idx1] = true;
                    considered_count++;
                    if (!pset->points[idx1]->is_same(pset->points[idx2]))
                    {
                        hyperplane *h = new hyperplane(pset->points[idx1], pset->points[idx2]);
                        if (R->check_relation(h) == 0)
                            choose = true;
                    }
                }
                if (considered_count >= (M * (M + 1)) / 2 && (!choose))
                {
                    orderd = true;
                    break;
                }
                if (ccount >= near_set.size() && !choose)
                {
                    //find the neearby point which can be the pixel
                    bool not_use = true;
                    for(int j = 0; j < near_set.size() && not_use; ++j)
                    {
                        not_use = false;
                        for (int q = 0; q < considered_id.size(); q++)
                        {
                            if (considered_id[q] == near_set[j]->id)
                            {
                                not_use = true;
                                break;
                            }
                        }
                        idx1 = j;
                    }
                    if(not_use) //cd + 1 >= near_set.size())
                    {
                        orderd = true;
                        break;
                    }
                    considered_id.push_back(pset->points[idx1]->id);
                    near_set.clear();
                    for (int q = 0; q < pset->points.size(); q++)
                    {
                        if (q != idx1)
                        {
                            pset->points[q]->value = 0;
                            for (int j = 0; j < dim; j++)
                            {
                                pset->points[q]->value += (pset->points[q]->attr[j] - pset->points[idx1]->attr[j]) *
                                                   (pset->points[q]->attr[j] - pset->points[idx1]->attr[j]);
                            }
                            pset->points[q]->value = sqrt(pset->points[q]->value);
                            near_set.push_back(pset->points[q]);
                        }
                    }
                    sort(near_set.begin(), near_set.end(), valueCmp());
                    ccount = 0;
                }
            }
        }
        for(int i=0; i< M;i++)
            delete []is_considered[i];
        delete []is_considered;
        if (orderd)
            break;

        point_t *p1 = pset->points[idx1];
        point_t *p2 = pset->points[idx2];
        double v1 = p1->dot_product(u);
        double v2 = p2->dot_product(u);
        int bestID;
        if (v1 > v2)
        {
            bestID = p1->id;
            hyperplane *h = new hyperplane(pset->points[idx2], pset->points[idx1]);
            R->hyperplanes.push_back(h);
        }
        else
        {
            bestID = p2->id;
            hyperplane *h = new hyperplane(pset->points[idx1], pset->points[idx2]);
            R->hyperplanes.push_back(h);
            idx1 = idx2;
        }
        R->set_ext_pts();

        //update candidate set
        R->find_possible_topK(pset, k);
        //find whether the better point exist
        last_exist = false;
        for (int i = 0; i < pset->points.size(); i++)
        {
            if (pset->points[i]->id == bestID)
            {
                last_exist = true;
                break;
            }
        }
    }


    pset->printResult("UH-Simplex", numOfQuestion, s, t1, 0, 0);

    point_set *resultSet = new point_set();
    pset->findTopk(u, k, resultSet);
    bool resultcheck = true;
    for(int i = 0; i < s; ++i)
    {
        if(!resultSet->checkExist(pset->points[i]))
        {
            resultcheck = false;
            break;
        }
    }
    std::cout << resultcheck <<"\n";


    return numOfQuestion;
}



/**
 * @brief The interactive algorithm UH-Simplex-Adapt.
 *        Find a points/return a point which satisfy the regret ratio
 *        Modified stopping condition, prune points which are not able to be top-k
 * @param original_P     the original dataset
 * @param u              the unkonwn utility vector
 * @param s              the question size
 * @param k              The threshold top-k
 */
int UHSimplex(point_set *originalSet, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    point_set *pset = new point_set(originalSet);
    for(int i = 0; i < pset->points.size(); ++i)
        pset->points[i]->place = i;
    int dim = pset->points[0]->dim, numOfQuestion = 0, ccount; //the number of questions asked

    bool orderd = false, last_exist = true;
    std::vector<point_t *> near_set;
    std::vector<int> considered_id;
    Partition *R = new Partition(dim);
    point_t *pidx1 = pset->points[pset->findBest(R->average_point())], *pidx2;

    double M = pset->points.size(), considered_count = 0;
    considered_id.clear();
    //considered_id.shrink_to_fit();
    bool **is_considered = new bool *[M];
    for (int i = 0; i < M; i++)
        is_considered[i] = new bool[M];
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
            is_considered[i][j] = false;
    }


    // interactively reduce the candidate set and shrink the candidate utility range
    while (pset->points.size() > k && (!orderd))
    {
        std::cout << pset->points.size() << "\n";
        numOfQuestion++;
        considered_count = 0;
        bool choose = false;

        considered_id.push_back(pidx1->id); //record the point id, we have used for pixel
        //build near_set
        near_set.clear();
        for (int q = 0; q < pset->points.size(); q++)
        {
            if (pset->points[q]->id != pidx1->id)
            {
                pset->points[q]->value = 0;
                for (int j = 0; j < dim; j++)
                {
                    pset->points[q]->value += (pset->points[q]->attr[j] - pidx1->attr[j]) *
                                              (pset->points[q]->attr[j] - pidx1->attr[j]);
                }
                pset->points[q]->value = sqrt(pset->points[q]->value);
                near_set.push_back(pset->points[q]);
            }
        }
        sort(near_set.begin(), near_set.end(), valueCmp());
        //initial
        ccount = -1;

        while (!choose)
        {
            ++ccount;
            pidx2 = near_set[ccount];
            if (!is_considered[pidx1->place][pidx2->place])
            {
                is_considered[pidx1->place][pidx2->place] = true;
                is_considered[pidx2->place][pidx1->place] = true;
                considered_count++;
                if (!pidx1->is_same(pidx2))
                {
                    hyperplane *h = new hyperplane(pidx1, pidx2);
                    if (R->check_relation(h) == 0)
                        choose = true;
                    delete h;
                }
            }
            double MMM = (M * (M + 1)) / 2;
            if (considered_count + 1 >= MMM && (!choose))
            {
                orderd = true;
                break;
            }
            if (ccount >= near_set.size() - 1 && !choose)
            {
                //find the neearby point which can be the pixel
                bool not_use = true;
                for (int j = 0; j < near_set.size() && not_use; ++j)
                {
                    not_use = false;
                    for (int q = 0; q < considered_id.size(); q++)
                    {
                        if (considered_id[q] == near_set[j]->id)
                        {
                            not_use = true;
                            break;
                        }
                    }
                    pidx1 = near_set[j];
                }
                if (not_use) //cd + 1 >= near_set.size())
                {
                    orderd = true;
                    break;
                }
                considered_id.push_back(pidx1->id);
                near_set.clear();
                for (int q = 0; q < pset->points.size(); q++)
                {
                    if (pset->points[q]->id != pidx1->id)
                    {
                        pset->points[q]->value = 0;
                        for (int j = 0; j < dim; j++)
                        {
                            pset->points[q]->value += (pset->points[q]->attr[j] - pidx1->attr[j]) *
                                                      (pset->points[q]->attr[j] - pidx1->attr[j]);
                        }
                        pset->points[q]->value = sqrt(pset->points[q]->value);
                        near_set.push_back(pset->points[q]);
                    }
                }
                sort(near_set.begin(), near_set.end(), valueCmp());
                ccount = -1;
            }
        }

        if (orderd)
            break;


        double v1 = pidx1->dot_product(u);
        double v2 = pidx2->dot_product(u);
        //std::cout << "P1:" << pidx1->id <<"  P2: " << pidx2->id << "\n";
        if (v1 > v2)
        {
            hyperplane *h = new hyperplane(pidx2, pidx1);
            R->hyperplanes.push_back(h);
        } else
        {
            hyperplane *h = new hyperplane(pidx1, pidx2);
            R->hyperplanes.push_back(h);
            pidx1 = pidx2;
        }
        R->set_ext_pts();
        //R->print();

        //update candidate set
        R->find_possible_topK(pset, k);

        //find whether the better point exist
        last_exist = false;
        for (int i = 0; i < pset->points.size(); i++)
        {
            if (pset->points[i]->id == pidx1->id)
            {
                last_exist = true;
                break;
            }
        }
        if(!last_exist)
            pidx1 = pset->points[pset->findBest(R->average_point())];
    }

    if(pset->points.size() > k)
    {
        point_set *ppt = new point_set();
        pset->findTopk(R->average_point(), k, ppt);
        pset = ppt;
    }
    pset->printResult("UH-Simplex", numOfQuestion, s, t1, 0, 0);
    point_set *resultSet = new point_set();
    originalSet->findTopk(u, k, resultSet);
    bool resultcheck = true;
    for(int i = 0; i < s; ++i)
    {
        if(!resultSet->checkExist(pset->points[i]))
        {
            resultcheck = false;
            break;
        }
    }
    std::cout << resultcheck <<"\n";


    return numOfQuestion;
}
























