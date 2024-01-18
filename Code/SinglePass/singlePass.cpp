#include "singlePass.h"


void singlePass(point_set *pset, point_t *u, int k, point_set *realResult, std::vector<std::string> &category,
                    std::vector<std::vector<double>> &bound, bool error)
{

    timeval t1;
    gettimeofday(&t1, 0);


    point_set *rankingSet = new point_set(), *resultSet = new point_set();;
    int numOfQuestion = 0;
    double theta = 5.0/8.0;
    std::vector<point_set *> S;
    point_set *P = new point_set(), *filter = new point_set();
    point_set *data = new point_set(pset);

    while(1)
    {
        double p_size = ceil(64 * log(2 * pset->points.size()));
        for (int i = 0; i < data->points.size(); i++)
        {
            //check if it is pruned
            bool is_pruned = false;
            for (int j = 0; j < S.size(); ++j)
            {
                if (S[j]->is_prune(data->points[i]))
                {
                    is_pruned = true;
                    break;
                }
            }
            if (is_pruned)
                continue;

            //fill in P
            if (P->points.size() < p_size)
            {
                P->points.push_back(data->points[i]);
                continue;
            }

            //add the point to the filter
            int left = 0, right = filter->points.size() - 1;
            while (left <= right)
            {
                ++numOfQuestion;
                int mid = (left + right) / 2;
                double v1 = filter->points[mid]->dot_product(u);
                double v2 = data->points[i]->dot_product(u);
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
                    left = mid + 1;
                else
                    right = mid - 1;
            }
            filter->points.insert(filter->points.begin() + left, data->points[i]);

            //check if there are any points in P that can be pruned by the filter
            point_set *PP = new point_set();
            for (int j = 0; j < P->points.size(); ++j)
            {
                if (filter->is_prune(P->points[j]))
                {
                    PP->points.push_back(P->points[j]);
                }
            }

            //filter->print();
            //PP->print();
            //update
            if (PP->points.size() >= theta * P->points.size())
            {
                S.push_back(filter);
                filter = new point_set();
                P->subtract(PP);
            }
        }

        if(filter->points.size() > 0)
            S.push_back(filter);
        point_t *best = P->points[0];
        if(P->points.size() <= 0)
        {
            best = new point_t(S[0]->points[0]->dim);
            for(int i = 0; i < S[0]->points[0]->dim; ++i)
                best->attr[i] = 0;
        }

        //find the best from P
        for (int i = 1; i < P->points.size(); ++i)
        {
            ++numOfQuestion;
            double v1 = best->dot_product(u);
            double v2 = P->points[i]->dot_product(u);
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
            if(v1 < v2 && sumr < prob || v1 > v2 && sumr > prob)
                best = P->points[i];
        }
        //find the best from S
        int s_index = -1;
        for (int i = 0; i < S.size(); ++i)
        {
            ++numOfQuestion;
            double v1 = best->dot_product(u);
            double v2 = S[i]->points[0]->dot_product(u);
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
            if(v1 < v2 && sumr < prob || v1 > v2 && sumr > prob)
            {
                best = S[i]->points[0];
                s_index = i;
            }
        }

        best->print();
        rankingSet->points.push_back(best);
        resultSet->points.clear();
        rankingSet->fairTopk(k, category, bound,resultSet);
        if(resultSet->points.size() >= k)
            break;

        pset->prunePt(best);
        data->points.clear();
        if (s_index > -1)
        {
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (S[s_index]->is_prune(pset->points[i]))
                    data->points.push_back(pset->points[i]);
            }
            if(s_index != S.size() - 1 || filter->points.size() <= 0)
            {
                S.erase(S.begin() + s_index);
            }
        }
        else
        {
            P->prunePt(best);
        }
        if(filter->points.size() > 0)
            S.pop_back();
    }

    resultSet->printResult("singlePass", numOfQuestion, realResult, u, t1, 0, 0);
}