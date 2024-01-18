#include "structure/data_utility.h"
#include "structure/data_struct.h"
#include "structure/point_set.h"
#include "structure/define.h"
#include <vector>
#include <ctime>
#include <sys/time.h>
#include "UH/UH.h"
#include "GroundTruth/Groundtruth.h"
#include "TDIA//PI.h"
#include "RH/RH.h"
#include "ActiveRanking/ActiveRanking.h"
#include "PreferenceLearning/preferenceLearning.h"
#include "UH/UH.h"
#include "UH/Median_Hull.h"
#include "HDIA/Piecewise_exhaustion.h"
#include "UH/UH.h"
#include "SinglePass/singlePass.h"

int main(int argc, char *argv[])
{
    //input all the configuration
    ifstream config("../config.txt");
    string alg_name, data_name;
    double k; int numCate;
    config >> alg_name >> data_name >> k;
    cout << "Algorithm: " << alg_name << "   Dataset: " << data_name << "   k: " << k << "\n";
    point_set *pSet = new point_set(data_name);

    //skyband
    vector<string> category; vector<vector<double>> bound;
    pSet->skyband_with_clusters(category, k, bound);
    pSet->write("../input/skyband.txt");

    /*bound
    vector<vector<double>> bound(numCate, vector<double>(2, 0));
    for(int i = 0; i < numCate; ++i)
        config >> bound[i][0] >> bound[i][1];
    double low_sum = 0, up_sum = 0;
    for(int i = 0; i < numCate; ++i)
    {
        low_sum += bound[i][0];
        up_sum += bound[i][1];
    }
    if(low_sum > k || up_sum < k)
    {
        printf("Error: k is not compatible with the bounds\n");
        return 0;
    }
    */

    //utility vector
    int dim = pSet->points[0]->dim;
    point_t *u = new point_t(dim);
    for(int i = 0; i < dim; ++i)
        config >> u->attr[i];
    cout << "utility vector: "; u->print();

    //ground truth
    printf("-----------------------------------------------------------------------------------\n");
    printf("|%15s |%15s |%15s |%15s |%10s |\n", "Algorithm", "# of Questions", "Time", "Accuracy", "Point #ID");
    printf("-----------------------------------------------------------------------------------\n");
    point_set *resultSet = ground_truth(pSet, u, k, category, bound);

    if(alg_name == "PreferenceLearning")
    {
        PreferenceLearning(pSet, u, k, resultSet, category, bound, 0);
    }
    else if(alg_name == "PreferenceLearningError")
    {
        PreferenceLearning(pSet, u, k, resultSet, category, bound, 1);
    }
    else if(alg_name == "ActiveRanking")
    {
        ActiveRanking(pSet, u, k, resultSet, category, bound);
    }
    else if(alg_name == "ActiveRankingRobust")
    {
        ActiveRanking_robust(pSet, u, k, resultSet, category, bound);
    }
    else if(alg_name == "RH")
    {
        RH(pSet, u, k, resultSet, category, bound, 0);
    }
    else if(alg_name == "RHError")
    {
        RH(pSet, u, k, resultSet, category, bound, 1);
    }
    else if(alg_name == "TDIA")
    {
        twoPI(pSet, u, k, resultSet, category, bound);
    }
    else if(alg_name == "HDIA")
    {
        Piecewise_exhaustion(pSet, u, k, resultSet, category, bound);
    }
    else if(alg_name == "HDIAR")
    {
        Piecewise_exhaustion_robust(pSet, u, k, resultSet, category, bound);
    }
    else if(alg_name == "UHSimplex")
    {
        max_utility(pSet, u, k, resultSet, category, bound, 2, SIMPLEX, 0);
    }
    else if(alg_name == "UHSimplexError")
    {
        max_utility(pSet, u, k, resultSet, category, bound, 2, SIMPLEX, 1);
    }
    else if(alg_name == "singlePass")
    {
        singlePass(pSet, u, k, resultSet, category, bound, 0);
    }
    else if(alg_name == "singlePassError")
    {
        singlePass(pSet, u, k, resultSet, category, bound, 1);
    }
    else
        printf("Error: algorithm name is not correct\n");


    delete pSet;
    return 0;
}
