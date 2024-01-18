#include "interaction_tree.h"


/**
 * @brief Constructor
 */
it_node::it_node()
{
    pt = NULL;
    ng = NULL;
    bestpt = NULL;
    bestng = NULL;
    bestdivideHyper = NULL;
}

/**
 * @brief Constructor
 */
it_node::it_node(int id, it_node *org)
{
    ID = id;
    origin = org;
    pt = NULL;
    ng = NULL;
}



/**
 * @brief Constructor
 */
it_node::it_node(int id, std::vector<point_t*> &ps)
{
    ID = id;
    ptSet = ps;
}


/**
 * @brief Insert a hyper-plane into a node
 * @param h
 */
int it_node::insert(hyperplane_set* hset, std::vector<point_t*> ptSet, int height, int partitionLeft)
{
    this->ptSet = ptSet;

    //if it reaches the leaf node
    if (ptSet.size() <= partitionLeft)
        return height;

    //if it reaches the internal node
    int minTree = 999999;
    std::vector<std::vector<point_t *>> upSet, downSet;
    std::vector<int> divideSizes;
    for (int t = 0; t < hset->hyperplanes.size(); ++t)
    {
        std::vector<point_t *> upSettest, downSettest;
        for (int j = 0; j < ptSet.size(); ++j)
        {
            if (hset->hyperplanes[t]->check_position(ptSet[j]) == 1)
                upSettest.push_back(ptSet[j]);
            else
                downSettest.push_back(ptSet[j]);
        }
        divideSizes.push_back(std::min(upSettest.size(), downSettest.size()));
        upSet.push_back(upSettest);
        downSet.push_back(downSettest);
    }

    //std::cout << hset->hyperplanes.size() << std::endl;
    for (int r = 0; r < 4 && r < hset->hyperplanes.size(); r++)
    {
        int IndHyper, maxPartition = -1;
        for (int t = 0; t < hset->hyperplanes.size(); ++t)
        {
           if(divideSizes[t] > maxPartition)
           {
                maxPartition = divideSizes[t];
                IndHyper = t;
           }
        }
        divideSizes[IndHyper] = -1;

        if(upSet[IndHyper].size() == 0 || downSet[IndHyper].size() == 0)
            continue;
        divideHyper = hset->hyperplanes[IndHyper];
        int maxH = std::max(ceil(log2(upSet[IndHyper].size())), ceil(log2(downSet[IndHyper].size())));
        if (maxH < minTree)
        {
            pt = new it_node(2 * ID, this);
            ng = new it_node(2 * ID + 1, this);
            hset->hyperplanes.erase(hset->hyperplanes.begin() + IndHyper);
            int m1 = pt->insert(hset, upSet[IndHyper], height + 1, partitionLeft);
            int m2 = ng->insert(hset, downSet[IndHyper], height + 1, partitionLeft);
            hset->hyperplanes.insert(hset->hyperplanes.begin() + IndHyper, divideHyper);
            maxH = std::max(m1, m2);
            if (maxH < minTree)
            {
                minTree = maxH;
                bestpt = pt;
                bestng = ng;
                bestdivideHyper = divideHyper;
            }
        }

    }
    return minTree + 1;
}


/**
 * @brief Print the information of the node
 */
void it_node::print()
{
    std::cout << "Node: " << ID <<"   ";
    if(pt == NULL && ng == NULL)
        std::cout << "Leave \n";
    else
        std::cout << "Internal Node\n";
    for (int i = 0; i < ptSet.size(); ++i)
    {
        ptSet[i]->print();
    }

}

/**
 * @brief Destructor
 */
it_node::~it_node()
{
    if(pt != NULL)
        delete pt;
    if(ng != NULL)
        delete ng;
}





/**
 * @Constrcutor
 */
it_tree::it_tree()
{
    root = new it_node();
    root->ID = 1;
    root->origin = NULL;
}


/**
 * @Constrcutor
 */
it_tree::it_tree(hyperplane_set* hset, std::vector<point_t*> ptSet, int partitionLeft)
{
    root = new it_node(1, ptSet);
    root->origin = NULL;
    root->insert(hset, ptSet, 0, partitionLeft);
}




/**
 * @brief Print the tree
 */
void it_tree::print()
{
    std::cout << "TREE: " << "\n";
    std::vector<it_node*> L;
    L.push_back(root);
    while(L.size() > 0)
    {
        L[0]->print();
        if(L[0]->bestpt != NULL && L[0]->bestng != NULL)
        {
            L.push_back(L[0]->bestpt);
            L.push_back(L[0]->bestng);
        }
        L.erase(L.begin());
    }
}

