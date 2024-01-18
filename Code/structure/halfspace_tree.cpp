#include "halfspace_tree.h"



/**
 * @brief Constructor
 */
hp_node::hp_node()
{
    cell = new Partition();
    pt = NULL;
    ng = NULL;
}


/**
 * @brief Constructor
 */
hp_node::hp_node(hp_node *org)
{
    cell = new Partition();
    origin = org;
    pt = NULL;
    ng = NULL;
}

/**
 * @brief Constructor
 */
hp_node::hp_node(int dim)
{
    cell = new Partition(dim);
    origin = NULL;
    pt = NULL;
    ng = NULL;
}


/**
 * @brief Constructor
 */
hp_node::hp_node(Partition *hset)
{
    cell = hset;
    pt = NULL;
    ng = NULL;
}


/**
 * @brief Insert a hyper-plane into a node
 * @param h
 */
void hp_node::insert(hyperplane *h)
{
    int relation = cell->check_relation(h);
    if (relation == 0)
    {
        if (pt == NULL && ng == NULL) //if it has child, there must be two children
        {
            pt = new hp_node(); pt->ID = 2 * ID;
            ng = new hp_node(); ng->ID = 2 * ID + 1;
            if(!cell->divide(h, pt->cell, ng->cell))
            {
                pt = NULL; ng = NULL;
                return;
            }
            divideHyper = h;
            //pt->cell->print();
            //ng->cell->print();
        }
        else
        {
            pt->insert(h);
            ng->insert(h);
        }
    }
}



/**
 * @brief Print the partitions that fulfill the requirements
 * @param k     Parameter k
 */
void hp_node::print_leaves()
{
    if(pt == NULL && ng == NULL)
    {
        std::cout << "Partition\n";
        cell->print();
    }
    else
    {
        if (pt != NULL)
            pt->print_leaves();
        if (ng != NULL)
            ng->print_leaves();
    }
}

/**
 * @brief Print the information of the node
 */
void hp_node::print()
{
    std::cout << "Node: " << ID <<"   ";
    if(pt == NULL && ng == NULL)
        std::cout << "Leave \n";
    else
        std::cout << "Internal Node\n";
    cell->print();
}


/**
 * @brief Find all the pts in the leaf nodes and the height of the tree
 * @param leafPt    The vector to store the leaf points
 * @param height    The current of the tree
 * @return          The height of the tree
 */
int hp_node::find_leafPt_height(std::vector<point_t*> &leafPt, int height)
{
    if(pt == NULL && ng == NULL)
    {
        leafPt.push_back(cell->center);
        cell->center->id = leafPt.size();
        //cell->print();
        //std::cout <<"\n";
        return height;
    }
    else
    {
        ++height;
        int h1 = pt->find_leafPt_height(leafPt, height);
        int h2 = ng->find_leafPt_height(leafPt, height);
        return std::max(h1, h2);
    }
}


/**
 * @brief Destructor
 */
hp_node::~hp_node()
{
    if(pt != NULL)
        delete pt;
    if(ng != NULL)
        delete ng;
}


/**
 * @Constrcutor
 */
halfspace_tree::halfspace_tree(int dim)
{
    root = new hp_node(dim);
    root->ID = 1;
    root->origin = NULL;
}


/**
 * @Constrcutor
 */
halfspace_tree::halfspace_tree(Partition *cell)
{
    root = new hp_node(cell);
    root->ID = 1;
    root->origin = NULL;
}


/**
 * @brief Insert a hyper-plane into the tree
 * @param h
 */
void halfspace_tree::insert(hyperplane *h)
{
    root->insert(h);
}

/**
 * @brief Print the partitions that fulfill the requirements
 * @param k     Parameter k
 */
void halfspace_tree::print_leaves()
{
    root->print_leaves();
}

/**
 * @brief Find the leaf node that contains the point and
 * @param pt
 * @return
 */
int halfspace_tree::find_leafPt_height(std::vector<point_t*> &leafPt)
{
    return root->find_leafPt_height(leafPt, 0);
}


/**
 * @brief Print the tree
 */
void halfspace_tree::print()
{
    std::cout << "TREE: " << "\n";
    std::vector<hp_node*> L;
    L.push_back(root);
    while(L.size() > 0)
    {
        L[0]->print();
        if(L[0]->pt != NULL && L[0]->ng != NULL)
        {
            L.push_back(L[0]->pt);
            L.push_back(L[0]->ng);
        }
        L.erase(L.begin());
    }
}





