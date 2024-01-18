#ifndef RUN_HALFSPACE_TREE_H
#define RUN_HALFSPACE_TREE_H
#include "hyperplane_set.h"
#include "Partition.h"



class hp_node
{
public:
    int ID;
    Partition *cell;
    hyperplane *divideHyper;
    hp_node *pt, *ng, *origin;

    hp_node();
    hp_node(hp_node *org);
    hp_node(int dim);
    hp_node(Partition *cell);
    ~hp_node();
    void insert(hyperplane *h);
    void updateThreshold(int k);
    int find_leafPt_height(std::vector<point_t*> &leafPt, int height);

    void print();
    void print_leaves();
};


class halfspace_tree
{
public:
    hp_node *root;

    halfspace_tree(int dim);
    halfspace_tree(Partition *cell);
    void insert(hyperplane *h);
    void print_leaves();
    void print();
    int find_leafPt_height(std::vector<point_t*> &leafPt);
};


#endif //RUN_HALFSPACE_TREE_H
