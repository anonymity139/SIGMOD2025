#ifndef RUN_INTERACTION_TREE_H
#define RUN_INTERACTION_TREE_H
#include "hyperplane_set.h"

class it_node
{
public:
    int ID;
    std::vector<point_t*> ptSet;
    hyperplane *divideHyper, *bestdivideHyper;
    it_node *pt, *ng, *origin;
    it_node *bestpt, *bestng, *bestorigin;

    it_node();
    it_node(int id, it_node *org);
    it_node(int id, std::vector<point_t*> &ps);
    ~it_node();
    int insert(hyperplane_set* hset, std::vector<point_t*> ptSet, int height, int partitionLeft);

    void print();
    void print_leaves();
};


class it_tree
{
public:
    it_node *root;

    it_tree();
    it_tree(hyperplane_set* hset, std::vector<point_t*> ptSet, int partitionLeft);
    void print();
};




















#endif //RUN_INTERACTION_TREE_H
