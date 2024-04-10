#ifndef __BaselineAllNode
#define __BaselineAllNode
#include <iostream>
#include <algorithm>
#include <queue>
#include <cmath>
#include <set>
#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"
using namespace std;

class BaselineAllNode : public AlgorithmBase{
private:
    vector<double> alpha;
    map<pair<int,int>, double> beta;
    vector<double> tau;                
    double epsilon = 0.2;
    double obj = 0;
    double delta;
    double M;

public:
    BaselineAllNode(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not);
    vector<vector<int>> Dijkstra_Tree(vector<int> &terminal);
    bool checkResource(vector<vector<int>> &tree);
    void path_assignment();
    void entangle();
    void swap();
    void send();
    void next_time_slot();

};

#endif
