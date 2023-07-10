#ifndef __MYALGO_H
#define __MYALGO_H

#include <iostream>
#include <algorithm>
#include <queue>
#include <limits>
#include <string>
#include <cmath>
#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Request/WholeRequest.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"

using namespace std;

class MyAlgo:public AlgorithmBase {

private:
    map<pair<int,int>, double> X;
    map<pair<int,int>, double> Y;
    vector<double> alpha;
    map<pair<int,int>, double> beta;
    double delta;
    double M;
    vector<double> tau;                
public: 
    map<int, int> num_of_path_count;
    map<int, int> path_length_encode;
    map<int, int> path_length_cnt;
    void reweight_graph();
    void separation_oracle(int req_no);
    vector<int> Dijkstra(int src, int dst);
    void initialize(int epsilon);
    void run();
    MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha);
    ~MyAlgo();
};

#endif