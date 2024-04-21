#ifndef __MyAlgo_H
#define __MyAlgo_H

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <queue>
#include <limits>
#include <string>
#include <cmath>
#include <tuple>
#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Request/WholeRequest.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"
using namespace std;

class Label {
public:
  double theta, obj_bar2;
  int i;

  Label(double t, double o, int i) : theta(t), obj_bar2(o), i(i){};

    bool operator>(const Label &b) const {
    if (theta == b.theta)
      return obj_bar2 > b.obj_bar2;
    return theta > b.theta;
  }
};

class MyAlgo:public AlgorithmBase {
private:
    //map<pair<int,int>, vector<double>> X;
    vector<vector<vector<map<pair<int,int>, double>>>> Y;
    vector<double> alpha;
    vector<vector<vector<int>>> x_i_t_tree;
    vector<double> x_i_t;
    map<pair<int,int>, double> beta;
    vector<vector<vector<int>>> all_source_target_path;
    int change_edge_num = 0; 
    int diff_num = 0;
    double epsilon = 0.2;
    double obj = 0;
    double delta;
    double M;
    vector<double> tau;                
    double X(int u, int v, int req_no, int path_id, int middle);
public: 
    map<int, int> num_of_path_count;
    map<int, int> path_length_encode;
    map<int, int> path_length_cnt;
    void separation_oracle(int src, int dst, int req_no, int path_id, vector<vector<vector<int>>> &cur_tree, vector<vector<vector<double>>> &cur_label, vector<vector<int>> &n_graph, int &U);
    vector<int> Dijkstra(int src, int dst, int req, int path_id, vector<pair<double,double>>&dist, vector<vector<int>> &n_graph);
    void path_assignment();
    void calculate();
    void entangle();
    void swap();
    void send();
    vector<vector<int>> allPathsSourceTarget(int src, int dst);
    void next_time_slot();
    vector<map<vector<vector<int>>, int>> Greedy_rounding();
    void find_bottleneck(vector<vector<int>> &tree, int req_no);
    void initialize(int mid);
    void init_dual();
    void create_graph(vector<vector<vector<vector<int>>>> &n_graph);
    double changing_obj();
    void find_violate();
    vector<int> Dijkstra_ori(int src, int dst, int req_no, int path_id);
    void separation_oracle_ori(int src,int dst,int req_no, int path_id);
    MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha);
    MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha, double epsilon);
    ~MyAlgo();
};

#endif
