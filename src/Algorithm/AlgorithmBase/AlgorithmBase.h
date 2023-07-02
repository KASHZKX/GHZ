#ifndef __ALGORITHMBASE_H
#define __ALGORITHMBASE_H

#include    <iostream>
#include    <vector>
#include    <queue>
#include    <algorithm>
#include    <string>
#include    <map>
#include    <sys/time.h>
#include "../../Network/Graph/Graph.h"
#include "../../Request/Request.h"
#include "../../config.h"
using namespace std;


class AlgorithmBase{
    string algorithm_name;
protected:
    int timeslot, waiting_time;
    int time_limit;
    double swap_prob;
    map<string, double> res;
    vector<Request> requests;
    //Request* generate_new_request();
public:
    Graph graph;
    AlgorithmBase(string filename, string algorithm_name, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha);
    virtual~AlgorithmBase();
    double get_swap_prob();
    string get_name();
    double get_res(string s);
    virtual void path_assignment()=0;
    int find_width(vector<int> path);
    vector<int> BFS(int source, int destination);
    void assign_resource(vector<int> path, int reqno);
    void base_entangle();
    void base_swap();
    void base_send();
    void base_next_time_slot();
    void base_test_active();
    //void request_received();
    void check_resource();
    virtual void entangle()=0;
    virtual void swap()=0;
    virtual void send()=0;
    virtual void next_time_slot()=0;

    void run();
    void add_new_request(Request new_request);
    int total_throughput();

    Path *find_swap_path(vector<int> path_nodes, map<pair<int, int>, vector<Channel*>> &remain_channels);
};

#endif
