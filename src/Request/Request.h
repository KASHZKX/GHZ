#ifndef __REQUEST_H
#define __REQUEST_H

#include<iostream>
#include<vector>
#include<algorithm>
#include"../Network/Node/Node.h"
#include"../Path/Path.h"
#include "../config.h"
using namespace std;

const int REQUEST_UNFINISHED = 0;
const int REQUEST_SUCC = 1;
const int REQUEST_FAIL = -1;


class Request{
protected:
    int node1, node2, node3;
    int cur_send = 0;
    int throughput = 0;
    int status;
    int tree_num = 0;
    vector<vector<Path *>> trees;                                       //休學-[哪個tree][tree的哪個邊]
    vector<vector<double>> tree_prob_vt;                     
public:
    Request(int node1, int node2, int node3);
    ~Request(void);
    void set_path(int path_id, int edge_id, Path *p);                        //should delete old path before set new path
    int get_node1();
    int get_node2();
    int get_node3();
    int get_throughput();
    int get_tree_num();
    int get_cur_send();
    void add_cur(int num);
    vector<vector<Path *>> get_trees();
    void clear_trees();
    void refresh_trees();
    void entangle();
    void swap();
    void send();
    bool is_finished();
    bool is_success();
    void next_timeslot();
    void operator+=(vector<Path *>tree);
};

#endif
