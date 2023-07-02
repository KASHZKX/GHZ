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
    int source, destination;
    int time_limit;
    int status;
    int send_path_length;
    double fidelity;
    vector<Path *> paths;                                       //休學
public:
    Request(int source, int destination, const int& time_limit);
    ~Request(void);
    void set_path(int path_id, Path *p);                        //should delete old path before set new path
    int get_time_limit();
    int get_waiting_time();
    int get_source();
    int get_destination();
    int get_throughput();
    int get_send_path_length();
    double get_fidelity();
    vector<Path *> get_paths();
    
    void clear_paths();
    void refresh_paths();
    void add_one_throughput();
    void entangle();
    void swap();
    void send();
    bool is_finished();
    bool is_success();
    void next_timeslot();
    void operator+=(Path *path);
    void print();
};

#endif