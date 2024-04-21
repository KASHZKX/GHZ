#include <iostream>
#include <queue>
#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

int main(void){
    ifstream ifs;
    ifs.open("src/4_19/input/round_1.input");
    int num_node, num_of_link;
    int x,y,memory;
    float swap_prob,fusion_prob;
    ifs >> num_node;
    vector<pair<double,double>> pos(num_node,{0,0});
    for(int i = 0; i < num_node; i++){
        ifs >> pos[i].first >> pos[i].second >> memory >> swap_prob >> fusion_prob;
    }

    ifs >>  num_of_link;
    double total_dist = 0;
        for(int i = 0; i < num_of_link; i++){
        int node1,node2,link,fideilty;
        ifs >> node1 >> node2 >> link >> fideilty;
        total_dist += sqrt(pow( pos[node1].first - pos[node2].first ,2) + pow( pos[node1].second - pos[node2].second ,2));
    }
    ifs.close();
    cout << total_dist / num_of_link;
    
}