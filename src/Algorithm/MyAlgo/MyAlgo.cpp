#include "MyAlgo.h"

MyAlgo::MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha)
    :AlgorithmBase(filename, "MyAlgo", request_time_limit, node_time_limit, swap_prob, entangle_alpha){
    if(DEBUG) cerr<<"new MyAlgo"<<endl;
}

void MyAlgo::reweight_graph(){

}

void MyAlgo::initialize(int epsilon){
    M = graph.get_size() + graph.get_num_of_edge() + requests.size(); //V+E+I
    delta = (1 + epsilon) * pow(((1 + epsilon) * M), (-1 / epsilon));

    for(int i = 0; i < graph.get_size(); i++){
        alpha[i] = delta / graph.Node_id2ptr(i)->get_memory_cnt();                //alpha_set
        vector<int> temp = graph.get_neighbors_id(i);                             //beta_set
        for(auto it: temp){
            if(i < it){
                beta[make_pair(i,it)] = delta / (graph.get_channel_size(i, it));
            }
            else{
                beta[make_pair(it,i)] = delta / (graph.get_channel_size(i, it));
            }
        }
    }
    for(int i = 0; i < requests.size(); i++){
        tau[i] = delta / requests[i].get_send_limit();
    }
    for(int i = 0; i < graph.get_size(); i++){
        vector<int> temp = graph.get_neighbors_id(i);                             
        for(auto it: temp){
            if(i < it){
                X[{i, it}] = alpha[i] + alpha[it] + beta[{i, it}];
                Y[{i, it}] = exp(graph.Node_id2ptr(i)->distance(*graph.Node_id2ptr(it))*(-graph.get_entangle_alpha()));
            }
            else{
                X[{it, i}] = alpha[it] + alpha[i] + beta[{it, i}];
                Y[{it, i}] = exp(graph.Node_id2ptr(i)->distance(*graph.Node_id2ptr(it))*(-graph.get_entangle_alpha()));
            }
        }
    }

}

vector<int> MyAlgo::Dijkstra(int src, int dst){ 
    const double INF = numeric_limits<double>::infinity();
    int n = graph.get_size();
    vector<double> distance(n, INF);
    vector<int> parent(n, -1);
    vector<bool> used(n, false);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    distance[dst] = 0;
    pq.push({0, dst});
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int neighbor : graph.get_neighbors_id(cur_node)) {
            if(distance[cur_node] + X[{cur_node, neighbor}] < distance[neighbor]) {
                distance[neighbor] = distance[cur_node] + X[{cur_node, neighbor}];
                parent[neighbor] = cur_node;
                pq.push({distance[neighbor], neighbor});
            }
        }
    }

    if(distance[src] >= INF) return{};
/*
    int cur_node = src;
    vector<int> path;
    while(cur_node != -1) {
        path.push_back(cur_node);
        cur_node = parent[cur_node];
    }

    reverse(path.begin(), path.end());
*/
    return parent;
    
}
    //     3 
    //     |
    //     |  
    // 2---1
    //  \  |
    //   \ |
    //     0

void MyAlgo::separation_oracle(int req_no){
    vector<vector<int>> SPT;
    int src = requests[req_no].get_source();
    int dst =  requests[req_no].get_destination();
    SPT[0] = Dijkstra(src, dst);      //the first SPT
    int SPTnum = 0;
    vector<int> path;
    while(U(x,y)<1){
        map<pair<int, int>, bool> used;
        pair<int,int> new_edge;
        double minimum = numeric_limits<double>::infinity();
        for(int i = 0; i < graph.get_size(); i++){
            vector<int> neighbors = graph.get_neighbors_id(i);  
            for(auto neighbor : neighbors){
                double temp1, temp2;
                if(SPT[0][i] == neighbor || SPT[0][neighbor] == i){  // calculate b / a
                    continue;
                }else{                                   // X[i, j] + X[p_j] - X[p_i] / Y[i, j] + Y[p_j] - Y[p_i]
                    temp1 = X[{i, neighbor}];
                    temp2 = Y[{i, neighbor}];
                    int prev_node;
                    int cur_node = i;
                    while(cur_node != dst){
                        temp1 += X[{cur_node, SPT[0][cur_node]}];
                        temp2 += Y[{cur_node, SPT[0][cur_node]}];
                        cur_node = SPT[0][cur_node];
                    } 
                    cur_node = neighbor;
                    while(cur_node != dst){
                        temp1 -= X[{cur_node, SPT[0][cur_node]}];
                        temp2 -= Y[{cur_node, SPT[0][cur_node]}];
                        cur_node = SPT[0][cur_node];
                    }       

                    if(temp2 < 0 && temp1 > 0){
                        if(minimum > -temp1 / temp2){
                            new_edge = {i, neighbor};
                            minimum = -temp1 / temp2;
                        }
                    }
                }
            }

        }       


            

        //give me a path vector
        
        double len;
        double c = 0;
        double r = 0;
        for(int i = 0; i < path.size() - 1; i++){
            if(path[i] < path[i+1]){
                c += X[{path[i], path[i+1]}];
                r += Y[{path[i], path[i+1]}];
            }
            else{
                c += X[{path[i+1], path[i]}];  
                r += Y[{path[i+1], path[i]}];           
            }
        }
        len = c / r;
                                                                        //此次最小
                                                                        //路線修改
                                                                        //新的spt產生
    }
                                                                       //很多SPT
                                                                       //最小的U=the shortest path
}

void MyAlgo::run(){
     
    double obj = M * delta; 
    while(obj < 1){
        int req_no = 0;
        separation_oracle(req_no);
        req_no++;

    }
}   

