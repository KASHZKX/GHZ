#include "BaselineDual.h"

BaselineDual::BaselineDual(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not)
    :AlgorithmBase(filename, "BaselineDual", request_time_limit, node_time_limit, swap_prob, entangle_alpha , limit_r_or_not){
    if(DEBUG) cerr<<"new BaselineDual"<<endl;
}

void BaselineDual::dual_initialize(){
    M = graph.get_size() + graph.get_num_of_edge() + requests.size();              // M = V + E + I         
    delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * M, 1.0 / epsilon));         

    for(int i = 0; i < graph.get_size(); i++){
        alpha.emplace_back(delta / graph.Node_id2ptr(i)->get_memory_cnt());       // alpha = delta / memory on node v

        vector<int> temp = graph.get_neighbors_id(i);                             // beta = delta /  channel on edge(u, v) 
        for(auto it: temp){
            beta[make_pair(i,it)] = delta / (graph.get_channel_size(i, it));
        }
    }
    for(unsigned int i = 0; i < requests.size(); i++){                               // tau = delta  
        tau.emplace_back(delta);                                            
    }
}

double BaselineDual::calculateDual(int u, int v, int req_no, int path_id,int middle){
	double weight = alpha[u] + alpha[v] + beta[{u, v}];
	//if(path_id == 0 && (u == requests[req_no].get_node1() || v == requests[req_no].get_node1())) weight += tau[req_no];                                                         //[Need fix!!!!!!!!!!!!!]
	if(u == middle || v == middle) weight += tau[req_no]/3;
    return weight;
}

vector<vector<int>> BaselineDual::Dijkstra_Tree(vector<int> &terminal, int req_no){
    const double INF = numeric_limits<double>::infinity();
    int node_num = graph.get_size();
    double min_weight = INF;
    vector<vector<int>> tree;
    tree.clear();
    for(int i = 0; i < node_num; i++){
        if(i == terminal[0] || i == terminal[1] || i == terminal[2]) continue; 
        vector<vector<int>> tmp_tree;
        tmp_tree.clear();
        double weight_sum = 0;
        int j = 0;
        for(j = 0; j < 3; j++){
            int destination = terminal[j];
            vector<double> dist(node_num + 1, INF); 
            vector<int> prev(node_num + 1, -1);
            priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq; 
            
            dist[i] = 0;
            pq.push({0, i});

            while (!pq.empty()) {
                double distance = pq.top().first;
                int u = pq.top().second;
                pq.pop();

                if (u == destination) break; 

                for(int v : graph.get_neighbors_id(u)) {
                    if (u == v) continue; 

                    double weight = calculateDual(u, v, req_no, j,i);
                    if (dist[v] > dist[u] + weight) { 
                        bool is1_repeater = true, is2_repeater = true;
                        if(u == i) is1_repeater = false;
                        if(v == destination) is2_repeater = false;
                        if(graph.remain_resource_cnt(u, v, is1_repeater, is2_repeater) > 0){
                            dist[v] = dist[u] + weight;
                            prev[v] = u;
                            pq.push({dist[v], v});
                        }
                    }
                }
            }
            if(prev[destination] == -1) break;
            weight_sum += dist[destination];
            vector<int> path;
            path.clear();
            for (int at = destination; at != -1; at = prev[at])
                path.push_back(at);
            reverse(path.begin(), path.end());

            tmp_tree.push_back(path);
        }
        
        if(min_weight > weight_sum && weight_sum != 0 && j == 3){
            min_weight = weight_sum;
            tree = tmp_tree;          
        }
    }
    if(tree.size() == 0)
        return {};
    
    return tree;
}

void BaselineDual::updateDual(vector<vector<int>> &tree, int req_no){
    double min_s_u = numeric_limits<double>::infinity();
    double min_s_uv = numeric_limits<double>::infinity();

    vector<double> s_u(graph.get_size() + 5);                                              
    vector<int> used_u(graph.get_size() + 5, 0);
    map<pair<int,int>,int> used_uv;
    map<pair<int,int>,double> s_uv;
    // Calculate  memory and channel used 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < tree[i].size() - 1; j++){
            if(j == 0)
                used_u[tree[i][j]]++;
            else
                used_u[tree[i][j]] += 2;

            if(used_uv.find({tree[i][j],tree[i][j+1]}) != used_uv.end()){
                used_uv[{tree[i][j],tree[i][j+1]}]++;
                used_uv[{tree[i][j+1],tree[i][j]}]++;
            }
            else{
                used_uv[{tree[i][j],tree[i][j+1]}] = 1;
                used_uv[{tree[i][j+1],tree[i][j]}] = 1;
            }
        }
        used_u[tree[i][tree[i].size()-1]]++;
    }

    for(int i = 0; i < graph.get_size(); i++){
        if(used_u[i] == 0) continue;
        s_u[i] = graph.Node_id2ptr(i)->get_memory_cnt() / used_u[i] ;        
        if(s_u[i] < min_s_u)
            min_s_u = s_u[i];
    }

    for(auto &it:used_uv){
        s_uv[{it.first.first, it.first.second}] = graph.get_channel_size(it.first.first, it.first.second) / it.second;                   //channel min
        if(s_uv[{it.first.first, it.first.second}] < min_s_uv)
            min_s_uv = s_uv[{it.first.first, it.first.second}];
    }

    int rate = 1;
    double s = min(min_s_u, min(min_s_uv, 1.0));                                    
    
    //alter dual object
    for(int i = 0; i < rate; i++){
        
        for(int i = 0; i < used_u.size(); i++){
            if(used_u[i] == 0) continue; 
            alpha[i] = alpha[i] * (1 + epsilon * s / s_u[i]);
        }

        for(auto &it : used_uv){
            beta[{it.first.first, it.first.second}] = beta[{it.first.first, it.first.second}] * (1 + epsilon * s / s_uv[{it.first.first, it.first.second}]);
            beta[{it.first.second, it.first.first}] = beta[{it.first.second, it.first.first}] * (1 + epsilon * s / s_uv[{it.first.second, it.first.first}]);
        }
        tau[req_no] = tau[req_no] * (1 + epsilon * s / 1);    
    }
}

bool BaselineDual::checkResource(vector<vector<int>> &tree){
    bool resource_enough = true;
    vector<int> memory_use(graph.get_size(), 0);
    map<pair<int,int>,int> channel_use;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < tree[i].size() - 1; j++){
            memory_use[tree[i][j]]++;
            memory_use[tree[i][j+1]]++;
            if(channel_use.find({tree[i][j],tree[i][j+1]}) != channel_use.end()){
                channel_use[{tree[i][j],tree[i][j+1]}]++;
                channel_use[{tree[i][j+1],tree[i][j]}]++;
            }
            else{
                channel_use[{tree[i][j],tree[i][j+1]}] = 1;
                channel_use[{tree[i][j+1],tree[i][j]}] = 1;
            }
        }
    }

    for(int i = 0; i < graph.get_size(); i++){
        if(graph.Node_id2ptr(i)->get_remain() < memory_use[i]){
            resource_enough = false;
            break;
        }
    }
    for(auto &it:channel_use){
        //cout<<"Channel_remain : "<<graph.remain_channel(it.first.first,it.first.second)<<" "<<it.second<<endl;
        if(graph.remain_channel(it.first.first,it.first.second) < it.second){
            resource_enough  = false;
            break;
        }
    }
    return resource_enough;   
}

void BaselineDual::path_assignment(){
    dual_initialize();
    bool flag = true;
    while(flag){
        flag = false;
        for(int i = 0; i < (int)requests.size(); i++){
            Request &request = requests[i];
            vector<int> terminal = {request.get_node1(), request.get_node2(), request.get_node3()};
            vector<vector<int>> tree = Dijkstra_Tree(terminal, i);
            if(tree.size() == 0) continue;
            if(checkResource(tree)){
                flag = true;
                // cout << "assign path" << endl;
                assign_resource(tree, 1, i);
                updateDual(tree,i);                //[added]
            }
            
        }
    }
   
     
}

void BaselineDual::entangle(){
    AlgorithmBase::base_entangle();
}
void BaselineDual::swap(){
    AlgorithmBase::base_swap();
}
void BaselineDual::send(){
    AlgorithmBase::base_send();
}
void BaselineDual::next_time_slot(){
    AlgorithmBase::base_next_time_slot();
}

