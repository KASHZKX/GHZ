#include "Greedy.h"

Greedy::Greedy(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not)
    :AlgorithmBase(filename, "Greedy", request_time_limit, node_time_limit, swap_prob, entangle_alpha , limit_r_or_not){
    if(DEBUG) cerr<<"new Greedy"<<endl;
}

void Greedy::cal_need(vector<int>path, vector<int>& need_memory, map<pair<int,int>, int> &need_channel){
    for(int i = 0; i < path.size()-1; i++){
        if(i == 0){
            need_memory[path[i]]++;
        }
        else{
            need_memory[path[i]]+=2;
        }
        if(path[i]<path[i+1]){
            if(need_channel.find({path[i],path[i+1]}) != need_channel.end()){
                need_channel[{path[i],path[i+1]}] ++;
            }
            else{
                need_channel[{path[i],path[i+1]}] = 1;
            }
        }
        else{
            if(need_channel.find({path[i+1],path[i]}) != need_channel.end()){
                need_channel[{path[i+1],path[i]}] ++;
            }
            else{
                need_channel[{path[i+1],path[i]}] = 1;
            }            
        }

    }
    need_memory[path[path.size()-1]]++;
    // for(auto it:need_memory){
    //     cout<<it<<" ";
    // }
    // cout<<endl;
    // for(auto it:need_channel){
    //     cout<<it.first.first<<"-"<<it.first.second<<" with "<<it.second<<endl;
    // }
}

vector<int> Greedy::Dijkstra(int src, int dst){                          
    double INF = numeric_limits<double>::infinity();
    vector<bool> used( graph.get_size(), false);
    vector<int> parent( graph.get_size(), -1);
    vector<double> dist(graph.get_size(), INF);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    pq.push({0, src});
    dist[src] = 0;
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(cur_node == dst) break;
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int &neigh : graph.get_neighbors_id(cur_node)) {
            double tmep_prob =-log(exp(graph.Node_id2ptr(cur_node)->distance(*graph.Node_id2ptr(neigh))*(-graph.get_entangle_alpha())));
            if(dist[cur_node] + tmep_prob < dist[neigh]) {
                dist[neigh] = dist[cur_node] + tmep_prob;
                parent[neigh] = cur_node;                                        
                pq.push({dist[neigh], neigh});
            }
        }
    }
    if(dist[dst] == INF) return{};
    int cur_node = dst;
    vector<int>path;
    while(cur_node != -1){
        //cout<<cur_node<<endl;
        path.push_back(cur_node);
        cur_node = parent[cur_node];
    }
    return path;
}     

void Greedy::path_assignment(){
    // base_test_active();
    vector<int>used_memory(graph.get_size(),0);
    map<pair<int,int>, int> used_channel;

    bool flag = true;
    while(flag){
        flag = false;
        for(int reqno = 0; reqno < (int)requests.size(); reqno++){
            vector<int>need_memory(graph.get_size(),0);
            map<pair<int,int>, int> need_channel;            
            Request &request = requests[reqno];
            double x1,x2,x3,xm,y1,y2,y3,ym;
            tie(x1,y1) = graph.Node_id2ptr(request.get_node1())->get_pos();
            tie(x2,y2) = graph.Node_id2ptr(request.get_node2())->get_pos();
            tie(x3,y3) = graph.Node_id2ptr(request.get_node3())->get_pos();
            double smallest_fermat = numeric_limits<double>::infinity();
            int smallest_middle = -1;
            for(int middle = 0; middle < graph.get_size(); middle ++){
                if(middle == request.get_node1() || middle == request.get_node2() || middle == request.get_node3() || (used_memory[middle]+3) > graph.Node_id2ptr(middle)->get_memory_cnt()){continue;}
                tie(xm,ym) = graph.Node_id2ptr(middle)->get_pos();
                double cur_fermat = pow((x1-xm),2) + pow((y1-ym),2) + pow((x2-xm),2) + pow((y2-ym),2) + pow((x3-xm),2) + pow((y3-ym),2);
                if( cur_fermat < smallest_fermat){
                    smallest_fermat = cur_fermat;
                    smallest_middle = middle;
                }
            }
            if(smallest_fermat == numeric_limits<double>::infinity()){continue;}    //單獨找只能確保一條路徑有足夠資源，不能保證整個tree有資源(但已經盡量讓它不會continue了)
            vector<int> path1 = Dijkstra(smallest_middle, request.get_node1());        
            vector<int> path2 = Dijkstra(smallest_middle, request.get_node2());
            vector<int> path3 = Dijkstra(smallest_middle, request.get_node3());
            if(path1.size() == 0 || path2.size() == 0 || path3.size() == 0 ){continue;}

            cal_need(path1,need_memory,need_channel); cal_need(path2,need_memory,need_channel); cal_need(path3,need_memory,need_channel);
            bool continue_flag = false;
            for(int i = 0; i < graph.get_size(); i++){
                if(need_memory[i] == 0){continue;}
                if( (need_memory[i] + used_memory[i]) > graph.Node_id2ptr(i)->get_memory_cnt()){
                    continue_flag = true;
                    break;
                }
            }
            if(continue_flag) continue;
            for(auto edge:need_channel){
                if(edge.first.first>edge.first.second){
                    if(used_channel.find(edge.first) != used_channel.end()){
                        if(used_channel[edge.first] + edge.second > graph.get_channel_size(edge.first.first,edge.first.second)){
                            continue_flag = true;
                            break;
                        }
                    }
                    else if(edge.second > graph.get_channel_size(edge.first.first,edge.first.second)){
                            continue_flag = true;
                            break;
                    }
                }
                else{
                    pair<int,int>rev_edge = {edge.first.second,edge.first.first};
                    if(used_channel.find(rev_edge) != used_channel.end()){
                        if(used_channel[rev_edge] + edge.second > graph.get_channel_size(edge.first.first,edge.first.second)){
                            continue_flag = true;
                            break;
                        }
                    }
                    else if(edge.second > graph.get_channel_size(edge.first.first,edge.first.second)){
                            continue_flag = true;
                            break;
                    }                    
                }
            }
            if(continue_flag) continue;
            vector<vector<int>> tree;
            tree.push_back(path1);tree.push_back(path2);tree.push_back(path3);
            assign_resource(tree, 1,reqno);

            for(int i = 0; i < graph.get_size(); i++){
                if(need_memory[i] == 0){continue;}
                used_memory[i] += need_memory[i];
            }

            for(auto edge:need_channel){
                if(edge.first.first > edge.first.second){
                    if(used_channel.find(edge.first) != used_channel.end()){
                        used_channel[edge.first] += edge.second;
                    }
                    else{
                        used_channel[edge.first] = edge.second;
                    }
                }
                else{
                    pair<int,int>rev_edge = {edge.first.second,edge.first.first};
                    if(used_channel.find(rev_edge) != used_channel.end()){
                        used_channel[rev_edge] += edge.second;
                    }
                    else{
                        used_channel[rev_edge] = edge.second;
                    }                    
                }
            }
            flag = true;
        }
    }
}

void Greedy::entangle(){
    AlgorithmBase::base_entangle();
}
void Greedy::swap(){
    AlgorithmBase::base_swap();
}
void Greedy::send(){
    AlgorithmBase::base_send();
}
void Greedy::next_time_slot(){
    AlgorithmBase::base_next_time_slot();
}

