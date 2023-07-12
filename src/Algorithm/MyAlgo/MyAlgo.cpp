#include "MyAlgo.h"

MyAlgo::MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha)
    :AlgorithmBase(filename, "MyAlgo", request_time_limit, node_time_limit, swap_prob, entangle_alpha){
    if(DEBUG) cerr<<"new MyAlgo"<<endl;
}

MyAlgo::~MyAlgo(){
    if(DEBUG) cerr << "delete MyAlgo" << endl;
}

void MyAlgo::entangle(){
    AlgorithmBase::base_entangle();
}

void MyAlgo::swap(){
     AlgorithmBase::base_swap();
}

void MyAlgo::send(){
     AlgorithmBase::base_send();
}

void MyAlgo::next_time_slot(){
     AlgorithmBase::base_next_time_slot();
}

void MyAlgo::initialize(){
    M = graph.get_size() + graph.get_num_of_edge() + requests.size(); //V+E+I
    delta = (1 + epsilon) * pow(((1 + epsilon) * M), (-1 / epsilon));
    cout<<"V E I="<< graph.get_size() <<" "<< graph.get_num_of_edge() <<" "<< requests.size()<<endl;
    cout<<"delta:"<<delta<<endl;
    for(int i = 0; i < graph.get_size(); i++){
        alpha.emplace_back(delta / graph.Node_id2ptr(i)->get_memory_cnt());       //alpha_set
        //cout<<"alpha id:"<<i<<" is"<<alpha[i]<<endl;
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
    ///test beta
    /*
    for(int i = 0; i < graph.get_size(); i++){
        vector<int> temp = graph.get_neighbors_id(i);                            
        for(auto it: temp){
            if(i < it){
                cout<<"beta for["<<i<<","<<it<<"]:"<<beta[{i,it}]<<endl;
            }
            else{
                cout<<"beta for["<<it<<","<<i<<"]:"<<beta[{it,i}]<<endl;
            }
        }
    }
    */

    for(int i = 0; i < requests.size(); i++){
        tau.emplace_back(delta / requests[i].get_send_limit());
    }
    Y.resize(requests.size() + 5);
    for(int i = 0; i < graph.get_size(); i++){
        vector<int> temp = graph.get_neighbors_id(i);                             
        for(auto it: temp){
            if(i < it){
                X[{i, it}] = alpha[i] + alpha[it] + beta[{i, it}];
            }
            else{
                X[{it, i}] = alpha[it] + alpha[i] + beta[{it, i}];
            }
            for(int j = 0;j < requests.size(); j++){
                int src = requests[j].get_source();
                int des = requests[j].get_destination();
                double ent_p = exp(graph.Node_id2ptr(i)->distance(*graph.Node_id2ptr(it))*(-graph.get_entangle_alpha()));
                if(i<it){
                    if(i != src && i != des && it != src && it != des){
                        Y[j][{i, it}] = -(log(ent_p)/log(exp(1))) - (log(sqrt(graph.Node_id2ptr(i)->get_swap_prob())/log(exp(1)))) - (log(sqrt(graph.Node_id2ptr(it)->get_swap_prob())/log(exp(1))));
                    }
                    else if((i == src && it != des) || (i == des && it != src)){
                        Y[j][{i, it}] = -(log(ent_p)/log(exp(1))) - (log(sqrt(graph.Node_id2ptr(it)->get_swap_prob())/log(exp(1))));
                    }
                    else if((i == src && it != des) || (i == des && it != src)){
                        Y[j][{i, it}] = -(log(ent_p)/log(exp(1))) - (log(sqrt(graph.Node_id2ptr(i)->get_swap_prob())/log(exp(1))));
                    }
                    else{
                        Y[j][{i, it}] = -(log(ent_p)/log(exp(1)));
                    }
                }
                else{
                    if(i != src && i != des && it != src && it != des){
                        Y[j][{it, i}] = -(log(ent_p)/log(exp(1))) - (log(sqrt(graph.Node_id2ptr(i)->get_swap_prob())/log(exp(1)))) - (log(sqrt(graph.Node_id2ptr(it)->get_swap_prob())/log(exp(1))));
                    }
                    else if((i == src && it != des) || (i == des && it != src)){
                        Y[j][{it, i}] = -(log(ent_p)/log(exp(1))) - (log(sqrt(graph.Node_id2ptr(it)->get_swap_prob())/log(exp(1))));
                    }
                    else if((i == src && it != des) || (i == des && it != src)){
                        Y[j][{it, i}] = -(log(ent_p)/log(exp(1))) - (log(sqrt(graph.Node_id2ptr(i)->get_swap_prob())/log(exp(1))));
                    }
                    else{
                        Y[j][{it, i}] = -(log(ent_p)/log(exp(1)));
                    }                  
                }

            }
                ////test X、Y
                if(i<it){
                    cout<<"edge["<<i<<","<<it<<"]: X:"<<X[{i,it}]<<"   Y:"<<Y[0][{i,it}]<<endl;
                }
                else{
                    cout<<"edge["<<it<<","<<i<<"]: X:"<<X[{it,i}]<<"   Y:"<<Y[0][{it,i}]<<endl;
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
    //     3---4
    //  /  |  /
    // /   | / 
    // 2---1/
    //  \  |
    //   \ |
    //     0

vector<int> MyAlgo::separation_oracle(int req_no, double &req_Us){  // running time 
    cout<<"sparation\n";
    vector<int> SPT(graph.get_size() + 5);
    vector<int> best_path(graph.get_size() + 5);
    double best_len; 
    int src = requests[req_no].get_source();
    int dst =  requests[req_no].get_destination();
    cout<<"Dijkstra\n";
    SPT = Dijkstra(src, dst);                               //the first SPT
    cout<<"Dijkstra end\n";
    
    int cur_node=src;                                       //counting the first path's U(x,y)
    double c;
    double r;
    while(cur_node != dst){

        if(cur_node < SPT[cur_node]){
            c += X[{cur_node,SPT[cur_node]}];
            r += Y[req_no][{cur_node,SPT[cur_node]}];
        }
        else{
            c += X[{SPT[cur_node],cur_node}];  
            r += Y[req_no][{SPT[cur_node],cur_node}];           
        }
        best_path.push_back(cur_node);
        cur_node = SPT[cur_node];
    } 
    best_path.push_back(dst);
    best_len=c/r;

    vector<int> new_path;
    double a;
    double b = 1;
    cout<<"midpoint1\n";
    while(b != 0){                                              //creating many new SPT
        map<pair<int, int>, bool> used;
        pair<int,int> new_edge;
        double minimum = numeric_limits<double>::infinity();
        for(int i = 0; i < graph.get_size(); i++){
            vector<int> neighbors = graph.get_neighbors_id(i);  
            for(auto neighbor : neighbors){
                double temp1=0, temp2=0;
                if(SPT[i] == neighbor || SPT[neighbor] == i){   // checking used edge or unused
                    continue;   
                }else{                                          // X[i, j] + X[p_j] - X[p_i] / Y[i, j] + Y[p_j] - Y[p_i]
                    temp1 = X[{i, neighbor}];
                    temp2 = Y[req_no][{i, neighbor}];
                    int prev_node;
                    int cur_node = i;
                    while(cur_node != dst){
                        temp1 += X[{cur_node, SPT[cur_node]}];
                        temp2 += Y[req_no][{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    } 
                    cur_node = neighbor;
                    while(cur_node != dst){
                        temp1 -= X[{cur_node, SPT[cur_node]}];
                        temp2 -= Y[req_no][{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    }       

                    if(temp2 < 0 && temp1 > 0){               // we need the smallest edge to change the SPT
                        if(minimum > - temp1 / temp2){
                            // revise edge
                            new_edge = {i, neighbor};
                            minimum = - temp1 / temp2;
                            a = minimum;
                            b = 1 / 1 + a;
                        }
                    }
                }
            }
        }        // 找到最小的 edge 
        cout<<"midpoint2\n";
        if(minimum == numeric_limits<double>::infinity()){
            break;
        }else{
            new_path.clear();
        }
        ////
        cur_node = new_edge.first;
        int cur_node2 = new_edge.second;
        while(cur_node != dst && cur_node2 != dst){
            cur_node=SPT[cur_node];
            cur_node2=SPT[cur_node2];
        }
        if(cur_node == dst){
            SPT[new_edge.second] = new_edge.first;
        }
        else{
            SPT[new_edge.first] = new_edge.second;
        }
        cout<<"midpoint3\n";
        cur_node = src;                                   //question:did you change the path? or just give me the old path?  
        while(cur_node != -1) {
            new_path.push_back(cur_node);
            cur_node = SPT[cur_node];
        }       

        //give me a path vector
        double new_len = 0;
        c = 0;
        r = 0;
        for(int i = 0; i < new_path.size() - 1; i++){
            if(new_path[i] < new_path[i+1]){
                c += X[{new_path[i], new_path[i+1]}];
                r += Y[req_no][{new_path[i], new_path[i+1]}];
            }
            else{
                c += X[{new_path[i+1], new_path[i]}];  
                r += Y[req_no][{new_path[i+1], new_path[i]}];           
            }
        }
        new_len = c / r;
        if(new_len < best_len){
            best_len = new_len;
            req_Us = best_len;
            best_path = new_path;                                            //路線修改,新的spt產生
        } 
    }

    for(auto p : best_path){
        cout << "Best path: " << p << endl;
    }
    return best_path;
                                                                       
}

void MyAlgo::find_bottleneck(vector<int> path, int req_no){
    
    double min_s_u = numeric_limits<double>::infinity();
    double min_s_uv = numeric_limits<double>::infinity();
    double s_i = requests[req_no].get_send_limit();
    vector<double> s_u(graph.get_size() + 5);
    vector<double> s_uv(graph.get_size() + 5);                                               

   
    for(int i = 0; i < path.size(); i++){
        if(i == 0 || path.size() - 1)
            s_u[path[i]] = graph.Node_id2ptr(path[i])->get_memory_cnt();// 是否要考慮 width
        else
            s_u[path[i]] = graph.Node_id2ptr(path[i])->get_memory_cnt() / 2;
        if(s_u[path[i]] < min_s_u)
            min_s_u = s_u[path[i]];
    }

    for(int i = 0; i < path.size() - 1; i++){
        s_uv[i] = graph.get_channel_size(path[i], path[i+1]);
        if(s_u[i] < min_s_uv)
            min_s_uv = s_uv[i];
    }
    
    double s = min(min_s_u, min(min_s_uv, s_i));
 
    x_i_p[path] += s;
    
    for(auto id : path){
        alpha[id] = alpha[id] * (1 + epsilon * s / s_u[id]);
    }

    for(int i = 0; i < path.size() - 1; i++){
        beta[{path[i], path[i+1]}] = beta[{path[i], path[i+1]}] * (1 + epsilon * s / s_uv[i]);
    }

    tau[req_no] = tau[req_no] * (1 + epsilon * s / s_i);    

}

double MyAlgo::changing_obj(){
    double temp_obj = 0.0;
    cout << "add alpha" << endl;
    for(int i = 0; i < alpha.size(); i++){
        temp_obj += alpha[i] * graph.Node_id2ptr(i)->get_memory_cnt();
    }
    
    cout << "add beta" << endl;
    for(auto it : beta){
        cout << "edge" << it.first.first << " " << it.first.second << endl;
        temp_obj += it.second * graph.get_channel_size(it.first.first, it.first.second);
    }

    cout << "add tau" << endl;
    for(int i = 0;i < requests.size(); i++){
        temp_obj += tau[i] * requests[i].get_send_limit();
    }
    return temp_obj;
}

void MyAlgo::path_assignment(){
    cout<< "run\n";
    initialize();
    double obj = M * delta;
    
    vector<int> best_path;
    vector<int> cur_path;
    double U;
    while(obj < 1){ // 是否在裡面做完 entangele, swap, send, separation_oracle 一個 SD 只找一條path?

        int req_no = 0;
        double smallest_U = numeric_limits<double>::infinity();
        for(int i = 0; i < requests.size(); i++){
            cur_path =  separation_oracle(i, U);
            if(U < smallest_U){
                smallest_U  = U;
                best_path = cur_path;
                req_no = i;
            }
            
        }
        
        // compare
        cout<<"find_bottle\n";
        find_bottleneck(best_path, req_no);
        cout<<"End find_bottle\n";
        obj = changing_obj();
        cout<<"changing_obj\n";
        //
    }
    
    for(auto x : x_i_p){
        vector<int> temp = x.first;
        cout << "PATH: ";
        for(auto it:temp){
            cout << it <<" ";
        }
        cout<<":";
        cout << x.second << endl;
    }
    // 

}   

