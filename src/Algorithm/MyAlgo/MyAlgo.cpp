#include "MyAlgo.h"

MyAlgo::MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha)
    :AlgorithmBase(filename, "MyAlgo", request_time_limit, node_time_limit, swap_prob, entangle_alpha , true /*let it be alway true*/){
    if(DEBUG) cerr<<"new MyAlgo"<<endl;
}
MyAlgo::MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha, double epsilon)
    :AlgorithmBase(filename, "MyAlgo_" + to_string(epsilon) , request_time_limit, node_time_limit, swap_prob, entangle_alpha , true /*let it be alway true*/), epsilon(epsilon){
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

double MyAlgo::X(int u, int v){
	double weight = alpha[u] + alpha[v] + beta[{u, v}];
	// if(requests[i].get_node1() == u || requests[i].get_node1() == v){            //[unsure] tau is bad!!
	// 	weight += tau[i];
	// }
	return weight;
}

void MyAlgo::init_dual(){
    M = graph.get_size() + graph.get_num_of_edge() + requests.size();              //M=V+E+I         
    delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * M, 1.0 / epsilon));         

    for(int i = 0; i < graph.get_size(); i++){
        alpha.emplace_back(delta / graph.Node_id2ptr(i)->get_memory_cnt());       //alpha=delta/Mu

        vector<int> temp = graph.get_neighbors_id(i);                             //beta=delta/Cuv
        for(auto it: temp){
            beta[make_pair(i,it)] = delta / (graph.get_channel_size(i, it));
        }
    }

    for(unsigned  i = 0; i < requests.size(); i++){                               //tau=delta/ri
        tau.emplace_back(delta / 1 );                                             //[unsure] r(i) = 1?
    }
}

void MyAlgo::initialize(int mid){
    for(unsigned  j = 0;j < requests.size(); j++){                       //Y=-ln(edge)-ln(repeater_1^(1/2))-ln(repeater_2^(1/2))
        int node1 = requests[j].get_node1();
        int node2 = requests[j].get_node2();
        int node3 = requests[j].get_node3();
        vector<int>three_node = {node1 , node2 , node3};
        for(int i = 0; i < graph.get_size(); i++){
        vector<int> temp = graph.get_neighbors_id(i);                             
            for(auto it: temp){
                double ent_p = exp(graph.Node_id2ptr(i)->distance(*graph.Node_id2ptr(it))*(-graph.get_entangle_alpha()));
                if( i == mid){                                                                                                                                      //if u=x,
                    if( find(three_node.begin() , three_node.end() ,it) != three_node.end() ){                                                                      //and if v={}                           
                        Y[mid][j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob()) / 3;              
                    }
                    else{                                                                                                                                           //else if v!={}       
                        Y[mid][j][{i, it}] = -log(ent_p) - (log(graph.Node_id2ptr(it)->get_swap_prob()) / 2) - (log(graph.Node_id2ptr(i)->get_swap_prob()) / 3);               
                    }
                }
                else if( it == mid ){                                                                                                                               //if v=x                                            
                    if( find(three_node.begin() , three_node.end() ,i) != three_node.end() ){                                                                       //and if u={}                          
                        Y[mid][j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 3;             
                    }
                    else{                                                                                                                                           //else if u!={}
                        Y[mid][j][{i, it}] = -log(ent_p) - (log(graph.Node_id2ptr(i)->get_swap_prob()) / 2) - (log(graph.Node_id2ptr(it)->get_swap_prob()) / 3);               
                    }
                }                
                else if( i != mid && find(three_node.begin() , three_node.end() ,it) != three_node.end()){                                                          //if u!=x and v={}                                           
                    Y[mid][j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob()) / 2;             
                } 
                else if( it != mid && find(three_node.begin() , three_node.end() ,i) != three_node.end()){                                                          //if v!=x and u={}                                                                                                                                                                                       
                    Y[mid][j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 2;             
                } 
                else{                                                                                                                                               //if u、v!=x and v、u!={}
                    Y[mid][j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 2 - log(graph.Node_id2ptr(i)->get_swap_prob()) / 2;             
                }
            }
        }
    }
}

vector<int> MyAlgo::Dijkstra(int src, int dst, int req, vector<pair<double,double>>&dist){                           
    double INF=numeric_limits<double>::infinity();
    vector<bool> used( graph.get_size(), false);
    vector<int> parent( graph.get_size(), -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    dist[src] = {0,0};                                                     //找src到所有點的最短距離(針對X-obj)
    pq.push({0, src});
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int neigh : graph.get_neighbors_id(cur_node)) {
            if(dist[cur_node].first + X(cur_node, neigh) < dist[neigh].first) {
                dist[neigh].first = dist[cur_node].first + X(cur_node, neigh);      //(1)store d
                dist[neigh].second = dist[cur_node].second + Y[dst][req][{cur_node, neigh}];
                parent[neigh] = cur_node;                                             //(1)store pred
                pq.push({dist[neigh].first, neigh});
            }
            if(dist[cur_node].first + X(cur_node, neigh) == dist[neigh].first && dist[neigh].second > (dist[cur_node].second + Y[dst][req][{cur_node, neigh}])) {
                dist[neigh].first = dist[cur_node].first + X(cur_node, neigh);      //(1)store d
                dist[neigh].second = dist[cur_node].second + Y[dst][req][{cur_node, neigh}];
                parent[neigh] = cur_node;                                             //(1)store pred
                pq.push({dist[neigh].first, neigh}); 
            }
        }
    }

    if(dist[dst].first >= INF) return{};
    return parent;
}       
   

void MyAlgo::separation_oracle(int src, int dst, int req_no, int path_id, vector<vector<vector<int>>> &cur_tree, vector<vector<vector<double>>> &cur_label){     
    double INF=numeric_limits<double>::infinity();
    vector<pair<double,double>>dist (graph.get_size(),{INF,INF});            
    vector<int> pred (graph.get_size(),-1);                                  
    vector<int> SDpath;
    vector<int> best_path;
    double smallest_U = INF;
    vector<double> theta_table(graph.get_size(),INF);                        
    vector<double> obj_bar1_table(graph.get_size(),INF);                     
    vector<double> obj_bar2_table(graph.get_size(),INF);
    vector<int> cpred_table(graph.get_size(),-1);                            
    priority_queue<Label,vector<Label>, greater<Label>> pq;                                          
    pred = Dijkstra(src, dst, req_no, dist);                                    //(1)

    double last_ratio = 0;                                              //(2)start
    double d1 = dist[dst].first,d2 = dist[dst].second;                  
    int current = dst;
    while(current != -1){
        SDpath.push_back(current);
        current = pred[current];
    }
    reverse(SDpath.begin(),SDpath.end());                               //(2)end

    if(d1 * exp(d2) < smallest_U){
        smallest_U = d1 * exp(d2);
        best_path = SDpath;
    }
    // cout<<"Begin path:";
    // for(auto it: SDpath){
    //     cout<<it<<" ";
    // }
    // cout<<endl;
    // cout<<"U = "<< d1 * exp(d2)<<endl;
    for(int i = 0; i < graph.get_size(); i++){                               //(4)compute lexmin
        if(i == src){continue;}
        for(auto &neigh:graph.get_neighbors_id(i)){                               
            double cur_obj_bar2 = dist[neigh].second + Y[dst][req_no][{neigh,i}] - dist[i].second;
            if(cur_obj_bar2 < 0){
                double cur_theta = -(dist[neigh].first + X(neigh,i) - dist[i].first) / cur_obj_bar2 ;
                if(cur_theta <= theta_table[i]){
                    if(cur_theta < theta_table[i]){
                        theta_table[i] = cur_theta;
                        obj_bar1_table[i] = dist[neigh].first + X(neigh,i) - dist[i].first;
                        obj_bar2_table[i] = cur_obj_bar2;
                        cpred_table[i] = neigh;
                    }
                    else if(cur_obj_bar2 < obj_bar2_table[i]){
                        theta_table[i] = cur_theta;
                        obj_bar1_table[i] = dist[neigh].first + X(neigh,i) - dist[i].first;
                        obj_bar2_table[i] = cur_obj_bar2;
                        cpred_table[i] = neigh;
                    }
                }
            }
        }
        if(cpred_table[i] != -1){
            pq.push({theta_table[i],obj_bar2_table[i], i});
            //cout<<"Push1: "<<theta_table[i]<<" "<<obj_bar2_table[i]<<" "<<i<<"\n";
        }
    }

    while(!pq.empty()){                                         //(6)
        double min_theta,min_obj_bar2;int min_i;
        Label temp = pq.top();                                   //(7)
        min_theta = temp.theta,min_obj_bar2=temp.obj_bar2,min_i=temp.i;
        pq.pop();                                                          //會有INF是這裡???
        if(min_obj_bar2 != obj_bar2_table[min_i]){continue;}
        //cout<<"pop: "<<min_theta<<" "<<min_obj_bar2<<" "<<min_i<<endl;
        dist[min_i].first += obj_bar1_table[min_i];
        dist[min_i].second += obj_bar2_table[min_i];
        pred[min_i] = cpred_table[min_i];
        if(min_i == dst){                                       //(9)
            if(min_theta > last_ratio){                         //(10)
                //cout << "[" << last_ratio << "," << min_theta << "] with " << d1 * exp(d2) << endl;
                vector<double>temp = {last_ratio,min_theta, (d1 * exp(d2))};
                last_ratio = min_theta;
                cur_tree[path_id].push_back(SDpath);
                cur_label[path_id].push_back(temp);
                // cout<<"Extremepath: ";
                // for(auto it:SDpath){
                //     cout << it <<" ";
                // }
                // cout<<endl<<endl;
            }
            SDpath.clear();
            current = dst;
            while(current != -1){
                SDpath.push_back(current);
                current = pred[current];
            }
            reverse(SDpath.begin(),SDpath.end());              //(12)
            d1 = dist[dst].first,d2 = dist[dst].second;

            if(d1 * exp(d2) < smallest_U){
                smallest_U = d1 * exp(d2);
                best_path = SDpath;
            }
        }

        theta_table[min_i] = INF,cpred_table[min_i] = -1,obj_bar1_table[min_i] = INF,obj_bar2_table[min_i] = INF;       //(13)
        for(auto neigh:graph.get_neighbors_id(min_i)){                             
            double cur_obj_bar2 = dist[neigh].second + Y[dst][req_no][{neigh,min_i}] - dist[min_i].second;
            if(cur_obj_bar2 < 0){
                double cur_theta = -(dist[neigh].first + X(neigh,min_i) - dist[min_i].first) / cur_obj_bar2 ;
                if(cur_theta <= theta_table[min_i]){
                    if(cur_theta < theta_table[min_i]){
                        theta_table[min_i] = cur_theta;
                        obj_bar1_table[min_i] = dist[neigh].first + X(neigh,min_i) - dist[min_i].first;
                        obj_bar2_table[min_i] = cur_obj_bar2;
                        cpred_table[min_i] = neigh;
                    }
                    else if(cur_obj_bar2 < obj_bar2_table[min_i]){
                        theta_table[min_i] = cur_theta;                                                                //no use just for read
                        obj_bar1_table[min_i] = dist[neigh].first + X(neigh,min_i) - dist[min_i].first;
                        obj_bar2_table[min_i] = cur_obj_bar2;
                        cpred_table[min_i] = neigh;
                    }
                }
            }
        }
        if(cpred_table[min_i] != -1){
            pq.push({theta_table[min_i],obj_bar2_table[min_i],min_i});
            //cout<<"Push2: "<<theta_table[min_i]<<" "<<obj_bar2_table[min_i]<<" "<<min_i<<"\n";
        }
        for(auto neigh:graph.get_neighbors_id(min_i)){
            double temp_obj_bar2 = dist[min_i].second + Y[dst][req_no][{min_i,neigh}] - dist[neigh].second;
            if(temp_obj_bar2 < 0){
                double temp_obj_bar1 = dist[min_i].first + X(min_i,neigh) - dist[neigh].first;
                if( -(temp_obj_bar1/temp_obj_bar2) < theta_table[neigh] || (-(temp_obj_bar1/temp_obj_bar2) == theta_table[neigh] && temp_obj_bar2 < obj_bar2_table[neigh])){
                    pq.push({-temp_obj_bar1/temp_obj_bar2,temp_obj_bar2,neigh});
                    theta_table[neigh] = -(temp_obj_bar1/temp_obj_bar2);
                    obj_bar1_table[neigh] = temp_obj_bar1;
                    obj_bar2_table[neigh] = temp_obj_bar2;
                    cpred_table[neigh] = min_i;
                    //cout<<"Push3: "<<theta_table[neigh]<<" "<<obj_bar2_table[neigh]<<" "<<neigh<<"\n";       
                }    
            }
        }
    }
    // cout<<"Last path :";
    // for(auto it:SDpath){
    //     cout<<it<<" ";
    // }
    vector<double>temp = {last_ratio,INF,(d1 * exp(d2))};
    cur_tree[path_id].push_back(SDpath);
    cur_label[path_id].push_back(temp);        
    //cout<<endl<<"U : "<< d1 * exp(d2) <<" and ratio: " <<last_ratio<<endl;
}

void MyAlgo::find_bottleneck(vector<vector<int>> &tree, int req_no){
    double min_s_u = numeric_limits<double>::infinity();
    double min_s_uv = numeric_limits<double>::infinity();
    //double s_i = 1;                                                         //send_limit = 1???
    vector<double> s_u(graph.get_size() + 5);                                              
    vector<int> used_u(graph.get_size() + 5, 0);
    map<pair<int,int>,int> used_uv;
    map<pair<int,int>,double> s_uv;
    // Calculate used 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < tree[i].size() - 1; j++){
            if(j == 0){
                used_u[tree[i][j]]++;
            }
            else{
                used_u[tree[i][j]]+=2;
            }
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
        if(used_u[i] == 0){continue;}
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
    double s = min(min_s_u, min_s_uv);                                      //request limit=1?
    for(int i = 0; i < rate; i++){
        bool find = false;
        for(int j = 0; j < x_i_t_tree.size(); j++){                                                   //add flow to path
            if(x_i_t_tree[j][0] == tree[0] && x_i_t_tree[j][1] == tree[1] && x_i_t_tree[j][2] == tree[2]){
                x_i_t[j] += s;
                find = true;
                break;
            }
        }                                        
        if(find == false){
            x_i_t.push_back(s);
            x_i_t_tree.push_back(tree);
        }    

        for(int i = 0; i < used_u.size(); i++){
            if(used_u[i] == 0){continue;}
            obj += (alpha[i] * (1 + epsilon * s / s_u[i]) - alpha[i]) * graph.Node_id2ptr(i)->get_memory_cnt();;            //alter alpha,beta,tau
            alpha[i] = alpha[i] * (1 + epsilon * s / s_u[i]);
        }

        for(auto &it : used_uv){
            obj += (beta[{it.first.first, it.first.second}] * (1 + epsilon * s / s_uv[{it.first.first, it.first.second}]) -  beta[{it.first.first, it.first.second}]) * graph.get_channel_size(it.first.first, it.first.second);;
            beta[{it.first.first, it.first.second}] = beta[{it.first.first, it.first.second}] * (1 + epsilon * s / s_uv[{it.first.first, it.first.second}]);
            beta[{it.first.second, it.first.first}] = beta[{it.first.second, it.first.first}] * (1 + epsilon * s / s_uv[{it.first.second, it.first.first}]);
        }
        obj += (tau[req_no] * (1 + epsilon * s / 1) - tau[req_no]) * 1;                            //request limit = 1?
        tau[req_no] = tau[req_no] * (1 + epsilon * s / 1);    
    }
}

void MyAlgo::find_violate(){
    vector<double> used_memory(graph.get_size());
    map<vector<int>, double> used_channel;
    map<tuple<int, int, int>, double> used_request;

    for(int it = 0; it < x_i_t_tree.size() ; it++){
        vector<vector<int>> tree = x_i_t_tree[it];
        int node1 = tree[0][0];
        int node2 = tree[1][0];
        int node3 = tree[2][0];
        if(used_request.find({node1, node2, node3}) != used_request.end())
            used_request[{node1, node2, node3}] += x_i_t[it];
        else
            used_request[{node1, node2, node3}] = x_i_t[it];
        
        for(int path_id = 0; path_id < 3; path_id ++){
            for(unsigned int i = 0; i < x_i_t_tree[it][path_id].size() - 1; i++){
                used_memory[x_i_t_tree[it][path_id][i]] += x_i_t[it];                         //memory add
                used_memory[x_i_t_tree[it][path_id][i+1]] += x_i_t[it];
                if(x_i_t_tree[it][path_id][i] < x_i_t_tree[it][path_id][i+1]){
                    auto iter = used_channel.find({x_i_t_tree[it][path_id][i], x_i_t_tree[it][path_id][i+1]});
                    if(iter != used_channel.end()){    //channel add
                        used_channel[{x_i_t_tree[it][path_id][i], x_i_t_tree[it][path_id][i+1]}] += x_i_t[it];
                    }
                    else{
                        used_channel[{x_i_t_tree[it][path_id][i], x_i_t_tree[it][path_id][i+1]}] = x_i_t[it];
                    }
                }
                else{
                    auto iter = used_channel.find({x_i_t_tree[it][path_id][i+1], x_i_t_tree[it][path_id][i]});
                    if(iter != used_channel.end()){
                        used_channel[{x_i_t_tree[it][path_id][i+1], x_i_t_tree[it][path_id][i]}] += x_i_t[it];
                    }
                    else{
                        used_channel[{x_i_t_tree[it][path_id][i+1], x_i_t_tree[it][path_id][i]}] = x_i_t[it];
                    }  
                }
            }
        }
    }


    double max_magni = 0.0;
    double cur_magni;

    for(auto &it : used_request){
        cur_magni = it.second / 1;
        if(cur_magni > max_magni){
            max_magni = cur_magni;
        }
    }

    for(auto &it : used_channel){
        cur_magni = it.second / graph.get_channel_size(it.first[0],it.first[1]);
        if(cur_magni > max_magni){
            max_magni = cur_magni;
        }
    }

    for(int i = 0; i < graph.get_size(); i++){
        cur_magni = used_memory[i] / graph.Node_id2ptr(i)->get_memory_cnt();
        if(cur_magni > max_magni){
            max_magni = cur_magni;
        }
    }


    

    cout << "Magnification:" << max_magni << endl;

    for(auto &x : x_i_t){
        x /= max_magni;
    }
}

vector<map<vector<int>, int>> MyAlgo::rounding(){
    // vector<map<vector<int>, double>> each_request(requests.size());
    // vector<map<vector<int>, int>> I_request(requests.size());
    // for(auto it : x_i_p){
    //     vector<int> path = it.first;
    //     int src = path[0];
    //     int dst = path.back();
    //     for(unsigned int i = 0; i < requests.size(); i++){
    //         if(src == requests[i].get_source() && dst == requests[i].get_destination()){
    //             each_request[i][path] = it.second;
    //             break;
    //         }
    //     }
    // }

    // for(unsigned int i = 0; i < requests.size(); i++){
    //     double total_prob = 0;
    //     double used_prob = 0;
    //     int used_I=0;
    //     int distri_I=0;
    //     vector<double>accumulate;
    //     accumulate.push_back(0.0);                                              // [0,a1,a2,...,0]
    //     for(auto it : each_request[i]){                    
    //         double frac_prob;

    //         int i_prob = it.second;                                             //每個path先取整數部分=>確定分配
    //         I_request[i][it.first] = i_prob;
    //         used_I+=i_prob;

    //         frac_prob = it.second - i_prob;                                     //total_prob代表random區間,丟進accumulate
    //         total_prob += frac_prob;
    //         accumulate.push_back(total_prob);
    //         used_prob += it.second;
    //     }
    //     used_I += (int)(requests[i].get_send_limit()- used_prob);               //unused_I=取底[ri - sum(request.I) - (unused.I)]
    //     distri_I=requests[i].get_send_limit()-used_I;

    //     accumulate.push_back(0.0);
    //     // cout<<"total_prob:"<<total_prob<<" distri_I:"<<distri_I<<endl;
    //     // cout<<"accumulate:";
    //     // for(auto it:accumulate){
    //     //     cout<<it<<" ";
    //     // }
    //     // cout<<endl;
    //     for(int j = 0; j < distri_I; j++){
    //         random_device rd;  
    //         mt19937 gen(rd()); 
    //         uniform_real_distribution<double> dis(0.0, total_prob);
    //         double temp = dis(gen);
    //         for(unsigned int k = 0; k < accumulate.size() - 1; k++){
    //             // cout<<"ramdom:"<<temp<<endl;
    //             if(temp > accumulate[k] && temp < accumulate[k+1]){
    //                 unsigned int index = 0;
    //                 for(auto it : each_request[i]){
    //                     if(index == k){
    //                         I_request[i][it.first]++;
    //                         //cout<<"random prob="<<temp<<" so "<<k<<" added!\n";
    //                         break;
    //                     }
    //                     index += 1;
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }

    // return I_request;
}

void MyAlgo::check_enough(vector<map<vector<int>, int>> &path){
    // vector<int> memory_used(graph.get_size());
    // map<vector<int>,int> channel_used; 
    // vector<int> over_memory(graph.get_size());
    // map<vector<int>,int> over_channel;
    // map<vector<int>,int>::iterator iter;
    // for(int i = 0; i <(int)path.size(); i++){
    //     for(auto it : path[i]){
    //         vector<int> cur_path = it.first;
    //         for(int j = 0; j < (int)cur_path.size() - 1; j++){
    //             memory_used[cur_path[j]] += it.second;
    //             memory_used[cur_path[j+1]] += it.second;
    //             iter = channel_used.find({cur_path[j],cur_path[j+1]});
    //             if(iter != channel_used.end()){
    //                 channel_used[{cur_path[j], cur_path[j+1]}] += it.second;
    //                 channel_used[{cur_path[j+1], cur_path[j]}] += it.second;
    //             }
    //             else{
    //                 channel_used[{cur_path[j], cur_path[j+1]}] = it.second;
    //                 channel_used[{cur_path[j+1], cur_path[j]}] = it.second;
    //             }
    //         }
    //     }
    // }

    // for(int i = 0; i < graph.get_size(); i++){
    //     over_memory[i] = memory_used[i] - graph.Node_id2ptr(i)->get_memory_cnt();
    //     for(auto it : graph.get_neighbors_id(i)){
    //         iter = over_channel.find({i, it});
    //         if(iter != over_channel.end()){
    //            over_channel[{i, it}] -= graph.get_channel_size(i, it) / 2;
    //            over_channel[{it, i}] -= graph.get_channel_size(i, it) / 2; 
    //         }
    //         else{
    //            over_channel[{i, it}] = -graph.get_channel_size(i, it) / 2;
    //            over_channel[{it, i}] = -graph.get_channel_size(i, it) / 2;
    //         }
    //     }
    // }

    // for(auto &it : over_channel){
    //     iter = channel_used.find(it.first);
    //     if(iter != channel_used.end()){
    //         it.second = channel_used[{it.first}] + it.second ; 
    //     }
    // }

    // bool flag;
    // while(1){
    //     flag = true;
    //     for(int i = 0; i < (int)over_memory.size(); i++){
    //         if(over_memory[i] > 0){
    //             flag = false;
    //         }
    //     }
    //     for(auto it : over_channel){
    //         if(it.second > 0){
    //             flag = false;
    //         }
    //     }
    //     if(flag == true){
    //         // cout << "before" << endl;
    //         // for(unsigned int i = 0; i < path.size(); i++){
    //         //     for(auto it:path[i]){
    //         //         vector<int>Final_path =it.first;
    //         //         for(auto it2:Final_path){
    //         //             cout<<it2<<" ";
    //         //         }
    //         //         cout<<"     Qubits:"<<it.second<<endl;
    //         //         requests[i].add_cur(it.second);
    //         //     }
    //         // }
    //         cout<<"--------------Reduece finish-------------\n";
    //         readd(path,over_memory,over_channel);  
    //         break;
    //     }
    //     int long_len = 0;
    //     int long_req = -1;
    //     vector<int> long_path;
    //     for(int i = 0; i < (int)path.size(); i++){
    //         for(auto it : path[i]){
    //             int associate_flag=false;
    //             /*
    //             for(auto temp:it.first){
    //                 cout<<temp<<" ";
    //             }
    //             cout<<"-----------"<<endl;
    //             */
    //             for(int j=0;j<(int)it.first.size()-1;j++){

    //                 //cout<<"memory check:"<<j<<"||"<<over_memory[it.first[j]]<<endl;
    //                 if(over_memory[it.first[j]]>0){
    //                     associate_flag=true;
    //                     break;
    //                 }
    //                 //cout<<"channel check:"<<j<<"/"<<j+1<<"||"<<over_channel[{it.first[j],it.first[j+1]}]<<endl;
    //                 iter = over_channel.find({it.first[j],it.first[j+1]});
    //                 if(iter!=over_channel.end() && over_channel[{it.first[j],it.first[j+1]}]>0){
    //                     associate_flag=true;
    //                     break;
    //                 }

    //             }
    //             if(over_memory[it.first[it.first.size()-1]]>0){
    //                 associate_flag=true;
    //             }

    //             if(associate_flag==true && (int)it.first.size() > long_len && it.second > 0){
    //                 long_len = it.first.size();
    //                 long_path = it.first;
    //                 long_req = i;
    //             }
    //         }
    //     }
    //     for(int i = 0; i < (int)long_path.size() - 1; i++){
    //         over_memory[long_path[i]]--;
    //         over_memory[long_path[i+1]]--;
    //         over_channel[{long_path[i], long_path[i+1]}]--;
    //         over_channel[{long_path[i+1], long_path[i]}]--;
    //     }
    //     path[long_req][long_path]--;
    // }  
}  

void MyAlgo::readd(vector<map<vector<int>, int>> &path,vector<int> &over_memory,map<vector<int>,int> &over_channel){
    // for(unsigned int i = 0; i < path.size(); i++){
    //     for(auto it : path[i]){
    //         requests[i].add_cur(it.second);
    //     }
    // }
    // vector<pair<vector<int>, int>> re;
    // int max = -1;
    // for(unsigned int i = 0; i < requests.size(); i++){
    //     for(auto it : path[i]){
    //         if(max < it.second){
    //             max = it.second;
    //         }
    //     }
    // }

    // for(int i = max; i >= 0; i--){
    //     for(unsigned int j = 0; j < requests.size(); j++){
    //         for(auto it : path[j]){
    //             if(i == it.second){
    //                 re.push_back({it.first, j});
    //             }
    //         }
    //     }
    // }
    // bool flag = true;
    // while(flag){
    //     flag = false;
    //     for(unsigned int i = 0; i < re.size(); i++){
    //         if(requests[re[i].second].get_send_limit() > requests[re[i].second].get_cur_send()){
    //             vector<int> each_path = re[i].first;
    //             bool assign = true;
    //             for(unsigned int j = 0; j < each_path.size() - 1; j++){
    //                 if(j == 0){
    //                     if(over_memory[each_path[j]] >= 0){
    //                         assign = false;
    //                     }
    //                 }
    //                 else{
    //                     if(over_memory[each_path[j]] >= -1){
    //                         assign = false;
    //                     }
    //                 }
    //                 if(over_channel[{each_path[j],each_path[j+1]}] >= 0){
    //                     assign = false;
    //                 }
    //             }
    //             if(over_memory[each_path[each_path.size()-1]] >= 0){
    //                 assign = false;
    //             }
    //             if(assign == true ){
    //                 requests[re[i].second].add_cur(1);
    //                 for(auto it : path[re[i].second]){
    //                     if(it.first == re[i].first){
    //                         path[re[i].second][it.first] += 1;
    //                         cout << "!!PATH +++" << endl;
    //                         flag = true;
    //                         for(int j = 0; j < (int)each_path.size() - 1; j++){
    //                             over_memory[each_path[j]]++;
    //                             over_memory[each_path[j+1]]++;
    //                             over_channel[{each_path[j], each_path[j+1]}]++;
    //                             over_channel[{each_path[j+1], each_path[j]}]++;
    //                         }
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // // for(auto x : over_memory){
    // //     cout << "node: " << x << endl;
    // // }
    // // for(auto x : over_channel){
    // //     cout<< "EDGE: ";
    // //     for(auto a : x.first ){
    // //         cout<< a << " ";
    // //     }
    // //     cout << x.second << endl;
    // // }

}

void MyAlgo::dfs(int src, int dst, vector<vector<int>> &ans, vector<int> &path, vector<bool> &visited){
        //base case
    visited[src] = true;
    path.push_back(src);
    if(src == dst){
        ans.push_back(path);
    } 
    else{
        for(auto i : graph.get_neighbors_id(src)){ 
            if(!visited[i]){
                dfs(i, dst, ans, path, visited);
            }
        }
    }
    visited[src] = false;
    path.pop_back();

}


void MyAlgo::calculate(){
    // double sum=0.0;

    // int t = 1;

    // for(auto it:x_i_p){
    //     double prob=1;
    //     vector<int>path=it.first;
    //     for(unsigned int i=0;i < it.first.size() - 1;i++){
    //         prob*=exp(graph.Node_id2ptr(path[i])->distance(*graph.Node_id2ptr(path[i+1]))*(-graph.get_entangle_alpha()));
    //     }
    //     for(unsigned int i=1;i<it.first.size()-1;i++){
    //         prob*=graph.Node_id2ptr(path[i])->get_swap_prob();  
    //     }
    //     sum+=it.second*prob;

    //     t++;
    // }
    // cerr << "sum = " << sum << endl;
    // res["primal"] = sum / (1 - epsilon) / (1 - epsilon);
}

vector<map<vector<vector<int>>, double>> MyAlgo::Greedy_rounding(){
    vector<map<vector<vector<int>>, double>> each_request(requests.size());
    vector<map<vector<vector<int>>, int>> I_request(requests.size());
    for(int i = 0; i < x_i_t_tree.size(); i++){
        vector<vector<int>> tree = x_i_t_tree[i];
        int node1 = tree[0][0];
        int node2 = tree[1][0];
        int node3 = tree[2][0];
        for(unsigned int j = 0; j < requests.size(); j++){
            if(node1 == requests[j].get_node1() && node2 == requests[j].get_node2() && node3 == requests[j].get_node3()){
                each_request[j][tree] = x_i_t[i];
                break;
            }
        }
    }
    for(auto &it:each_request){
        cout<<"Request ------"<<endl;
        for(auto &it2:it){
            cout<<"Value : "<<it2.second<<endl;;
            for(auto &it3:it2.first){
                for(auto &it4:it3){
                    cout<<it4<<" ";
                }
                cout<<endl;
            }
        }
    }
/*
    for(unsigned int i = 0; i < each_request.size(); i++){
        for(auto it:each_request[i]){
            vector<int>undistri_path =it.first;
        }
    }
	vector<int> used_I(requests.size());										//第 i 個 request 目前用了幾調 path
	vector< tuple<double, int, vector<int>> > fractional_xip;	
    for(unsigned int i = 0; i < requests.size(); i++){
        used_I[i] = 0;
		for(auto it : each_request[i]){                    
            double frac_prob;

            int i_prob = it.second;                                             //每個path先取整數部分=>確定分配
            I_request[i][it.first] = i_prob;
            used_I[i] += i_prob;
			assign_resource(it.first, i_prob, i);
            frac_prob = it.second - i_prob;                                     //total_prob代表random區間,丟進accumulate
			fractional_xip.emplace_back(frac_prob*graph.find_success_probability(it.first), i, it.first);
        }                                                   //unused_I=取底[ri - sum(request.I) - (unused.I)]
    }

	// 對 x^i_p 由小到大排序
	sort(fractional_xip.begin(), fractional_xip.end());
	reverse(fractional_xip.begin(), fractional_xip.end());

	// 如果資源足夠讓 x 變成 1 ，就直接讓 x 變成 1 
	for(auto it:fractional_xip){
		vector<int> extra_path;
		double x_prob;
		int request_id;
		tie(x_prob, request_id, extra_path) = it;
		if(find_width(extra_path) >= 1 && used_I[request_id] < requests[request_id].get_send_limit()){
			assign_resource(extra_path, 1, request_id);
			used_I[request_id] += 1;
            I_request[request_id][extra_path]++;
		}
	}
	// 如果還有剩下資源的話，盡量塞爆
	for(auto it:fractional_xip){
		vector<int> extra_path;
		double x_prob;
		int request_id;
		tie(x_prob, request_id, extra_path) = it;
		int width = 0;
		int extra_send_limit = requests[request_id].get_send_limit() - used_I[request_id];
		width = min(find_width(extra_path), extra_send_limit);
		if(width >= 1){
			assign_resource(extra_path, width, request_id);
            used_I[request_id] += width;
			I_request[request_id][extra_path]++;
		}
	}
	for(int request_id=0;request_id<(int)requests.size();request_id++){
		// cerr<<"GG: It work?"<<endl;
		while(requests[request_id].get_send_limit() - used_I[request_id] > 0){
			vector<int> extra_path = BFS(requests[request_id].get_source(), requests[request_id].get_destination());
			int width = 0;
			if(extra_path.size() != 0){
				width = min(find_width(extra_path), requests[request_id].get_send_limit() - used_I[request_id]);
				assign_resource(extra_path, width, request_id);
				used_I[request_id] += width;
			}
			if(width == 0){
				break;
			}
		}
	}
*/
    return each_request;

}

void MyAlgo::path_assignment(){
    double INF=numeric_limits<double>::infinity();
    Y.resize(graph.get_size());
    for(int i = 0; i < graph.get_size(); i++){
        int mem = graph.Node_id2ptr(i)->get_memory_cnt();
        graph.Node_id2ptr(i)->revise(mem);
        Y[i].resize(requests.size() + 1);
    }
    init_dual();
    for(int i = 0; i < graph.get_size(); i++){
        initialize(i);
    }
    // for(auto it1:Y){
    //     cout<<"------------------------------"<<endl;
    //     for(auto it2:it1){
    //         cout<<"next request"<<endl;
    //         for(auto it3:it2){
    //             cout<<it3.first.first<<"-"<<it3.first.second<<":"<<it3.second<<endl;
    //         }
            
    //     }
    // }

    // for(unsigned int i = 0; i < requests.size(); i++){
    //     int src = requests[i].get_source();
    //     int dst = requests[i].get_destination();
    //     all_source_target_path.push_back(allPathsSourceTarget(src, dst));
    // }

    obj = M * delta;
    while(obj < 1){
        int req_no = 0;
        vector<vector<int>> best_tree(3,vector<int>());
        double global_U = numeric_limits<double>::infinity();
        for(int middle = 0; middle < graph.get_size(); middle++){
            double smallest_U = numeric_limits<double>::infinity();
            int smallest_index = -1;
            vector<vector<vector<int>>> cur_tree(3,vector<vector<int>>());                //[path_id][extreme_tree_id][eles]
            vector<vector<vector<double>>> cur_label(3,vector<vector<double>>());         //[path_id][extreme_tree_id][last_ratio, new_ratio, U]
            for(unsigned int i = 0; i < requests.size(); i++){
                if(middle == requests[i].get_node1() || middle == requests[i].get_node2() || middle == requests[i].get_node3()){continue;}
                separation_oracle(requests[i].get_node1(), middle, i, 0, cur_tree, cur_label);
                separation_oracle(requests[i].get_node2(), middle, i, 1, cur_tree, cur_label);
                separation_oracle(requests[i].get_node3(), middle, i, 2, cur_tree, cur_label);
                // for(int fir = 0; fir < cur_tree.size(); fir++){
                //     for(int sec = 0; sec < cur_tree[fir].size(); sec++){
                //         cout<<"["<<cur_label[fir][sec][0]<<" ~ "<<cur_label[fir][sec][1]<<"] : "<<cur_label[fir][sec][2]<<endl;
                //         for(auto it:cur_tree[fir][sec]){
                //             cout<<it<<" ";
                //         }
                //         cout<<"\n-----------------------"<<endl;
                //     }
                //     cout<<"\nNEXT\n";
                // }
                vector<double>label_area;
                label_area.push_back(0);
                label_area.push_back(INF);
                for(int fir = 0; fir < cur_label.size(); fir++){
                    for(int sec = 0; sec < cur_label[fir].size()-1; sec++){
                        label_area.push_back(cur_label[fir][sec][1]);
                    }
                }
                sort(label_area.begin(),label_area.end());
                vector<double>total_U(label_area.size()-1,0);
                for(int fir = 0; fir < cur_label.size(); fir++){
                    for(int sec = 0; sec < cur_label[fir].size(); sec++){
                        for(int third = 1; third < label_area.size(); third++){
                            if(cur_label[fir][sec][0] <= label_area[third-1] && cur_label[fir][sec][1] >= label_area[third]){
                                total_U[third-1] += cur_label[fir][sec][2];
                            }
                        }
                    }
                }
                // for(auto it:label_area){
                //     cout<<it<<" ";
                // }
                // cout<<endl;
                // for(auto it:total_U){
                //     cout<<it<<" ";
                // }
                // cout<<endl;              
                for(int j = 0; j < total_U.size(); j++){
                    if(smallest_U > total_U[j]){
                        smallest_index = j;
                        smallest_U = total_U[j];
                    }
                }
                if(smallest_U < global_U){
                    global_U = smallest_U;
                    for(int fir = 0; fir < cur_label.size(); fir++){
                        for(int sec = 0; sec < cur_label[fir].size(); sec++){
                            if(cur_label[fir][sec][0] <= label_area[smallest_index] &&  label_area[smallest_index+1] <= cur_label[fir][sec][1]){
                                best_tree[fir] = cur_tree[fir][sec];
                                break;
                            }
                        }
                    } 
                    // cout<<"Request "<<i<<" with : "<<smallest_U <<" \n";
                    // for(auto it:best_tree){
                    //     for(auto it1:it){
                    //         cout<<it1<<" ";
                    //     }
                    //     cout<<endl;
                    // }
                    // cout<<"--------------------------"<<endl;
                }

            }
        }
        // cout<<"THIS ROUND ANS: "<<global_U<<endl;
        // for(auto it:best_tree){
        //     for(auto it2:it){
        //         cout<<it2<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        find_bottleneck(best_tree, req_no);
        //cout <<"OBJ : "<< obj << endl;
    }
    find_violate();
    for(int i = 0; i < x_i_t_tree.size(); i++){
        cout<<"\nTree with "<<x_i_t[i]<<endl;
        for(auto &it2:x_i_t_tree[i]){
            for(auto &it3:it2){
                cout<<it3<<" ";
            }
            cout<<endl;
        }
    }
    vector<map<vector<vector<int>>, double>>path = Greedy_rounding();
}
    
//     calculate();

//     res["change_edge_num"] = change_edge_num;
//     res["diff_edge_num"] = diff_num;
// }   


// //                   Our Old Friend                                //
// vector<int> MyAlgo::Dijkstra_ori(int src, int dst, int req_no){ 
//     const double INF = numeric_limits<double>::infinity();
//     //const double INF = 100;
//     int n = graph.get_size();
//     vector<double> distance(n, INF);
//     vector<int> parent(n, -1);
//     vector<bool> used(n, false);
//     priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

//     distance[dst] = 0;
//     pq.push({0, dst});
//     while(!pq.empty()) {
//         int cur_node = pq.top().second;
//         pq.pop();
//         if(used[cur_node]) continue;
//         used[cur_node] = true;
//         for(int neighbor : graph.get_neighbors_id(cur_node)) {
//             if(distance[cur_node] + X(cur_node, neighbor) < distance[neighbor]) {
//                 distance[neighbor] = distance[cur_node] + X(cur_node, neighbor);
//                 parent[neighbor] = cur_node;
//                 pq.push({distance[neighbor], neighbor});
//             }
//         }
//     }

//     if(distance[src] >= INF) return{};
// /*
//     int cur_node = src;
//     vector<int> path;
//     while(cur_node != -1) {
//         path.push_back(cur_node);
//         cur_node = parent[cur_node];
//     }

//     reverse(path.begin(), path.end());
// */
//     return parent;
// }

// void MyAlgo::separation_oracle_ori(int src,int dst,int req_no){     
//     vector<int> SPT;                  //nodes' parent in the spanning tree
//     vector<int> best_path;
//     double best_len; 
//     SPT = Dijkstra_ori(src, dst, req_no);                               //the first SPT is get by dijkstra
//     int cur_node = src;                                     //counting the first path's U(X,Y)=c* e^r
//     double c = 0;                                           //c=SUM[u,v]:alpha(u)+alpha(v)+beta(u,v)==X[u,v]
//     double r = 0;                                           //r=SUM[u,v]:-ln(Pr(u,v))==Y[u,v]
//     while(cur_node != dst){
//         // if(cur_node < SPT[cur_node]){                       //[can alter]no need if else
//         //     c += X[{cur_node,SPT[cur_node]}][req_no];               
//         //     r += Y[req_no][{cur_node,SPT[cur_node]}];
//         // }
//         // else{
//         //     c += X[{SPT[cur_node],cur_node}][req_no];  
//         //     r += Y[req_no][{SPT[cur_node],cur_node}];           
//         // }
//         c += X(SPT[cur_node], cur_node);  
//         r += Y[dst][req_no][{SPT[cur_node],cur_node}];  
//         best_path.push_back(cur_node);
//         cur_node = SPT[cur_node];
//     } 
//     best_path.push_back(dst);
//     best_len = c * exp(r);

//     cout << "\nBegin path: ";
//     for(auto p : best_path){
//             cout << p << "->";
//     }
//     cout << endl;
//     cout << "U = "<< best_len<<endl;
//     map<pair<int, int>, bool> used_edge;
//     vector<int> new_path;   
//     pair<int,int> new_edge;

//     for(unsigned int i = 0; i < SPT.size(); i++){
//         int cur_node=i;
//         while(cur_node!=dst){
//             if(used_edge.find({cur_node,SPT[cur_node]})!=used_edge.end() && used_edge[{cur_node,SPT[cur_node]}] != false){
//                 break;
//             }
//             used_edge[{cur_node,SPT[cur_node]}] = true;
//             used_edge[{SPT[cur_node],cur_node}] = true;
//             cur_node=SPT[cur_node];
//         }
//     }

//     while(1){
//         double minimum = numeric_limits<double>::infinity();
//         for(int i = 0; i < graph.get_size(); i++){                 //creating many new SPT
//             vector<int> neighbors = graph.get_neighbors_id(i);  
//             for(auto neighbor : neighbors){
//                 double temp1 = 0, temp2 = 0;
//                 if(SPT[i] == neighbor || SPT[neighbor] == i){      // checking used edge or unused
//                     continue;   
//                 }else{                                             // if unused
//                     temp1 = X(i, neighbor);
//                     temp2 = Y[dst][req_no][{i, neighbor}];
//                     int cur_node = i;
//                     while(cur_node != dst){
//                         temp1 += X(cur_node, SPT[cur_node]);
//                         temp2 += Y[dst][req_no][{cur_node, SPT[cur_node]}];
//                         cur_node = SPT[cur_node];
//                     } 

//                     cur_node = neighbor;
//                     while(cur_node != dst){
//                         temp1 -= X(cur_node, SPT[cur_node]);
//                         temp2 -= Y[dst][req_no][{cur_node, SPT[cur_node]}];
//                         cur_node = SPT[cur_node];
//                     }       
//                     if(temp2 < 0 && temp1 > 0 ) {               // we need the smallest edge to change the SPT
//                         if(i<neighbor){
//                             if(used_edge.find({i, neighbor}) != used_edge.end() && used_edge[{i, neighbor}] != false){
//                                 continue;
//                             }
//                         }
//                         else{
//                             if(used_edge.find({neighbor ,i}) != used_edge.end() && used_edge[{neighbor ,i}] != false){
//                                 continue;
//                             }
//                         }
//                         if(minimum > - temp1 / temp2){
//                             new_edge = {i, neighbor};
//                             minimum = - temp1 / temp2;
//                         }
//                     }
//                 }
//             }
//         }        // 找到最小的 edge 

//         if(minimum == numeric_limits<double>::infinity()){   //原本設計是有break,但之後用不到
//             break;
//         }else{
//             new_path.clear();
//         }
//         used_edge[{new_edge.second,SPT[new_edge.second]}] = false;
//         used_edge[{SPT[new_edge.second],new_edge.second}] = false;
//         used_edge[new_edge] = true;
//         used_edge[{new_edge.second,new_edge.first}] = true;
//         change_edge_num++;
//         SPT[new_edge.second]=new_edge.first;

//         cur_node = src;                                   
//         while(cur_node != dst) {
//             new_path.push_back(cur_node);
//             cur_node = SPT[cur_node];
//         }       
//         new_path.push_back(dst);
        
//         double new_len = 0;                                         //counting the new path's U(X,Y)=c* e^r
//         c = 0;
//         r = 0;
//         for(unsigned int i = 0; i < new_path.size() - 1; i++){
//             if(new_path[i] < new_path[i+1]){                        //[can alter]no need if else
//                 c += X(new_path[i], new_path[i+1]);
//                 r += Y[dst][req_no][{new_path[i], new_path[i+1]}];
//             }
//             else{
//                 c += X(new_path[i+1], new_path[i]);  
//                 r += Y[dst][req_no][{new_path[i+1], new_path[i]}];           
//             }
//             //cout<<"PATH:["<<new_path[i]<<" "<<new_path[i+1]<<"] with tot "<< c <<" / " << r <<endl;
//         }
//         new_len =  c * exp(r);
//         if(new_len < best_len){
//             cout<<"find new\n";
//             best_len = new_len;
//             best_path = new_path;                                            //路線修改,新的spt產生
//         } 
        
       
//     }
//     cout << "Best path: ";
//     for(auto p : best_path){
//         cout <<p << " ";
//     }
//     cout << endl;
//     cout << "U: " << best_len << endl;                                         
// }