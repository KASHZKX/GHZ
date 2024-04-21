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

double MyAlgo::X(int u, int v, int req_no, int path_id, int middle){
    int node_num = graph.get_size();
    if((u == middle && v == (requests[req_no].get_node2() + node_num)) || (u == middle + node_num && v == (requests[req_no].get_node3() + node_num * 2)))
        return 0;
	double weight = alpha[u % node_num] + alpha[v % node_num] + beta[{u % node_num, v % node_num}];
	//if(path_id == 0 && (u == requests[req_no].get_node1() || v == requests[req_no].get_node1())) weight += tau[req_no];    //[check]                                                     //[Need fix!!!!!!!!!!!!!]
    if(u % node_num == middle || v % node_num == middle) weight += tau[req_no] / 3;
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
    int node_num = graph.get_size();
    for(unsigned  j = 0; j < requests.size(); j++){                      
        for(int k = 0; k < 3; k++){
            int node1 = requests[j].get_node1();
            int node2 = requests[j].get_node2();
            int node3 = requests[j].get_node3();
            vector<int> three_node;                                               //在轉圖中，G1只視端點1為起始點，G2只視端點2，...
            if(k == 0){  three_node.push_back(node1);}
            else if(k == 1){ three_node.push_back(node2);}
            else if(k == 2){ three_node.push_back(node3);}

            Y[mid][j][k][{mid, node_num + node2}] = 1;
            Y[mid][j][k][{mid + node_num, node_num * 2  + node3}] = 1;
            for(int i = 0; i < node_num; i++){                           //找所有邊[i,it]
                vector<int> temp = graph.get_neighbors_id(i);
                for(int l = 0; l < 3; l++){                             
                    for(auto it: temp){
                        double ent_p = exp(graph.Node_id2ptr(i)->distance(*graph.Node_id2ptr(it))*(-graph.get_entangle_alpha()));
                        if( i == mid && find(three_node.begin() , three_node.end() ,it) != three_node.end()){           //v=x,u=s                                                                          
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_fusion_prob()) / 3;              
                        }
                        else if( it == mid && find(three_node.begin() , three_node.end() ,i) != three_node.end()){      //u=x,v=s                                                                                                                                       
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_fusion_prob()) / 3;                   
                        }
                        else if( i != mid && find(three_node.begin() , three_node.end() ,it) != three_node.end()){      //v!=x,u=s                                                                                                                          
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob()) / 2;
                        }
                        else if( it != mid && find(three_node.begin() , three_node.end() ,i) != three_node.end()){      //u!=x,v=s                                                                                                                                                                                                                         
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 2;             
                        } 
                        else if( i == mid && find(three_node.begin() , three_node.end() ,it) == three_node.end()){      //v=x,u!=s、x
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - (log(graph.Node_id2ptr(it)->get_swap_prob()) / 2) - (log(graph.Node_id2ptr(i)->get_fusion_prob()) / 3);                                                                                        
                        }                                         

                        else if( it == mid && find(three_node.begin() , three_node.end() ,i) == three_node.end()){      //u=x,v!=s、x                                                                                                                                 
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - (log(graph.Node_id2ptr(i)->get_swap_prob()) / 2) - (log(graph.Node_id2ptr(it)->get_fusion_prob()) / 3);               
                        }                
                        else{                                                                                           //else                                                 
                            Y[mid][j][k][{i + l * node_num, it + l * node_num}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 2 - log(graph.Node_id2ptr(i)->get_swap_prob()) / 2;             
                        }
                    }
                }
            }
        }
    }

}

vector<int> MyAlgo::Dijkstra(int src, int dst, int req_no, int path_id, vector<pair<double,double>>&dist, vector<vector<int>> &n_graph){                           
    double INF=numeric_limits<double>::infinity();
    int node_num = graph.get_size();
    vector<bool> used(node_num * 3 + 1, false);
    vector<int> parent(node_num * 3 + 1, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    dist[src] = {0,0};                                                                  //找src到所有點的最短距離(針對C1-obj)
    pq.push({0, src});
    while(!pq.empty()){
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int neigh = 0; neigh < node_num * 3; neigh++){
            if(n_graph[cur_node][neigh] == 0) continue;
            if(dist[cur_node].first + X(cur_node, neigh, req_no, path_id, dst) < dist[neigh].first) {
                dist[neigh].first = dist[cur_node].first + X(cur_node, neigh, req_no, path_id, dst);      //(1)store d
                dist[neigh].second = dist[cur_node].second + Y[dst][req_no][path_id][{cur_node, neigh}];
                parent[neigh] = cur_node;                                             //(1)store pred
                pq.push({dist[neigh].first, neigh});
            }
            if(dist[cur_node].first + X(cur_node, neigh, req_no, path_id, dst) == dist[neigh].first && dist[neigh].second > (dist[cur_node].second + Y[dst][req_no][path_id][{cur_node, neigh}])) {
                dist[neigh].first = dist[cur_node].first + X(cur_node, neigh,req_no,path_id,dst);      //(1)store d
                dist[neigh].second = dist[cur_node].second + Y[dst][req_no][path_id][{cur_node, neigh}];
                parent[neigh] = cur_node;                                             //(1)store pred
                pq.push({dist[neigh].first, neigh}); 
            }

        }
    }

    if(dist[dst].first >= INF) return{};
    return parent;
}       
   

void MyAlgo::separation_oracle(int src, int dst, int req_no, int path_id, vector<vector<vector<int>>> &cur_tree, vector<vector<vector<double>>> &cur_label, vector<vector<int>> &n_graph, int &U){     
    double INF = numeric_limits<double>::infinity();                    //以下號碼根據[2015 Raith] A Dijkstra-like method算法
    int node_num = graph.get_size();
    vector<pair<double,double>>dist (node_num * 3 + 1, {INF,INF});       //存放d^1_t、d^2_t  
    vector<int> pred (node_num * 3 + 1,-1);                             //存放pred     
    vector<int> SDpath;
    vector<int> best_path;
    double smallest_U = INF;
    vector<double> theta_table(node_num * 3 + 1, INF);                   //存放theta
    vector<double> obj_bar1_table(node_num * 3 + 1, INF);                //存放oba_bar1
    vector<double> obj_bar2_table(node_num * 3 + 1, INF);                //存放oba_bar2
    vector<int> cpred_table(node_num * 3 + 1, -1);                       //存放cpred
    priority_queue<Label,vector<Label>, greater<Label>> pq;             //Label <theta,oba_bar2,cpred>                         
    pred = Dijkstra(src, dst, req_no, path_id, dist, n_graph);                   //(1)找C1最小的path

    double last_ratio = 0;                                              
    double d1 = dist[dst].first, d2 = dist[dst].second;                  
    int current = dst;
    while(current != -1){
        SDpath.push_back(current);
        current = pred[current];
    }
    reverse(SDpath.begin(),SDpath.end());                               

    if(d1 * exp(d2) < smallest_U){                                      //[測試]與old separation比較時使用
        smallest_U = d1 * exp(d2);
        best_path = SDpath;

    }
    // cout<<"Begin path:";                                             
    // for(auto it: SDpath){
    //     cout<<it<<" ";
    // }
    // cout<<endl;
    // cout<<"U = "<< d1 * exp(d2)<<endl;
    for(int i = 0; i < node_num * 3; i++){                         //(4)compute lexmin
        if(i == src) continue;
        for(int neigh = 0; neigh < node_num * 3; neigh++){
            if(n_graph[i][neigh] == 0) continue;
            double cur_obj_bar2 = dist[neigh].second + Y[dst][req_no][path_id][{neigh, i}] - dist[i].second;
            if(cur_obj_bar2 < 0){                                      //確認oba_bar2有沒有小於0
                double cur_theta = -(dist[neigh].first + X(neigh, i, req_no, path_id, dst) - dist[i].first) / cur_obj_bar2 ;
                if(cur_theta <= theta_table[i]){                       //計算對於i來說 theta最小的predecessor是誰
                    if(cur_theta < theta_table[i]){
                        theta_table[i] = cur_theta;
                        obj_bar1_table[i] = dist[neigh].first + X(neigh, i, req_no, path_id, dst) - dist[i].first;
                        obj_bar2_table[i] = cur_obj_bar2;
                        cpred_table[i] = neigh;
                    }
                    else if(cur_obj_bar2 < obj_bar2_table[i]){         //如果theta相同，則oba_bar2比較小的作為代表
                        theta_table[i] = cur_theta;
                        obj_bar1_table[i] = dist[neigh].first + X(neigh, i, req_no, path_id, dst) - dist[i].first;
                        obj_bar2_table[i] = cur_obj_bar2;
                        cpred_table[i] = neigh;
                    }
                }
            }
        }
        if(cpred_table[i] != -1){                                       //如果有找到比現在好的，丟到priority_queue
            pq.push({theta_table[i],obj_bar2_table[i], i});
            //cout<<"Push1: "<<theta_table[i]<<" "<<obj_bar2_table[i]<<" "<<i<<"\n";
        }
    }

    while(!pq.empty()){                                                 //(6)若priority_queue不為空
        double min_theta, min_obj_bar2;
        int min_i;
        Label temp = pq.top();                                          //(7)把theta最小的抓出來
        min_theta = temp.theta, min_obj_bar2 = temp.obj_bar2, min_i = temp.i;
        pq.pop();                                                       //根據換新的邊，做dist修改與parent修改
        if(min_obj_bar2 != obj_bar2_table[min_i]){continue;}
        //cout<<"pop: "<<min_theta<<" "<<min_obj_bar2<<" "<<min_i<<endl;
        dist[min_i].first += obj_bar1_table[min_i];
        dist[min_i].second += obj_bar2_table[min_i];
        pred[min_i] = cpred_table[min_i];
        if(min_i == dst){                                               //(9)如果換到終點，則找到extreme path
            if(min_theta > last_ratio){                                 //(10)
                //cout << "[" << last_ratio << "," << min_theta << "] with " << d1 * exp(d2) << endl;
                vector<double>temp = {last_ratio,min_theta, (d1 * exp(d2))};
                last_ratio = min_theta;
                cur_tree[path_id].push_back(SDpath);                    //extreme path丟入
                cur_label[path_id].push_back(temp);                     //丟入該extreme tree的[las_ratio , cur_ratio]
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
            reverse(SDpath.begin(),SDpath.end());                      //(12)修改當前的路線
            d1 = dist[dst].first,d2 = dist[dst].second;

            if(d1 * exp(d2) < smallest_U){
                smallest_U = d1 * exp(d2);
                best_path = SDpath;
                U = smallest_U;
            }
        }

        theta_table[min_i] = INF, cpred_table[min_i] = -1, obj_bar1_table[min_i] = INF, obj_bar2_table[min_i] = INF;       //(13)reset 剛剛換邊的那個點的所有資訊
        
        for(int neigh = 0; neigh < node_num * 3; neigh++){
            if(n_graph[min_i][neigh] == 0) continue;
            double cur_obj_bar2 = dist[neigh].second + Y[dst][req_no][path_id][{neigh,min_i}] - dist[min_i].second;
            if(cur_obj_bar2 < 0){
                double cur_theta = -(dist[neigh].first + X(neigh,min_i,req_no,path_id,dst) - dist[min_i].first) / cur_obj_bar2 ;
                if(cur_theta <= theta_table[min_i]){
                    if(cur_theta < theta_table[min_i]){
                        theta_table[min_i] = cur_theta;
                        obj_bar1_table[min_i] = dist[neigh].first + X(neigh,min_i,req_no,path_id,dst) - dist[min_i].first;
                        obj_bar2_table[min_i] = cur_obj_bar2;
                        cpred_table[min_i] = neigh;
                    }
                    else if(cur_obj_bar2 < obj_bar2_table[min_i]){
                        theta_table[min_i] = cur_theta;                                                                //no use just for read
                        obj_bar1_table[min_i] = dist[neigh].first + X(neigh,min_i,req_no,path_id,dst) - dist[min_i].first;
                        obj_bar2_table[min_i] = cur_obj_bar2;
                        cpred_table[min_i] = neigh;
                    }
                }
            }
        }

        if(cpred_table[min_i] != -1){                                                                                   //有找到就push進priority_queue
            pq.push({theta_table[min_i],obj_bar2_table[min_i],min_i});
            //cout<<"Push2: "<<theta_table[min_i]<<" "<<obj_bar2_table[min_i]<<" "<<min_i<<"\n";
        }
        for(int neigh = 0; neigh < node_num * 3; neigh++){ //(14)~(19)對successor找最小的替換邊
            if(n_graph[min_i][neigh] == 0) continue;
            double temp_obj_bar2 = dist[min_i].second + Y[dst][req_no][path_id][{min_i,neigh}] - dist[neigh].second;
            if(temp_obj_bar2 < 0){
                double temp_obj_bar1 = dist[min_i].first + X(min_i,neigh,req_no,path_id,dst) - dist[neigh].first;
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
    // cout<<"Best path :";
    // for(auto it:best_path){
    //     cout<<it<<" ";
    // }
    vector<double>temp = {last_ratio,INF,(d1 * exp(d2))};                          //(20)
    cur_tree[path_id].push_back(SDpath);
    cur_label[path_id].push_back(temp);        
    //cout<<endl<<"U : "<< smallest_U <<endl;
}

void MyAlgo::find_bottleneck(vector<vector<int>> &tree, int req_no){
    double min_s_u = numeric_limits<double>::infinity();                    //memory違反最多幾倍
    double min_s_uv = numeric_limits<double>::infinity();                   //channel違反最多幾倍
    double s_i = 1;                                                         //send-limit違反最多幾倍
    vector<double> s_u(graph.get_size() + 5);                               
    vector<int> used_u(graph.get_size() + 5, 0);
    map<pair<int,int>,int> used_uv;
    map<pair<int,int>,double> s_uv;
    // Calculate used 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < tree[i].size() - 1; j++){
            if(j == 0){
                used_u[tree[i][j]]++;                                       //端點memory花一個
            }
            else{
                used_u[tree[i][j]]+=2;                                      //中繼點都花兩個
            }
            if(used_uv.find({tree[i][j],tree[i][j+1]}) != used_uv.end()){   //邊花一條
                used_uv[{tree[i][j],tree[i][j+1]}]++;
                used_uv[{tree[i][j+1],tree[i][j]}]++;
            }
            else{
                used_uv[{tree[i][j],tree[i][j+1]}] = 1;
                used_uv[{tree[i][j+1],tree[i][j]}] = 1;
            }
        }
        used_u[tree[i][tree[i].size()-1]]++;                               //中繼點 貢獻給一條路線一個memory，有三條路，共3個memories
    }

    for(int i = 0; i < graph.get_size(); i++){
        if(used_u[i] == 0){continue;}
        s_u[i] = graph.Node_id2ptr(i)->get_memory_cnt() / used_u[i] ;      //計算memory超過倍數    
        if(s_u[i] < min_s_u)
            min_s_u = s_u[i];
    }

    for(auto &it:used_uv){
        s_uv[{it.first.first, it.first.second}] = graph.get_channel_size(it.first.first, it.first.second) / it.second;    //計算channel超過倍數
        if(s_uv[{it.first.first, it.first.second}] < min_s_uv)
            min_s_uv = s_uv[{it.first.first, it.first.second}];
    }

    int rate = 10;
    double s = min(min_s_u, min(min_s_uv, s_i));                          //違反最多的倍數
    for(int i = 0; i < rate; i++){
        bool find = false;
        for(int j = 0; j < x_i_t_tree.size(); j++){                       //add flow to path
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
            obj += (beta[{it.first.first, it.first.second}] * (1 + epsilon * s / s_uv[{it.first.first, it.first.second}]) -  beta[{it.first.first, it.first.second}]) * graph.get_channel_size(it.first.first, it.first.second);
            beta[{it.first.first, it.first.second}] = beta[{it.first.first, it.first.second}] * (1 + epsilon * s / s_uv[{it.first.first, it.first.second}]);
            beta[{it.first.second, it.first.first}] = beta[{it.first.second, it.first.first}] * (1 + epsilon * s / s_uv[{it.first.second, it.first.first}]);
        }
        obj += (tau[req_no] * (1 + epsilon * s ) - tau[req_no]);                            //request limit = 1?
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


    

    //cout << "Magnification:" << max_magni << endl;

    for(auto &x : x_i_t){
        x /= max_magni;
    }
}

void MyAlgo::calculate(){
    double UB = 0;
    double entangle_alpha = graph.get_entangle_alpha();
    for(int i = 0; i < x_i_t_tree.size(); i++){
        double prob = 1;
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < x_i_t_tree[i][j].size() - 1; k++){
                Node* node1=graph.Node_id2ptr(x_i_t_tree[i][j][k]);
                Node* node2=graph.Node_id2ptr(x_i_t_tree[i][j][k+1]);
                prob *= exp(-entangle_alpha * (node1->distance(*node2)));
            }
            for(int k = 1; k < x_i_t_tree[i][j].size() - 1; k++){
                Node* node=graph.Node_id2ptr(x_i_t_tree[i][j][k]);
                prob *= node->get_swap_prob();
            }
            prob *= pow(graph.Node_id2ptr(x_i_t_tree[i][j][x_i_t_tree[i][j].size() - 1])->get_fusion_prob(),(1/3));
        }
        
        bool check_flag = false;
        for(int j = 0; j < requests.size(); j++){
            if(requests[j].get_node1() == x_i_t_tree[i][0][0] && requests[j].get_node2() == x_i_t_tree[i][1][0] && requests[j].get_node3() == x_i_t_tree[i][2][0]){
                check_flag = true;
                UB += requests[j].get_value() * prob * x_i_t[i] * pow(0.8,-2);
            }
        }
        if(check_flag == false){
            cout<<"No matching request......."<<endl;
            exit(1);
        }
    }
    res["UB"] = UB;
}

vector<map<vector<vector<int>>, int>> MyAlgo::Greedy_rounding(){
    vector<vector<int>> temp;
    vector<tuple<double, int, vector<vector<int>>>> each_request(requests.size(),{-1,-1,temp});
    vector<map<vector<vector<int>>, int>> I_request(requests.size());
    vector<int>num(requests.size(),0);
    for(int i = 0; i < x_i_t_tree.size(); i++){
        vector<vector<int>> tree = x_i_t_tree[i];
        int node1 = tree[0][0];
        int node2 = tree[1][0];
        int node3 = tree[2][0];
        for(unsigned int j = 0; j < requests.size(); j++){
            
            if(node1 == requests[j].get_node1() && node2 == requests[j].get_node2() && node3 == requests[j].get_node3()){
                num[j]++;
                double prob; int req;
                tie(prob,req,temp) =  each_request[j];
                if(x_i_t[i] > prob){
                    each_request[j] = tie(x_i_t[i], j, tree);
                }
                break;
            }
        }
    }
    sort(each_request.begin(),each_request.end(),greater<tuple<double, int, vector<vector<int>>>>());
    // for(auto &it:each_request){
    //     double prob; int req;
    //     tie(prob,req,temp) = it;
    //     cout<<"Request :"<<req<<" ----------------------- "<<prob<<endl;
    //     for(auto &it2:temp){
    //         for(auto &it3:it2){
    //             cout<<it3<<" ";
    //         }
    //         cout<<endl;
    //     }
    //     cout<<"-------------"<<endl;
    // }

	// 如果資源足夠讓 x 變成 1 ，就直接讓 x 變成 1 
    bool global_enough;
    do{
        global_enough = false;
        for(auto &it:each_request){
            double prob;
            int request_id;
            vector<vector<int>> tree;
            tie(prob, request_id, tree) = it;

            if(request_id == -1){continue;}                                                //代表有請求沒有tree QQ
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
            //is_assinable get_remain
            bool enough=true;
            for(int i = 0; i < graph.get_size(); i++){
                if(graph.Node_id2ptr(i)->get_remain() < memory_use[i]){
                    enough = false;
                    break;
                }
            }
            for(auto &it:channel_use){
                //cout<<"Channel_remain : "<<graph.remain_channel(it.first.first,it.first.second)<<" "<<it.second<<endl;
                if(graph.remain_channel(it.first.first,it.first.second) < it.second){
                    enough = false;
                    break;
                }
            }
            if(enough){
                //cout<<"hi"<<endl;
                global_enough = true;
                assign_resource(tree, 1, request_id);
                if(I_request[request_id].find(tree) != I_request[request_id].end()){
                    I_request[request_id][tree]++;
                }
                else{
                    I_request[request_id][tree] = 1;
                }
                
            }
        }

    }while(global_enough);

    // heuristic algorithm
    do{
        global_enough = false;
        for(auto &tree:x_i_t_tree){
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
            //is_assinable get_remain
            bool enough=true;
            for(int i = 0; i < graph.get_size(); i++){
                if(graph.Node_id2ptr(i)->get_remain() < memory_use[i]){
                    enough = false;
                    break;
                }
            }
            for(auto &it:channel_use){
                //cout<<"Channel_remain : "<<graph.remain_channel(it.first.first,it.first.second)<<" "<<it.second<<endl;
                if(graph.remain_channel(it.first.first,it.first.second) < it.second){
                    enough = false;
                    break;
                }
            }
            if(enough){
                //cout<<"hi"<<endl;
                global_enough = true;
                int node1 = tree[0][0];
                int node2 = tree[1][0];
                int node3 = tree[2][0];
                int request_id;
                for(int i = 0; i < requests.size(); i++){
                    if(node1 == requests[i].get_node1() && node2 == requests[i].get_node2() && node3 == requests[i].get_node3()){
                        request_id = i;
                    }
                }
                assign_resource(tree, 1, request_id);
                if(I_request[request_id].find(tree) != I_request[request_id].end()){
                    I_request[request_id][tree]++;
                }
                else{
                    I_request[request_id][tree] = 1;
                }
            }
        }
	}while(global_enough);
    return I_request;
}

void MyAlgo::create_graph(vector<vector<vector<vector<int>>>> &n_graph){
    int node_num = graph.get_size();
    for(int i = 0; i < node_num; i++){
        for(int j = 0; j < requests.size(); j++){
            int node1 = requests[j].get_node1();
            int node2 = requests[j].get_node2();
            int node3 = requests[j].get_node3();
            if(i == node1 || i == node2 || i == node3){continue;}
            n_graph[i][j][i][node_num + node2] = 1;
            n_graph[i][j][i + node_num][node_num * 2 + node3] = 1;
            for(int k = 0; k < 2; k++){
                for(int l = 0; l < node_num; l++){
                    vector<int> tmp = graph.get_neighbors_id(l);
                    for(auto neighbor : tmp){
                        n_graph[i][j][l + node_num * k][neighbor + node_num * k] = 1;
                    }
                }
            }

        }
    }
}

void MyAlgo::path_assignment(){
    double INF = numeric_limits<double>::infinity();
    int node_num = graph.get_size();
    vector<vector<vector<vector<int>>>> n_graph(node_num, vector<vector<vector<int>>>(requests.size(), vector<vector<int>>(node_num * 3 + 1, vector<int>(node_num * 3 + 1, 0))));
    Y.resize(node_num);                                     // Y結構為[middle_node][request][3] 計算以midddle node為X、request三個端點為端點的Y圖
    for(int i = 0; i < node_num; i++){
        Y[i].resize(requests.size() + 1);
        for(int j = 0; j < requests.size() + 1; j++){
            Y[i][j].resize(3);
        }
    }    

    init_dual();                                                  
    for(int i = 0; i < node_num; i++){
        initialize(i);                                              // 計算Y
    }


    obj = M * delta;
    while(obj < 1){
        int req_no = -1;
        vector<vector<int>> best_tree(3, vector<int>());                                         //存放在所有middle、request計算中最好的tree
        double global_U = numeric_limits<double>::infinity();
        for(int middle = 0; middle < node_num; middle++){
            for(unsigned int i = 0; i < requests.size(); i++){
                double smallest_U = numeric_limits<double>::infinity();
                int smallest_index = -1;
                vector<vector<vector<int>>> cur_tree(3,vector<vector<int>>());                  //[path_id][extreme_tree_id][eles]
                vector<vector<vector<double>>> cur_label(3,vector<vector<double>>());           //[path_id][extreme_tree_id][last_ratio, new_ratio, U]
                if(middle == requests[i].get_node1() || middle == requests[i].get_node2() || middle == requests[i].get_node3()){continue;}      //middle_node不可以是三個端點
                separation_oracle(requests[i].get_node1(), request[i].get_node3() + node_num * 2, i, 0, cur_tree, cur_label, n_graph[middle][i], smallest_U);  //對端點一到X做溜滑梯

                      
               
                smallest_U /= requests[i].get_value();   //added
                if(smallest_U < global_U){
                    global_U = smallest_U;
                    req_no = i;

                    // cout<<"Request "<<i<<" with : "<<smallest_U <<" \n";
                    // for(auto it:best_tree){
                    //     for(auto it1:it){
                    //         cout<<it1<<" ";
                    //     }
                    //     cout<<endl;
                    // }
                    // cout<<endl;
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
    calculate();
    vector<map<vector<vector<int>>, int>>path = Greedy_rounding();
}


//                   Our Old Friend                                //
vector<int> MyAlgo::Dijkstra_ori(int src, int dst, int req_no, int path_id){ 
    const double INF = numeric_limits<double>::infinity();
    //const double INF = 100;
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
            if(distance[cur_node] + X(cur_node, neighbor,req_no,path_id,dst) < distance[neighbor]) {
                distance[neighbor] = distance[cur_node] + X(cur_node, neighbor,req_no,path_id,dst);
                parent[neighbor] = cur_node;
                pq.push({distance[neighbor], neighbor});
            }
        }
    }

    if(distance[src] >= INF) return{};

    return parent;
}

void MyAlgo::separation_oracle_ori(int src,int dst,int req_no, int path_id){     
    vector<int> SPT;                  //nodes' parent in the spanning tree
    vector<int> best_path;
    double best_len; 
    SPT = Dijkstra_ori(src, dst, req_no, path_id);                               //the first SPT is get by dijkstra
    int cur_node = src;                                     //counting the first path's U(X,Y)=c* e^r
    double c = 0;                                           //c=SUM[u,v]:alpha(u)+alpha(v)+beta(u,v)==X[u,v]
    double r = 0;                                           //r=SUM[u,v]:-ln(Pr(u,v))==Y[u,v]
    while(cur_node != dst){
        c += X(SPT[cur_node], cur_node,req_no,path_id,dst);  
        r += Y[dst][req_no][path_id][{SPT[cur_node],cur_node}];  
        best_path.push_back(cur_node);
        cur_node = SPT[cur_node];
    } 
    best_path.push_back(dst);
    best_len = c * exp(r);

    // cout << "\nBegin path[2]: ";
    // for(auto p : best_path){
    //         cout << p << "->";
    // }
    // cout << endl;
    // cout << "U = "<< best_len<<endl;
    map<pair<int, int>, bool> used_edge;
    vector<int> new_path;   
    pair<int,int> new_edge;

    for(unsigned int i = 0; i < SPT.size(); i++){
        int cur_node=i;
        while(cur_node!=dst){
            if(used_edge.find({cur_node,SPT[cur_node]})!=used_edge.end() && used_edge[{cur_node,SPT[cur_node]}] != false){
                break;
            }
            used_edge[{cur_node,SPT[cur_node]}] = true;
            used_edge[{SPT[cur_node],cur_node}] = true;
            cur_node=SPT[cur_node];
        }
    }

    while(1){
        double minimum = numeric_limits<double>::infinity();
        for(int i = 0; i < graph.get_size(); i++){                 //creating many new SPT
            vector<int> neighbors = graph.get_neighbors_id(i);  
            for(auto neighbor : neighbors){
                double temp1 = 0, temp2 = 0;
                if(SPT[i] == neighbor || SPT[neighbor] == i){      // checking used edge or unused
                    continue;   
                }else{                                             // if unused
                    temp1 = X(i, neighbor,req_no,path_id,dst);
                    temp2 = Y[dst][req_no][path_id][{i, neighbor}];
                    int cur_node = i;
                    while(cur_node != dst){
                        temp1 += X(cur_node, SPT[cur_node],req_no,path_id,dst);
                        temp2 += Y[dst][req_no][path_id][{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    } 

                    cur_node = neighbor;
                    while(cur_node != dst){
                        temp1 -= X(cur_node, SPT[cur_node],req_no,path_id,dst);
                        temp2 -= Y[dst][req_no][path_id][{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    }       
                    if(temp2 < 0 && temp1 > 0 ) {               // we need the smallest edge to change the SPT
                        if(i<neighbor){
                            if(used_edge.find({i, neighbor}) != used_edge.end() && used_edge[{i, neighbor}] != false){
                                continue;
                            }
                        }
                        else{
                            if(used_edge.find({neighbor ,i}) != used_edge.end() && used_edge[{neighbor ,i}] != false){
                                continue;
                            }
                        }
                        if(minimum > - temp1 / temp2){
                            new_edge = {i, neighbor};
                            minimum = - temp1 / temp2;
                        }
                    }
                }
            }
        }        // 找到最小的 edge 

        if(minimum == numeric_limits<double>::infinity()){   //原本設計是有break,但之後用不到
            break;
        }else{
            new_path.clear();
        }
        used_edge[{new_edge.second,SPT[new_edge.second]}] = false;
        used_edge[{SPT[new_edge.second],new_edge.second}] = false;
        used_edge[new_edge] = true;
        used_edge[{new_edge.second,new_edge.first}] = true;
        change_edge_num++;
        SPT[new_edge.second]=new_edge.first;

        cur_node = src;                                   
        while(cur_node != dst) {
            new_path.push_back(cur_node);
            cur_node = SPT[cur_node];
        }       
        new_path.push_back(dst);

        double new_len = 0;                                         //counting the new path's U(X,Y)=c* e^r
        c = 0;
        r = 0;
        for(unsigned int i = 0; i < new_path.size() - 1; i++){
            if(new_path[i] < new_path[i+1]){                        //[can alter]no need if else
                c += X(new_path[i], new_path[i+1],req_no,path_id,dst);
                r += Y[dst][req_no][path_id][{new_path[i], new_path[i+1]}];
            }
            else{
                c += X(new_path[i+1], new_path[i],req_no,path_id,dst);  
                r += Y[dst][req_no][path_id][{new_path[i+1], new_path[i]}];           
            }
            //cout<<"PATH:["<<new_path[i]<<" "<<new_path[i+1]<<"] with tot "<< c <<" / " << r <<endl;
        }
        new_len =  c * exp(r);
        if(new_len < best_len){
            //cout<<"find new\n";
            best_len = new_len;
            best_path = new_path;                                            //路線修改,新的spt產生
        } 


    }
    cout << "Best path[2]: ";
    for(auto p : best_path){
        cout <<p << " ";
    }
    cout << endl;
    cout << "U: " << best_len << endl<<endl;                                         
}
