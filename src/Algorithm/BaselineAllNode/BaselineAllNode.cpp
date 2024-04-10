#include "BaselineAllNode.h"

BaselineAllNode::BaselineAllNode(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not)
    :AlgorithmBase(filename, "BaselineAllNode", request_time_limit, node_time_limit, swap_prob, entangle_alpha , limit_r_or_not){
    if(DEBUG) cerr<<"new BaselineAllNode"<<endl;
}

vector<vector<int>> BaselineAllNode::Dijkstra_Tree(vector<int> &terminal){
    const double INF = numeric_limits<double>::infinity();
    int node_num = graph.get_size();

    double min_weight = INF;
    double max_prob = 0;

    vector<vector<int>> tree;
    tree.clear();
    for(int i = 0; i < node_num; i++){
        if(i == terminal[0] || i == terminal[1] || i == terminal[2]) continue;
        vector<vector<int>> tmp_tree;
        tmp_tree.clear();

        double prob_sum = 1; //double weight_sum = 0;

        int j = 0;
        for(j = 0; j < 3; j++){ 
            int destination = terminal[j];
            vector<double> prob(node_num, 0); //vector<double> dist(node_num, INF);
            vector<int> prev(node_num, -1);
            priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

            prob[i] = 1;    //dist[i] = 0;
            pq.push({1, i});    //pq.push({0, i});

            while (!pq.empty()) {
                //double porbability = pq.top().first;  //double distance = pq.top().first;
                int u = pq.top().second;
                pq.pop();

                if (u == destination) break;

                for(int v : graph.get_neighbors_id(u)) {
                    if (u == v) continue;

                    Node *node1 = graph.Node_id2ptr(u);
                    //prob *= node->get_swap_prob();

                    double probability = node1->get_swap_prob(); //get swap prob of u   //double weight = alpha[u] + beta[{u, v}] + tau[v];
                    Node* node2 = graph.Node_id2ptr(v);
                    probability *= exp(-graph.get_entangle_alpha() * (node1->distance(*node2)));   //weight = (swap prob of u)*(entangle prob of u and v)

                    if (prob[v] < prob[u] * probability) {   //if (dist[v] > dist[u] + weight) {
                        bool is1_repeater = true, is2_repeater = true;
                        if(u == i) is1_repeater = false;
                        if(v == destination) is2_repeater = false;

                        if(graph.remain_resource_cnt(u, v, is1_repeater, is2_repeater) > 0){
                            prob[v] = prob[u]*probability;   //dist[u] + weight;
                            prev[v] = u;
                            pq.push({prob[v], v});
                        }
                    }
                }
            }
            if(prev[destination] == -1) break;

            prob_sum *= prob[destination];  //weight_sum += dist[destination];

            vector<int> path;
            path.clear();
            for (int at = destination; at != -1; at = prev[at])
                path.push_back(at);
            reverse(path.begin(), path.end());

            tmp_tree.push_back(path);
        }

        bool flag = true;

        for(int k = 0; k < 3 && j == 3; k++){
            for(int l = 1; l < tmp_tree[k].size(); l++){
                bool is1_repeater = true, is2_repeater = true;
                if(l == 1) is1_repeater = false;
                if(l ==  tmp_tree[k].size()-1) is2_repeater = false;
                if(graph.remain_resource_cnt( tmp_tree[k][l-1],  tmp_tree[k][l], is1_repeater, is2_repeater) <= 0){
                    flag = false;
                }
                if(!flag)
                    break;
            }
            if(!flag)
                break;
        }

        //if(min_weight > weight_sum && weight_sum != 0 && j == 3 && flag){
        if(max_prob < prob_sum && prob_sum != 1 && j == 3 && flag){
            max_prob = prob_sum;    //min_weight = weight_sum;
            tree = tmp_tree;
        }
    }
    if(tree.size() == 0)
        return {};

    return tree;
}

bool BaselineAllNode::checkResource(vector<vector<int>> &tree){
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

void BaselineAllNode::path_assignment(){
    bool flag = true;
    while(flag){
        flag = false;
        for(int i = 0; i < (int)requests.size(); i++){
            Request &request = requests[i];
            vector<int> terminal = {request.get_node1(), request.get_node2(), request.get_node3()}; //J: get the 3 users' nodes
            vector<vector<int>> tree = Dijkstra_Tree(terminal);
            if(tree.size() == 0 ) continue;
            if(checkResource(tree)){
                flag = true;
                assign_resource(tree, 1, i);
            }
        }
    }
}

void BaselineAllNode::entangle(){
    AlgorithmBase::base_entangle();
}
void BaselineAllNode::swap(){
    AlgorithmBase::base_swap();
}
void BaselineAllNode::send(){
    AlgorithmBase::base_send();
}
void BaselineAllNode::next_time_slot(){
    AlgorithmBase::base_next_time_slot();
}

