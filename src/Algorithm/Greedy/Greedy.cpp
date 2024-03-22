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
        if(need_channel.find({path[i],path[i+1]}) != need_channel.end()){
            need_channel[{path[i],path[i+1]}] ++;
        }
        else{
            need_channel[{path[i],path[i+1]}] = 1;
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
            vector<int> path1 = BFS(request.get_node1(), smallest_middle);          //狠一點 別給它dij
            vector<int> path2 = BFS(request.get_node2(), smallest_middle);
            vector<int> path3 = BFS(request.get_node3(), smallest_middle);
            cal_need(path1,need_memory,need_channel); cal_need(path2,need_memory,need_channel); cal_need(path3,need_memory,need_channel);
            for(auto it :path1){
                cout << it << " ";
            }
            cout <<endl;
            for(auto it :path2){
                cout << it << " ";
            }
            cout <<endl;
            for(auto it :path3){
                cout << it << " ";
            }
            cout <<endl;
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
            if(continue_flag) continue;
            vector<vector<int>> tree;
            tree.push_back(path1);tree.push_back(path2);tree.push_back(path3);
            assign_resource(tree, 1,reqno);

            for(int i = 0; i < graph.get_size(); i++){
                if(need_memory[i] == 0){continue;}
                used_memory[i] += need_memory[i];
            }

            for(auto edge:need_channel){
                if(used_channel.find(edge.first) != used_channel.end()){
                    used_channel[edge.first] += edge.second;
                }
                else{
                    used_channel[edge.first] = edge.second;
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

