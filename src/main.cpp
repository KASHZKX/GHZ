//             #memory version
#include <iostream>
#include <queue>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/Greedy/Greedy.h"
#include "Algorithm/MyAlgo/MyAlgo.h"
#include "Algorithm/BaselineDual/BaselineDual.h"
#include "Algorithm/BaselineAllNode/BaselineAllNode.h"
// #include "Algorithm/MyGreedyAlgo/MyGreedyAlgo.h"

using namespace std;

Request generate_new_request(int num_of_node, int time_limit){
    //亂數引擎 
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif1(0, num_of_node-1);
    int node1 = unif1(generator), node2 = unif1(generator) , node3 = unif1(generator);
    while(node1 == node2) node2 = unif1(generator);
    while(node3 == node2 || node3 == node1) node3 = unif1(generator);

    return Request(node1, node2, node3, 1.0);
}

Request generate_fix_request(int node1, int node2, int node3, int time_limit, int request){//demo
    return Request(node1, node2, node3 , 1.0);
}

void create_dir_if_not_exists(const std::string &path) {
	struct stat info;
	if (stat(path.c_str(), &info) == 0 && info.st_mode & S_IFDIR) {
		return;
	}
	mkdir(path.c_str(), 0777);
	return;
}



int main(int argc, char *argv[]){
    string file_path = "../data/";
    if(argc > 1){
	cout<<argv[1]<<endl;
	file_path = string(argv[1]) + '/';
    	cerr<<"the data is store in "<<file_path<<endl;
	create_dir_if_not_exists(file_path);
	create_dir_if_not_exists(file_path+"ans");
	create_dir_if_not_exists(file_path+"input");
	create_dir_if_not_exists(file_path+"log");

    }
    map<string, double> default_setting;
    default_setting["num_of_node"] = 50;
    default_setting["area_alpha"] = 0.1;
    default_setting["memory_cnt_avg"] = 12;
    default_setting["channel_cnt_avg"] = 5;
    default_setting["resource_ratio"] = 1;
    
    default_setting["swap_prob"] = 0.875;
    default_setting["entangle_alpha"] = 0.0002;
    default_setting["new_request_cnt"] = 70;
    default_setting["total_time_slot"] = 1;
    default_setting["epsilon"] = 0.2;    
    default_setting["fusion_prob"] =0.5;

    // not used in this paper
    default_setting["node_time_limit"] = 1;
    default_setting["social_density"] = 0.5;
    default_setting["min_fidelity"] = 0.7;
    default_setting["max_fidelity"] = 0.95;
    default_setting["request_time_limit"] = 1;
    default_setting["service_time"] = 100;
    
    map<string, vector<double>> change_parameter;
    change_parameter["swap_prob"] = {0.75,0.8,0.85,0.9,0.95};

    vector<string> X_names =  { "swap_prob"}; 
    vector<string> Y_names =  { /*"max_over_ratio",*/"total_earn"};
    vector<string> algo_names = { "MyAlgo", "Greedy", "BaselineDual", "BaselineAllNode" }; 

    // init result
    for(string X_name : X_names) {
        for(string Y_name : Y_names){
            string filename = "ans/" + X_name + "_" + Y_name + ".ans";
            fstream file( file_path + filename, ios::out );
        }
    }

    int round = 60;
    map<string, double> input_parameter = default_setting;
    vector<vector<Request>>request_list(round);
    vector<map<string, map<string, double>>> result(round);
    map<string,map<string,vector<double>>> sum_vt;    
    
    int num_of_node = input_parameter["num_of_node"];
    // double social_density = input_parameter["social_density"];
    double area_alpha = input_parameter["area_alpha"];
    double resource_ratio = input_parameter["resource_ratio"];
    int min_memory_cnt = input_parameter["memory_cnt_avg"] * resource_ratio - 2;
    int max_memory_cnt = input_parameter["memory_cnt_avg"] * resource_ratio + 2;
    int min_channel_cnt = input_parameter["channel_cnt_avg"] * resource_ratio - 2;
    int max_channel_cnt = input_parameter["channel_cnt_avg"] * resource_ratio + 2;
    double min_fidelity = input_parameter["min_fidelity"];
    double max_fidelity = input_parameter["max_fidelity"];

    double swap_prob = input_parameter["swap_prob"], entangle_alpha = input_parameter["entangle_alpha"];
    double min_swap_prob = input_parameter["swap_prob"] - 0.05;
    double max_swap_prob = input_parameter["swap_prob"] + 0.05;
    double fusion_prob = input_parameter["fusion_prob"];
    double min_fusion_prob = input_parameter["fusion_prob"] - 0.5;
    double max_fusion_prob = input_parameter["fusion_prob"] + 0.5;

    int node_time_limit = input_parameter["node_time_limit"];
    int new_request_cnt = input_parameter["new_request_cnt"];
    int service_time = input_parameter["service_time"];
    int request_time_limit = input_parameter["request_time_limit"];
    int total_time_slot = input_parameter["total_time_slot"];

    random_device rd2;
    default_random_engine generator2 = default_random_engine(rd2());
    std::normal_distribution<double> normal_distribution(1.0,1.5);

    // make T張 graph
    #pragma omp parallel for
    for(int T = 0; T < round; T++){
        string round_str = to_string(T);
        string filename = file_path + "input/round_" + round_str + ".input";
        string command = "python3 main.py ";
        string parameter = to_string(num_of_node) + " " + to_string(min_channel_cnt) + " " + to_string(max_channel_cnt) + " " + to_string(min_memory_cnt) + " " + to_string(max_memory_cnt) + " " + to_string(min_fidelity) + " " + to_string(max_fidelity) + " " + to_string(area_alpha) + " " + to_string(min_swap_prob) + " " +  to_string(max_swap_prob) + " " + to_string(min_fusion_prob) + " " + to_string(max_fusion_prob);
        //cout<<command + filename + " " + parameter<<endl;
        if(system((command + filename + " " + parameter).c_str()) != 0){
            cerr<<"error:\tsystem proccess python error"<<endl;
            exit(1);
        }
        for(int q = 0; q < new_request_cnt; q++){                                               //同時把所有requests都放到request_list
            bool check_no_repeat;
            do{
                check_no_repeat=true;
                Request new_request = generate_new_request(num_of_node, request_time_limit);
                for(auto &it:request_list[T]){
                    if(it.get_node1() == new_request.get_node1() && it.get_node2() == new_request.get_node2() == it.get_node3() == new_request.get_node3()){
                        check_no_repeat = false;
                        break;
                    }
                }
                if(check_no_repeat==true){
                    int node1 = new_request.get_node1(), node2 = new_request.get_node2(), node3 = new_request.get_node3();

                    //------------------setting value to requests------------------
                    double willness=0;
                    while(willness<1.0 || willness>3.0){
                        willness = normal_distribution(generator2);
                    }
                    //new_request.set_willness(willness);                            //[Wait1]
                    //cout<<"willness:"<<new_request.get_willness()<<endl;

                    //--------------------read graph-------------------------
                    ifstream graph_input;
                    int num_of_node_copy;
                    int num_of_edge_copy;
                    graph_input.open (filename);
                    graph_input >> num_of_node_copy;
                    vector<vector<int>> neighbor;
                    neighbor.resize(num_of_node_copy);
                    vector<pair<double,double>> pos(num_of_node_copy);
                    double pos_x, pos_y, new_swap_prob, new_fusion_prob;
                    int memory_cnt;
                    for(int i = 0; i < num_of_node; i++){
                        graph_input >> pos_x >> pos_y >> memory_cnt >> new_swap_prob >> new_fusion_prob;
                        pos[i] = {pos_x, pos_y};
                    }
                    int node_id1, node_id2; int channel_cnt; double fidelity; double dis_sum = 0; double prob_sum = 0;
                    graph_input >> num_of_edge_copy;
                    for(int i = 0;i < num_of_edge_copy; i++){
                        graph_input >> node_id1 >> node_id2 >> channel_cnt >> fidelity;
                        neighbor[node_id1].emplace_back(node_id2);
                        neighbor[node_id2].emplace_back(node_id1);
                    }
                    //--------------Dijastra by hops---------------------
                    double x1,x2,x3,xm,y1,y2,y3,ym;
                    tie(x1,y1) = pos[new_request.get_node1()];
                    tie(x2,y2) = pos[new_request.get_node2()];
                    tie(x3,y3) = pos[new_request.get_node3()];
                    double smallest_fermat = numeric_limits<double>::infinity();
                    int smallest_middle = -1;
                    for(int middle = 0; middle < num_of_node_copy; middle ++){
                        if(middle == new_request.get_node1() || middle == new_request.get_node2() || middle == new_request.get_node3()){continue;}
                        tie(xm,ym) = pos[middle];
                        double cur_fermat = pow((x1-xm),2) + pow((y1-ym),2) + pow((x2-xm),2) + pow((y2-ym),2) + pow((x3-xm),2) + pow((y3-ym),2);
                        if( cur_fermat < smallest_fermat){
                            smallest_fermat = cur_fermat;
                            smallest_middle = middle;
                        }
                    }
                    vector<bool> used( num_of_node_copy, false);
                    vector<int> parent( num_of_node_copy, -1);
                    vector<int> dist(num_of_node_copy, num_of_node_copy+1);
                    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
                    pq.push({0, smallest_middle});
                    dist[smallest_middle] = 0;
                    while(!pq.empty()) {
                        int cur_node = pq.top().second;
                        pq.pop();
                        if(used[cur_node]) continue;
                        used[cur_node] = true;
                        if(used[node1] && used[node2] && used[node3]) break;
                        for(int &neigh : neighbor[cur_node]) {
                            if(dist[cur_node] + 1 < dist[neigh]) {
                                dist[neigh] = dist[cur_node] + 1;
                                parent[neigh] = cur_node;                                        
                                pq.push({dist[neigh], neigh});
                            }
                        }
                    }
                    if(dist[node1] == (num_of_node_copy+1) || dist[node2] == (num_of_node_copy+1) || dist[node3] == (num_of_node_copy+1) ){
                        cout<<"No way to give value\n";
                        break;
                    }
                    int total_hops = dist[node1]+dist[node2]+dist[node3];
                    new_request.set_value((total_hops + total_hops * 1/*swap_ratio*/) * willness);
                    request_list[T].push_back(new_request);
                    graph_input.close();
                }   
            }while(check_no_repeat==false);
        }
    }

    for(int change_index = 0;change_index < change_parameter["swap_prob"].size();change_index++) {     //對於不同的參數數值，跑同round的圖與request_list
        time_t now;
        char* dt;
        double swap_prob = change_parameter["swap_prob"][change_index];
        min_swap_prob = swap_prob - 0.05;
        max_swap_prob = swap_prob + 0.05;
        #pragma omp parallel for
        for(int T = 0; T < round; T++){
            string round_str = to_string(T);
            string filename = file_path + "input/round_" + round_str + ".input";
            string command = "python3 alter_graph.py ";
            string parameter = to_string(min_swap_prob) + " " +  to_string(max_swap_prob);
            if(system((command + filename + " " + parameter).c_str()) != 0){
                cerr<<"error:\tsystem proccess python error"<<endl;
                exit(1);
            }  
        }
      
        now = time(0);
        dt = ctime(&now);
        #pragma omp parallel for
        for(int T = 0; T < round; T++){
            ofstream ofs;
            string round_str = to_string(T);
            cout<<"swap_prob:"<<swap_prob<<" round:"<<round_str<<"-------------------------"<<endl;
            ofs.open(file_path + "log/" + "swap_prob" + "_in_" + to_string(swap_prob) + "_Round_" + round_str + ".log");
            cerr  << "時間 " << dt << endl << endl; 
            ofs  << "時間 " << dt << endl << endl; 
            
            string filename = file_path + "input/round_" + round_str + ".input";
            vector<AlgorithmBase*> algorithms;
            algorithms.emplace_back(new MyAlgo(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha));
            algorithms.emplace_back(new Greedy(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0));
            algorithms.emplace_back(new BaselineDual(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0));
            algorithms.emplace_back(new BaselineAllNode(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0));

            ofs<<"---------------in round " <<T<<" -------------" <<endl;
                for(int t = 0; t < total_time_slot; t++){
                    ofs<<"---------------in timeslot " <<t<<" -------------" <<endl;
                    cout<< "---------generating requests in main.cpp----------" << endl;
                    for(auto &req:request_list[T]){
                        if(t >= service_time){break;}
                        Request new_request = req;
                        //cout<<"[REQUEST] : "<<new_request.get_node1()<<" "<<new_request.get_node2()<<" "<<new_request.get_node3()<<endl;
                        for(auto &algo:algorithms){
                            result[T][algo->get_name()]["total_request"]++; 
                            algo->add_new_request(new_request);
                        }
                    }
                    cout<< "---------generating requests in main.cpp----------end" << endl;
                    
                    #pragma omp parallel for 
                    for(int i = 0; i < (int)algorithms.size(); i++){
                        auto &algo = algorithms[i];
                        ofs<<"-----------run "<< algo->get_name() << " ---------"<<endl;
                        algo->run();
                        ofs<<"total_earn : "<<algo->get_res("total_earn")<<endl;
                        ofs<<"-----------run "<<algo->get_name() << " ---------end"<<endl;
                    }
                    
                }
                ofs<<"---------------in round " <<T<<" -------------end" <<endl;
                ofs << endl;
                for(auto &algo:algorithms){
                    ofs<<"("<<algo->get_name()<<")total earn = "<<algo->get_res("total_earn")<<endl;
                }
                cout<<"---------------in round " <<T<<" -------------end" <<endl;
                cout << endl;
                for(auto &algo:algorithms){
                    cout<<"("<<algo->get_name()<<")total earn = "<<algo->get_res("total_earn")<<endl;
                }
                
                for(auto &algo:algorithms){
                    for(string Y_name : Y_names) {
                        result[T][algo->get_name()][Y_name] = algo->get_res(Y_name);
                        if(Y_name == "total_earn" && (algo->get_name() == "MyAlgo") ){
                            result[T][algo->get_name()]["UB"] = algo->get_res("UB");
                            cout<<"total_earn : "<<algo->get_res("UB")<<endl;
                        }
                    }
                }

                now = time(0);
                dt = ctime(&now);
                cerr  << "時間 " << dt << endl << endl; 
                ofs  << "時間 " << dt << endl << endl; 
                ofs.close();
            
                for(auto &algo:algorithms){
                    delete algo;
                }
                algorithms.clear();
        }
        
        map<string, map<string, double>> sum_res;
            for(string algo_name : algo_names){
                for(int T = 0; T < round; T++){
                // result[T][algo_name]["waiting_time"] /= result[T][algo_name]["total_request"];
                // result[T][algo_name]["encode_ratio"] = result[T][algo_name]["encode_cnt"] / (result[T][algo_name]["encode_cnt"] + result[T][algo_name]["unencode_cnt"]);
                // result[T][algo_name]["succ-finished_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["finished_throughputs"];
                // result[T][algo_name]["fail-finished_ratio"] = 1 - result[T][algo_name]["succ-finished_ratio"];
                // result[T][algo_name]["path_length"] = result[T][algo_name]["path_length"] / result[T][algo_name]["finished_throughputs"];
                // result[T][algo_name]["divide_cnt"] = result[T][algo_name]["divide_cnt"] / result[T][algo_name]["finished_throughputs"];
                //  result[T][algo_name]["use_memory_ratio"] = result[T][algo_name]["use_memory"] / result[T][algo_name]["total_memory"];
                //  result[T][algo_name]["use_channel_ratio"] = result[T][algo_name]["use_channel"] / result[T][algo_name]["total_channel"];
                //  result[T][algo_name]["throughput_memory_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["use_memory"];
                //  result[T][algo_name]["throughput_channel_ratio"] = result[T][algo_name]["throughputs"] /  result[T][algo_name]["use_channel"];
                }
            }
        for(int T = 0; T < round; T++){
            // result[T]["MyAlgo"]["diff_rate"] = result[T]["MyAlgo"]["change_edge_num"] / result[T]["MyAlgo"]["diff_edge_num"];
            // result[T]["MyAlgo"]["edge_difference"] = result[T]["MyAlgo"]["change_edge_num"] - result[T]["MyAlgo"]["change_edge_num"];
        }

        double min_UB;
        for(int T = 0; T < round; T++){
            cout<<result[T]["MyAlgo"]["UB"]<<endl;
            min_UB=result[T]["MyAlgo"]["UB"];
            sum_res["MyAlgo"]["UB"] += min_UB;
        }
        
        for(string Y_name : Y_names) {
            string filename = "ans/swap_prob_" + Y_name + ".ans";
            ofstream ofs;
            ofs.open(file_path + filename, ios::app);
            ofs << change_parameter["swap_prob"][change_index] << ' ';
            
            for(string algo_name : algo_names){
                for(int T = 0; T < round; T++){
                    sum_res[algo_name][Y_name] += result[T][algo_name][Y_name];
                }
                ofs << sum_res[algo_name][Y_name] / round << ' ';
            }
            if(Y_name == "total_earn"){
                ofs << sum_res["MyAlgo"]["UB"] / round << " ";
            }
            ofs << endl;
            ofs.close();
        }

        
        // for(string Y_name : Y_names){
            
        //     string filename = "ans/" + X_name + "_" + Y_name + "_before_ent_path_prob_vt.ans";
        //     ofstream ofs;
        //     ofs.open(file_path + filename, ios::app);
        //     ofs << change_value << endl;
            
        //     for(string algo_name : algo_names){
        //         ofs<<algo_name<<endl;
        //         sort(sum_vt[algo_name][Y_name].begin(),sum_vt[algo_name][Y_name].end());
        //         for(auto it:sum_vt[algo_name][Y_name]){
        //             ofs << it << " ";
        //         }
        //         ofs << endl;
        //     }
        //     ofs << endl;
        //     ofs.close();
        // }
        
    }
    return 0;
}