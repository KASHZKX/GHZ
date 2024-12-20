#include <iostream>
#include <algorithm>
#include <omp.h>
#include <queue>
#include <vector>
#include <map>
#include <limits>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <random>
using namespace std;

int num_of_node;
vector<vector<int>> neightbor;

class Label {
public:
  double theta, obj_bar2;
  int i;

  Label(double t, double o, int i) : theta(t), obj_bar2(o), i(i){};

    bool operator>(const Label &b) const {
    if (theta == b.theta)
      return obj_bar2 > b.obj_bar2;
    return theta > b.theta;
  }
};

vector<int> Dijkstra(int src, int dst, map<pair<int,int>,double>&X){ 
    const double INF = numeric_limits<double>::infinity();
    //const double INF = 100;
    int n = num_of_node;
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
        for(int neighbor : neightbor[cur_node]) {
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

vector<int> separation_oracle1(int src,int dst , map<pair<int,int>,double>&X, map<pair<int,int>,double>&Y){     
    vector<int> SPT;                  //nodes' parent in the spanning tree
    vector<int> best_path;
    double best_len; 
    
    SPT = Dijkstra(src, dst,X);                               //the first SPT is get by dijkstra
    int cur_node = src;                                     //counting the first path's U(X,Y)=c* e^r
    double c = 0;                                           //c=SUM[u,v]:alpha(u)+alpha(v)+beta(u,v)==X[u,v]
    double r = 0;                                           //r=SUM[u,v]:-ln(Pr(u,v))==Y[u,v]
    while(cur_node != dst){
        // if(cur_node < SPT[cur_node]){                       //[can alter]no need if else
        //     c += X[{cur_node,SPT[cur_node]}][req_no];               
        //     r += Y[req_no][{cur_node,SPT[cur_node]}];
        // }
        // else{
        //     c += X[{SPT[cur_node],cur_node}][req_no];  
        //     r += Y[req_no][{SPT[cur_node],cur_node}];           
        // }
        c += X[{SPT[cur_node], cur_node}];  
        r += Y[{SPT[cur_node],cur_node}];  
        best_path.push_back(cur_node);
        cur_node = SPT[cur_node];
    } 
    best_path.push_back(dst);
    best_len = c * exp(r);


    cout << "V1_First_Path: ";
    for(auto p : best_path){
            cout << p << " ";
    }
    cout <<"with "<<best_len<< endl;

    map<pair<int, int>, bool> used_edge;
    vector<int> new_path;   
    pair<int,int> new_edge;

    for(unsigned int i = 0; i < SPT.size(); i++){
        int cur_node=i;
        while(cur_node!=dst){
            if(used_edge.find({cur_node,SPT[cur_node]})!=used_edge.end()){
                break;
            }
            used_edge[{cur_node,SPT[cur_node]}] = true;
            used_edge[{SPT[cur_node],cur_node}] = true;
            cur_node=SPT[cur_node];
        }
    }

    while(1){
        double minimum = numeric_limits<double>::infinity();
        for(int i = 0; i < num_of_node; i++){                 //creating many new SPT
            vector<int> neighbors = neightbor[i];  
            for(auto neighbor : neighbors){
                double temp1 = 0, temp2 = 0;
                if(SPT[i] == neighbor || SPT[neighbor] == i){      // checking used edge or unused
                    continue;   
                }else{                                             // if unused
                    temp1 = X[{i, neighbor}];
                    temp2 = Y[{i, neighbor}];
                    int cur_node = i;
                    while(cur_node != dst){
                        temp1 += X[{cur_node, SPT[cur_node]}];
                        temp2 += Y[{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    } 

                    cur_node = neighbor;
                    while(cur_node != dst){
                        temp1 -= X[{cur_node, SPT[cur_node]}];
                        temp2 -= Y[{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    }       
                    if(temp2 < 0 && temp1 > 0 ) {               // we need the smallest edge to change the SPT
                        if(i<neighbor){
                            if(used_edge.find({i, neighbor}) != used_edge.end()  && used_edge[{i, neighbor}]!=false){
                                continue;
                            }
                        }
                        else{
                            if(used_edge.find({neighbor ,i}) != used_edge.end() && used_edge[{neighbor,i}]!=false){
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
        used_edge[new_edge] = true;
        used_edge[{new_edge.second,new_edge.first}] = true;
        used_edge[{new_edge.second,SPT[new_edge.second]}] = false;
        used_edge[{SPT[new_edge.second],new_edge.second}] = false;

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
                c += X[{new_path[i], new_path[i+1]}];
                r += Y[{new_path[i], new_path[i+1]}];
            }
            else{
                c += X[{new_path[i+1], new_path[i]}];  
                r += Y[{new_path[i+1], new_path[i]}];           
            }
            //cout<<"PATH:["<<new_path[i]<<" "<<new_path[i+1]<<"] with tot "<< c <<" / " << r <<endl;
        }
        new_len =  c * exp(r);
        if(new_len < best_len){
            best_len = new_len;
            best_path = new_path;                                            //路線修改,新的spt產生
        } 
    }
    for(auto it:best_path){
        cout<<it<<" ";
    }
    cout<<"with "<<best_len;
    return best_path;                                             
}

//----------------------------------------------------------------
vector<int> Dijkstra2(int src, int dst, map<pair<int,int>,double>&X, map<pair<int,int>,double>&Y, vector<pair<double,double>>&dist){                           
    double INF=numeric_limits<double>::infinity();
    vector<bool> used(num_of_node, false);
    vector<int> parent(num_of_node, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    dist[src] = {0,0};                                                     //找src到所有點的最短距離(針對X-obj)
    pq.push({0, src});
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int neigh : neightbor[cur_node]) {
            if(dist[cur_node].first + X[{cur_node, neigh}] < dist[neigh].first) {
                dist[neigh].first = dist[cur_node].first + X[{cur_node, neigh}];      //(1)store d
                dist[neigh].second = dist[cur_node].second + Y[{cur_node, neigh}];
                parent[neigh] = cur_node;                                             //(1)store pred
                pq.push({dist[neigh].first, neigh});
            }
            if(dist[cur_node].first + X[{cur_node, neigh}] == dist[neigh].first && dist[neigh].second > (dist[cur_node].second + Y[{cur_node, neigh}])) {
                dist[neigh].first = dist[cur_node].first + X[{cur_node, neigh}];      //(1)store d
                dist[neigh].second = dist[cur_node].second + Y[{cur_node, neigh}];
                parent[neigh] = cur_node;                                             //(1)store pred
                pq.push({dist[neigh].first, neigh}); 
            }
        }
    }

    if(dist[dst].first >= INF) return{};
    return parent;
}   

vector<int> separation_oracle2(int src, int dst, map<pair<int,int>,double>X, map<pair<int,int>,double>Y){
    double INF=numeric_limits<double>::infinity();
    vector<pair<double,double>>dist (num_of_node,{INF,INF});            
    vector<int> pred (num_of_node,-1);                                  
    vector<int> SDpath;
    vector<int> best_path;
    double smallest_U = INF;
    vector<double> theta_table(num_of_node,INF);                        
    vector<double> obj_bar1_table(num_of_node,INF);                     
    vector<double> obj_bar2_table(num_of_node,INF);
    vector<int> cpred_table(num_of_node,-1);                            
    priority_queue<Label,vector<Label>, greater<Label>> pq;                                          
    pred = Dijkstra2(src, dst, X, Y, dist);                             //(1)

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
    for(int i = 0; i < num_of_node; i++){                               //(4)compute lexmin
        if(i == src){continue;}
        for(auto &neigh:neightbor[i]){                                    
            double cur_obj_bar2 = dist[neigh].second + Y[{neigh,i}] - dist[i].second;
            if(cur_obj_bar2 < 0){
                double cur_theta = -(dist[neigh].first + X[{neigh,i}] - dist[i].first) / cur_obj_bar2 ;
                if(cur_theta <= theta_table[i]){
                    if(cur_theta < theta_table[i]){
                        theta_table[i] = cur_theta;
                        obj_bar1_table[i] = dist[neigh].first + X[{neigh,i}] - dist[i].first;
                        obj_bar2_table[i] = cur_obj_bar2;
                        cpred_table[i] = neigh;
                    }
                    else if(cur_obj_bar2 < obj_bar2_table[i]){
                        theta_table[i] = cur_theta;
                        obj_bar1_table[i] = dist[neigh].first + X[{neigh,i}] - dist[i].first;
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

        // for(auto i1:dist){
        //     cout<<"Dist:"<<i1.first<<","<<i1.second<<"   ";
        // }
        // cout<<endl;

        // for(int i=1;i<num_of_node;i++){
        //     cout<<"Cpred,pred:"<<cpred_table[i]<<","<<pred[i]<<"   ";
        // }
        // cout<<endl;

        // for(auto i2:theta_table){
        //     cout<<"theta:"<<i2<<"   ";
        // }
        // cout<<endl<<endl;


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
                cout << "[" << last_ratio << "," << min_theta << "] with " << d1 * exp(d2) << endl;
                last_ratio = min_theta;
                cout<<"Extremepath: ";
                for(auto it:SDpath){
                    cout << it <<" ";
                }
                cout<<endl<<endl;
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
        for(auto neigh:neightbor[min_i]){                             
            double cur_obj_bar2 = dist[neigh].second + Y[{neigh,min_i}] - dist[min_i].second;
            if(cur_obj_bar2 < 0){
                double cur_theta = -(dist[neigh].first + X[{neigh,min_i}] - dist[min_i].first) / cur_obj_bar2 ;
                if(cur_theta <= theta_table[min_i]){
                    if(cur_theta < theta_table[min_i]){
                        theta_table[min_i] = cur_theta;
                        obj_bar1_table[min_i] = dist[neigh].first + X[{neigh,min_i}] - dist[min_i].first;
                        obj_bar2_table[min_i] = cur_obj_bar2;
                        cpred_table[min_i] = neigh;
                    }
                    else if(cur_obj_bar2 < obj_bar2_table[min_i]){
                        theta_table[min_i] = cur_theta;                                                                //no use just for read
                        obj_bar1_table[min_i] = dist[neigh].first + X[{neigh,min_i}] - dist[min_i].first;
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
        for(auto neigh:neightbor[min_i]){
            double temp_obj_bar2 = dist[min_i].second + Y[{min_i,neigh}] - dist[neigh].second;
            if(temp_obj_bar2 < 0){
                double temp_obj_bar1 = dist[min_i].first + X[{min_i,neigh}] - dist[neigh].first;
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
    for(auto it:SDpath){
        cout<<it<<" ";
    }
    cout<<"with "<< d1 * exp(d2) <<endl;
    return best_path;
}

int main(){
    ofstream ofs;
    ofs.open("Quantum/Exp_output.txt",ios::app);
    ofs<<"\nNew round\n";   
    ofs.close(); 
    ifstream ifs;
    ifs.open("Quantum/data/input/round_7.input");
    double edge_num,g1,g2,g3,g4;
    ifs >> num_of_node;
    neightbor.resize(num_of_node);
    for(int i=0;i<num_of_node;i++){
        ifs>>g1>>g2>>g3>>g4;
        
    }
    ifs >> edge_num;
    for(int i=0;i<edge_num;i++){
        ifs>>g1>>g2>>g3>>g4;
        neightbor[g1].push_back(g2);
        neightbor[g2].push_back(g1);
    }
    ifs.close();
    map<pair<int,int>,double>X;
    map<pair<int,int>,double>Y;

    random_device rd;  
    mt19937 gen(rd()); 
    uniform_real_distribution<double> dis1(0.1, 5.0);
    uniform_real_distribution<double> dis2(0.1, 1.0);
    for(int i=0;i<num_of_node;i++){
        for(auto it:neightbor[i]){
            if(i<it){
                X[{i,it}] = dis1(gen);
                X[{it,i}] = X[{i,it}];
                Y[{i,it}] = dis2(gen);
                Y[{it,i}] = Y[{i,it}];
            }
        }
    }
    // for(auto it:X){
    //     cout<<it.first.first<<"-"<<it.first.second<<":"<<it.second<<" | "<<Y[it.first]<<endl;
    // }
    #pragma omp parallel for
    for(int round = 0;round <1;round ++){
        default_random_engine generator = default_random_engine(rd());
        uniform_int_distribution<int> unif(0, num_of_node-1);
        int src=unif(generator);
        int dst;
        do{
            dst=unif(generator);
        }while(src == dst);
        cout<<"\n\n S:"<<src<<">>>>>>>>>>>>>>>>>>>>"<<"D:"<<dst<<endl;
        vector<int>path1=separation_oracle1(src,dst,X,Y);
        cout<<"\n------------------------------------\n";
        vector<int>path2=separation_oracle2(src,dst,X,Y);
        cout<<"Ans : ";
        for(auto it:path2){
            cout<<it<<" ";
        }
        if(path1 != path2){
            ofstream ofs;
            ofs.open("Quantum/Exp_output.txt",ios::app);
            ofs<<"1";
            cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
            ofs.close();
        }
    }
    


}