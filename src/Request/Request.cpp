#include"Request.h"

Request::Request(int node1, int node2, int node3):node1(node1),
    node2(node2), node3(node3), status(REQUEST_UNFINISHED){
    if(DEBUG)cerr<<"new Request"<<endl;
}

Request::~Request(void){
    if(DEBUG)cerr<<"delete Request"<<endl;
    for(int i=0;i<(int)trees.size();i++){
        for(int j=0; j<3; j++){
            trees[i][j]->release();
            delete trees[i][j];
            trees[i][j] = nullptr;
        }
        trees[i].clear();
    }
    trees.clear();
}

void Request::set_path(int path_id, int edge_id,Path *p){
    if(path_id >= (int)trees.size() || edge_id >= 3 ){
        cerr<<"can't set this path!"<<endl;
        exit(1);
    }
    trees[path_id][edge_id] = p;
}

int Request::get_node1(){ 
    return node1;
}
int Request::get_node2(){
    return node2;
}
int Request::get_node3(){
    return node3;
}

int Request::get_throughput(){
    return this->throughput;
}

int Request::get_tree_num(){
    return tree_num;
}

int Request::get_cur_send(){
    return cur_send;
}

void Request::add_cur(int num){
    cur_send += num;
}

vector<vector<Path *>> Request::get_trees(){
    return trees;
}

void Request::clear_trees(){
    for(int i=0;i < (int)trees.size();i++){
        for(int j=0; j<3;j ++){
            if(trees[i][j] != nullptr){
                delete trees[i][j];
                trees[i][j] = nullptr;
            }   
        }
        trees[i].clear();
    }
    trees.clear();
}

void Request::refresh_trees(){
    for(auto &tree : trees){
        for(Path *path : tree){
            path->refresh();
        }
    }
    status = REQUEST_UNFINISHED;
}

void Request::entangle(){
    for(auto &tree:trees){
        for(auto &path:tree){
            path->entangle(); 
        }
    }
}

void Request::swap(){
    for(auto &tree:trees){
        throughput ++;
    }
}

void Request::send(){

}

bool Request::is_finished(){
    return status != REQUEST_UNFINISHED;
}
bool Request::is_success(){
    if(status == REQUEST_UNFINISHED){
        cerr<<"the request is unfinshed!"<<endl;
        exit(1);
    }
    return status == REQUEST_SUCC;
}

void Request::next_timeslot(){
    for(auto &tree:trees){
        for(auto &path:tree){
            if(path != nullptr){
                path->release();
                delete path;
            }
        }
        tree.clear();
    }
    trees.clear();
}

void Request::operator+=(vector<Path *>tree){
    trees.emplace_back(tree);
}



