#include <iostream>
#include <queue>
#include <map>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

int main(void){
    map<vector<vector<int>>,double> test;
    vector<vector<int>>t1 = {{0,1},{2,3},{4,5}};
    test[t1] = 10;
    for(auto it:test){
        cout<<it.second<<endl;
        for(auto it2:it.first){
            for(auto it3:it2){
                cout<<it3<<" ";
            }
            cout<<endl;
        }
    }
}