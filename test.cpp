#include <iostream>
#include <queue>
#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

int main(void){
    ifstream ifs;
    ifs.open("src/test/ans/num_of_node_fusion_vt.ans");
    if (!ifs) {
        cerr << "File open error!" << endl;
        return 1;
    }
    vector<double> area = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    vector<vector<double>>array;
    for (int algo = 0; algo < 4; algo++) {
        int total_num = 0;
        vector<float> count(11, 0);
        double value;
        while (ifs >> value && value != -1) {
            cout << value<<endl;
            for (int i = 0; i < area.size(); i++) {
                if (value <= area[i]) {
                    count[i] += 1;
                    break;
                }
            }
            total_num++;
        }

        for (int i = 0; i < count.size(); i++) {
            count[i] = count[i] / total_num;
            //cout<<count[i]<<"-";
        }
        //cout<<endl;
        vector<double>temp;
        temp.push_back(count[0]);
        for (int i = 1; i < count.size(); i++) {
            count[i] += count[i - 1];
            temp.push_back(count[i]);
        }
        array.push_back(temp);
    }
    for(int i = 0; i < 11; i++){
        cout<<area[i]<<" ";
        for(int j = 0; j < 4; j++){
            cout<<array[j][i]*100<<" ";
        }
        cout<<endl;
    }
    
    ifs.close();
}