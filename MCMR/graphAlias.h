//
// Created by 林丹丹 on 2018/5/20.
//

#ifndef GRAPHALIAS_H
#define GRAPHALIAS_H

#include <random>
#include <algorithm>
#include <stack>
#include <iostream>
#include "Random.h"
using namespace std;

class graphAlias{
public:
    double** p;
    int** h;
    int** map1;
    //int** map2;
    int* outdegree;
    int totalNodes;
    graphAlias(){
        totalNodes =0;
    }
    void inputGraphAlias(vector<vector<pair<int, double> >> outAdjLists, int numNodes){
        totalNodes = numNodes;
        p = new double*[numNodes];
        h = new int*[numNodes];
        map1 = new int*[numNodes];
        outdegree = new int[numNodes];

        for(int i = 0; i < numNodes; i++){
            //cout << "Processing node" << i<<endl;
            vector<pair<int, double> > outAdjList = outAdjLists[i];

            double sum =0;
            int size = outAdjList.size();
            outdegree[i] = size;
            //cout << "alias size" << outdegree[i]<< endl;

            stack<int> small;
            stack<int> big;

            p[i] = new double[size];
            h[i] = new int[size];
            map1[i] = new int[size];
            for(int j =0; j < size; j++){
                sum += outAdjList[j].second;
                map1[i][j] = outAdjList[j].first;
            }
            //cout << "sum of weights: " << sum << endl;

            for(int j =0; j < size; j++){
                p[i][j] = outAdjList[j].second * size / sum;
                if(p[i][j] > 1){
                    big.push(j);
                }
                else{
                    small.push(j);
                }
            }

            while(!small.empty() && ! big.empty()){
                int smallId = small.top();
                small.pop();
                int bigId = big.top();
                h[i][smallId] = bigId;
                p[i][bigId] -= (1-p[i][smallId]);
                if(p[i][bigId] < 1){
                    small.push(bigId);
                    big.pop();
                }
            }

            if(big.size()>0){
                while(!big.empty()){
                    int bigId = big.top();
                    h[i][bigId] = bigId;
                    big.pop();
                }
            }
            if(small.size()>0){
                while(!small.empty()){
                    int smallId = small.top();
                    h[i][smallId] = smallId;
                    small.pop();
                }
            }
/*
            for(int j=0; j < size; j++){
                cout<< "node " <<j<<" : " << p[i][j] <<endl;
                cout<< "mapped node from " << j << ": " << map1[i][j] << endl;
                cout<< "after hashing from " << j << ": " << h[i][j]<< endl;
                cout<< "mapped node after hashing from "<< j << ": " << map1[i][h[i][j]] << endl;
            }
*/
        }


    }

    ~graphAlias(){
        for(int i =0; i < totalNodes; i++){
            delete [] p[i];
            delete[] h[i];
            delete[] map1[i];
        }
        delete[] p;
        delete[] h;
        delete[] map1;
        delete[] outdegree;
    }
    int generateRandom(int vert, Random& R){
        //cout<< "Generating Random Node" <<endl;
        int firstId = R.drand() * outdegree[vert];
        //cout << "first ID: " << firstId << " p[firstID]" << p[vert][firstId]<< endl;
        int secondId = R.drand() < p[vert][firstId] ? map1[vert][firstId] : map1[vert][h[vert][firstId]];

        return secondId;
    }
    int generateRandom_t(int vert, Random& R){
        int firstId = R.drand_t() * outdegree[vert];
        int secondId = R.drand_t() < p[vert][firstId] ? map1[vert][firstId] : map1[vert][h[vert][firstId]];
        return secondId;
    }
};
#endif //TOPPPR_GRAPHALIAS_H
