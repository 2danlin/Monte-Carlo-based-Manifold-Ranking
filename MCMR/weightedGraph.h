//
// Created by 林丹丹 on 2018/5/14.
//

#ifndef WEIGHTEDGRAPH_H
#define WEIGHTEDGRAPH_H

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <utility>
#include <iostream>
#include <fstream>
#include "alias.h"
#include "graphAlias.h"
#include "Random.h"
using namespace std;

class weightedGraph{

public:
    int n; // # of nodes
    long m; // # of edges

    vector<vector<pair<int, double> >> inAdjList; // int: the in neighbor, double: the weight
    vector<vector<pair<int, double> >> outAdjList; // int: the out neighbor, double: the weight
    vector<vector<int>> outNeiList;

    graphAlias vertex_alias; //vertex_alias[0]: the alias of the first node for picking up the next node during the random walk

    double* sumWeights; //sumWeights[0]: the sum of weights of the first node

    int* indegree; // indegree[0]: the in-degree of the first node
    int* outdegree; // outdegree[0]: the out-degree of the first node
    double* maxWeight; // maxweight[0]: the max value of the weights of all the out-neighbors of the first node
    double weight_scalar;

    weightedGraph(){
        n=m=0;
    }

    ~weightedGraph(){

        delete sumWeights;
        delete indegree;
        delete outdegree;
        delete maxWeight;

        //vector<Alias>().swap(vertex_alias);
       // vertex_alias.~graphAlias();

    }

    void inputWeightedGraph(string filename){
        ifstream infile(filename.c_str());
        infile >> n;
        cout << "n= " << n << endl;

        m=0;

        cout << " Initializing the data containers:" <<endl;
        inAdjList = vector<vector<pair<int, double> >>(n, vector<pair<int, double>>()); //initiallize the vector
        //cout << " InAdjList finished. " <<endl;
        outAdjList = vector<vector<pair<int, double> >>(n, vector<pair<int, double>>());
        outNeiList = vector<vector<int> >(n,vector<int>());
        //cout << " InAdjList finished. " <<endl;
        //vertex_alias = vector<Alias>(n, Alias(vector<pair<int, double> >()));

        sumWeights = new double[n];
        maxWeight = new double[n];
        indegree = new int[n];
        outdegree = new int[n];
        for(int i=0; i < n; i++){
            sumWeights[i] =0.0;
            maxWeight[i] = 0.0;
            indegree[i] = 0;
            outdegree[i] = 0;
        }
        weight_scalar =  1.0;

        int fromNode, toNode;
        double weight;
        clock_t t0 = clock();


        while(infile.good()){
            infile >> fromNode >> toNode >> weight;
            //cout << "(" << fromNode << ", " << toNode << "): " << weight <<endl;
            weight = weight / weight_scalar;
            if(fromNode>= 0 && fromNode < n && fromNode != toNode && toNode >=0 && toNode <n) {

                //if this edge haven't been inserted

                if ((find(outNeiList[fromNode].begin(), outNeiList[fromNode].end(), toNode) == outNeiList[fromNode].end())) {
                    //cout<<"here"<<endl;
                    outAdjList[fromNode].push_back(make_pair(toNode, weight));
                    outAdjList[toNode].push_back(make_pair(fromNode, weight));
                    outNeiList[fromNode].push_back(toNode);
                    outNeiList[toNode].push_back(fromNode);
                    sumWeights[fromNode] += weight;
                    sumWeights[toNode] += weight;
                    outdegree[fromNode]++;
                    outdegree[toNode]++;
                    indegree[toNode]++;
                    indegree[fromNode]++;
                    if (weight > maxWeight[fromNode])
                        maxWeight[fromNode] = weight;
                    if (weight > maxWeight[toNode])
                        maxWeight[toNode] = weight;

                    m = m + 2;

                }
            }

        }
        infile.close();

        vertex_alias.inputGraphAlias(outAdjList, n);

        clock_t t1 = clock();
        cout << "m = :" << m << endl;
        cout << "read file time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;

    }

    int getInSize(int vert){
        return indegree[vert];
    }

    int getOutSize(int vert)
    {
        return outdegree[vert];
    }
    int getOutVert(int vert, int pos){
        /*if(vert == 0)
            return outEdge[pos];
        else
            return outEdge[outCount[vert-1] + pos];*/
        return outAdjList[vert][pos].first;
    }
    double getMaxWeight(int vert){
        return maxWeight[vert];
    }

    double getSumWeights(int vert){
        return sumWeights[vert];
    }
    void saveSumWeight(string filename){
        stringstream ss;
        ss << "data_preprocessing/sum_weight_"<<filename<<".txt";
        string outputFile = ss.str();
        //cout << "file: " << outputFile << endl;
        ofstream fout(outputFile);

        for(int i =0; i < n; i++){
            fout<< sumWeights[i]<<"\n";
        }

        fout.close();
        cout<< "done sum weight file ... " <<endl;
    }

    void saveNormarlizedWeight(string filename){
        stringstream ss;
        ss << "data_preprocessing/normarlized_"<<filename<<".txt";
        string outputFile = ss.str();
        //cout << "file: " << outputFile << endl;
        ofstream fout(outputFile);

        for(int i =0; i < n; i++){
            int size_i = outAdjList[i].size();
            double sum_weight_i = sumWeights[i];
            for(int j=0; j<size_i; j++){
                int neighbor = outAdjList[i][j].first;
                double weight = outAdjList[i][j].second;
                double nor_weight = weight / sqrt(sum_weight_i) / sqrt(sumWeights[neighbor]);
                fout<< i << " " << neighbor << " " << nor_weight <<endl;
            }
        }

        fout.close();
        cout<< "done normalized weight file ... " <<endl;
    }

    double getMaxSumWeights(){
        double maxSum = 0.0;
        for(int i =0; i < n; i++){
            if(maxSum < sumWeights[i]){
                maxSum = sumWeights[i];
            }
        }
        return maxSum;
    }

    int getNextNode(int vert, Random& R){

        int nextNode = vertex_alias.generateRandom(vert, R);
        return nextNode;
    }
private:
};


#endif //TOPPPR_WEIGHTEDGRAPH_H
