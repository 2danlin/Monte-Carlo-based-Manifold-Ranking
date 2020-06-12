//
// Created by 林丹丹 on 2018/5/17.
//

#ifndef MANIFOLDRANKING_H
#define MANIFOLDRANKING_H

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <fstream>
#include <future>
#include <string>
#include <sstream>
#include "weightedGraph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <sys/time.h>
#include <time.h>
#include <unordered_map>
#include "algo.h"
//#include "SparseMatrix.h"
#include <cstdlib>
#include <Eigen/Sparse>
using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<double> SparseMtx;
typedef Eigen::Triplet<double> T;

bool maxScoreCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

bool maxIdCmp(int a, int b){
    return a>b;
}
class MR{
    friend void mr_t_PowerMethod(MR* mr, vector<int> nodeList, int iterations);
public:
    double avg_pre;
    int avg_walk_num;
    double avg_time;
    int k;
    double forward_c;
    weightedGraph g;
    Random R;
    int vert;
    double alpha;
    string target_filename; //dataset name
    string graph_name;
    double* vert_count;
    int value_count;
    int* value_verts;
    unsigned NUM_CORES;
    int** multiVertCount;
    double* resultList;
    double back_walk_num;
    double avg_L1_error;
    double avg_max_error;
    double max_max_error;
    double avg_avg_error;
    //int* r_hash_arr;
    //bool* r_item_arr;
    //double* r_max_c;
    //double avg_forward_time;
    //double avg_rsum;
    //vector<vector<pair<int, double> > > r_vec;

    void PowerMethodMulti(int iterations, int node_count, int num_thread);
    //const static int NUMTHREAD = 20;
    Random* Rs;
    vector<stack<pair<int, double> >> random_walk_storation;
    vector<stack<pair<int, double> >> random_walk_storation_one_iteration;
   //vector<vector<pair<int, double> > > random_walk_storation_one_iteration;


    MR(string full_filename, string graphname, double forward_para_c, int k_num, double input_alpha) {
        k = k_num;
        forward_c = forward_para_c;
        avg_L1_error = 0;
        avg_max_error = 0;
        max_max_error = 0;
        avg_avg_error = 0;
        avg_time = 0;
        avg_pre = 0 ;
        avg_walk_num =0;
        back_walk_num = 0;
        target_filename = full_filename;
        graph_name = graphname;
        string filename = "dataset/" + full_filename + ".txt";
        if(!file_exists_test(filename)){
            cout<<"File does not exist: " << filename<<endl;
            exit(0);
        }
        g.inputWeightedGraph(filename);
        cout << "edge num: " << g.m << endl;
        vert = g.n;

        alpha = input_alpha;
        srand(unsigned(time(0)));
        R = Random(unsigned(rand()));
        vert_count = new double[vert];
        resultList = new double[vert];
        value_count = 0;
        value_verts = new int[vert];
        random_walk_storation = vector<stack<pair<int, double> > >(vert);
        random_walk_storation_one_iteration = vector<stack<pair<int, double>> >(vert);

        for(int i =0 ; i < vert; i++){
            resultList[i] = 0;
            vert_count[i] = 0;
            value_verts[i] = -1;
            //vector<pair<int, double> > temp_walk_storation = vector<pair<int, double> >(vert);
            //random_walk_storation_one_iteration.push_back(temp_walk_storation);
            //vector<stack<double>> random_walk_storation_i = vector<stack<double>>(vert);
            //random_walk_storation[i].push_back(random_walk_storation_i);
        }
        NUM_CORES = std::thread::hardware_concurrency();
        assert(NUM_CORES >= 2);
        cout << "thread core: " << NUM_CORES << endl;
        multiVertCount = new int*[NUM_CORES];
        Rs = new Random[NUM_CORES];
        for(int i = 0; i < NUM_CORES; i++){
            Rs[i] = Random(unsigned(rand()));
            multiVertCount[i] = new int[vert];
            for(int j = 0; j < vert; j++){
                multiVertCount[i][j] = 0;
            }
        }

        cout << "init done! " << endl;
    }

    ~MR() {
        for(int i = 0; i < NUM_CORES; i++){
            delete[] multiVertCount[i];

        }

        delete[] multiVertCount;

        delete[] vert_count;
        delete[] value_verts;
        delete[] Rs;

        //vector<vector<pair<int, double> >>().swap(random_walk_storation_one_iteration);
        delete[] resultList;

        //g.~weightedGraph();
        //cout<<"here"<<endl;
    }

    //generate random query node
    void generateQueryNode(int nodeNum, ofstream& fout){
        for(int i = 0; i < nodeNum; i++){
            int tempNode = R.generateRandom() % vert;
            if(g.getOutSize(tempNode) == 0){
                i--;
                continue;
            }
            fout << tempNode << endl;
        }
    }

    //取s点的groundtruth
    vector<int> getRealTopK(int s, int k){
        stringstream ss;
        ss << "real_mr/" << graph_name << "/" << s << ".txt";
        string infile = ss.str();
        ifstream real(infile);
        vector<int> realList;
        vector<double> simList;
        for(int i = 0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            if(i >= k){
                double dif = abs(tempSim - simList[k-1]);
                //cout<<"dif " << dif <<endl;
                if(dif != 0.0){
                   // cout<< "ignore this node: " << tempSim << " " <<simList[k-1] <<endl;
                    break;
                }

            }
            realList.push_back(tempId);
            simList.push_back(tempSim);

            //cout<< tempId <<" " <<tempSim <<endl;
        }
        real.close();
        return realList;
    }

    unordered_map<int, double> getRealTopKMap(int s, int k){
        unordered_map<int, double> answer_map;
        stringstream ss;
        ss << "real_mr/" << graph_name << "/" << s << ".txt";
        string infile = ss.str();
        ifstream real(infile);
        double k_Sim = 0;
        for(int i = 0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            if(i == k - 1){
                k_Sim = tempSim;
            }
            if(i >= k && tempSim < k_Sim){
                break;
            }
            answer_map[tempId] = tempSim;
        }
        real.close();
        return answer_map;
    }

    //MCMR+: forward + random, for top-k query
    void MCMRPlus(int s, double r_threshold, int k){
        //r_threshold = 0.5/log(k);
        cout << "node: " << s << endl;

        vector<int> realList = getRealTopK(s, k);
        int* forwardCandidate = new int[vert];
        double* pi_f; //reserve
        double* residue_f; //residue
        bool* isInQueue;
        double* betas; //Berstein Inequality: the bound
        pi_f = new double[vert];
        residue_f = new double[vert];
        betas = new double[vert];
        isInQueue = new bool[vert]; //the queue for forward pushing
        queue<int> r_queue;

        //the knn graph does not exist this kind of problem
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            delete[] pi_f;
            delete[] residue_f;
            delete[] betas;
            delete[] isInQueue;
            delete[] forwardCandidate;
            //return resultList;
        }

        //vector<pair<int, double> > random_walk_storation_one_iteration = vector<pair<int, double> >(vert, vector<stack<double>>(vert));

        double total_t1 = 0, total_t2 = 0, total_t3 = 0, total_t4 = 0;
        int upgrad_count = 0;
        int walk_num = 0;
        double r_sum = 0;

        int candidate_size = 0;
        int fora_candidate = 0;
        double upperBound = 0, lowerBound = 0;
        int fwd_walk_iteration = 0;
        while(true){
            candidate_size = 0;
            fora_candidate = 0;
            for(int i = 0; i < vert; i++){
                pi_f[i] = 0;
                residue_f[i] = 0;
                resultList[i] = 0;
                forwardCandidate[i] = 0;
                betas[i] = 0;
                isInQueue[i] = false;
                vert_count[i] = 0;
                value_verts[i] = 0;
            }
            residue_f[s] = 1;

            if(1 * g.getMaxWeight(s) > r_threshold){
                r_queue.push(s);
                isInQueue[s] = true;
            }

            clock_t t0 = clock();
            unsigned int total_walk_num =0;
            //the forward search process
            while(r_queue.size() > 0){
                int tempNode = r_queue.front();
                r_queue.pop();
                //cout << "current processing node: " << tempNode <<endl;
                isInQueue[tempNode] = false;
                int tempOutSize = g.getOutSize(tempNode); // it must be larger than or equal to k_nn
                double maxWeight = g.getMaxWeight(tempNode);
                double tempR = residue_f[tempNode]; //the residue of tempNode
                residue_f[tempNode] = 0; //the existing residue is set to be zero
                pi_f[tempNode] +=(1-alpha) * tempR; // MR setting: convert a proportion of residue into reserve of tempNode
                double remaining_residue = tempR * (alpha); //MR setting: transfer the remaining residues to its neighbours

                if(tempOutSize == 0){ //this case won't appear in MR due to the knn graph
                    residue_f[s] += remaining_residue;
                    if(g.getOutSize(s) > 0 && !isInQueue[s] && residue_f[s] * g.getMaxWeight(s)> r_threshold){
                        isInQueue[s] = true;
                        r_queue.push(s);
                    }
                }
                else{
                    vector<pair<int, double> > out_neighbors  = g.outAdjList[tempNode];
                    for(int i = 0; i < tempOutSize; i++){
                        int newNode = out_neighbors[i].first;
                        double weight_newNode = out_neighbors[i].second;
                        residue_f[newNode] += remaining_residue * weight_newNode;
                        if(!isInQueue[newNode] && residue_f[newNode] * g.getMaxWeight(newNode) > r_threshold){
                            isInQueue[newNode] = true;
                            r_queue.push(newNode);
                        }
                    }
                }
                if(resultList[tempNode] == 0){
                    forwardCandidate[fora_candidate++] = tempNode;
                }
                resultList[tempNode] = pi_f[tempNode];
            }

            vector<pair<int, double> > aliasP;
            r_sum = 0;
            for(int i = 0; i < fora_candidate; i++){
                int tempNode = forwardCandidate[i];
                if(residue_f[tempNode] > 0){
                    r_sum += residue_f[tempNode];
                    aliasP.push_back(pair<int, double>(tempNode, residue_f[tempNode]));
                }
            }
            cout << "rsum: " << r_sum << endl;

            Alias alias = Alias(aliasP); //based on the remaining residual, chosen to be source node of random walk
            clock_t ta = clock();
            walk_num = max(1000, (int) ((ta - t0) / (double) CLOCKS_PER_SEC * 5000 * log(k) * log(log(g.m))));
            cout << "walk num: " << walk_num << endl;
            total_walk_num += walk_num;
            double increment = r_sum / (double) walk_num; // increment is different for each random walk depends on the length of random walk
            double temp_count = 0;
            value_count = 0;
            //cout << "Starting Random Walks..." << endl;
            for(int i = 0; i < walk_num; i++){
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;//chosen to be source node of a random walk
                if(fwd_walk_iteration < 4){
                    //cout<< "First Time doing Random walk ..." << endl;
                    double scalar_prob = 1;
                    while(R.drand() < alpha){
                        // cout << "Sum of weights of this source node ... " << scalar_prob << endl;
                        tempNode = g.getNextNode(tempNode, R); // the next node be chosen
                        //tempNode = vertex_alias[tempNode].generateRandom(R);
                        scalar_prob = g.getSumWeights(tempNode);
                        //cout<< "Select the next node: " << tempNode <<endl;
                    }

                    if(vert_count[tempNode] == 0){
                        value_verts[value_count++] = tempNode;
                    }
                    if(resultList[tempNode] == 0){
                        forwardCandidate[fora_candidate++] = tempNode;
                    }
                    resultList[tempNode] += increment * scalar_prob; // due to the random generator
                    //cout << "Updated MR value of node " << tempNode << ": " << resultList[tempNode] <<endl;
                    vert_count[tempNode]++;
                    //cout << "Ending Walking ..." << endl;
                    pair<int, double> temp_pair_prob = make_pair(tempNode, scalar_prob);
                    random_walk_storation[sourceNode].push(temp_pair_prob);
                }
                else{
                    //first checking whether there is an existing random walk
                    //cout<< "Checking Wheter there is an existing random walk..." << endl;
                    //bool find_existing_walk = false;
                    if(!random_walk_storation[sourceNode].empty()){
                        pair<int, double> temp_pair = random_walk_storation[sourceNode].top();
                        tempNode = temp_pair.first;
                        double current_scalar_prob = temp_pair.second;
                        if(vert_count[tempNode] == 0){
                            value_verts[value_count++] = tempNode;
                        }
                        if(resultList[tempNode] == 0){
                            forwardCandidate[fora_candidate++] = tempNode;
                        }
                        resultList[tempNode] += increment * current_scalar_prob;
                        vert_count[tempNode]++;
                        random_walk_storation_one_iteration[sourceNode].push(temp_pair);
                        random_walk_storation[sourceNode].pop();

                    }
                    else{

                        double scalar_prob = 1;
                        while(R.drand() < alpha){

                            // cout << "Sum of weights of this source node ... " << scalar_prob << endl;
                            tempNode = g.getNextNode(tempNode, R); // the next node be chosen
                            //tempNode = vertex_alias[tempNode].generateRandom(R);
                            scalar_prob = g.getSumWeights(tempNode);
                            //cout<< "Select the next node: " << tempNode <<endl;
                        }

                        if(vert_count[tempNode] == 0){
                            value_verts[value_count++] = tempNode;
                        }
                        if(resultList[tempNode] == 0){
                            forwardCandidate[fora_candidate++] = tempNode;
                        }
                        resultList[tempNode] += increment * scalar_prob; // due to the random generator
                        //cout << "Updated MR value of node " << tempNode << ": " << resultList[tempNode] <<endl;
                        vert_count[tempNode]++;
                        //cout << "Ending Walking ..." << endl;
                        pair<int, double> temp_pair_prob = make_pair(tempNode, scalar_prob);
                        random_walk_storation_one_iteration[sourceNode].push(temp_pair_prob);
                    }

                }

            }

            //cout << "Ending Random Walk..." <<endl;
            clock_t t1 = clock();
            double gap1 = forward_c * sqrt(log(2*vert)) / (double) walk_num;
            double gap2 = forward_c * alpha * 7 / (double)3 * r_sum * log(vert) / (double) walk_num;
            double multiCoff = r_sum * r_sum * alpha * alpha / increment;
            clock_t t2 = clock();


            pq upper_pq(pqcompare(true));
            pq lower_pq(pqcompare(true));

            double temp_upper_min = 0, temp_lower_min = 0;
            for(int i = 0; i < fora_candidate; i++){
                int tempNode = forwardCandidate[i];
                betas[tempNode] = gap1 * sqrt(multiCoff * (resultList[tempNode] - pi_f[tempNode])) + gap2;
                double tempUpperBound = resultList[tempNode] + betas[tempNode];
                double tempLowerBound = resultList[tempNode] - betas[tempNode];
                if(i < k){
                    lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                    upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                }
                else if(i == k){
                    upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                    lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                    lower_pq.pop();
                    temp_lower_min = lower_pq.top().second;
                }
                else if(i == k + 1){
                    if(tempLowerBound > temp_lower_min){
                        lower_pq.pop();
                        lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                        temp_lower_min = lower_pq.top().second;
                    }
                    upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                    upper_pq.pop();
                    temp_upper_min = upper_pq.top().second;
                }
                else{
                    if(tempLowerBound > temp_lower_min){
                        lower_pq.pop();
                        lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                        temp_lower_min = lower_pq.top().second;
                    }
                    if(tempUpperBound > temp_upper_min){
                        upper_pq.pop();
                        upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                        temp_upper_min = upper_pq.top().second;
                    }
                }
            }
            upperBound = upper_pq.top().second;
            lowerBound = lower_pq.top().second;
            int leftCandidateSize = 0;

            for(int i = 0; i < fora_candidate; i++){
                int tempNode = forwardCandidate[i];
                double lowBound = resultList[tempNode] - betas[tempNode];
                double upBound = resultList[tempNode] + betas[tempNode];
                if(lowBound > upperBound){
                    leftCandidateSize++;
                }
                else if(lowBound <= upperBound && upBound >= lowerBound){
                    candidate_size++;
                }
            }
            clock_t t3 = clock();
            total_t1 += ta - t0;
            total_t2 += t1 - ta;
            total_t3 += t2 - t1;
            total_t4 += t3 - t2;

            if(leftCandidateSize >= k){
                clock_t t3_1 = clock();
                vector<pair<int, double> > leftVec, middleVec;
                for(int i = 0; i < fora_candidate; i++){
                    int tempNode = forwardCandidate[i];
                    double lowBound = resultList[tempNode] - betas[tempNode];
                    double upBound = resultList[tempNode] + betas[tempNode];

                    if(lowBound > upperBound){
                        leftVec.push_back(pair<int, double>(tempNode, resultList[tempNode]));
                    }
                    else if(upBound >= lowerBound){
                        middleVec.push_back(pair<int, double>(tempNode, resultList[tempNode]));
                    }
                }

                clock_t t4 = clock();


               vector<int> candidate_left;
                for(int i =0 ; i < leftVec.size(); i++){
                    candidate_left.push_back(leftVec[i].first);
                }


                cout << "left size: " << candidate_left.size() << endl;
                clock_t t5 = clock();
                double fora_precision = calPrecision(candidate_left, realList, candidate_left.size());
                cout << "fora precision: " << fora_precision << endl;
                cout << "forward time: " << total_t1 / (double) CLOCKS_PER_SEC << endl;
                cout << "random walk num: " << total_walk_num << endl;
                cout << "random walk time: " << total_t2 / (double) CLOCKS_PER_SEC << endl;
                cout << "constant time: " << total_t3 / (double) CLOCKS_PER_SEC << endl;
                cout << "sort time: " << total_t4 / (double) CLOCKS_PER_SEC << endl;
                cout << "bound time: " << (t4 - t3_1) / (double) CLOCKS_PER_SEC << endl;
                //cout << "useless time: " << (t5 - t4) / (double) CLOCKS_PER_SEC << endl;
                //cout << "back time: " << 0 << endl;
                double totalTime =
                        (total_t1 + total_t2 + total_t3 + total_t4 + t4 - t3_1) / (double) CLOCKS_PER_SEC;
                cout << "total time: " << totalTime << endl;
                avg_pre += fora_precision;
                avg_time += totalTime;
                avg_walk_num += total_walk_num;
                break;
            }

            for(int i =0; i < vert; i++){
                while(!random_walk_storation_one_iteration[i].empty()){
                    pair<int, double> new_temp_pair = random_walk_storation_one_iteration[i].top();
                    random_walk_storation[i].push(new_temp_pair);
                    random_walk_storation_one_iteration[i].pop();
                }
            }

            fwd_walk_iteration += 1;
            r_threshold *= 0.5;
        }

        delete[] pi_f;
        delete[] residue_f;
        delete[] isInQueue;
        delete[] betas;
        delete[] forwardCandidate;

        for(int i = 0; i < vert; i++){
            while(!random_walk_storation_one_iteration[i].empty()){
                pair<int, double> new_temp_pair = random_walk_storation_one_iteration[i].top();
                //random_walk_storation[i].push(new_temp_pair);
                random_walk_storation_one_iteration[i].pop();
            }
        }
        //vector<vector<stack<double> >>().swap(random_walk_storation_one_iteration);
       // return resultList;
    }

    //MCMR: only random walks
    void MonteCarlo(int u, unsigned int walk_num){
        vector<int> realList = getRealTopK(u, k);

        //for(int i=0; i< g.getOutSize(u); i++){
          //  cout<< g.outAdjList[u][i].first<<endl;
        //}

        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;

        clock_t t0 = clock();
        //double total_scalar_probs = 0.0;
        for(double i = 0; i < walk_num; i++){
            int tempNode = u;
            double scalar_probs = 1;
            while(R.drand() > (1 - alpha)* g.weight_scalar){
                //int length = g.getOutSize(tempNode);
                double weight = g.getSumWeights(tempNode);

                tempNode = g.getNextNode(tempNode, R);
                scalar_probs = scalar_probs * weight;
            }

            //cout<<"scalar of this walk: " << scalar_probs <<endl;
            resultList[tempNode] += scalar_probs;
        }

        clock_t t1 = clock();

        vector<pair<int, double>> mrs;

        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) walk_num;
            mrs.push_back(make_pair(i, resultList[i]));
            //cout<<"Source Node:" << u << ", Target Node: " << i << ", MR value: " << resultList[i] << endl;
        }

        sort(mrs.begin(), mrs.end(), maxScoreCmp);
        vector<int> topks;
        for(int i =0; i < k; i++){
            topks.push_back(mrs[i].first);
            //cout << "Top " << (i+1) << ": " << topks[i] << " " << mrs[i].second <<endl;
        }

        double precision_mr = calPrecision(topks, realList, k);
        cout << "precision: " << precision_mr << endl;
        avg_pre += precision_mr;
        avg_time +=  (t1 - t0) / (double) CLOCKS_PER_SEC ;

        delete [] resultList;

    }

    //MCMR_prune: random walks with the prune strategy (which is the same as MCMRPlus)
    void MonteCarlo_with_prune(int u, int k) {
        vector<int> realList = getRealTopK(u, k);

        int candidate =0;
        int* mc_candidate =  new int[vert];
        double *resultList = new double[vert];
        double * betas =  new double[vert];
        for (int i = 0; i < vert; i++){
            resultList[i] = 0;
            mc_candidate[i] = 0;
            betas[i] = 0;
        }

        clock_t t0 = clock();
        int walk_num = 100;
        int total_num  = walk_num;
        double upperBound = 0, lowerBound = 0;

        while(true) {
            cout<<"current walk num: " << walk_num << endl;
            //start the random walk with this new walk num
            for (double i = 0; i < walk_num; i++) {
                int tempNode = u;
                double scalar_probs = 1;
                while (R.drand() > 1 - alpha) {
                    //int length = g.getOutSize(tempNode);
                    tempNode = g.getNextNode(tempNode, R);
                    scalar_probs = scalar_probs * g.getSumWeights(tempNode);
                }
                resultList[tempNode] += scalar_probs;
                if(mc_candidate[tempNode] == 0){
                    mc_candidate[candidate++] = tempNode;
                }
            }

            vector<pair<int, double>> mrs;
            for (int i = 0; i < candidate; i++) {
                int tempNode = mc_candidate[i];
                resultList[tempNode] /= (double) total_num;
                mrs.push_back(make_pair(tempNode, resultList[tempNode]));
                //cout<<"Source Node:" << u << ", Target Node: " << i << ", MR value: " << resultList[i] << endl;
            }

            double gap1 = forward_c * sqrt(log(2 * vert)) / (double) total_num;
            double gap2 = forward_c * alpha * 7 / (double) 3 * log(vert) / (double) total_num;
            double multiCoff = alpha * alpha;

            typedef priority_queue<pair<int, double>, vector<pair<int, double> >, pqcompare> pq;
            pq upper_pq(pqcompare(true));
            pq lower_pq(pqcompare(false));

            double temp_upper_min = 0, temp_lower_min = 0;
            for (int i = 0; i < candidate; i++) {
                int tempNode = mc_candidate[i];
                betas[tempNode] = gap1 * sqrt(multiCoff * (resultList[tempNode])) + gap2;
                double tempUpperBound = resultList[tempNode] + betas[tempNode];
                double tempLowerBound = resultList[tempNode] - betas[tempNode];
                if (i < k) {
                    lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                    upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                } else if (i == k) {
                    upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                    lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                    lower_pq.pop();
                    temp_lower_min = lower_pq.top().second;
                } else if (i == k + 1) {
                    if (tempLowerBound > temp_lower_min) {
                        lower_pq.pop();
                        lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                        temp_lower_min = lower_pq.top().second;
                    }
                    upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                    upper_pq.pop();
                    temp_upper_min = upper_pq.top().second;
                } else {
                    if (tempLowerBound > temp_lower_min) {
                        lower_pq.pop();
                        lower_pq.push(pair<int, double>(tempNode, tempLowerBound));
                        temp_lower_min = lower_pq.top().second;
                    }
                    if (tempUpperBound > temp_upper_min) {
                        upper_pq.pop();
                        upper_pq.push(pair<int, double>(tempNode, tempUpperBound));
                        temp_upper_min = upper_pq.top().second;
                    }
                }
            }
            upperBound = upper_pq.top().second;
            lowerBound = lower_pq.top().second;
            int leftCandidateSize = 0;

            for (int i = 0; i < candidate; i++) {
                int tempNode = mc_candidate[i];
                double lowBound = resultList[tempNode] - betas[tempNode];
                double upBound = resultList[tempNode] + betas[tempNode];
                if (lowBound > upperBound) {
                    leftCandidateSize++;
                }
            }

            if (leftCandidateSize >= k) {
                clock_t t3_1 = clock();
                vector<pair<int, double> > leftVec;
                for (int i = 0; i < candidate; i++) {
                    int tempNode = mc_candidate[i];
                    double lowBound = resultList[tempNode] - betas[tempNode];
                    double upBound = resultList[tempNode] + betas[tempNode];

                    if (lowBound > upperBound) {
                        leftVec.push_back(pair<int, double>(tempNode, resultList[tempNode]));
                    }
                }

                clock_t t1 = clock();

                vector<int> candidate_left;
                for (int i = 0; i < leftVec.size(); i++) {
                    candidate_left.push_back(leftVec[i].first);
                }

                cout << "left size: " << candidate_left.size() << endl;
                clock_t t5 = clock();
                double cur_precision = calPrecision(candidate_left, realList, candidate_left.size());
                cout << "num walk: " << walk_num<< endl;
                cout << "precision: " << cur_precision << endl;

                //double precision_mr = calPrecision(topks, realList, k);
                //cout << "precision: " << precision_mr << endl;

                avg_pre += cur_precision;
                avg_time += (t1 - t0) / (double) CLOCKS_PER_SEC;

                break;
            }

            walk_num *=2;
            total_num += walk_num;
        }

        delete[] betas;
        delete[] mc_candidate;
        delete[] resultList;
    }

    void saveSumWeight(){
        g.saveSumWeight(target_filename);
    }

    void saveNormWeight(){
        g.saveNormarlizedWeight(target_filename);
    }

    /*double* MonteCarloFullPath(int u, unsigned int walk_num){
        vector<int> realList = getRealTopK(u, k);
        unordered_map<int, double> realMap = getRealTopKMap(u, k);
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;

        clock_t t0 = clock();
        double total_scalar_probs = 0.0;
        for(double i = 0; i < walk_num; i++){
            int tempNode = u;
            double scalar_probs = g.getSumWeights(u);
            vector<pair<int, double>> path;
            double current_weight = scalar_probs;
            while(R.drand() < alpha){
                //int length = g.getOutSize(tempNode);
                tempNode = g.getNextNode(tempNode, R);
                path.push_back(make_pair(tempNode, current_weight));
                scalar_probs = scalar_probs * g.getSumWeights(tempNode);
                current_weight = g.getSumWeights(tempNode);
            }
            for(auto it = path.begin(); it != path.end(); it++){
                int node_in_path = it->first;
                double weight_to_add = it->second;
                resultList[node_in_path] += (weight_to_add / scalar_probs);
            }
            //resultList[tempNode] += scalar_probs;
            total_scalar_probs += scalar_probs;
        }

        cout<<"total_scalar_probs" << total_scalar_probs<<endl;
        clock_t t1 = clock();
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;

        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) total_scalar_probs;
            //cout<<"Source Node:" << u << ", Target Node: " << i << ", MR value: " << resultList[i] << endl;
        }
        for(int i = 0; i<k; i++){
            int tempNode = realList[i];
            avg_top_err[i] += abs(realMap[tempNode] - resultList[tempNode]);
        }

        vector<int> newLeftVec;
        typedef priority_queue<pair<int, double>, vector<pair<int,double> >, pqcompare> pq;
        pq upper_pq(pqcompare(true));

        double UpperBound = 0;
        for(int i = 0; i < vert; i++){
            if(i < k){
                upper_pq.push(pair<int, double>(i, resultList[i]));
                if(i == k - 1)
                    UpperBound = upper_pq.top().second;
            }
            else{
                if(resultList[i] > UpperBound){
                    upper_pq.pop();
                    upper_pq.push(pair<int, double>(i, resultList[i]));
                    UpperBound = upper_pq.top().second;
                }
            }
        }
        for(int i = 0; i < k; i++){
            newLeftVec.push_back(upper_pq.top().first);
            upper_pq.pop();
        }
        for(int i =0; i < k; i ++){
            cout << "Top " << (i+1) << ": " << newLeftVec[i]<<endl;
        }

        cout << "precision: " << calPrecision(newLeftVec, realList, k) << endl;
        avg_pre += calPrecision(newLeftVec, realList, k);
        avg_time +=  (t1 - t0) / (double) CLOCKS_PER_SEC ;
        return resultList;
    }
*/

    void Power(int iterations, int u, int k){
        //unordered_map<int, double> map_residual;
        //map_residual.clear();
        //map_residual[u] = 1.0;
        clock_t t0 = clock();
        int num_iter=0;
        double* map_mr = new double[vert];
        for(int i = 0; i < vert; i++){
            map_mr[i] = 0.0;
        }
        map_mr[u] = 1.0;

        //double abs_error = 1.0;

        while(num_iter <= iterations){
            //abs_error =0.0;
            cout << u << ": iter " << num_iter << endl;
            num_iter++;
            for(int i =0; i< vert; i++){
                vector<pair<int, double> > weights_i = g.outAdjList[i];
                int size_i = weights_i.size();
                //cout<<"size" <<size_i<<endl;
                double new_mr = 0.0;
                for(int j=0; j<size_i; j++){
                    int neighbor = weights_i[j].first;
                    if(map_mr[neighbor] != 0.0){
                        double weight = weights_i[j].second;
                        new_mr += alpha * weight * map_mr[neighbor];
                    }
                }
                //old_ppr = map_ppr[i];
                map_mr[i] = new_mr;

                if(i == u){
                    map_mr[i] += (1-alpha);
                }

                //abs_error += abs(map_ppr[i]-old_ppr);
            }
            //cout<< "abs_error: "<< abs_error << endl;
        }



        vector<pair<int, double> > pprs;
        for(int j = 0; j < vert; j++){
            pprs.push_back(pair<int, double>(j, map_mr[j]));
        }
        sort(pprs.begin(), pprs.end(), maxScoreCmp);
        vector<int> topks;
        for(int i =0; i < k ; i++){
             topks.push_back(pprs[i].first);
         }
        clock_t t1 = clock();

         vector<int> realList = getRealTopK(u, k);
         //unordered_map<int, double> realMap = getRealTopKMap(u, k);

         cout << "precision: " << calPrecision(topks, realList, k) << endl;
         avg_pre += calPrecision(topks, realList, k);
         avg_time +=  (t1 - t0) / (double) CLOCKS_PER_SEC ;

        delete[] map_mr;
    }

    void PowerMethodK(int iterations, string outputFile, int u, int k){
        //unordered_map<int, double> map_residual;
        //map_residual.clear();
        //map_residual[u] = 1.0;
        clock_t t0 = clock();
        int num_iter=0;
        vector<double> map_ppr(vert);

        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0;
        }
        map_ppr[u] = 1;

        //double abs_error = 1.0;

        while(num_iter <= iterations){
            //abs_error =0.0;
            cout << u << ": iter " << num_iter << endl;
            num_iter++;
            vector<double> old_pprs = map_ppr;

            for(int i=0; i<g.n; i++){
                map_ppr[i] = 0;
            }
            map_ppr[u] = 1-alpha;

            for(int i =0; i<g.n; i++){
                vector<pair<int, double> > weights_i = g.outAdjList[i];
                int size_i = weights_i.size();
                //cout<<"size" <<size_i<<endl;
                if(size_i != 0){
                    for(int j=0; j<size_i; j++){
                        int neighbor = weights_i[j].first;
                        double weight = weights_i[j].second;
                        map_ppr[neighbor] += alpha * weight * old_pprs[i];
                    }
                }
                else{
                    map_ppr[u] += alpha * old_pprs[i];
                }
            }

        }

        clock_t t1 = clock();
        ofstream fout(outputFile);
        vector<pair<int, double> > pprs;
        for(int j = 0; j < vert; j++){
            pprs.push_back(pair<int, double>(j, map_ppr[j]));
        }
        sort(pprs.begin(), pprs.end(), maxScoreCmp);
        /* vector<int> topks;
         for(int i =0; i < k ; i++){
             topks.push_back(pprs[i].first);
         }
         */
        //double sum = 0.0;
        //int nodes_to_output = min(k, vert);
         for(int j = 0; j < vert; j++) {
             if (pprs[j].second >= 0) {
                 fout << pprs[j].first << " " << pprs[j].second << "\n";
                 //sum += pprs[j].second;
             }
         }

         //fout << "sum: " << sum;
         fout.close();


       /* vector<int> realList = getRealTopK(u, k);
        unordered_map<int, double> realMap = getRealTopKMap(u, k);

        cout << "precision: " << calPrecision(topks, realList, k) << endl;
        avg_pre += calPrecision(topks, realList, k);
 */
        avg_time +=  (t1 - t0) / (double) CLOCKS_PER_SEC ;

        //delete[] map_ppr;
    }

    void t_PowerMethod(vector<int> nodeList, int iterations){
        for(int i = 0; i < nodeList.size(); i++){
            int tempNode = nodeList[i];
            stringstream ss;
            ss << "real_mr/" << graph_name << "/" << tempNode << ".txt";
            string outputFile = ss.str();
            cout << "file: " << outputFile << endl;
            PowerMethodK(iterations, outputFile, tempNode, 20);
            cout << outputFile << "done!"  << endl;
        }
    }


    double calPrecision(vector<int> topK1, vector<int> realList, int k, bool isShowMissing = false){
        cout<<"computing precision ... "<<endl;
        unordered_set<int> set_to_be_searched;
        int hitCount = 0;
        int denominator = 0;
        if(realList.size() >=  topK1.size()){
            //for each node in topk list, check whether it is in the real list
            for(int i=0; i<realList.size(); i++){
                set_to_be_searched.insert(realList[i]);
            }

            for(int i=0; i< topK1.size(); i++){
                if(set_to_be_searched.count(topK1[i]) > 0)
                    hitCount++;
                //else
                    //cout<<"useless node: " << topK1[i] <<endl;
            }

            denominator = topK1.size();
        }
        else{
            //for each node in the real list, check whether it is in the top k list
            for(int i=0; i<topK1.size(); i++){
                set_to_be_searched.insert(topK1[i]);
            }

            for(int i=0; i< realList.size(); i++){
                if(set_to_be_searched.count(realList[i]) > 0){
                    hitCount++;
                }
                else
                    cout<<"useless node: " << topK1[i] <<endl;
            }

            denominator = realList.size();
        }

        double result = hitCount / (double) denominator;
        return result;
    }

    vector<int> getTopKList(double* resultList, int u, int k){
        vector<int> topKList;
        vector<pair<int,double> > sim;
        for(int i = 0; i < vert; i++){
            sim.push_back(pair<int, double>(i, resultList[i]));
        }
        sort(sim.begin(), sim.end(), maxScoreCmp);
        for(int i = 0; i < k; i++){
            topKList.push_back(sim[i].first);
        }
        return topKList;
    }

    //在一个大的pair<int, double>数组中取TopK大个（按照double的值大小）
    vector<pair<int, double> > getTopK(vector<pair<int, double> > target, int k){
        typedef priority_queue<pair<int, double>, vector<pair<int,double> >, pqcompare> pq;
        pq upper_pq(pqcompare(true));

        double UpperBound = 0;
        for(int i = 0; i < target.size(); i++){
            if(i < k){
                upper_pq.push(target[i]);
                if(i == k - 1)
                    UpperBound = upper_pq.top().second;
            }
            else{
                if(target[i].second > UpperBound){
                    upper_pq.pop();
                    upper_pq.push(target[i]);
                    UpperBound = upper_pq.top().second;
                }
            }
        }

        vector<pair<int, double> > answer;
        for(int i = 0; i < k; i++){
            answer.push_back(upper_pq.top());
            upper_pq.pop();
        }
        return answer;
    }

    double compute_y_mogul(int real_i, double q_vec_i, double L_ii, unordered_map<int, double> L_in_i, const vector<int> non_zero_ids_y, double* y_vec){
        //in this domain, all the node ids used are the new ids after permutation
        double y_i =0.0;
        if(non_zero_ids_y.size() == 0){
            if(q_vec_i == 0.0){
                y_i = 0.0;
            }
            else{
                //real_i is the query node
                y_i = q_vec_i / L_ii;
            }
        }
        else {
            double temp_y = 0.0;
            if (L_in_i.size() > non_zero_ids_y.size()) {
                for (int j = 0; j < non_zero_ids_y.size(); j++) {
                    int real_j = non_zero_ids_y[j]; //the new node id after permutation
                    if ((real_j < real_i) && L_in_i.count(real_j) > 0) {
                        //cout<<"L_" << real_i<<" " <<real_j<<": " << L_in_i[real_j] <<" * " << y_vec[real_j]<<endl;
                        temp_y += L_in_i[real_j] * y_vec[real_j];
                        //cout<<temp_y<<endl;
                    } else if (real_j >= real_i) {
                        break;
                    }
                }

            } else {
                for (auto L_iter = L_in_i.begin(); L_iter != L_in_i.end(); ++L_iter) {
                    pair<int, double> temp_pair = *L_iter;
                    int real_j = temp_pair.first;
                    if (real_j < real_i && real_j >= 0) {
                        //cout<<"L_" << real_i<<" " <<real_j<<": " << temp_pair.second <<" * " << y_vec[real_j]<<endl;
                        temp_y += temp_pair.second * y_vec[real_j];
                        //cout<<temp_y<<endl;
                    } else {
                        break;
                    }
                }
            }
            y_i = (q_vec_i - temp_y) / L_ii;
        }
        return y_i;
    }

    double compute_x_mogul(int real_i, double y_vec_i, double U_ii, unordered_map<int, double> U_in_i, const vector<int> non_zero_ids_x, VectorXd x_vec ){
        //in this domain, all the node ids used are the new ids after permutation
        double x_i;
        double temp_x = 0.0;
        if(U_in_i.size() > non_zero_ids_x.size()){
            for(int j=0; j<non_zero_ids_x.size(); j++){
                int real_j = non_zero_ids_x[j];
                if(real_j > real_i && U_in_i.count(real_j) > 0){
                    temp_x += U_in_i[real_j] * x_vec(real_j);
                } else if(real_j <= real_i){
                    break;
                }
            }

        }
        else{
            for(auto U_iter = U_in_i.begin(); U_iter != U_in_i.end(); ++U_iter){
                pair<int, double> temp_pair = * U_iter;
                int real_j = temp_pair.first;
                if(real_j > real_i && x_vec(real_j) != 0.0){
                    temp_x += temp_pair.second * x_vec(real_j);
                }else if(real_j > g.n -1){
                    break;
                }
            }
        }
        x_i = (y_vec_i - temp_x) / U_ii;

        return x_i;
    }

    double compute_U_bar_cluster_i(vector<int> nodes_in_C_i, SparseMtx U ){
        //all the nodes in this domain, are after permutation
        double U_bar_cluster_i = 0.0;

        vector<double> U_values;
        for(int j=0; j<nodes_in_C_i.size(); j++){
            for(int k=0; k < nodes_in_C_i.size(); k++){
                if(nodes_in_C_i[j] != nodes_in_C_i[k]){
                    if(U_bar_cluster_i < abs(U.coeff(nodes_in_C_i[j], nodes_in_C_i[k]))){
                        U_bar_cluster_i = abs(U.coeff(nodes_in_C_i[j], nodes_in_C_i[k]));
                    }
                }
            }
        }
        return U_bar_cluster_i;
    }

    double compute_U_bar_cluster_i_j(vector<int> nodes_in_C_i, int j, SparseMtx U){
        //in this domain, all the node ids are after permutation
        double U_bar_cluster_i_j = 0.0;
        vector<double> U_values;

        for(int i=0; i<nodes_in_C_i.size(); i++){
            if( U_bar_cluster_i_j < abs(U.coeff(nodes_in_C_i[i], j))){
                U_bar_cluster_i_j = abs(U.coeff(nodes_in_C_i[i], j));
            }
        }

        return U_bar_cluster_i_j;
    }

    void mogul_in_query_phase(vector<int> sourceNodes, string algo){
        //read  the pre-computed files
        SparseMtx per_mtrx(g.n, g.n);
        vector<int> new_node_ids(g.n);
        vector<int> pre_node_ids(g.n);
        for(int i=0; i<g.n; i++){
            new_node_ids[i] = -1;
            pre_node_ids[i] = -1;
        }
        cout<<"reading permutation file... "<<endl;
        get_permute_node_id(per_mtrx, new_node_ids, pre_node_ids);

        SparseMtx chol_L(g.n, g.n);
        SparseMtx chol_D(g.n, g.n);

        cout<<"reading decomposed matrices ..."<<endl;
        if(algo == "MOGUL"){
            string chol_L_filename = "./index/" + target_filename + "/incom_chol_L.txt";
            get_chol_L_matrix(chol_L, chol_L_filename);
            string chol_D_filename = "./index/" + target_filename + "/incom_chol_D.txt";
            get_chol_D_matrix(chol_D, chol_D_filename);
        }
        else if(algo == "MOGUL_E"){
            string chol_L_filename = "./index/" + target_filename + "/chol_L.txt";
            get_chol_L_matrix(chol_L, chol_L_filename);
            string chol_D_filename = "./index/" + target_filename + "/chol_D.txt";
            get_chol_D_matrix(chol_D, chol_D_filename);
        }

/*
        for(int i=0; i<g.n; i++){
            cout<<i<< ": ";
            for(int j=0; j<=i; j++){
                //if(chol_L.coeff(i,j) != 0.0){
                    cout << chol_L.coeff(i, j) <<" ";
                //}
            }
            cout<<endl;
        }

        for(int i = 0; i < g.n; i++){
            cout<<i<< ": ";
            for(int j=0; j<=i; j++)
                cout<< chol_D.coeff(i, j) <<" ";

            cout<<endl;
        }
*/
        cout<<"reading clusters ..."<<endl;
        string cluster_filename = "./index/" + target_filename + "/reordered_cluster.txt";
        string cluster_before_reorder_filename = "./index/" + target_filename + "/cluster.txt";
        vector<int> node_cluster(g.n);
        for(int i=0; i< g.n; i++){
            node_cluster[i] = -1;
        }
        bool totally_discon = false;
        vector<vector<int>> cluster_nodes = get_clusters(node_cluster, cluster_filename, cluster_before_reorder_filename, totally_discon);


        cout<<"starting query phase ... "<<endl;
        for(int i=0; i<sourceNodes.size(); i++){
            int source = sourceNodes[i];
            cout <<i<< "-th node: " << source << endl;

            vector<int> realList = getRealTopK(source, k);
            vector<int> topk_list = mogul(source, k, totally_discon, node_cluster, cluster_nodes, new_node_ids, pre_node_ids, chol_L, chol_D, per_mtrx);

            cout<<"yes mogul successful "<<endl;
            double cur_precision = calPrecision(topk_list, realList, realList.size());
            cout<<"precision: " << cur_precision<<endl;
            avg_pre += cur_precision;

            cout<< "------------"<<endl;
        }

    }



    vector<int> mogul(int source, int k, bool totally_discon, const vector<int> node_cluster, const vector<vector<int>> cluster_nodes, const vector<int> new_node_ids, const vector<int> pre_node_ids, const SparseMtx chol_L, const SparseMtx chol_D, const SparseMtx per_mtx){

        BinaryHeap<double, less<double> > topk_mrs(g.n, less<double>());

        VectorXd q_vec (g.n);
        VectorXd x_vec (g.n);
        double* y_vec = new double[g.n];
        for(int i=0; i < g.n; i++) {
            if(i == source)
                q_vec(i) = 1-alpha; //(1-alpha)Pq
            else
                q_vec(i) = 0.0;

            x_vec(i) = 0.0;
            y_vec[i] = 0.0;
        }

        clock_t start_time = clock();

        q_vec = per_mtx * q_vec; //new query vector

        int real_source = new_node_ids[source];
        cout<<"real source id: " <<real_source <<endl;
        double delta = 0; //the k-th largest value


        SparseMtx new_L(g.n, g.n);
        new_L = (chol_L*chol_D).pruned();
        SparseMtx U(g.n, g.n);
        U = chol_L.transpose();

        /*
        for(int i=0; i<g.n; i++){
            cout<<i<< ": ";
            for(int j=0; j<=i; j++){
                //if(chol_L.coeff(i,j) != 0.0){
                cout << new_L.coeff(i, j) <<" ";
                //}
            }
            cout<<endl;
        }

        for(int i = 0; i < g.n; i++){
            cout<<i<< ": ";

            for(int j=0; j<g.n; j++){
                if(j < i){
                    cout << 0 <<" ";
                }
                else{
                    cout<< U.coeff(i, j) <<" ";
                }
            }

            cout<<endl;
        }*/

        vector<vector<int> > non_zero_ids_new_L =  vector<vector<int>>(g.n, vector<int>());
        vector<vector<int> > non_zero_ids_U = vector<vector<int>>(g.n, vector<int>());
        for(int i=0; i<g.n; i++){
            vector<int> temp_vec;
            //temp_vec.reserve(g.n);
            non_zero_ids_new_L[i] = temp_vec;
            vector<int> temp_U_vec;
            //temp_U_vec.reserve(g.n);
            non_zero_ids_U[i] = temp_U_vec;
        }

        for(int i=0; i<new_L.outerSize(); i++){
            for (SparseMatrix<double>::InnerIterator it(new_L,i); it; ++it)
            {
                if(it.row() != it.col())
                    non_zero_ids_new_L[it.row()].push_back(it.col());
            }
        }
        for(int i=0; i<U.outerSize(); i++){
            for (SparseMatrix<double>::InnerIterator it(U,i); it; ++it)
            {
                if(it.row() != it.col())
                    non_zero_ids_U[it.row()].push_back(it.col());
            }
        }
/*
        for(int i = 0; i < g.n; i++){
            cout<<i<< ": ";

            for(int j=0; j<non_zero_ids_new_L[i].size(); j++){

                cout<< non_zero_ids_new_L[i][j]<<" ";
            }

            cout<<endl;
        }

        for(int i = 0; i < g.n; i++){
            cout<<i<< ": ";

            for(int j=0; j<non_zero_ids_U[i].size(); j++){

                cout<< non_zero_ids_U[i][j]<<" ";
            }

            cout<<endl;
        }
*/

        vector<int> non_zero_ids_y; //the new node ids after permutation
        vector<int> non_zero_ids_x; //the new node ids after permutation



        if(totally_discon || node_cluster[source] == cluster_nodes.size()-1) {//no C_n exists where C_n is the cluster that contains the nodes with cross cluster edges
            //or query node is in C_n
            //in this case, we don't need to compute the upper bound for each cluster
            cout<<"C_q and C_N are the same cluster"<<endl;
            //firstly, compute y values for each node in C_Q
            int C_q = node_cluster[source];
            //cout<<C_q<<endl;
            vector<int> new_node_ids_in_C_q;
            for(int i=0; i<cluster_nodes[C_q].size(); i++){
                //cout<<cluster_nodes[C_q][i] <<" ";
                new_node_ids_in_C_q.push_back(new_node_ids[cluster_nodes[C_q][i]]);
            }

            //cout<<endl;
            sort(new_node_ids_in_C_q.begin(), new_node_ids_in_C_q.end());
            cout<<"computing y for nodes in C_q" <<endl;
            for (int i = 0; i < new_node_ids_in_C_q.size(); i++) {
                //int temp_node = cluster_nodes[C_q][i]; //node id in the graph

                int temp_new_node_id = new_node_ids_in_C_q[i]; //new node id after the permutation
                //cout<<temp_new_node_id<<endl;

                unordered_map<int, double> L_in_i;//put the non zero id into this cluster
                for(int j = 0; j<non_zero_ids_new_L[temp_new_node_id].size(); j++){
                    int temp_j_id = non_zero_ids_new_L[temp_new_node_id][j]; //the new id after permuting
                    L_in_i.insert(make_pair(temp_j_id, new_L.coeff(temp_new_node_id, temp_j_id)));
                }

                y_vec[temp_new_node_id] = compute_y_mogul(temp_new_node_id, q_vec(temp_new_node_id),new_L.coeff(temp_new_node_id, temp_new_node_id),L_in_i,non_zero_ids_y, y_vec);

                if(y_vec[temp_new_node_id] !=0.0){
                    non_zero_ids_y.push_back(temp_new_node_id);
                }
                L_in_i.clear();

                //cout<<temp_new_node_id << " " << y_vec[temp_new_node_id]<<endl;
            }


            //secondly, compute the x values for the last node and compute the top k list
            int counter = 0;// count the number of nodes in the heap
            cout<<"computing x for each node ... " <<endl;
            for(int i= 0; i < g.n; i++){
                int real_i = g.n - 1 -i; //backward node id
                //cout<<real_i<<endl;


                if(real_i == g.n -1){
                    x_vec(real_i) = y_vec[real_i] / U.coeff(real_i, real_i);
                }
                else{
                    unordered_map<int, double> U_in_i;//put the non zero id of this row
                    for(int j = 0; j<non_zero_ids_U[real_i].size(); j++){
                        int temp_j_id = non_zero_ids_U[real_i][j]; //the new id after permuting
                        if(temp_j_id < 0 || temp_j_id >= g.n){
                            cout<<"incorrect: "<< real_i <<" " <<temp_j_id<<endl;
                            cout<<non_zero_ids_U[real_i].size()<<endl;
                        }else{
                            U_in_i.insert(make_pair(temp_j_id, U.coeff(real_i, temp_j_id)));
                        }
                    }
                    x_vec(real_i) = compute_x_mogul(real_i, y_vec[real_i], U.coeff(real_i, real_i), U_in_i, non_zero_ids_x, x_vec);
                }

                if(x_vec(real_i) != 0.0)
                    non_zero_ids_x.push_back(real_i);

                if(counter < k){//adding this node to the heap
                    topk_mrs.insert(real_i, x_vec(real_i));
                    //topk_mrs.display();
                    counter++; // the # nodes in this heap increases by 1
                }
                else{// we need to remove a node from the topk heap and insert a node from this
                    //cout<<"more than k nodes in the topk list"<<endl;
                    //cout<<counter<<endl;
                    auto top = topk_mrs.extract_top();
                    delta = top.first; //the value
                    //cout<< top.second <<" " << top.first<<endl;
                    if(x_vec(real_i) > delta){//add it into min heap
                        topk_mrs.delete_top();
                        topk_mrs.insert(real_i, x_vec(real_i));
                        //topk_mrs.display();
                    }
                }


                //cout<<real_i <<" " << x_vec(real_i)<<endl;
            }

        }
        else{//there exists C_n
            cout<<"C_q and C_N are different ..." <<endl;
            vector<int> nodes_in_C_q_C_N; //the node id after permutations
            int C_N = cluster_nodes.size() - 1;
            int C_q = node_cluster[source];
            int c_N = g.n-1; //the minimum node id in the last cluster
            for(int i =0; i<cluster_nodes[C_q].size(); i++){
                int temp_node = cluster_nodes[C_q][i]; //node id in the graph
                int temp_new_node_id = new_node_ids[temp_node]; //new node id after the permutation
                nodes_in_C_q_C_N.push_back(temp_new_node_id);
            }
            for(int i =0; i<cluster_nodes[C_N].size(); i++){
                int temp_node = cluster_nodes[C_N][i]; //node id in the graph
                int temp_new_node_id = new_node_ids[temp_node]; //new node id after the permutation
                nodes_in_C_q_C_N.push_back(temp_new_node_id);
                if(temp_new_node_id < c_N)
                    c_N = temp_new_node_id;
            }
            cout<<"# of nodes in C_q and C_N: " << nodes_in_C_q_C_N.size()<<endl;
            sort(nodes_in_C_q_C_N.begin(), nodes_in_C_q_C_N.end()); //in the increasing order of new node ids (after permutation)

            for (int i = 0; i < nodes_in_C_q_C_N.size(); i++) {
                int temp_new_node_id = nodes_in_C_q_C_N[i]; //node id after permutation

                unordered_map<int, double> L_in_i;//put the non zero id into this cluster
                for(int j = 0; j<non_zero_ids_new_L[temp_new_node_id].size(); j++){
                    int temp_j_id = non_zero_ids_new_L[temp_new_node_id][j]; //the new id after permuting
                    L_in_i.insert(make_pair(temp_j_id, new_L.coeff(temp_new_node_id, temp_j_id)));
                }
                y_vec[temp_new_node_id] = compute_y_mogul(temp_new_node_id, q_vec(temp_new_node_id), new_L.coeff(temp_new_node_id, temp_new_node_id), L_in_i, non_zero_ids_y, y_vec);

                if(y_vec[temp_new_node_id] !=0.0){
                    non_zero_ids_y.push_back(temp_new_node_id);
                }
            }

            //compute x values for nodes in C_Q nad C_N
            int counter = 0;// count the number of nodes in the heap
            for(int i= nodes_in_C_q_C_N.size()-1; i>=0; i--){
                int real_i = nodes_in_C_q_C_N[i]; //backward node id
                unordered_map<int, double> U_in_i;//put the non zero id of this row
                for(int j = 0; j<non_zero_ids_U[real_i].size(); j++){
                    int temp_j_id = non_zero_ids_U[real_i][j]; //the new id after permuting
                    U_in_i.insert(make_pair(temp_j_id, U.coeff(real_i, temp_j_id)));
                }

                if(real_i == g.n -1){
                    x_vec(real_i) = y_vec[real_i] / U.coeff(real_i, real_i);
                }
                else{
                    x_vec(real_i) = compute_x_mogul(real_i, y_vec[real_i], U.coeff(real_i, real_i), U_in_i, non_zero_ids_x, x_vec);
                }

                if(x_vec(real_i) != 0.0)
                    non_zero_ids_x.push_back(real_i);

                if(counter < k){//adding this node to the heap
                    topk_mrs.insert(real_i, x_vec(real_i));
                    counter++; // the # nodes in this heap increases by 1
                }
                else{// we need to remove a node from the topk heap and insert a node from this
                    auto top = topk_mrs.extract_top();
                    delta = top.first; //the value
                    if(x_vec(real_i) > delta){//add it into min heap
                        topk_mrs.delete_top();
                        topk_mrs.insert(real_i, x_vec(real_i));
                    }
                }
            }

            //for each other cluster, compute their bound estimation and then go into the computation of mrs.
            unordered_map<int, double> cluster_upper_bound; // containing the upper bound for each cluster
            for(int i=0; i<cluster_nodes.size(); i++){
                if(i != C_q && i != C_N){
                    //put the nodes in C_i into a new cluster
                    vector<int> nodes_in_C_i;
                    for(int k =0; i<cluster_nodes[k].size(); k++){
                        int temp_node = cluster_nodes[i][k]; //node id in the graph
                        int temp_new_node_id = new_node_ids[temp_node]; //new node id after the permutation
                        nodes_in_C_i.push_back(temp_new_node_id);
                    }

                    //compute the U_bar_C_i
                    double U_bar_C_i = compute_U_bar_cluster_i(nodes_in_C_i, U);

                    //compute X_bar_C_i
                    double X_bar_C_i =0.0;
                    for(int j = c_N; j< g.n; j++){
                        if(x_vec(j)!=0.0){
                            double U_bar_C_i_j = compute_U_bar_cluster_i_j(nodes_in_C_i, j, U);
                            X_bar_C_i += U_bar_C_i_j * abs(x_vec(j));
                        }
                    }

                    double upper_bound_C_i = X_bar_C_i * pow(1+U_bar_C_i, nodes_in_C_i.size() -1);

                    auto top = topk_mrs.extract_top();
                    delta = top.first;
                    if(upper_bound_C_i >= delta){
                        //compute the x values for the nodes in this
                        sort(nodes_in_C_i.begin(), nodes_in_C_i.end());
                        for(int i= nodes_in_C_i.size()-1; i>=0; i--){
                            int real_i = nodes_in_C_i[i]; //backward node id
                            unordered_map<int, double> U_in_i;//put the non zero id of this row
                            for(int j = 0; j<non_zero_ids_U[real_i].size(); j++){
                                int temp_j_id = non_zero_ids_U[real_i][j]; //the new id after permuting
                                U_in_i.insert(make_pair(temp_j_id, U.coeff(real_i, temp_j_id)));
                            }

                            if(real_i == g.n -1){
                                x_vec(real_i) = y_vec[real_i] / U.coeff(real_i, real_i);
                            }
                            else{
                                x_vec(real_i) = compute_x_mogul(real_i, y_vec[real_i], U.coeff(real_i, real_i), U_in_i, non_zero_ids_x, x_vec);
                            }

                            if(x_vec(real_i) != 0.0)
                                non_zero_ids_x.push_back(real_i);

                            if(counter < k){//adding this node to the heap
                                topk_mrs.insert(real_i, x_vec(real_i));
                                counter++; // the # nodes in this heap increases by 1
                            }
                            else{// we need to remove a node from the topk heap and insert a node from this
                                auto top = topk_mrs.extract_top();
                                delta = top.first; //the value
                                if(x_vec(real_i) > delta){//add it into min heap
                                    topk_mrs.delete_top();
                                    topk_mrs.insert(real_i, x_vec(real_i));
                                }
                            }
                        }

                    }
                    //otherwise do nothing
                }
            }

        }

        //for each node in the topk heap find out its true node id
        cout<<"permuting the nodes into the previous node id "<<endl;
        vector<int> topk_list; //containing the real node id
        vector<pair<int, double>> sorted_topk_list;
        //cout<<"node  mr"<<endl;


        int return_size = min(k, topk_mrs.size());
        cout<<"return size: " << return_size<<endl;
        for(int i = 0; i<return_size; i++){
            //cout<<"i: " << i<<endl;
            auto top = topk_mrs.extract_top();
            //convert into the previous id
            //cout<<top.second << " " << top.first<<endl;
            int pre_node_id = pre_node_ids[top.second];
            //cout<<pre_node_id <<" " <<top.first<<endl;
            sorted_topk_list.push_back(make_pair(pre_node_id, top.first));
            topk_mrs.delete_top();
            //topk_mrs.display();
            //cout<<topk_mrs.size()<<endl;
        }
        sort(sorted_topk_list.begin(), sorted_topk_list.end(), maxScoreCmp);
        for(int i=0; i<sorted_topk_list.size(); i++){
            topk_list.push_back(sorted_topk_list[i].first);
            cout<<sorted_topk_list[i].first<<" " <<sorted_topk_list[i].second<<endl;
        }

        clock_t end_time = clock();
        double query_time = (end_time - start_time) /(double) CLOCKS_PER_SEC;
        cout<<"query time: " << query_time<<endl;
        avg_time += query_time;

        delete [] y_vec;

        return topk_list;
    }

    void get_unweighted_graph(){
        string full_filename = "./dataset/unweighted/" + target_filename +".txt";
        ofstream graph_file(full_filename);
        graph_file<<vert;

        for(int i=0; i<vert; i++){
            int temp_deg = g.getOutSize(i);
            for(int j=0; j< temp_deg; j++){
                int temp_neigh  = g.outNeiList[i][j];
                graph_file<<endl<<i<<" " << temp_neigh;
            }
        }
        graph_file.close();
    }

    void preprocessing_mogul(){
        string full_filename= "./dataset/unweighted/" + target_filename +".txt";
        string cluster_file = "./index/" + target_filename +"/cluster.txt";
        string time_file = "./index/" + target_filename + "/preprocessing.time";

        string reor_cluster_file = "./index/" + target_filename + "/reordered_cluster.txt";

        string chol_L_filename = "./index/" + target_filename + "/incom_chol_L.txt";
        string chol_D_filename = "./index/" + target_filename + "/incom_chol_D.txt";

        Graph ug;
        ug.inputGraph(full_filename);

        ofstream t_file(time_file);
        double cluster_time = graph_clustering(ug, cluster_file);
        t_file <<"clustering time: " << cluster_time<<endl;

        double reorder_time = get_final_clusters(ug, cluster_file);

        t_file<<"reorder time: " << reorder_time<<endl;

        double permute_time = reorder_nodes(ug, reor_cluster_file);

        t_file<<"permute time: " << permute_time<<endl;

        double decomp_time = pre_in_chol_decomp(chol_L_filename, chol_D_filename);

        t_file <<"decomposition time: " << decomp_time<<endl;
        t_file<<"total preprocessing time: " << reorder_time + cluster_time + permute_time + decomp_time<<endl;
        t_file.close();

        cout<<"total preprocessing time: " << reorder_time + cluster_time + permute_time + decomp_time<<endl;

    }

    void preprocessing_mogul_exact(){
        string full_filename= "./dataset/unweighted/" + target_filename +".txt";
        string cluster_file = "./index/" + target_filename +"/cluster.txt";
        string time_file = "./index/" + target_filename + "/preprocessing_mogul_E.time";

        string reor_cluster_file = "./index/" + target_filename + "/reordered_cluster.txt";
        string chol_L_filename = "./index/" + target_filename + "/chol_L.txt";
        string chol_D_filename = "./index/" + target_filename + "/chol_D.txt";

        Graph ug;
        ug.inputGraph(full_filename);

        ofstream t_file(time_file);
        double cluster_time = graph_clustering(ug, cluster_file);
        t_file <<"clustering time: " << cluster_time<<endl;

        double reorder_time = get_final_clusters(ug, cluster_file);

        t_file<<"reorder time: " << reorder_time<<endl;

        double permute_time = reorder_nodes(ug, reor_cluster_file);

        t_file<<"permute time: " << permute_time<<endl;

        double decomp_time = pre_chol_decomp(chol_L_filename, chol_D_filename);

        t_file <<"decomposition time: " << decomp_time<<endl;
        t_file<<"total preprocessing time: " << reorder_time + cluster_time + permute_time + decomp_time<<endl;
        t_file.close();

        cout<<"total preprocessing time: " << reorder_time + cluster_time + permute_time + decomp_time<<endl;

    }

    //put the nodes with inter-edges into the N-th cluster
    double get_final_clusters(Graph& ug, string cluster_filename) {

        int num_clusters;
        ifstream cluster_file(cluster_filename);
        cluster_file >> num_clusters;

        cout << num_clusters << endl;
        vector<vector<int>> cluster_nodes; //cluster_id, node_id
        final_cluster.clear(); // node_id, cluster_id
        final_cluster.reserve(ug.n);

        cluster_nodes = vector<vector<int> >(num_clusters+1, vector<int>());

        bool* has_cross_edge = new bool[ug.n];


        for (int i = 0; i < ug.n; i++) {
            int temp_node;
            int temp_cluster;
            cluster_file >> temp_cluster >> temp_node;
            //cout<< temp_cluster <<endl;
            cluster_nodes[temp_cluster].push_back(temp_node);
            final_cluster[temp_node] = temp_cluster;

            has_cross_edge[i] = false;
        }
        cluster_file.close();
        cout << "reading clusters done ... " << endl;

        clock_t t0 = clock();
        //put the nodes who has the cross-cluster edges into the last cluster
        bool changed = false;
        for (int i = 0; i < num_clusters; i++) {
            int num_nodes = cluster_nodes[i].size();
            //cout<<"num_nodes: " << num_nodes<<endl;
            int num_deleted = 0; //total # of nodes that have been deleted.
            for (int j = 0; j < num_nodes; j++) {
                int temp_node = cluster_nodes[i][j-num_deleted];
                //
                if(has_cross_edge[temp_node] != true){//currently, temp_node has not been checked.
                    //cout<<temp_node<<endl;
                    for (int k = 0; k < ug.getOutSize(temp_node); k++) {
                       // int neighbor = ug.outNeiList[temp_node][k];
                        int neighbor = ug.outAdjList[temp_node][k];
                        //cout<<k<< " neighbor: " << neighbor <<endl;
                        if (final_cluster[neighbor] != i) { // temp_node has an cross edge
                            //put temp_node and neighbor to the last cluster
                            //cout<<"putting temp node to the last cluster .."<<endl;
                            cluster_nodes[num_clusters].push_back(temp_node);
                            cluster_nodes[i].erase(cluster_nodes[i].begin() + j - num_deleted); //remove temp_node from c_i
                            num_deleted++;
                            //cout<<"num that have been deleted: " << num_deleted <<endl;
                            final_cluster[temp_node] = num_clusters;
                            has_cross_edge[temp_node] = true;
                            changed = true;
                            break;
                        }
                    }

                }


            }
        }

        clock_t t1 = clock();
        double total_time = (t1 - t0) / (double) CLOCKS_PER_SEC;

        cout << "clustering settle down ... " << endl;
        stringstream ss2;
        ss2 << "index/" << target_filename << "/reordered_cluster.txt";
        ofstream final_cluster_file(ss2.str());
        //if(cluster_nodes[num_clusters].size() != g.n){
        if(changed){
            final_cluster_file << num_clusters + 1 << endl;
            for (int i = 0; i < num_clusters + 1; i++) {
                int num_nodes = cluster_nodes[i].size();
                for (int j = 0; j < num_nodes; j++) {
                    final_cluster_file << i <<" " <<cluster_nodes[i][j] << endl;
                }
            }
            final_cluster_file.close();
        }
        else{
            final_cluster_file << num_clusters << endl;
            for (int i = 0; i < num_clusters; i++) {
                int num_nodes = cluster_nodes[i].size();
                for (int j = 0; j < num_nodes; j++) {
                    final_cluster_file << i <<" " <<cluster_nodes[i][j] << endl;
                }
            }
            final_cluster_file.close();
        }

        delete [] has_cross_edge;

        return total_time;
    }

    //get the permutation matrix and return time cost
    double reorder_nodes(Graph& ug, string cluster_filename){
        //read the nodes in each cluster
        int num_clusters;
        ifstream cluster_file(cluster_filename);
        cluster_file >> num_clusters;

        //cout<< num_clusters << endl;
        vector<vector<int>> cluster_nodes;
        cluster_nodes = vector<vector<int> >(num_clusters, vector<int>());

        for(int i=0; i < ug.n; i++){
            int temp_node;
            int temp_cluster;
            cluster_file >> temp_cluster >> temp_node;
            cluster_nodes[temp_cluster].push_back(temp_node);
        }
        cluster_file.close();

        cout<<"computing permutation ..." <<endl;

        string permute_filename = "./index/" + target_filename + "/permutation.txt";
        ofstream permute_file(permute_filename);
        clock_t t0 = clock();
        int k=0;
        static BinaryHeap<int, less<long long> > min_heap_edges(ug.n, less<long long>());
        for(int i=0; i<num_clusters; i++){
            int num_nodes = cluster_nodes[i].size();
           // cout<<i<<endl;
            //cout<<num_nodes<<endl;

            min_heap_edges.clear();
            //cout<<min_heap_edges.size()<<endl;
            //vector<bool> processed(num_nodes);
            for(int j=0; j<num_nodes; j++){
                //processed[j] = false;
                min_heap_edges.insert(cluster_nodes[i][j], ug.getOutSize(cluster_nodes[i][j]));
            }


            while(min_heap_edges.size() > 0){
                auto top = min_heap_edges.extract_top();
                int node_to_permute = top.second;
                min_heap_edges.delete_top();

                //permute
                permute_file << k <<" " << node_to_permute << " " << 1 <<endl;
                k++;
            }

        }
        permute_file.close();
        cout << "permuting done .. " <<endl;

        clock_t t1 = clock();
        double time  = (t1 - t0) / (double) CLOCKS_PER_SEC;

        return time;
    }

    //incomplete cholesky decomp
    void incom_cholesky_decomp(vector<vector<int>> non_zero_idx, SparseMtx& W_mtrx, SparseMtx& chol_L, SparseMtx& chol_D){
        //store the idx whose W entry is non zero.
        //for each row of this matrix

        vector<vector<pair<int, double> > > temp_L_values = vector<vector<pair<int, double>>> (g.n, vector<pair<int, double>>());//row id i: j, L_ij which is non zero
        double * temp_D_values =  new double[g.n]; // i, i

        vector<T> L_tripletList; //chol_L
        vector<T> D_tripletList; //chol_D
        L_tripletList.reserve(g.m * 2);
        D_tripletList.reserve(g.n);

        //start decomposing from each row
        for(int i=0; i<g.n; i++){
            if(i == 0){
                double W_ii = W_mtrx.coeff(i,i);
                cout<< W_ii <<endl;
                temp_L_values[i].push_back(make_pair(i, 1.0));
                L_tripletList.push_back(T(i,i, 1));

                temp_D_values[i] = W_ii;
                D_tripletList.push_back(T(i,i, W_ii));
                //chol_D(i,i) = W_ii; //D_11 = W_11 ,here 1 is 0
                //chol_L(i,i) = 1;
            }
            else{
                for(int j=0; j<non_zero_idx[i].size(); j++){
                    int real_j = non_zero_idx[i][j];
                    if(real_j > i){
                        break;
                    }
                    else if(real_j == i){//L_ii and D_ii
                        //chol_L(i, real_j) = 1;
                        temp_L_values[i].push_back(make_pair(real_j, 1));
                        L_tripletList.push_back(T(i, real_j, 1));

                        double W_ij = W_mtrx.coeff(i, real_j);
                        double D_ij = 0;


                        for(int k = 0; k<temp_L_values[real_j].size(); k++){
                            pair<int, double> temp_pair = temp_L_values[real_j][k];
                            int real_k = temp_pair.first;
                            if(real_k < real_j){
                                D_ij += temp_pair.second * temp_pair.second * temp_D_values[real_k];
                            }
                            else{
                                break;
                            }
                        }

                        temp_D_values[real_j] = W_ij - D_ij;
                        D_tripletList.push_back(T(i, real_j, W_ij-D_ij));
                        //chol_D(i, real_j) = W_mtrx(i, real_j) - D_ij;

                    }
                    else{// real_j < i: compute L_ij only
                        double W_ij = W_mtrx.coeff(i, real_j);
                        double L_ij = 0;

                        unordered_map<int, double> L_in_i; //the values of L_ik where k<j-1
                        for(int k=0; k<temp_L_values[i].size(); k++){
                            pair<int, double> temp_pair = temp_L_values[i][k];
                            int real_k = temp_pair.first;
                            if(real_k < real_j){
                                L_in_i.insert(make_pair(real_k, temp_pair.second));
                            }
                        }

                        for(int k=0; k<temp_L_values[real_j].size(); k++){
                            pair<int, double> temp_pair = temp_L_values[real_j][k];
                            int real_k = temp_pair.first;
                            if((real_k < real_j) && (L_in_i.count(real_k)>0)){
                                L_ij += L_in_i[real_k] * temp_pair.second * temp_D_values[real_k];
                            }
                        }

                        L_ij = (W_ij - L_ij) / temp_D_values[real_j];

                        temp_L_values[i].push_back(make_pair(real_j, L_ij));
                        L_tripletList.push_back(T(i, real_j, L_ij));

                    }
                }
            }

        }

        chol_L.setFromTriplets(L_tripletList.begin(), L_tripletList.end());
        chol_D.setFromTriplets(D_tripletList.begin(), D_tripletList.end());

        delete [] temp_D_values;
    }

    double pre_in_chol_decomp(string chol_L_filename, string chol_D_filename){

        SparseMtx adj_mtrx(g.n, g.n);
        SparseMtx diag_mtrx(g.n, g.n);
        convert_graph_to_matrix_eigen(adj_mtrx, diag_mtrx);
        SparseMtx per_mtrx(g.n, g.n);
        get_permutation_matrix(per_mtrx);
        SparseMtx per_trans = per_mtrx.transpose();

        clock_t permute_start = clock();
        //doing permutation on the matrices
        adj_mtrx = (per_mtrx * adj_mtrx).pruned();
        adj_mtrx = (adj_mtrx * per_trans).pruned();

        diag_mtrx = (per_mtrx * diag_mtrx).pruned();
        diag_mtrx = (diag_mtrx * per_trans).pruned();

        SparseMtx W_mtrx(g.n, g.n);
        W_mtrx = (diag_mtrx * adj_mtrx).pruned();
        W_mtrx = (W_mtrx * diag_mtrx).pruned();

        SparseMtx I (g.n, g.n);
        I.setIdentity();
        W_mtrx = (I - alpha * W_mtrx).pruned();

        /*
        cout<<"W: " <<endl;
        for(int i=0; i<g.n; i++){
            cout<<i<<": ";
            for(int j=0; j<=i; j++){
                cout<<W_mtrx.coeff(i, j)<<" ";
            }
            cout<<endl;
        }

        vector<vector<int>> non_zero_ids_W;
        for(int i=0; i<g.n; i++){
            vector<int> temp_vec;
            non_zero_ids_W.push_back(temp_vec);
        }

        for(int i=0; i< W_mtrx.outerSize(); i++){
            for (SparseMatrix<double>::InnerIterator it(W_mtrx,i); it; ++it)
            {
                non_zero_ids_W[it.row()].push_back(it.col());
            }
        }

        for(int i=0; i<g.n; i++){
            cout<<i <<": ";
            for(int j=0; j<non_zero_ids_W[i].size(); j++){
                cout<<non_zero_ids_W[i][j]<<" ";
            }
            cout<<endl;
        }
         */

        clock_t  permute_end = clock();
        double permute_time = (permute_end - permute_start) /(double) CLOCKS_PER_SEC;
        cout<<"permute time: " << permute_time <<endl;

        cout<<"decomposing ... " <<endl;

        clock_t decom_start = clock();
        SparseMtx chol_L(g.n, g.n);
        SparseMtx chol_D(g.n, g.n);

        //get the non-zero idxes
        vector<vector<int>> non_zero_idx = vector<vector<int>> (g.n, vector<int>());//row id: a list of col id whose entry is non-zero

        for(int i=0; i<g.n; i++){
            vector<int> temp_row;
            //cout<<i <<" ";
            for(int j=0; j<g.n; j++){
                if(W_mtrx.coeff(i, j) != 0){
                    temp_row.push_back(j);
                    //cout<< j <<" ";
                }
            }
            //cout<<endl;
            non_zero_idx[i] = temp_row;
        }

        incom_cholesky_decomp(non_zero_idx, W_mtrx,chol_L, chol_D);

        clock_t decom_end = clock();
        double decom_time = (decom_end - decom_start) / (double) CLOCKS_PER_SEC;
        cout<<"decomposing time: " << decom_time <<endl;

        /*
        cout<<"L: " <<endl;
        for(int i=0; i<g.n; i++){
            cout<<i<<": ";
            for(int j=0; j<=i; j++){
                cout<<chol_L.coeff(i, j)<<" ";
            }
            cout<<endl;
        }

        cout<<"D: " <<endl;
        for(int i=0; i<g.n; i++){
            cout<<i<<": ";
            for(int j=0; j<=i; j++){
                cout<<chol_D.coeff(i, j)<<" ";
            }
            cout<<endl;
        }
        */
        //saving the indices;
        ofstream chol_L_file(chol_L_filename);

        for(int i=0; i<g.n; i++){
            for(int j=0; j<=i; j++){
                if(chol_L.coeff(i,j) != 0.0){
                    if(i == g.n -1 && j == i)
                        chol_L_file<< i<<" " <<j <<" " << chol_L.coeff(i, j);
                    else
                        chol_L_file<< i<<" " <<j <<" " << chol_L.coeff(i, j)<<endl;
                }
            }
        }
        chol_L_file.close();

        ofstream chol_D_file(chol_D_filename);

        for(int i = 0; i < g.n; i++){
            if(i == g.n - 1)
                chol_D_file<<chol_D.coeff(i, i);
            else
                chol_D_file << chol_D.coeff(i, i) <<endl;
        }
        chol_D_file.close();

        double total_time = permute_time + decom_time;
        return total_time;

    }

    //complete cholesky decomp
    void cholesky_decomp(SparseMtx& W_mtrx, SparseMtx& chol_L, SparseMtx& chol_D){
        //store the idx whose W entry is non zero.
        //for each row of this matrix

        vector<vector<pair<int, double> > > temp_L_values;//row id i: j, L_ij which is non zero

        double * temp_D_values =  new double[g.n]; // i, i
        for(int i=0; i<g.n ; i++){
            temp_D_values[i] = 0.0;
            vector<pair<int, double>> temp_i;
            temp_L_values.push_back(temp_i);
        }

        vector<T> L_tripletList; //chol_L
        vector<T> D_tripletList; //chol_D
        L_tripletList.reserve(g.n*g.n);
        D_tripletList.reserve(g.n);

        //start decomposing from each row
        for(int i=0; i<g.n; i++){
            if(i == 0){
                double W_ii = W_mtrx.coeff(i,i);
                cout<< W_ii <<endl;
                temp_L_values[i].push_back(make_pair(i, 1.0));
                L_tripletList.push_back(T(i,i, 1));

                temp_D_values[i] = W_ii;
                D_tripletList.push_back(T(i,i, W_ii));
                //chol_D(i,i) = W_ii; //D_11 = W_11 ,here 1 is 0
                //chol_L(i,i) = 1;
            }
            else{
                for(int j=0; j<= i; j++){
                    int real_j = j;
                    if(real_j == i){//L_ii and D_ii
                        //chol_L(i, real_j) = 1;
                        temp_L_values[i].push_back(make_pair(real_j, 1.0));
                        L_tripletList.push_back(T(i, real_j, 1.0));

                        double W_ij = W_mtrx.coeff(i, real_j);
                        //cout<< W_ij <<endl;
                        double D_ij = 0;


                        for(int k = 0; k<temp_L_values[real_j].size(); k++){
                            pair<int, double> temp_pair = temp_L_values[real_j][k];
                            int real_k = temp_pair.first;
                            if(real_k < real_j){
                                D_ij += temp_pair.second * temp_pair.second * temp_D_values[real_k];
                            }
                            else{
                                break;
                            }
                        }

                        temp_D_values[real_j] = W_ij - D_ij;
                        D_tripletList.push_back(T(i, real_j, W_ij-D_ij));
                        //cout<<temp_D_values[real_j]<<endl;
                        //chol_D(i, real_j) = W_mtrx(i, real_j) - D_ij;

                    }
                    else{// real_j < i: compute L_ij only
                        double W_ij = W_mtrx.coeff(i, real_j);
                        //cout<<i<< " "<<j <<" " << W_ij <<endl;
                        double L_ij = 0;

                        unordered_map<int, double> L_in_i; //the values of L_ik where k<j-1
                        for(int k=0; k<temp_L_values[i].size(); k++){
                            pair<int, double> temp_pair = temp_L_values[i][k];
                            int real_k = temp_pair.first;
                            if(real_k < real_j){
                                L_in_i.insert(make_pair(real_k, temp_pair.second));
                            }
                        }

                        for(int k=0; k<temp_L_values[real_j].size(); k++){
                            pair<int, double> temp_pair = temp_L_values[real_j][k];
                            int real_k = temp_pair.first;
                            if((real_k < real_j) && (L_in_i.count(real_k)>0)){
                                L_ij += L_in_i[real_k] * temp_pair.second * temp_D_values[real_k];
                            }
                        }

                        L_ij = (W_ij - L_ij) / temp_D_values[real_j];

                        if(L_ij != 0.0){
                            temp_L_values[i].push_back(make_pair(real_j, L_ij));
                        }
                        L_tripletList.push_back(T(i, real_j, L_ij));

                        //cout<<L_ij<<endl;
                    }
                }
            }

        }

        chol_L.setFromTriplets(L_tripletList.begin(), L_tripletList.end());
        chol_D.setFromTriplets(D_tripletList.begin(), D_tripletList.end());

        delete [] temp_D_values;
    }

    double pre_chol_decomp(string chol_L_filename, string chol_D_filename){

        SparseMtx adj_mtrx(g.n, g.n);
        SparseMtx diag_mtrx(g.n, g.n);
        convert_graph_to_matrix_eigen(adj_mtrx, diag_mtrx);
        SparseMtx per_mtrx(g.n, g.n);
        get_permutation_matrix(per_mtrx);
        SparseMtx per_trans = per_mtrx.transpose();

        clock_t permute_start = clock();
        //doing permutation on the matrices
        adj_mtrx = (per_mtrx * adj_mtrx).pruned();
        adj_mtrx = (adj_mtrx * per_trans).pruned();

        diag_mtrx = (per_mtrx * diag_mtrx).pruned();
        diag_mtrx = (diag_mtrx * per_trans).pruned();

        SparseMtx W_mtrx(g.n, g.n);
        W_mtrx = (diag_mtrx * adj_mtrx).pruned();
        W_mtrx = (W_mtrx * diag_mtrx).pruned();

        SparseMtx I (g.n, g.n);
        I.setIdentity();
        W_mtrx = (I - alpha * W_mtrx).pruned();

        /*cout<<"W: " <<endl;
        for(int i=0; i<g.n; i++){
            cout<<i<<": ";
            for(int j=0; j<=i; j++){
                cout<<W_mtrx.coeff(i, j)<<" ";
            }
            cout<<endl;
        }*/


        clock_t  permute_end = clock();
        double permute_time = (permute_end - permute_start) /(double) CLOCKS_PER_SEC;
        cout<<"permute time: " << permute_time <<endl;

        cout<<"decomposing ... " <<endl;

        clock_t decom_start = clock();
        SparseMtx chol_L(g.n, g.n);
        SparseMtx chol_D(g.n, g.n);


        //using complete cholesky decomposition
        cholesky_decomp(W_mtrx,chol_L, chol_D);

        clock_t decom_end = clock();
        double decom_time = (decom_end - decom_start) / (double) CLOCKS_PER_SEC;
        cout<<"decomposing time: " << decom_time <<endl;

        /*cout<<"L: " <<endl;
        for(int i=0; i<g.n; i++){
            cout<<i<<": ";
            for(int j=0; j<=i; j++){
                cout<<chol_L.coeff(i, j)<<" ";
            }
            cout<<endl;
        }

        cout<<"D: " <<endl;
        for(int i=0; i<g.n; i++){
            cout<<i<<": ";
            for(int j=0; j<=i; j++){
                cout<<chol_D.coeff(i, j)<<" ";
            }
            cout<<endl;
        }*/

        //saving the indices;
        ofstream chol_L_file(chol_L_filename);

        for(int i=0; i<g.n; i++){
            for(int j=0; j<=i; j++){
                if(chol_L.coeff(i,j) != 0.0){
                    chol_L_file<< i<<" " <<j <<" " << chol_L.coeff(i, j)<<endl;
                }
            }
        }
        chol_L_file.close();

        ofstream chol_D_file(chol_D_filename);

        for(int i = 0; i < g.n; i++){
            chol_D_file << chol_D.coeff(i, i) <<endl;
        }
        chol_D_file.close();

        double total_time = permute_time + decom_time;
        return total_time;

    }
//the matrix-vector multiplications by using eigen lib
    void convert_graph_to_matrix_eigen(SparseMtx& adj_mtrx, SparseMtx& diag_mtrx){
        //store the graph in SparseMatrix

        vector<T> tripletList;
        vector<T> diag_tripletList;
        tripletList.reserve(g.m);
        diag_tripletList.reserve(g.n);

        for (int i = 0; i < g.n; i++) {
            int outdeg = g.getOutSize(i);
            double sum_weight = 0.0;
            //cout<< i<<endl;
            if (outdeg != 0) {
                //double trans_prob = 1 / (double) outdeg;
                for (int j = 0; j < outdeg; j++) {
                    pair<int, double> neighbor = g.outAdjList[i][j];
                    int neigh_id = neighbor.first;
                    double weight = neighbor.second;
                    //double W_ij = (-1) * alpha * weight;
                    tripletList.push_back(T(i, neigh_id, weight));
                    sum_weight += weight;
                }
            }
            diag_tripletList.push_back(T(i, i, 1 / sqrt(sum_weight)));
        }
        adj_mtrx.setFromTriplets(tripletList.begin(), tripletList.end());
        diag_mtrx.setFromTriplets(diag_tripletList.begin(), diag_tripletList.end());

    }

    void get_permutation_matrix(SparseMtx& per_mtrx){
        vector<T> tripletList;
        tripletList.reserve(g.n);

        string permute_filename = "./index/" + target_filename + "/permutation.txt";
        ifstream permute_file(permute_filename);

        for(int i=0; i<g.n; i++){
            int new_id;
            int previous_id;
            int value;

            permute_file>> new_id >> previous_id >> value;
            tripletList.push_back(T(new_id, previous_id, value));
        }
        permute_file.close();
        per_mtrx.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void get_permute_node_id(SparseMtx& per_mtrx, vector<int>& new_node_ids, vector<int>& pre_node_ids){
        //unordered_map<int, int> new_node_ids;// node_id, new_id

        vector<T> tripletList;
        tripletList.reserve(g.n);

        string permute_filename = "./index/" + target_filename + "/permutation.txt";
        ifstream permute_file(permute_filename);

        for(int i=0; i<g.n; i++){
            int new_id;
            int previous_id;
            int value;

            permute_file>> new_id >> previous_id >> value;
            tripletList.push_back(T(new_id, previous_id, value));
            new_node_ids[previous_id] = new_id;
            pre_node_ids[new_id] = previous_id;
           // new_node_ids.insert(make_pair(previous_id, new_id));
        }
        permute_file.close();
        per_mtrx.setFromTriplets(tripletList.begin(), tripletList.end());

        //return new_node_ids;
    }

    void get_chol_L_matrix(SparseMtx& chol_L, string filename){
        vector<T> tripletList;
        tripletList.reserve(g.n*g.n);

        ifstream infile(filename);

        while(infile.good()){
            int row_id;
            int col_id;
            double value;

            infile >> row_id >> col_id >> value;
            tripletList.push_back(T(row_id, col_id, value));
        }
        infile.close();
        chol_L.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void get_chol_D_matrix(SparseMtx& chol_D, string filename){
        vector<T> tripletList;
        tripletList.reserve(g.n);

        ifstream infile(filename);

        for(int i=0; i<g.n; i++){
            double value;
            infile>>value;
            tripletList.push_back(T(i, i, value));
        }
        infile.close();
        chol_D.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    vector<vector<int>> get_clusters(vector<int>& node_cluster, string filename, string general_filename, bool& totally_discon){
        int num_clusters;
        ifstream reorder_cluster_file(filename);
        reorder_cluster_file >> num_clusters;

        //cout<< num_clusters << endl;
        vector<vector<int>> cluster_nodes;
        cluster_nodes = vector<vector<int> >(num_clusters, vector<int>());

        for(int i=0; i < g.n; i++){
            int temp_node;
            int temp_cluster;
            reorder_cluster_file >> temp_cluster >> temp_node;
            cluster_nodes[temp_cluster].push_back(temp_node);
            node_cluster[temp_node] = temp_cluster;
        }
        reorder_cluster_file.close();

        ifstream cluster_file(general_filename);
        int num_clusters_before_reorder;
        cluster_file >> num_clusters_before_reorder;
        cluster_file.close();

        if(num_clusters == num_clusters_before_reorder){
            totally_discon = true;
        }

        return cluster_nodes;
    }
    void symmetric_normalization(){
        SparseMtx adj_mtrx(g.n, g.n);
        SparseMtx diag_mtrx(g.n, g.n);
        convert_graph_to_matrix_eigen(adj_mtrx, diag_mtrx);

        SparseMtx norm_mtrx(g.n, g.n);
        norm_mtrx = (diag_mtrx * adj_mtrx).pruned();
        norm_mtrx = (norm_mtrx * diag_mtrx).pruned();
        stringstream ss;
        ss << "./dataset/sym_" << target_filename << ".txt";
        string filename = ss.str();
        ofstream outfile(filename);
        outfile << g.n <<endl;
        for (int i=0; i<norm_mtrx.outerSize(); i++)
            for (SparseMatrix<double>::InnerIterator it(norm_mtrx,i); it; ++it)
            {
                outfile << it.row() << " " << it.col() <<" "<< it.value() <<endl;
            }

        outfile.close();
    }

    void comparison_lu_chol_eigen(int source){

        SparseMtx adj_mtrx(g.n, g.n);
        SparseMtx diag_mtrx(g.n, g.n);
        convert_graph_to_matrix_eigen(adj_mtrx, diag_mtrx);
        SparseMtx per_mtrx(g.n, g.n);
        get_permutation_matrix(per_mtrx);
        SparseMtx per_trans = per_mtrx.transpose();

        clock_t permute_start = clock();
        //doing permutation on the matrices
        adj_mtrx = (per_mtrx * adj_mtrx).pruned();
        adj_mtrx = (adj_mtrx * per_trans).pruned();

        diag_mtrx = (per_mtrx * diag_mtrx).pruned();
        diag_mtrx = (diag_mtrx * per_trans).pruned();

        SparseMtx norm_mtrx(g.n, g.n);
        norm_mtrx = (diag_mtrx * adj_mtrx).pruned();
        norm_mtrx = (norm_mtrx * diag_mtrx).pruned();
        clock_t  permute_end = clock();
        double permute_time = (permute_end - permute_start) /(double) CLOCKS_PER_SEC;


        cout<< adj_mtrx.nonZeros() <<endl;
        cout<< diag_mtrx.nonZeros() <<endl;
        cout<< norm_mtrx.nonZeros() <<endl;

        clock_t lu_start_time = clock();
        SparseMtx I (g.n, g.n);
        I.setIdentity();
        norm_mtrx = (I - alpha * norm_mtrx).pruned();
        SparseLU<SparseMtx> ldlt(norm_mtrx);
        SparseMtx inv = ldlt.solve(I);
         inv = (1-alpha) * inv;
        cout<< inv.nonZeros()<<endl;

        VectorXd x_vec(g.n);

        for (int i=0; i<inv.outerSize(); i++)
            for (SparseMatrix<double>::InnerIterator it(inv,i); it; ++it)
            {
                if(it.col() == source) {
                    x_vec(it.row()) = it.value();
                    cout << it.row() << " " << it.value() << endl;
                }
            }
        clock_t lu_end_time = clock();
        double lu_time = (lu_end_time - lu_start_time) / (double) CLOCKS_PER_SEC;
        cout<< "lu time: " << lu_time <<endl;

        clock_t chol_start_time = clock();
        VectorXd q_vec (g.n);
        for(int i=0; i< g.n; i++){
            if(i == source)
                q_vec(i) = 1;
            else
                q_vec(i) = 0.0;
        }
        SimplicialCholesky<SparseMtx> chol(norm_mtrx);
        VectorXd x = chol.solve(q_vec);
        x = (1-alpha) * x;
        clock_t chol_end_time = clock();
        double chol_time = (chol_end_time - chol_start_time) / (double) CLOCKS_PER_SEC;
        cout<<"chol time:" << chol_time << endl;

        double abs_error = 0.0;
        for(int i=0; i< g.n; i++){
            cout<< i << " " << x(i) <<endl;
            if(abs(x(i) - x_vec(i)) > abs_error)
                abs_error = abs(x(i) - x_vec(i));
        }
        cout << "avg error " << abs_error <<endl;

    }

    void inverse_eigen(int source){

        SparseMtx adj_mtrx(g.n, g.n);
        SparseMtx diag_mtrx(g.n, g.n);
        convert_graph_to_matrix_eigen(adj_mtrx, diag_mtrx);
       // SparseMtx per_mtrx(g.n, g.n);
       // get_permutation_matrix(per_mtrx);
        //SparseMtx per_trans = per_mtrx.transpose();

        /*clock_t permute_start = clock();
        //doing permutation on the matrices
        adj_mtrx = (per_mtrx * adj_mtrx).pruned();
        adj_mtrx = (adj_mtrx * per_trans).pruned();

        diag_mtrx = (per_mtrx * diag_mtrx).pruned();
        diag_mtrx = (diag_mtrx * per_trans).pruned();
*/
        SparseMtx norm_mtrx(g.n, g.n);
        norm_mtrx = (diag_mtrx * adj_mtrx).pruned();
        norm_mtrx = (norm_mtrx * diag_mtrx).pruned();
        //clock_t  permute_end = clock();
        //double permute_time = (permute_end - permute_start) /(double) CLOCKS_PER_SEC;


        cout<< adj_mtrx.nonZeros() <<endl;
        cout<< diag_mtrx.nonZeros() <<endl;
        cout<< norm_mtrx.nonZeros() <<endl;

        clock_t lu_start_time = clock();
        SparseMtx I (g.n, g.n);
        I.setIdentity();
        norm_mtrx = (I - alpha * norm_mtrx).pruned();
        SparseLU<SparseMtx> ldlt(norm_mtrx);
        SparseMtx inv = ldlt.solve(I);
        inv = (1-alpha) * inv;
        cout<< inv.nonZeros()<<endl;

        VectorXd q_vec (g.n);
        for(int i=0; i< g.n; i++){
            if(i == source)
                q_vec(i) = 1;
            else
                q_vec(i) = 0.0;
        }

        VectorXd x_vec(g.n);

        for (int i=0; i<inv.outerSize(); i++)
            for (SparseMatrix<double>::InnerIterator it(inv,i); it; ++it)
            {
                if(it.col() == source) {
                    x_vec(it.row()) = it.value();
                    cout << it.row() << " " << it.value() << endl;
                }
            }
        clock_t lu_end_time = clock();
        double lu_time = (lu_end_time - lu_start_time) / (double) CLOCKS_PER_SEC;
        cout<< "inv time with help of lu: " << lu_time <<endl;

    }

};

void mr_t_PowerMethod(MR* mr, vector<int> nodeList, int iterations){
    return mr->t_PowerMethod(nodeList, iterations);
}

void MR::PowerMethodMulti(int iterations, int node_count, int num_thread){
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
    long start = ((long)t_start.tv_sec)*1000+(long)t_start.tv_usec/1000;
    string inputFile = "./dataset/" + graph_name + ".query";
    ifstream node_file(inputFile);
    vector<int> nodes;
    for(int i = 0; i < node_count; i++){
        int temp_node;
        node_file >> temp_node;
        if(g.getOutSize(temp_node) == 0){
            i--;
            cout << "illegal : " << temp_node << endl;
            continue;
        }
        nodes.push_back(temp_node);
    }
    node_file.close();

    vector<thread> threads;
    for(int i = 0; i < num_thread-1; i++){
        vector<int> t_nodes;
        for(int j = 0; j < node_count / num_thread; j++){
            t_nodes.push_back(nodes[i * node_count / num_thread + j]);
        }
        threads.push_back(thread(mr_t_PowerMethod, this, t_nodes, iterations));
    }
    vector<int> t_nodes;
    for(int j = 0; j < node_count / num_thread; j++){
        t_nodes.push_back(nodes[(num_thread-1) * node_count / num_thread + j]);
    }
    t_PowerMethod(t_nodes, iterations);
    for (int i = 0; i < num_thread - 1; i++){
        threads[i].join();
    }
    gettimeofday(&t_end, NULL);
    long end = ((long)t_end.tv_sec)*1000+(long)t_end.tv_usec/1000;
    int cost_time = end - start;

    cout << "cost: " << cost_time / (double) 1000 << endl;
}
#endif //TOPPPR_MANIFOLDRANKING_H
