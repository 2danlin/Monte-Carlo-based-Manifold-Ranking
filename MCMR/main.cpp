#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "manifoldRanking.h"
#include "algo.h"
#include "Graph.h"
//#include <stdio.h>
#include <cstdlib>
#include <unordered_set>    
#include <cstdlib>
#include <math.h>
#include <cstring>
//#include "SparseMatrix.h"

void usage() {
    cerr << "MCMR [-d <dataset>] [-e epsilon] [-k topk (default 1)] [-n node_count (default 1)] [-a alpha (default 0.99)] [-t power_t (default 1000)] [-c mcmr_c (default 100)]" << endl;
}

int check_inc(int i, int max) {
    if (i == max) {
        usage();
        exit(1);
    }
    return i + 1;
}

int main(int argc, char *argv[]){
    int i = 1;
    char *endptr;
    string filename;
    string graph_name;
    double eps=0.1;                   //relative error
    int k = 1;                      //parameter k
    double alpha = 0.99;            //propagation probability
    int node_count = 1;            //query node size
    double forward_c = 1;
    int power_t = 1000;             //number of iterations for power method / ground truth
    int mcmr_c = 100;                //parameter c for mcmr to decide the number of walks
    double gap_k = 1/(double) mcmr_c;           //1/c
    string algo = "MCMR";
    if(argc < 4){
        usage();
        exit(1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-d")) {
            i = check_inc(i, argc);
            filename = argv[i];
        } 
        else if (!strcmp(argv[i], "-algo")) {
            i = check_inc(i, argc);
            algo = argv[i];
        }
        else if (!strcmp(argv[i], "-e")) {
            i = check_inc(i, argc);
            eps = strtod(argv[i], &endptr);
            if ((eps == 0 || eps > 1) && endptr) {
                cerr << "Invalid eps argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-k")) {
            i = check_inc(i, argc);
            k = strtod(argv[i], &endptr);
            if ((k < 0) && endptr) {
                cerr << "Invalid k argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-n")) {
            i = check_inc(i, argc);
            node_count = strtod(argv[i], &endptr);
            if ((node_count < 0) && endptr) {
                cerr << "Invalid node_count argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-t")) {
            i = check_inc(i, argc);
            power_t= strtod(argv[i], &endptr);
            if ((power_t < 0) && endptr) {
                cerr << "Invalid  power_t argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-a")) {
            i = check_inc(i, argc);
            alpha = strtod(argv[i], &endptr);
            if (((alpha < 0) || (alpha > 1)) && endptr) {
                cerr << "Invalid alpha argument" << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-c")) {
            i = check_inc(i, argc);
            mcmr_c = strtod(argv[i], &endptr);
            if ((k < 0) && endptr) {
                cerr << "Invalid c argument" << endl;
                exit(1);
            }
        }
        else {
            usage();
            exit(1);
        }
        i++;
    }

    //for the algorithms that to build the unnorm knn graph by exploiting the pre-computed large_knn

    if(algo == "BUILD_KNN_GRAPH"){
        cout<<"building the knn graph" <<endl;

        construct_knn_graph(filename, k);


        stringstream ss;
        ss<< filename <<"_" << k;
        string new_filename = ss.str();
        graph_name = new_filename;
        MR mr = MR(new_filename, graph_name, forward_c, k, alpha);
        mr.symmetric_normalization();
    }
    else{
        //for the algorithms that the knn graphs are already ready.
        graph_name = filename;
        cout<<"graph name: " << graph_name <<endl;
        if(algo != "MOGUL" && algo !="MOGUL_E" && algo != "PRE_MOGUL" && algo != "PRE_MOGUL_E" && algo != "SYM_NORM"){
            string folder_name = "sym_";
            folder_name += filename;
            filename = folder_name;
            cout<<"current filename:" << filename<<endl;
            //filename = strcpy("unnorm_graphs/", filename.c_str());
        }

        MR mr = MR(filename, graph_name, forward_c, k, alpha);
        if(algo == "GEN_QUERY"){
            ofstream outFile("dataset/" + graph_name + ".query");
            //ofstream outFile2("dataset/unnorm_graphs/" + filename + ".query");
            mr.generateQueryNode(node_count, outFile);
            outFile.close();
            //outFile2.close();
        }
        else if(algo == "GEN_GROUND_TRUTH"){
            cout<<"starting generating ground truth ... "<<endl;
            mr.PowerMethodMulti(power_t, node_count, 10);
            /*  多线程PowerMethparameter: iteration loops, node size, thread num */

            stringstream ss;
            ss << "estimated_mr/" << graph_name << "/POWER_GRT.txt";
            string time_filename = ss.str();
            ofstream time_file(time_filename);
            time_file <<"No. of iterations: " << power_t <<endl;
            time_file <<"avg query time: " << mr.avg_time / (double) node_count * 10  << endl;
            time_file <<"avg precision: " << 1.0 << endl;
            time_file.close();
        }
        else if(algo =="POWER") {
            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            while (!nodes_file.eof()) {
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
            }
            cout << "read done!" << endl;
            for (int para = 0; para < 1; para++) {
                if (k > mr.g.n) {
                    break;
                }
                int realCount = 0;
                for (int t = 0; t < node_count; t++) {
                    int test_node = test_nodes[realCount++];
                    cout << t + 1 << "-th node: " << test_node << " " << eps << " " << k << endl;
                    //power
                    mr.Power(power_t, test_node, k);
                }

                cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/POWER_" << k << ".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file << "k: " << k << endl;
                time_file << "avg query time: " << mr.avg_time / (double) node_count * 10 << endl;
                time_file << "avg precision: " << 1.0 << endl;
                time_file.close();
            }
        }
        else if(algo == "DATA_PREPROCESS"){
            mr.saveSumWeight();
            mr.saveNormWeight();
        }
        else if(algo == "MCMRPlus") {
            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            while (!nodes_file.eof()) {
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
            }
            cout << "read done!" << endl;
            for (int para = 0; para < 1; para++) {
                if (k > mr.g.n) {
                    break;
                }
                int realCount = 0;
                for (int t = 0; t < node_count; t++) {
                    int test_node = test_nodes[realCount++];
                    cout<< t+1 << "-th node: " << test_node << " " << eps << " " << k << endl;

                    double real_eps =
                            eps * sqrt(1 / (double) mr.g.m / (double) mr.g.n / log(mr.g.n)) / log(k);
                    cout << "r_max: " << real_eps << endl;
                    //MCMR+
                    mr.MCMRPlus(test_node, real_eps, k);
                }
                //cout << "avg random walk num: " << mr.avg_walk_num / (double) node_count << endl;
                cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;
                cout <<"avg walk num: " << mr.avg_walk_num/ (double) node_count <<endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/MCMRPlus_"<<k<<".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file <<"k: " << k <<endl;
                time_file <<"avg query time: " << mr.avg_time / (double) node_count  << endl;
                time_file <<"avg precision: " << mr.avg_pre / (double) node_count  << endl;
                time_file <<"avg walk num: " << mr.avg_walk_num/ (double) node_count <<endl;
                time_file.close();

                mr.avg_time = 0;
                mr.avg_pre = 0;
                mr.avg_walk_num =0;
            }
            for (int i = 0; i < argc; i++) {
                cout << argv[i] << " ";
            }
            cout << endl;
        }else if(algo == "GET_MAX_SUM"){
            double max = mr.g.getMaxSumWeights();
            cout<< "max sum weight: " << max << endl;
            cout<< "alpha * max = " << max*alpha << endl;
        }
        else if(algo == "MCMR"){
            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            while(!nodes_file.eof()){
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
            }
            cout << "read done!" << endl;

            int c[1] = {10000000};

            for(int para = 0; para < 1; para++){
                if(k > mr.g.n){
                    break;
                }
                //mcmr_c = c[para];
                cout<<"current c: " << mcmr_c <<endl;

                int realCount = 0;
                for(int t = 0; t < node_count; t++){
                    int test_node = test_nodes[realCount++];
                    cout << t+1 << "-th node: " << test_node << " " << eps << " " << k << endl;
                    //double relative_error = 0.5;
                    //double gap = 0.0001;
                    gap_k = 1/(double) mcmr_c;
                    unsigned int walk_num = 10 * log(mr.g.n)/ 3.0 /gap_k;
                    //unsigned int maximum_walk = 3 * mr.g.n * log(2* mr.g.n) / eps / eps;
                    //if ( walk_num > maximum_walk){
                    //    walk_num = maximum_walk;
                    //}
                    cout<< "num walks: " << walk_num << endl;
                    cout<< walk_num << endl;
                    mr.MonteCarlo(test_node, walk_num);
                }
                cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/MCMR_"<<k << "_" << mcmr_c <<".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file <<"k: " << k <<endl;
                time_file <<"avg query time: " << mr.avg_time / (double) node_count  << endl;
                time_file <<"avg precision: " << mr.avg_pre / (double) node_count  << endl;
                time_file.close();
                mr.avg_time = 0;
                mr.avg_pre = 0;
            }
            for(int i = 0; i < argc; i++){
                cout << argv[i] << " ";
            }
            cout << endl;
        }
        else if(algo == "MCMR_TOPK"){
            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            while(!nodes_file.eof()){
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
            }
            cout << "read done!" << endl;
            for(int para = 0; para < 1; para++){
                if(k > mr.g.n){
                    break;
                }
                int realCount = 0;
                for(int t = 0; t < node_count; t++){
                    int test_node = test_nodes[realCount++];
                    cout << t+1 << "-th node: " << test_node << " " << eps << " " << k << endl;
                    //double relative_error = 0.5;
                    //double gap = 0.0001;
                    gap_k = 1/(double) mcmr_c;
                    //unsigned int walk_num = 10*log(mr.g.n)/ 3.0 /gap_k;
                    //cout<< "num walks: " << walk_num << endl;
                    //cout<< walk_num << endl;
                    mr.MonteCarlo_with_prune(test_node, k);
                }
                cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/MCMR_TOPK_"<<k <<".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file <<"k: " << k <<endl;
                time_file <<"avg query time: " << mr.avg_time / (double) node_count  << endl;
                time_file <<"avg precision: " << mr.avg_pre / (double) node_count  << endl;
                time_file.close();
                mr.avg_time = 0;
                mr.avg_pre = 0;
            }
            for(int i = 0; i < argc; i++){
                cout << argv[i] << " ";
            }
            cout << endl;
        }
        else if(algo =="SYM_NORM"){

            mr.symmetric_normalization();
        } else if(algo == "INV"){

            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            while(!nodes_file.eof()){
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
            }
            cout << "read done!" << endl;
            cout << test_nodes.size() << endl;
            for(int para = 0; para < 1; para++){
                if(k > mr.g.n){
                    break;
                }
                //int realCount = 0;
                for(int t = 0; t < node_count; t++){
                    int test_node = test_nodes[t];
                    cout << t+1 << "-th node: " << test_node << " " << k << endl;

                    mr.inverse_eigen(test_node);

                }
                //cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/inverse_"<<k<<".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file <<"k: " << k <<endl;
                time_file <<"avg query time: " << mr.avg_time / (double) node_count  << endl;
                //time_file <<"avg precision: " << mr.avg_pre / (double) node_count  << endl;
                time_file.close();
                mr.avg_time = 0;
                //mr.avg_pre = 0;
            }
            for(int i = 0; i < argc; i++){
                cout << argv[i] << " ";
            }
            cout << endl;

        }
        else if(algo == "PRE_MOGUL"){
            mr.get_unweighted_graph();
            mr.preprocessing_mogul();
        }
        else if(algo == "PRE_MOGUL_E"){
            mr.get_unweighted_graph();
            mr.preprocessing_mogul_exact();
        }
        else if(algo == "MOGUL") {
            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            int counter = 0;
            while (!nodes_file.eof() && counter < node_count) {
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
                counter++;
            }
            cout << "read done!" << endl;
            for (int para = 0; para < 1; para++) {
                if (k > mr.g.n) {
                    break;
                }

                mr.mogul_in_query_phase(test_nodes, algo);
                cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout<<"total pre: " << mr.avg_pre <<endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/mogul_"<<k<<".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file <<"k: " << k <<endl;
                time_file <<"avg query time: " << mr.avg_time / (double) node_count  << endl;
                time_file <<"avg precision: " << mr.avg_pre / (double) node_count  << endl;
                time_file.close();

                mr.avg_time = 0;
                mr.avg_pre = 0;
            }
            for (int i = 0; i < argc; i++) {
                cout << argv[i] << " ";
            }
            cout << endl;
        }
        else if(algo == "MOGUL_E") {
            ifstream nodes_file("dataset/" + graph_name + ".query");
            vector<int> test_nodes;
            int counter = 0;
            while (!nodes_file.eof() && counter < node_count) {
                int temp_node;
                nodes_file >> temp_node;
                test_nodes.push_back(temp_node);
                counter++;
            }
            cout << "read done!" << endl;
            for (int para = 0; para < 1; para++) {
                if (k > mr.g.n) {
                    break;
                }

                mr.mogul_in_query_phase(test_nodes, algo);
                cout << "avg precision: " << mr.avg_pre / (double) node_count << endl;
                cout<<"total pre: " << mr.avg_pre <<endl;
                cout << "avg time: " << mr.avg_time / (double) node_count << endl;

                stringstream ss;
                ss << "estimated_mr/" << graph_name << "/mogul_E_"<<k<<".txt";
                string time_filename = ss.str();
                ofstream time_file(time_filename);
                time_file <<"k: " << k <<endl;
                time_file <<"avg query time: " << mr.avg_time / (double) node_count  << endl;
                time_file <<"avg precision: " << mr.avg_pre / (double) node_count  << endl;
                time_file.close();

                mr.avg_time = 0;
                mr.avg_pre = 0;
            }
            for (int i = 0; i < argc; i++) {
                cout << argv[i] << " ";
            }
            cout << endl;
        }
    }



    return 0;
};
