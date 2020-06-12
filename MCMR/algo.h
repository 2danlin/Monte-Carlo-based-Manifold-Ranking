//
// Created by 林丹丹 on 2020/6/2.
//

#ifndef TOPPPR_ALGO_H
#define TOPPPR_ALGO_H


#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <queue>
#include "weightedGraph.h"
#include "Graph.h"
#include <unordered_map>
#include <unordered_set>
#include <time.h>
#include <limits>
#include "heap.h"
#include <functional>

using namespace std;

typedef struct cluster cluster;
struct cluster{
    int id; //the cluster id for this cluster
    int self_loop_edge_weight; //the weight of the self-loop edge of this cluster

    int outer_edge_weight; //the weight of edges that incident to its neighbors

    //vector<int> nodes_in_cluster;//the nodes in this cluster
    unordered_set<int> nodes_in_cluster;

   // cluster* adj_clusters; //the clusters that incidents this cluster
    //vector<int> cluster_incident; //the id of clusters that incidents this cluster
    unordered_map<int, int> cluster_incident; //the id of clusters that incident this cluster, and the weight between them.
};

class pqcompare
{
    bool reverse;
public:
    pqcompare(const bool& revparam=false)
    {reverse=revparam;}
    bool operator() (const pair<int, double>& lhs, const pair<int, double>&rhs) const
    {
        if (reverse) return (lhs.second > rhs.second);
        else return (lhs.second < rhs.second);
    }
};

typedef priority_queue<pair<int, double>, vector<pair<int,double> >, pqcompare> pq;

vector<int> final_cluster; //node_id, cluster_id
vector<cluster> target_set;//containing the clusters whose delta Q are to be computed
//int infinity = numeric_limits<int>::max();

cluster merge_two_clusters(cluster c1, cluster c2){
    // update the self-loop weight by checking how many the number of edges between this two clusters
    int weight_c1_c2 = c1.cluster_incident[c2.id];
    //cout<<"weight between two cluster: " <<weight_c1_c2 <<endl;
    //cout<<"previous self loop: " << c1.self_loop_edge_weight <<" " << c2.self_loop_edge_weight<<endl;
    c1.self_loop_edge_weight = c1.self_loop_edge_weight + c2.self_loop_edge_weight + weight_c1_c2;

    //cout<<"new self loop: " << c1.self_loop_edge_weight<<endl;
    //udpdate the outer-weight
   // cout<<"previous outer edge: " << c1.outer_edge_weight <<" " << c2.outer_edge_weight<<endl;
    c1.outer_edge_weight = c1.outer_edge_weight + c2.outer_edge_weight - 2 * weight_c1_c2;

    //cout<<"new outer edge: " << c1.outer_edge_weight <<endl;
    //put the nodes in c2 into c1
    for(auto c2_iter = c2.nodes_in_cluster.begin(); c2_iter != c2.nodes_in_cluster.end(); ++c2_iter){
        int temp_node = * c2_iter;
        c1.nodes_in_cluster.insert(temp_node);//only insert if there is no the same node in this cluster.
        final_cluster[temp_node] = c1.id; //update it to the latest cluster id
    }


    //put the incident clusters into c1;
    //cout<<"cluster in c2"<<endl;
    for(auto c2_iter = c2.cluster_incident.begin(); c2_iter != c2.cluster_incident.end(); ++c2_iter){
        pair<int, int> temp_pair = * c2_iter;
        int temp_cluster = temp_pair.first;
        int temp_weight = temp_pair.second;

       // cout<<temp_cluster <<" " <<temp_weight <<endl;

        if(c1.cluster_incident.count(temp_cluster) >0){
            //cluster already incident to c1
            c1.cluster_incident[temp_cluster] += temp_weight;
            target_set[temp_cluster].cluster_incident[c1.id] += temp_weight;

        }
        else if (temp_cluster != c1.id){
            c1.cluster_incident.insert(temp_pair);
            target_set[temp_cluster].cluster_incident.insert(make_pair(c1.id, temp_weight));
        }
        target_set[temp_cluster].cluster_incident.erase(c2.id);
    }
    c1.cluster_incident.erase(c2.id);

    /*
    int sum_out_weight = 0;
    cout<<"cluster in new c1:" <<endl;
    for(auto c1_iter = c1.cluster_incident.begin(); c1_iter != c1.cluster_incident.end(); ++c1_iter){
        pair<int, int> temp_pair = * c1_iter;
        sum_out_weight += temp_pair.second;
        cout<<temp_pair.first <<" " << temp_pair.second<<endl;
    }

    if(sum_out_weight != c1.outer_edge_weight){
        cout<<"error error error" <<endl;
        exit(1);
    }*/

    return c1;
}

//compute delta Q: e_uv, a_u, a_v, # of edges in the graph
double compute_detal_Q(int weight_c1_c2, int a_c1, int a_c2, int m){
    double delta_Q = weight_c1_c2 / (double) m - a_c1 / (double) m * a_c2 / (double) m / 2.0;

    return delta_Q;
}
//doing the clustering and return the time cost
double graph_clustering(Graph& ug, string cluster_file){
    //here g is an undirected graph and unweighted graph
    final_cluster.clear();
    final_cluster.resize(ug.n);
    target_set.clear();
    target_set.reserve(ug.n);

    static BinaryHeap<int, less<long long> > min_heap_outer_weight(ug.n, less<long long>());//find out the min degree in the current target set
    //initialize 2 types of sets of nodes.
    //vector<cluster> prunable_set; //the clusters that could be ignored
    min_heap_outer_weight.clear();


    bool* in_target_set = new bool[ug.n];
    int cur_num_clusters_in_target_set = ug.n;

    clock_t start_time = clock();
    //initially, put each node into this cluster
    for(int i = 0; i < ug.n; i++){
        cluster temp_cluster;
        temp_cluster.id = i;
        final_cluster[i] = i;
        temp_cluster.self_loop_edge_weight = 0;
        temp_cluster.outer_edge_weight = ug.getOutSize(i);
        //cout<<temp_cluster.outer_edge_weight<<endl;
        temp_cluster.nodes_in_cluster.insert(i);
        //cout<<"i: " << i<<" ";
        //if(ug.getOutSize(i)>0){
        for(int j=0; j < ug.getOutSize(i); j++){
            int temp_neighbor = ug.outAdjList[i][j];
            //cout<<temp_neighbor<<" ";
            temp_cluster.cluster_incident.insert(make_pair(temp_neighbor, 1));
        }
            //cout<<endl;
            // temp_cluster.adj_clusters =NULL;
        target_set.push_back(temp_cluster);
        if(ug.getOutSize(i) > 0){
            in_target_set[i] = true;
            // min_outer_weight.push(make_pair(temp_cluster.id, temp_cluster.outer_edge_weight));
            min_heap_outer_weight.insert(i, ug.getOutSize(i));
        }
        else{
            in_target_set[i] = false;
            cur_num_clusters_in_target_set--;
        }
        //}

    }
    //min_heap_outer_weight.display();

    int iter=0;
    cout<<"clustering ..." <<endl;

    while(cur_num_clusters_in_target_set > 0){

        //cout<<cur_num_clusters_in_target_set<<endl;
        iter++;
        //cout<<"iter: " << iter<<endl;

        //find out the prunable clusters and prune it from the target set
        //cout<<"setting prunable set ..." <<endl;
        for(int i=0; i < target_set.size(); i++){
            //cout<<"i " <<i<<endl;
            cluster temp_cluster =  target_set[i];
            //if this cluster is not removed
            //and if this cluster's outer weight is equal to 1
            if(in_target_set[i] && temp_cluster.cluster_incident.size() == 1){
                //cout<<"prune "<<i<<endl;
                //find out its neighbor and merge it into its neighbor
                //cout<<"outer edges of this cluster "<< temp_cluster.cluster_incident.size()<<endl;
                pair<int, int> neigh_id  = *temp_cluster.cluster_incident.begin();
                if(neigh_id.first <= i){
                    cluster temp_neighbor = target_set[neigh_id.first];
                    cluster new_cluster = merge_two_clusters(temp_neighbor, temp_cluster);//put this cluster into its neighbor
                    in_target_set[i] = false; //prune it from the target set
                    cur_num_clusters_in_target_set--;
                    //cout<< "merge cluster " << i << " into cluster " << neigh_id.first<<endl;
                    target_set[neigh_id.first] = new_cluster;

                    //min_heap_outer_weight.modify(i, infinity);
                    min_heap_outer_weight.delete_a_node(i);
                    //cout<< target_set[neigh_id.first].outer_edge_weight<<endl;
                    if(target_set[neigh_id.first].outer_edge_weight == 0){
                        in_target_set[neigh_id.first] = false;
                        cur_num_clusters_in_target_set--;
                       // min_heap_outer_weight.modify(neigh_id.first, infinity);
                        min_heap_outer_weight.delete_a_node(neigh_id.first);
                    }
                    else{
                        min_heap_outer_weight.modify(neigh_id.first, new_cluster.outer_edge_weight);
                    }
                }
                else{
                    cluster temp_neighbor = target_set[neigh_id.first];
                    cluster new_cluster = merge_two_clusters(temp_cluster, temp_neighbor);
                    in_target_set[neigh_id.first] = false; //prune it from the target set
                    cur_num_clusters_in_target_set--;
                    //cout<< "merge cluster " << neigh_id.first << " into cluster " << i<<endl;
                    target_set[i] = new_cluster;

                    //min_heap_outer_weight.modify(neigh_id.first, infinity);
                    min_heap_outer_weight.delete_a_node(neigh_id.first);
                    //cout<< target_set[neigh_id.first].outer_edge_weight<<endl;
                    if(target_set[i].outer_edge_weight == 0){
                        in_target_set[i] = false;
                        cur_num_clusters_in_target_set--;
                        //min_heap_outer_weight.modify(i, infinity);
                        min_heap_outer_weight.delete_a_node(i);
                    }
                    else{
                        min_heap_outer_weight.modify(i, new_cluster.outer_edge_weight);
                    }
                }

                //min_heap_outer_weight.modify(neigh_id.first, target_set[neigh_id.first].outer_edge_weight);
                //min_outer_weight.push(make_pair(temp_neighbor.id, temp_neighbor.outer_edge_weight));//since neighbor has been updated
            }
        }

        //min_heap_outer_weight.display();

        if(cur_num_clusters_in_target_set > 1 && min_heap_outer_weight.size() > 1){

            //select a cluster u from the current target set that has the smallest outer edge
            //cout<<"starting target pruning" <<endl;
            //int cluster_u_id = min_outer_weight.top().first;
            auto top = min_heap_outer_weight.extract_top();
            int cluster_u_id = top.second;

            //cout<<"u: " << cluster_u_id <<" " << top.first<<endl;
            cluster c_u = target_set[cluster_u_id];
            int a_u = c_u.self_loop_edge_weight + c_u.outer_edge_weight;
            pq max_delta_Q(pqcompare(false));

            for(auto u_neigh_iter = c_u.cluster_incident.begin(); u_neigh_iter != c_u.cluster_incident.end(); ++ u_neigh_iter){
                pair<int, int> u_neigh = * u_neigh_iter;
                //cout<<"u neighbor: " << u_neigh.first <<" " << u_neigh.second<<endl;
                if(in_target_set[u_neigh.first]){
                    cluster c_u_neigh = target_set[u_neigh.first];
                    int a_u_neigh = c_u_neigh.self_loop_edge_weight + c_u_neigh.outer_edge_weight;
                    double temp_delta_Q = compute_detal_Q(u_neigh.second, a_u, a_u_neigh, ug.m);
                    //cout<<"temp delta Q: " << u_neigh.first <<" " << temp_delta_Q<<endl;
                    max_delta_Q.push(make_pair(u_neigh.first, temp_delta_Q));
                }
            }

            int cluster_v_id = max_delta_Q.top().first;
            //cout<<"v: " << cluster_v_id <<endl;
            double delta_Q = max_delta_Q.top().second;

            if(delta_Q > 0){
                if(cluster_u_id <= cluster_v_id){
                    cluster new_c_u = merge_two_clusters(c_u, target_set[cluster_v_id]);
                    target_set[cluster_u_id] = new_c_u;
                    in_target_set[cluster_v_id] = false;
                    cur_num_clusters_in_target_set--;
                    //cout<< "merge cluster " << cluster_v_id << " into cluster " << cluster_u_id<<endl;
                    //min_heap_outer_weight.modify(cluster_v_id, infinity);
                    min_heap_outer_weight.delete_a_node(cluster_v_id);
                    //cout<< target_set[cluster_u_id].outer_edge_weight<<endl;
                    if(new_c_u.outer_edge_weight == 0){
                        in_target_set[new_c_u.id] = false;
                        cur_num_clusters_in_target_set--;
                        //min_heap_outer_weight.modify(cluster_u_id, infinity);
                        min_heap_outer_weight.delete_a_node(cluster_u_id);
                    }
                    else{
                        min_heap_outer_weight.modify(cluster_u_id, new_c_u.outer_edge_weight);
                    }

                }
                else {
                    cluster new_c_v = merge_two_clusters(target_set[cluster_v_id], c_u);
                    target_set[cluster_v_id] = new_c_v;
                    in_target_set[cluster_u_id] = false;
                    cur_num_clusters_in_target_set--;
                    //cout << "merge cluster " << cluster_u_id << " into cluster " << cluster_v_id << endl;

                    //min_heap_outer_weight.modify(cluster_u_id, infinity);
                    //cout<< target_set[cluster_v_id].outer_edge_weight<<endl;
                    min_heap_outer_weight.delete_a_node(cluster_u_id);
                    if (new_c_v.outer_edge_weight == 0) {
                        in_target_set[new_c_v.id] = false;
                        cur_num_clusters_in_target_set--;
                        //min_heap_outer_weight.modify(cluster_v_id, infinity);
                        min_heap_outer_weight.delete_a_node(cluster_v_id);
                    }
                    else {
                        min_heap_outer_weight.modify(cluster_v_id, new_c_v.outer_edge_weight);

                    }
                }
            }
            else{
                in_target_set[cluster_u_id] = false;
                cur_num_clusters_in_target_set--;
            }


        }
        else{//only one node in this graph
            for(int i=0; i < target_set.size(); i++){
                if(in_target_set[i]){
                    in_target_set[i] = false;
                    cur_num_clusters_in_target_set--;
                    min_heap_outer_weight.delete_a_node(i);
                }
            }
        }

        //min_heap_outer_weight.display();
    }

    //print out the final id for each node
    int max_cluster_id = 0;
    //cout<<"node_id      cluster_id "<<endl;
    for(int i=0; i<ug.n; i++){
        //cout<<i <<"    " <<final_cluster[i]<<endl;
        if(final_cluster[i] > max_cluster_id)
            max_cluster_id = final_cluster[i];
    }
    cout<<"max cluster id: " << max_cluster_id<<endl;

    //convert to the normal idx
    bool* assigned = new bool[max_cluster_id];
    for(int i=0; i<max_cluster_id; i++){
        assigned[i] = false;
    }
    unordered_map<int, int> cluster_new_idx;
    int pointer = 0;
    for(int i=0; i<ug.n; i++){
        int cluster_id = final_cluster[i];
        if(assigned[cluster_id] != true){
            //cout<<"sss " ;
            cluster_new_idx.insert(make_pair(cluster_id, pointer));
            pointer++;
            assigned[cluster_id] = true;
        }
    }
    //cout<<pointer<<endl;
    vector<vector<int> > cluster_node = vector<vector<int> > (pointer, vector<int>());
    for(int i=0; i< ug.n; i++){
        int cluster_id  = final_cluster[i];
        int new_id = cluster_new_idx[cluster_id];
        cluster_node[new_id].push_back(i);
    }

    clock_t end_time = clock();

    double time_cost = (end_time - start_time) / (double) CLOCKS_PER_SEC;

    ofstream c_file(cluster_file);
    c_file << pointer<<endl;
    for(int i=0; i<pointer; i++){
        for(int j=0; j<cluster_node[i].size(); j++){
            c_file<< i <<" " << cluster_node[i][j] <<endl;
        }
    }
    c_file.close();


    delete [] assigned;
    delete [] in_target_set;

    return time_cost;

}

void construct_knn_graph(string dataset_name, int k){

    string index_filename = "./dataset/large_knn/" + dataset_name + ".txt";
    ifstream index_file (index_filename);

    int num_nodes;
    index_file>>num_nodes;

    vector<vector<pair<int, double>> > knn_set =  vector<vector<pair<int, double>>>(num_nodes, vector<pair<int, double>>());
    vector<int> knn_counter(num_nodes);
    for(int i=0; i<num_nodes; i++){
        knn_counter[i] = 0;
    }

    while(index_file.good()){
        int temp_node;
        int temp_neighbor;
        double weight;
        index_file >> temp_node >> temp_neighbor >> weight;
        if(knn_counter[temp_node] < k){
            knn_set[temp_node].push_back(make_pair(temp_neighbor, weight));
            knn_counter[temp_node]++;
            //cout<<temp_node <<" has " << knn_counter[temp_node] << " neighbors already ..." <<endl;
        }
    }
    index_file.close();

    stringstream ss;
    ss<<"./dataset/" << dataset_name << "_" << k << ".txt";
    string knn_filename = ss.str();
    ofstream knn_file(knn_filename);

    knn_file<<num_nodes;

    for(int i=0; i<num_nodes; i++){
        for(int j =0; j<knn_set[i].size(); j++){
            pair<int, double> temp_pair = knn_set[i][j];
            knn_file<<endl << i <<" " << temp_pair.first <<" " << temp_pair.second;
        }
    }

    knn_file.close();
    cout<<"construction finished " << knn_filename<<endl;

}

inline bool file_exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}



#endif //TOPPPR_ALGO_H
