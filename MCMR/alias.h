#ifndef ALIAS_H
#define ALIAS_H

#include <random>
#include <algorithm>
#include <stack>
#include "Random.h"
using namespace std;

class Alias{
public:
	double* p;
	int* h;
	int* map1;
	int* map2;
	int n;
	Alias(){
		cout << " empty alias" <<endl;
	}
	Alias(vector<pair<int, double> > pi){
		double sum = 0;
		n = pi.size();
		cout << "alias size" << n <<endl;
		stack<int> small;
		stack<int> big;
		p = new double[n];
		h = new int[n];
		map1 = new int[n];
		for(int i = 0; i < n; i++){
			sum += pi[i].second;
			map1[i] = pi[i].first;
		}
		//cout << "sum of weights: " << sum << endl;
		for(int i = 0; i < n; i++){
			p[i] = pi[i].second * n / sum;
			if(p[i] > 1)
				big.push(i);
			else
				small.push(i);
		}
		while(!small.empty() && !big.empty()){
			int smallId = small.top();
			small.pop();
			int bigId = big.top();
			h[smallId] = bigId;
			p[bigId] -= (1-p[smallId]);
			if(p[bigId] < 1 ){
				small.push(bigId);
				big.pop();
			}
		}
		if(big.size()>0){
			while(!big.empty()){
				int bigId = big.top();
				h[bigId] = bigId;
				big.pop();
			}
		}
		if(small.size()>0){
			while(!small.empty()){
				int smallId = small.top();
				h[smallId] = smallId;
				small.pop();
			}
		}
		/*
		for(int i=0; i < n; i++){
			cout<< "node " <<i<<" : " << p[i] <<endl;
			cout<< "mapped node from " << i << ": " << map1[i] << endl;
			cout<< "after hashing from " << i << ": " << h[i]<< endl;
			cout<< "mapped node after hashing from "<< i << ": " << map1[h[i]] << endl;
		}*/

	}

	~Alias(){
		delete[] p;
		delete[] h;
		delete[] map1;
	}
	int generateRandom(Random& R){
		//cout<< "Generating Random Node" <<endl;
		int firstId = R.drand() * n;
		//cout << "first ID: " << firstId << " p[firstID]" << p[firstId]<< endl;
		int secondId = R.drand() < p[firstId] ? map1[firstId] : map1[h[firstId]];

		return secondId;
	}
	int generateRandom_t(Random& R){
		int firstId = R.drand_t() * n;
		int secondId = R.drand_t() < p[firstId] ? map1[firstId] : map1[h[firstId]];
		return secondId;
	}
};

#endif