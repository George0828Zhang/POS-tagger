#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
#include <limits>
#include <cmath>
#include "Array.h"

// #define TAGNUM 60
#define MAXTAGNUM 200 
#define MAXWDLEN 200
#define MAXTAGLEN 10

using namespace std;


unordered_map<string, int> dict_map;
unordered_map<string, int> tag_map;
vector<int> EqvClass;
char tagname[MAXTAGNUM][MAXTAGLEN]={};
double initial[MAXTAGNUM];
double transition[MAXTAGNUM][MAXTAGNUM];
double trigram[MAXTAGNUM][MAXTAGNUM][MAXTAGNUM];
double* emission[MAXTAGNUM];

int TAGNUM;
const double inf = std::numeric_limits<double>::infinity();

void load_model(char* name);
int wordIndex(string const& word);
void POStag3(vector<string> const& sentence, vector<int>& tag);

int main(int argc, char** argv){
	assert(argc==2);
	load_model(argv[1]);
	string sentence;
	while(getline(cin, sentence)){
		if(sentence.size()<2) continue;
		vector<string> tokenized;
		for(int i = 0, ls = 0; i < sentence.size(); i++){
			if (i==ls) continue;
			if(sentence[i]==' ' || sentence[i]==0){
				tokenized.push_back(sentence.substr(ls, i-ls));
				ls = i+1;
			}else if (sentence[i]=='.' || sentence[i]==',' || sentence[i]=='!' || sentence[i]=='?'){
				tokenized.push_back(sentence.substr(ls, i-ls));
				tokenized.push_back(string(1, sentence[i]));
				ls = i+1;
			}
		}
		vector<int> tags(tokenized.size());
		POStag3(tokenized, tags);
		for(int i = 0; i < tags.size(); i++){
			cout << tokenized[i] << "(" << tagname[tags[i]] << ") ";
		}
		cout << endl;
	}
}


int wordIndex(string const& word){
	auto iter = dict_map.find(word);
	return iter==dict_map.end() ? -1 : EqvClass[iter->second];
}

Array<double> delta;
Array<int> phi;
void POStag3(vector<string> const& sentence, vector<int>& tag){
	int slen = sentence.size();
	delta.reshape({TAGNUM, TAGNUM, slen});
	phi.reshape({TAGNUM, TAGNUM, slen});

	// t = 0
	int index = wordIndex(sentence[0]);
	for(int i = 0; i < TAGNUM; i++){
		double emi_i = index==-1 ? (1.0/TAGNUM) : emission[i][index];		
		for(int j = 0; j < TAGNUM; j++)
			delta[{i,j,0}] = emi_i==0 ? -inf : ( log(initial[i]) + log(emi_i) );
	}

	// t = 1
	if(slen > 1){
		index = wordIndex(sentence[1]);
		for(int i = 0; i < TAGNUM; i++){
			for(int j = 0; j < TAGNUM; j++){
				double emi_i = index==-1 ? (1.0/TAGNUM) : emission[i][index];
				if(emi_i==0 || transition[j][i]==0)
					delta[{i,j,1}] = -inf;
				else
					delta[{i,j,1}] = delta[{j,0,0}] + log(transition[j][i]) + log(emi_i);
			}
		}

		// t > 0
		for(int t = 2; t < slen;t++){
			index = wordIndex(sentence[t]);

			for(int i = 0; i < TAGNUM; i++){
				for(int j = 0; j < TAGNUM; j++){
					double emi_i = index==-1 ? (1.0/TAGNUM) : emission[i][index];
					phi[{i,j,t}] = -1;
					if(emi_i==0){
						delta[{i,j,t}] = -inf;
					}
					else{
						double logprob = -inf;
						for(int k = 0; k < TAGNUM; k++){
							if(trigram[k][j][i]<=0){ cerr << "bad trigram." <<endl; exit(1);};
							double tmp = delta[{j,k,t-1}] + log(trigram[k][j][i]);
							if(tmp > logprob){
								logprob = tmp;
								phi[{i,j,t}] = k;
							}
						}
						delta[{i,j,t}] = logprob + log(emi_i);
					}
				}
			}		
		}
	}

	// find best tail
	int best_tail_i = -1, best_tail_j = -1;
	double best_prob = -inf;
	for(int i = 0; i < TAGNUM; i++){			
		for(int j = 0; j < TAGNUM; j++){
			// cout << "p=" << delta[i][j][sentence.size()-1] << " i=" << i << " j=" << j << endl;
			if(delta[{i,j,slen-1}]>best_prob){
				best_tail_i = i;
				best_tail_j = j;
				best_prob = delta[{i,j,slen-1}];
			}
		}
	}

	// back tracking	
	for(int t = slen-1; t >= 1; t--){
		tag[t] = best_tail_i;
		tag[t-1] = best_tail_j;
		int tmp = phi[{best_tail_i,best_tail_j,t}];
		best_tail_i = best_tail_j;
		best_tail_j = tmp;
	}
}


void load_model(char* name){
	char buffer[MAXWDLEN];
	double p;
	int d;
	FILE* fp = fopen(name, "rt");

	fscanf(fp, "%s", buffer);
	assert(strncmp(buffer, "#initial", 8)==0);
	fscanf(fp, "%d", &TAGNUM);
	// assert(d==TAGNUM);
	for(int i = 0; i < TAGNUM; i++){
		memset(buffer, 0, MAXWDLEN);
		fscanf(fp, "%s", buffer);
		strncpy(tagname[i], buffer, MAXTAGLEN);
		string tmp = buffer;
		tag_map[tmp] = i;
	}
	for(int i = 0; i < TAGNUM; i++){		
		fscanf(fp, "%lf", &initial[i]);
	}
	// fprintf(stderr, "[DEBUG] inital done.\n");

	fscanf(fp, "%s", buffer);
	assert(strncmp(buffer, "#transition", 11)==0);

	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			fscanf(fp, "%lf", &transition[i][j]);
		}
	}
	// fprintf(stderr, "[DEBUG] transition done.\n");

	fscanf(fp, "%s", buffer);
	assert(strncmp(buffer, "#trigram", 8)==0);
	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			for(int k = 0; k < TAGNUM; k++){		
				fscanf(fp, "%lf", &trigram[i][j][k]);
			}
		}
	}

	fscanf(fp, "%s", buffer);
	assert(strncmp(buffer, "#vocab", 6)==0);
	fscanf(fp, "%s", buffer);
	while(strncmp(buffer, "#emission", 9)!=0){
		fscanf(fp, "%d", &d);
		string vocab(buffer);
		// dictionary.push_back(vocab);
		EqvClass.push_back(d);
		dict_map[vocab] = EqvClass.size()-1;
		fscanf(fp, "%s", buffer);
	}
	cout << "[Info] Dictionary Size:\t" << EqvClass.size() << endl;

	int EQCLASS;
	fscanf(fp, "%d", &EQCLASS);
	for(int i = 0; i < TAGNUM; i++){
		emission[i] = new double[EQCLASS]();		
		for(int j = 0; j < EQCLASS; j++){
			fscanf(fp, "%lf", &emission[i][j]);
		}
	}
}


void POStag(vector<string> const& sentence, vector<int>& tag){
	double* delta[TAGNUM];
	int* phi[TAGNUM];
	for(int i = 0; i < TAGNUM; i++){
		delta[i] = new double[sentence.size()];
		phi[i] = new int[sentence.size()];
	}

	// t = 0
	int index = wordIndex(sentence[0]);	
	for(int i = 0; i < TAGNUM; i++){
		if(index==-1){
			// assume uniform distribution
			delta[i][0] = log(initial[i]) + log(1.0/TAGNUM);
		}
		else{
			delta[i][0] = emission[i][index]==0 ? -inf : ( log(initial[i]) + log(emission[i][index]) );
		}
	}

	// t > 0
	for(int t = 1; t < sentence.size();t++){
		index = wordIndex(sentence[t]);
		if(index==-1){
			cerr << sentence[t] << " not in dictionary." << endl;
		}
		for(int i = 0; i < TAGNUM; i++){			
			phi[i][t] = -1;
			double logprob = -inf;
			if(index==-1 || emission[i][index]){					
				for(int j = 0; j < TAGNUM; j++){
					if(transition[j][i]==0) continue;
					double tmp = delta[j][t-1] + log(transition[j][i]);
					if(tmp > logprob){
						logprob = tmp;
						phi[i][t] = j;
					}
				}
				logprob += index==-1 ? log(1.0/TAGNUM) : log(emission[i][index]);
			}		
			delta[i][t] = logprob;
		}
	}

	// find best tail
	int best_tail = -1;
	double best_prob = -inf;
	for(int i = 0; i < TAGNUM; i++){			
		if(delta[i][sentence.size()-1]>best_prob){
			best_tail = i;
			best_prob = delta[i][sentence.size()-1];
		}
	}

	// back tracking	
	for(int i = sentence.size()-1; i >= 0; i--){
		tag[i] = best_tail;
		best_tail = phi[best_tail][i];
	}


	for(int i = 0; i < TAGNUM; i++){
		delete[] delta[i];
		delete[] phi[i];
	}
}
