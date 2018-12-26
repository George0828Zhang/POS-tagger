#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include "Array.h"

// #define TAGNUM 60
// #define MAXTAGNUM 200 
// #define MAXWDLEN 200
// #define MAXTAGLEN 10
// #define RARETHRES 10
// #define MAXSUFFIX 2
constexpr int MAXTAGNUM = 62; 
constexpr int MAXWDLEN = 200;
constexpr int MAXTAGLEN = 10;
constexpr int RARETHRES = 10;
constexpr int MAXSUFFIX = 5;
constexpr int BEAMSIZE = 100;


using namespace std;


unordered_map<string, int> dict_map;
// unordered_map<string, int> tag_map;
unordered_map<string, int> suffix_map;
Array<double> suf_emission;
// vector<int> EqvClass;
char tagname[MAXTAGNUM][MAXTAGLEN]={};
double initial[MAXTAGNUM];
double transition[MAXTAGNUM][MAXTAGNUM];
double trigram[MAXTAGNUM][MAXTAGNUM][MAXTAGNUM];
double* emission[MAXTAGNUM];

int TAGNUM;
int EQCLASS;
const double inf = std::numeric_limits<double>::infinity();

void load_model(char* name);
int wordIndex(string const& word);
void POStag3(vector<string> const& sentence, vector<int>& tag);
void tokenizer(string const& sentence, vector<string>& tokenized);

int main(int argc, char** argv){
	assert(argc==2);
	load_model(argv[1]);
	string sentence;
	while(getline(cin, sentence)){
		// if(sentence.size()<2) continue;
		vector<string> tokenized;
		tokenizer(sentence, tokenized);
		// Capitalization
		string Cap(tokenized[0]);
		if(tokenized[0][0]>='A' && tokenized[0][0]<='Z') 
			tokenized[0][0] = tolower(tokenized[0][0]);
		// POS tagging
		vector<int> tags(tokenized.size());
		POStag3(tokenized, tags);
		tokenized[0] = Cap;
		for(int i = 0; i < tags.size(); i++){
			// cout << tokenized[i] << "(" << tagname[tags[i]] << ") ";
			cout << tokenized[i] << "/" << tagname[tags[i]] << " ";
		}
		cout << endl;
	}
}


int wordIndex(string const& word){
	auto iter = dict_map.find(word);
	return iter==dict_map.end() ? -1 : iter->second;
}
int suffIndex(string const& word){
	auto iter = suffix_map.find(word);
	return iter==suffix_map.end() ? -1 : iter->second;
}

double oovEmission(string const& word, int tg){
	int wlen = word.size();
	for(int l = min(MAXSUFFIX, wlen); l > 0; l--){
		int index = suffIndex(word.substr(wlen-l, l));
		if(index!=-1){
			return suf_emission[{tg, index}];
		}
	}
	return (1.0/TAGNUM);
}

Array<double> delta;
Array<int> phi;
Kmax<double, int> beamspace(BEAMSIZE);
vector<int> beam;
void POStag3(vector<string> const& sentence, vector<int>& tag){
	int slen = sentence.size();
	delta.reshape({TAGNUM, TAGNUM, slen});
	phi.reshape({TAGNUM, TAGNUM, slen});

	delta.clear(-inf);
	phi.clear(-1);

	// t = 0
	int index = wordIndex(sentence[0]);
	for(int i = 0; i < TAGNUM; i++){
		double emi_i = index==-1 ? oovEmission(sentence[0], i) : emission[i][index];	
		if(emi_i>0){
			for(int j = 0; j < TAGNUM; j++){
				delta[{i,j,0}] = ( log(initial[i]) + log(emi_i) );
			}
		}
	}

	// t = 1
	if(slen > 1){
		beamspace.clear();
		index = wordIndex(sentence[1]);
		for(int i = 0; i < TAGNUM; i++){
			double emi_i = index==-1 ? oovEmission(sentence[1], i) : emission[i][index];
			for(int j = 0; j < TAGNUM; j++){				
				if(emi_i>0 && transition[j][i]>0){
					delta[{i,j,1}] = delta[{j,0,0}] + log(transition[j][i]) + log(emi_i);	
				}
				beamspace.insert(delta[{i,j,1}], i*TAGNUM+j);
			}
		}

		// t > 1
		for(int t = 2; t < slen;t++){
			beamspace.extract(beam);
			beamspace.clear();

			index = wordIndex(sentence[t]);

			for(int i = 0; i < TAGNUM; i++){
				double emi_i = index==-1 ? oovEmission(sentence[t], i) : emission[i][index];
				for(auto& pair : beam){
					int j = pair/TAGNUM, k = pair%TAGNUM;					

					if(emi_i > 0 && trigram[k][j][i] > 0){
						double logprob = delta[{j,k,t-1}] + log(trigram[k][j][i]) + log(emi_i);
						if(logprob > delta[{i,j,t}]){
							delta[{i,j,t}] = logprob;
							phi[{i,j,t}] = k;							
						}
					}
					beamspace.insert(delta[{i,j,t}], i*TAGNUM+j);
				}
			}
		}
	}

	// find best tail
	int best_tail_i = -1, best_tail_j = -1;
	double best_prob = -inf;
	// for(int i = 0; i < TAGNUM; i++){			
	// 	for(int j = 0; j < TAGNUM; j++){
	// 		if(delta[{i,j,slen-1}]>best_prob){
	// 			best_tail_i = i;
	// 			best_tail_j = j;
	// 			best_prob = delta[{i,j,slen-1}];
	// 		}
	// 	}
	// }
	beamspace.extract(beam);
	for(auto& pair : beam){
		int i = pair/TAGNUM, j = pair%TAGNUM;					
		if(delta[{i,j,slen-1}]>best_prob){
			best_tail_i = i;
			best_tail_j = j;
			best_prob = delta[{i,j,slen-1}];
		}
	}


	// back tracking	
	for(int t = slen-1; t >= 1; t--){
		tag[t] = best_tail_i;
		tag[t-1] = best_tail_j;
		if (t==1) break;
		int tmp = phi[{best_tail_i,best_tail_j,t}];
		best_tail_i = best_tail_j;
		best_tail_j = tmp;
	}
}

void tokenizer2(string const& sentence, vector<string>& tokenized){
	const string common_ch = "~!@#$%%^&*()_+=`\\/.,<>?\":;";
	// "-" might be hyphenized text, and ' and . might be abreviation
	int slen = sentence.size();
	int tlen = 0;
	int c_at;
	for(int i = 0; i < slen; i++){
		char c = sentence[i];
		if(common_ch.find(c) != string::npos){
			if(tlen){
				// process sub-word tokenizing TODO
				string word = sentence.substr(i-tlen, tlen);
				c_at = word.find('\'');
				if(c_at != string::npos && c_at != 0){
					tokenized.push_back(word.substr(0, c_at));
					word.erase(0, c_at);
				}				
				tokenized.push_back(word);
			}
			tokenized.push_back(string(1, c));
			tlen = 0;
		}else if(tlen && (c == ' ' || c == '\n' || c == '\t' || c == '\r')){
			// process sub-word tokenizing TODO
			string word = sentence.substr(i-tlen, tlen);
			c_at = word.find('\'');
			if(c_at != string::npos && c_at != 0){
				tokenized.push_back(word.substr(0, c_at));
				word.erase(0, c_at);
			}				
			tokenized.push_back(word);
			tlen = 0;
		}else if(i == slen - 1){
			// process sub-word tokenizing TODO
			string word = sentence.substr(i-tlen, 1+tlen);
			c_at = word.find('\'');
			if(c_at != string::npos && c_at != 0){
				tokenized.push_back(word.substr(0, c_at));
				word.erase(0, c_at);
			}				
			tokenized.push_back(word);
			tlen = 0;
		}
		else
			tlen++;
	}
}

void tokenizer(string const& sentence, vector<string>& tokenized){

	int slen = sentence.size();
	int tlen = 0;
	int c_at;
	for(int i = 0; i < slen; i++){
		char c = sentence[i];
		if(tlen && (c == ' ' || c == '\n' || c == '\t' || c == '\r')){
			string word = sentence.substr(i-tlen, tlen);							
			tokenized.push_back(word);
			tlen = 0;
		}else if(i == slen - 1){
			string word = sentence.substr(i-tlen, 1+tlen);							
			tokenized.push_back(word);
			tlen = 0;
		}
		else
			tlen++;
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
	cerr << "[Info] Tagset Size: " << TAGNUM << endl;

	for(int i = 0; i < TAGNUM; i++){
		memset(buffer, 0, MAXWDLEN);
		fscanf(fp, "%s", buffer);
		strncpy(tagname[i], buffer, MAXTAGLEN);
		string tmp = buffer;
		// tag_map[tmp] = i;
	}
	for(int i = 0; i < TAGNUM; i++){		
		fscanf(fp, "%lf", &initial[i]);
	}

	fscanf(fp, "%s", buffer);
	assert(strncmp(buffer, "#transition", 11)==0);

	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			fscanf(fp, "%lf", &transition[i][j]);
		}
	}

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
	assert(strncmp(buffer, "#emission", 9)==0);
	
	// int EQCLASS;
	fscanf(fp, "%d", &EQCLASS);
	for(int i = 0; i < TAGNUM; i++){
		emission[i] = new double[EQCLASS]();		
		for(int j = 0; j < EQCLASS; j++){
			fscanf(fp, "%lf", &emission[i][j]);
		}
	}
	cerr << "[Info] Equivalent Classes: " << EQCLASS << endl;

	fscanf(fp, "%s", buffer);
	assert(strncmp(buffer, "#vocab_freq_Eqv", 15)==0);
	int freq;

	vector<vector<int> > f_S_Eqv;
	int sum_of_Suffix[EQCLASS];
	memset(sum_of_Suffix, 0, EQCLASS*sizeof(int));
	while(fscanf(fp, "%s %d %d", buffer, &freq, &d)!=EOF){
		string vocab(buffer);
		dict_map[vocab] = d;

		if(freq < RARETHRES){
			for(int i = 0; i < MAXSUFFIX && i < vocab.size(); i++){
				string suf = vocab.substr(vocab.size()-1-i,i+1);
				int suf_index = suffIndex(suf);
				if(suf_index==-1){					
					f_S_Eqv.push_back(vector<int>(EQCLASS, 0));
					suf_index = f_S_Eqv.size() - 1;
					suffix_map[suf] = suf_index;
				}
				f_S_Eqv[suf_index][d] += freq;
				sum_of_Suffix[d] += freq;
			}
		}
	}
	cerr << "[Info] Dictionary Size: " << dict_map.size() << endl;
	cerr << "[Info] Suffix Map Size: " << suffix_map.size() << endl;

	// TODO: remove suffixes that appeared only once?
	// suffix smoothing	
	suf_emission.reshape({TAGNUM, (int)suffix_map.size()});

	for(int tg = 0; tg < TAGNUM; tg++){
		for(int s = 0; s < suffix_map.size(); s++){
			double prob = 0.;
			for(int q = 0; q < EQCLASS; q++){
				double P_S_Eqv = sum_of_Suffix[q] ? ((double)f_S_Eqv[s][q])/(double)sum_of_Suffix[q] : 0.;
				prob += emission[tg][q] * P_S_Eqv;
			}
			suf_emission[{tg, s}] = prob;
		}		
	}
}

void POStag3_old(vector<string> const& sentence, vector<int>& tag){
	int slen = sentence.size();
	delta.reshape({TAGNUM, TAGNUM, slen});
	phi.reshape({TAGNUM, TAGNUM, slen});

	// t = 0
	int index = wordIndex(sentence[0]);
	for(int i = 0; i < TAGNUM; i++){
		double emi_i = index==-1 ? oovEmission(sentence[0], i) : emission[i][index];	
		for(int j = 0; j < TAGNUM; j++)
			delta[{i,j,0}] = emi_i==0 ? -inf : ( log(initial[i]) + log(emi_i) );
	}

	// t = 1
	if(slen > 1){
		index = wordIndex(sentence[1]);
		for(int i = 0; i < TAGNUM; i++){
			for(int j = 0; j < TAGNUM; j++){
				double emi_i = index==-1 ? oovEmission(sentence[1], i) : emission[i][index];
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
					double emi_i = index==-1 ? oovEmission(sentence[t], i) : emission[i][index];
					
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
