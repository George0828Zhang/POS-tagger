#include <iostream>
#include <fstream>
#include <sstream>//parsing line
// #include <unordered_map>
#include <string>
#include <cassert>
#include <array>
#include <vector>
#include <limits>// for inf
#include <cmath>// for log
#include <iomanip>// for set precision

#include "Tree.hpp"
#include "Array.hpp"
#include "tag.hpp"

constexpr int MAXTAGNUM = 62;
constexpr int RARETHRES = 10;
constexpr int MAXSUFFIX = 10;
constexpr int BEAMSIZE = 100;
constexpr double inf = std::numeric_limits<double>::infinity();

using Pair = std::array<int, 2>;

// std::unordered_map<std::string, int> dict_map;
lexiTree dict_map;
// std::unordered_map<std::string, int> suffix_map;
suffixTree suffix_map;
double initial[MAXTAGNUM];
double transition[MAXTAGNUM][MAXTAGNUM];
double trigram[MAXTAGNUM][MAXTAGNUM][MAXTAGNUM];
double* emission[MAXTAGNUM];
DyArray<double> suf_emission;

int TAGNUM;
int EQCLASS;
int CORPUS_N;

int wordIndex(std::string const& word){
	// auto iter = dict_map.find(word);
	// return iter==dict_map.end() ? -1 : iter->second;
	return dict_map.wordIndex(word);
}
int suffIndex(std::string const& word){
	// auto iter = suffix_map.find(word);
	// return iter==suffix_map.end() ? -1 : iter->second;
	return suffix_map.suffixIndex(word);
}

// double oovEmission(std::string const& word, int tg){
// 	int wlen = word.size();
// 	for(int l = std::min(MAXSUFFIX, wlen); l > 0; l--){
// 		int index = suffIndex(word.substr(wlen-l, l));
// 		if(index!=-1){
// 			return suf_emission[{tg, index}];
// 		}
// 	}
// 	return (1.0/TAGNUM);
// }

double oovEmission(std::string const& word, int tg){
	int index = suffix_map.suffixIndex(word);
	if(index!=-1){
		return suf_emission[{tg, index}];
	}else
	return (1.0/TAGNUM);
}


DyArray<double> delta;
DyArray<int> phi;
Kmax<double, Pair> beamspace(BEAMSIZE);
std::vector<Pair> beam;
void POStag3(std::vector<std::string> const& sentence, std::vector<int>& tag){
	int slen = sentence.size();
	delta.reshape({slen, TAGNUM, TAGNUM});
	phi.reshape({slen, TAGNUM, TAGNUM});

	delta.clear(-inf);
	phi.clear(-1);

	// t = 0
	int index = wordIndex(sentence[0]);
	for(int i = 0; i < TAGNUM; i++){
		double emi_i = index==-1 ? oovEmission(sentence[0], i) : emission[i][index];
		if(emi_i>0){
			for(int j = 0; j < TAGNUM; j++){
				delta[{0, i,j}] = ( log(initial[i]) + log(emi_i) );
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
					delta[{1, i,j}] = delta[{0, j,0}] + log(transition[j][i]) + log(emi_i);	
				}
				beamspace.insert(delta[{1, i,j}], {i,j});
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
					int j = pair[0], k = pair[1];					

					if(emi_i > 0 && trigram[k][j][i] > 0){
						double logprob = delta[{t-1, j,k}] + log(trigram[k][j][i]) + log(emi_i);
						if(logprob > delta[{t, i,j}]){
							delta[{t, i,j}] = logprob;
							phi[{t, i,j}] = k;							
						}
					}
					beamspace.insert(delta[{t, i,j}], {i,j});
				}
			}
		}
	}

	// find best tail
	int best_tail_i = -1, best_tail_j = -1;
	double best_prob = -inf;
	beamspace.extract(beam);
	for(auto& pair : beam){
		int i = pair[0], j = pair[1];					
		if(delta[{slen-1,i,j}]>best_prob){
			best_tail_i = i;
			best_tail_j = j;
			best_prob = delta[{slen-1,i,j}];
		}
	}


	// back tracking	
	for(int t = slen-1; t >= 1; t--){
		tag[t] = best_tail_i;
		tag[t-1] = best_tail_j;
		if (t==1) break;
		int tmp = phi[{t, best_tail_i,best_tail_j}];
		best_tail_i = best_tail_j;
		best_tail_j = tmp;
	}
}
// void POStag3(std::vector<std::string> const& sentence, std::vector<int>& tag){
// 	int slen = sentence.size();
// 	delta.reshape({TAGNUM, TAGNUM, slen});
// 	phi.reshape({TAGNUM, TAGNUM, slen});

// 	delta.clear(-inf);
// 	phi.clear(-1);

// 	// t = 0
// 	int index = wordIndex(sentence[0]);
// 	for(int i = 0; i < TAGNUM; i++){
// 		double emi_i = index==-1 ? oovEmission(sentence[0], i) : emission[i][index];
// 		if(emi_i>0){
// 			for(int j = 0; j < TAGNUM; j++){
// 				delta[{i,j,0}] = ( log(initial[i]) + log(emi_i) );
// 			}
// 		}
// 	}

// 	// t = 1
// 	if(slen > 1){
// 		beamspace.clear();
// 		index = wordIndex(sentence[1]);
// 		for(int i = 0; i < TAGNUM; i++){
// 			double emi_i = index==-1 ? oovEmission(sentence[1], i) : emission[i][index];
// 			for(int j = 0; j < TAGNUM; j++){				
// 				if(emi_i>0 && transition[j][i]>0){
// 					delta[{i,j,1}] = delta[{j,0,0}] + log(transition[j][i]) + log(emi_i);	
// 				}
// 				beamspace.insert(delta[{i,j,1}], {i,j});
// 			}
// 		}

// 		// t > 1
// 		for(int t = 2; t < slen;t++){
// 			beamspace.extract(beam);
// 			beamspace.clear();

// 			index = wordIndex(sentence[t]);

// 			for(int i = 0; i < TAGNUM; i++){
// 				double emi_i = index==-1 ? oovEmission(sentence[t], i) : emission[i][index];
// 				for(auto& pair : beam){
// 					int j = pair[0], k = pair[1];					

// 					if(emi_i > 0 && trigram[k][j][i] > 0){
// 						double logprob = delta[{j,k,t-1}] + log(trigram[k][j][i]) + log(emi_i);
// 						if(logprob > delta[{i,j,t}]){
// 							delta[{i,j,t}] = logprob;
// 							phi[{i,j,t}] = k;							
// 						}
// 					}
// 					beamspace.insert(delta[{i,j,t}], {i,j});
// 				}
// 			}
// 		}
// 	}

// 	// find best tail
// 	int best_tail_i = -1, best_tail_j = -1;
// 	double best_prob = -inf;
// 	beamspace.extract(beam);
// 	for(auto& pair : beam){
// 		int i = pair[0], j = pair[1];					
// 		if(delta[{i,j,slen-1}]>best_prob){
// 			best_tail_i = i;
// 			best_tail_j = j;
// 			best_prob = delta[{i,j,slen-1}];
// 		}
// 	}


// 	// back tracking	
// 	for(int t = slen-1; t >= 1; t--){
// 		tag[t] = best_tail_i;
// 		tag[t-1] = best_tail_j;
// 		if (t==1) break;
// 		int tmp = phi[{best_tail_i,best_tail_j,t}];
// 		best_tail_i = best_tail_j;
// 		best_tail_j = tmp;
// 	}
// }

void POStag(std::vector<std::string> const& sentence, std::vector<int>& tag){
	int slen = sentence.size();
	delta.reshape({slen, TAGNUM});
	phi.reshape({slen, TAGNUM});

	delta.clear(-1);
	phi.clear(-1);

	// t = 0
	int index = dict_map.wordIndex(sentence[0]);
	for(int i = 0; i < TAGNUM; i++){
		double emi_i = index==-1 ? oovEmission(sentence[0], i) : emission[i][index];
		delta[{0,i}] = initial[i] * emi_i;
	}

	// t >= 1
	for(int t = 1; t < slen;t++){
		index = dict_map.wordIndex(sentence[t]);

		for(int i = 0; i < TAGNUM; i++){
			double emi_i = index==-1 ? oovEmission(sentence[t], i) : emission[i][index];
			for(int j = 0; j < TAGNUM; j++){					
				double prob = delta[{t-1,j}] * transition[j][i] * emi_i;
				if(prob > delta[{t,i}]){
					delta[{t,i}] = prob;
					phi[{t,i}] = j;							
				}
			}
		}
	}

	// find best tail
	int best_tail_i = -1;
	double best_prob = -1;
	for(int i = 0; i < TAGNUM; i++){
		if(delta[{slen-1,i}]>best_prob){
			best_tail_i = i;
			best_prob = delta[{slen-1,i}];
		}
	}


	// back tracking	
	for(int t = slen-1; t >= 0; t--){
		tag[t] = best_tail_i;
		if (t==0) break;
		best_tail_i = phi[{t,best_tail_i}];
	}
}





DyArray<float> alpha;
DyArray<float> beta;
DyArray<float> Gamma;
DyArray<float> epsilon;
DyArray<float> epsilon_3;

void initEmission(std::vector<std::string> const& sentence, DyArray<float>& ehat);
void makeAlpha(std::vector<std::string> const& sentence, DyArray<float> & ehat);
void makeBeta(std::vector<std::string> const& sentence, DyArray<float>& ehat);
void makeGamma(int T);
void makeEpsilons(int T, DyArray<float>& ehat);

void refineHMM(std::string const& name){
	std::vector<std::vector<std::string>> sentences;
	std::ifstream source(name, std::ifstream::in);
	std::string buffer;
	while(getline(source, buffer)){
		std::vector<std::string> tokenized;
		tokenizer(buffer, tokenized);
		// Capitalization
		std::string Cap(tokenized[0]);
		if(tokenized[0][0]>='A' && tokenized[0][0]<='Z') 
			tokenized[0][0] = tolower(tokenized[0][0]);
		sentences.push_back(tokenized);
	}


	int samples = sentences.size();
	std::cerr << "[info] trainning size: " << samples << std::endl;
	DyArray<float> emit_hat;
	DyArray<float> initial_cummu({TAGNUM});
	DyArray<float> transition_cummu({TAGNUM, TAGNUM});
	DyArray<float> trigram_cummu({TAGNUM, TAGNUM, TAGNUM});
	DyArray<float> visit_cummu({TAGNUM});
	DyArray<float> visit_pair_cummu({TAGNUM, TAGNUM});
	DyArray<float> observ_cummu({TAGNUM, EQCLASS});

	initial_cummu.clear(0);
	transition_cummu.clear(0);
	trigram_cummu.clear(0);
	observ_cummu.clear(0);
	visit_cummu.clear(0);
	visit_pair_cummu.clear(0);

	// for(auto const& sent : sentences ){
	for(int s = 0; s < samples; s++ ){
		auto const& sent = sentences[s];

		int T = sent.size();
		initEmission(sent, emit_hat);
		makeAlpha(sent, emit_hat);
		makeBeta(sent, emit_hat);
		makeGamma(T);
		makeEpsilons(T, emit_hat);

		// for(int t = 0; t < T - 2; t++){
		for(int t = 0; t < T; t++){
			// cummulate for initial
			for(int i = 0; i < TAGNUM; i++){
				initial_cummu[{i}] += Gamma[{0, i}];
			}

			// cummulate for transition
			for(int i = 0; i < TAGNUM; i++){
				for(int j = 0; t<T-1 && j < TAGNUM; j++){
					transition_cummu[{i, j}] += epsilon[{t, i, j}];
				}
				visit_cummu[{i}] += Gamma[{t, i}];
			}

			// cummulate for trigram
			// for(int i = 0; i < TAGNUM; i++){
			// 	for(int j = 0; j < TAGNUM; j++){
			// 		for(int k = 0; k < TAGNUM; k++){
			// 			trigram_cummu[{i, j, k}] += epsilon_3[{t, i, j, k}];
			// 		}
			// 	}
			// }

			// cummulate for emission
			int index = dict_map.wordIndex(sent[t]);
			if(index != -1){
				for(int i = 0; i < TAGNUM; i++){
					observ_cummu[{i, index}] += Gamma[{t, i}];
				}
			}
		}

		std::string bar(30, '.');
		std::fill(bar.begin(), bar.begin() + (int)std::lround(30 * s / samples),'=');
		std::cerr << "(" << s+1 << "/" << samples << ") sentences done. Progress [" << bar << "]\r";
	}

	// reestimate model
	float sum = 0.;
	for(int i = 0; i < TAGNUM; i++){
		// initial
		initial[i] = initial_cummu[{i}] / samples;
		// sum += initial_cummu[{i}];
	}
	// for(int i = 0; i < TAGNUM; i++){
	// 	initial[i] = initial_cummu[{i}] / sum;
	// }

	for(int i = 0; i < TAGNUM; i++){
		// transition
		for(int j = 0; visit_cummu[{i}] > 0 && j < TAGNUM; j++){
			transition[i][j] = transition_cummu[{i, j}] / visit_cummu[{i}];
			// for(int k = 0; transition_cummu[{i, j}] > 0 && k < TAGNUM; k++){
			// 	trigram[i][j][k] = trigram_cummu[{i, j, k}] / transition_cummu[{i, j}];
			// }
		}

		for(int j = 0; visit_cummu[{i}] > 0 && j < EQCLASS; j++){
			emission[i][j] = observ_cummu[{i, j}] / visit_cummu[{i}];
		}
	}


}
void initEmission(
	std::vector<std::string> const& sentence, 
	DyArray<float>& ehat)
{
	int T = sentence.size();
	ehat.reshape({T, TAGNUM});
	for(int t = 0; t < T; t++){
		int index = dict_map.wordIndex(sentence[t]);
		for(int i = 0; i < TAGNUM; i++){
			ehat[{t, i}] = index==-1 ? (float)oovEmission(sentence[t], i) : (float)emission[i][index];
		}
	}
}
// void makeAlpha(
// 	std::vector<std::string> const& sentence, 
// 	DyArray<float>& ehat)
// {
// 	int T = sentence.size();

// 	alpha.reshape({T, TAGNUM});
// 	alpha.clear(0);

// 	// t = 0
// 	for(int i = 0; i < TAGNUM; i++){
// 		float emi_i = ehat[{0, i}];
// 		if(emi_i>0){
// 			alpha[{0, i}] = initial[i] * emi_i;
// 		}
// 	}

// 	if(T > 1){
// 		// t = 1
// 		// index = indices[1];
// 		for(int i = 0; i < TAGNUM; i++){
// 			float emi_i = ehat[{1, i}];
// 			for(int j = 0; j < TAGNUM; j++){				
// 				if(emi_i>0 && transition[j][i]>0){
// 					alpha[{1, i}] += alpha[{0, j}] * transition[j][i] * emi_i;	
// 				}
// 			}
// 		}

// 		// t > 1
// 		for(int t = 2; t < T;t++){
// 			// index = indices[t];
// 			for(int i = 0; i < TAGNUM; i++){
// 				float emi_i = ehat[{t, i}];
// 				if(emi_i>0){
// 					// int index2 = indices[t-1];
// 					for(int j = 0; j < TAGNUM; j++){
// 						float emi_j = ehat[{t-1, j}];
// 						if(emi_j>0){
// 							float prior = 0.;
// 							for(int k = 0; k < TAGNUM; k++){
// 								if(trigram[k][j][i] > 0){
// 									prior += alpha[{t-2, j, k}] * trigram[k][j][i]; //alpha buggggggggggggggg
// 								}
// 							}
// 							alpha[{t, i}] = prior * emi_i * emi_j;
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}	
// }
// void makeBeta(
// 	std::vector<std::string> const& sentence, 
// 	DyArray<float>& ehat)
// {
// 	int T = sentence.size();

// 	beta.reshape({T, TAGNUM});
// 	beta.clear(0);

// 	// t = T - 1
// 	for(int i = 0; i < TAGNUM; i++){
// 		beta[{T-1, i}] = 1.;
// 	}

// 	if(T > 1){
// 		// t = T - 2
// 		for(int j = 0; j < TAGNUM; j++){			
// 			for(int i = 0; i < TAGNUM; i++){	
// 				float emi_i = ehat[{T-1, i}];			
// 				if(emi_i>0 && transition[j][i]>0){
// 					beta[{T-2, j}] += beta[{T-1, i}] * transition[j][i] * emi_i;	
// 				}
// 			}
// 		}

// 		// t > 1
// 		for(int t = T-3; t >= 0;t--){
// 			for(int k = 0; k < TAGNUM; k++){
// 				for(int j = 0; j < TAGNUM; j++){	
// 					float emi_j = ehat[{t+1, j}];	
// 					if(emi_j>0){	
// 						for(int i = 0; i < TAGNUM; i++){	
// 							float emi_i = ehat[{t+2, i}];			
// 							if(emi_i>0){
// 								beta[{t, k}] += beta[{t+2, i}] * trigram[k][j][i] * emi_i * emi_j;	
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}	
// }
void makeGamma(int T){
	Gamma.reshape({T, TAGNUM});
	
	for(int t = 0; t < T; t++){
		float sum = 0.;
		for(int i = 0; i < TAGNUM; i++){
			Gamma[{t, i}] = alpha[{t, i}] * beta[{t, i}];
			sum += Gamma[{t, i}];
		}
		if(sum > 0){
			for(int i = 0; i < TAGNUM; i++){
				Gamma[{t, i}] /= sum;
			}
		}
	}
}
void makeEpsilons(int T, DyArray<float>& ehat){
	epsilon.reshape({T, TAGNUM, TAGNUM});
	
	for(int t = 0; t < T - 1; t++){
		float sum = 0;
		for(int j = 0; j < TAGNUM; j++){
			for(int i = 0; i < TAGNUM; i++){
				epsilon[{t, j, i}] = alpha[{t, j}] * transition[j][i] * ehat[{t+1, i}] * beta[{t+1, i}];
				sum += epsilon[{t, j, i}];
			}
		}
		if(sum > 0){
			for(int j = 0; j < TAGNUM; j++){
				for(int i = 0; i < TAGNUM; i++){
					epsilon[{t, j, i}] /= sum;
				}
			}
		}
	}

	// epsilon_3.reshape({T, TAGNUM, TAGNUM, TAGNUM});
	
	// for(int t = 0; t < T - 2; t++){
	// 	float sum = 0;
	// 	for(int k = 0; k < TAGNUM; k++){
	// 		for(int j = 0; j < TAGNUM; j++){
	// 			for(int i = 0; i < TAGNUM; i++){
	// 				epsilon_3[{t, k, j, i}] = alpha[{t, k}] * trigram[k][j][i] * ehat[{t+1, j}] * ehat[{t+2, i}] * beta[{t+2, i}];
	// 				sum += epsilon_3[{t, k, j, i}];
	// 			}
	// 		}
	// 	}
	// 	if(sum > 0){
	// 		for(int k = 0; k < TAGNUM; k++){
	// 			for(int j = 0; j < TAGNUM; j++){
	// 				for(int i = 0; i < TAGNUM; i++){
	// 					epsilon_3[{t, k, j, i}] /= sum;
	// 				}
	// 			}
	// 		}
	// 	}
	// }
}
void makeAlpha(
	std::vector<std::string> const& sentence, 
	DyArray<float>& ehat)
{
	int T = sentence.size();

	alpha.reshape({T, TAGNUM});
	alpha.clear(0);

	// t = 0
	for(int i = 0; i < TAGNUM; i++){
		alpha[{0, i}] = initial[i] * ehat[{0, i}];
	}

	// t > 1
	for(int t = 1; t < T;t++){
		for(int i = 0; i < TAGNUM; i++){
			float emi_i = ehat[{t, i}];
			// alpha[{t, i}] = 0;
			for(int j = 0; j < TAGNUM; j++){
				alpha[{t, i}] += alpha[{t-1, j}] * transition[j][i];
			}
			alpha[{t, i}] *= emi_i;
		}
	}
}
void makeBeta(
	std::vector<std::string> const& sentence, 
	DyArray<float>& ehat)
{
	int T = sentence.size();

	beta.reshape({T, TAGNUM});
	beta.clear(0);

	// t = T - 1
	for(int i = 0; i < TAGNUM; i++){
		beta[{T-1, i}] = 1.;
	}

	// t > 1
	for(int t = T-2; t >= 0;t--){
		for(int j = 0; j < TAGNUM; j++){
			for(int i = 0; i < TAGNUM; i++){	
				float emi_i = ehat[{t+1, i}];			
				beta[{t, j}] += beta[{t+1, i}] * transition[j][i] * emi_i;	
			}
		}
	}
}

void smooth_transition(){
// TODO: implement Kneser–Ney smoothing
	double lbd[] = {0, 0};
	for(int j = 0; j < TAGNUM; j++){
		for(int i = 0; i < TAGNUM; i++){
			if(transition[j][i]==0) continue;
			double visit_j = initial[j] * CORPUS_N;
			double visit_i = initial[i] * CORPUS_N;
			double fq = transition[j][i] * visit_j;

			double p1 = visit_j > 1 && fq >= 1 ? (fq-1)/(visit_j-1) : 0;
			double p2 = visit_i >= 1 ? (visit_i-1)/(CORPUS_N-1) : 0;
			if(p1 >= p2){
				lbd[1] += fq;
			}else{
				lbd[0] += fq;
			}
		}
	}

	assert(lbd[0]>=0 && lbd[1] >= 0 && lbd[0] + lbd[1] > 0);

	for(int j = 0; j < TAGNUM; j++){
		for(int i = 0; i < TAGNUM; i++){
			transition[j][i] = (lbd[0]*initial[i] + lbd[1]*transition[j][i]) / (lbd[0] + lbd[1]);
		}
	}
}
void smooth_trigram(){
// TODO: implement Kneser–Ney smoothing
	double lbd[] = {0, 0, 0};
	for(int k = 0; k < TAGNUM; k++){
		for(int j = 0; j < TAGNUM; j++){
			for(int i = 0; i < TAGNUM; i++){
				if(trigram[k][j][i]==0) continue;
				double visit_j = initial[j] * CORPUS_N;
				double visit_i = initial[i] * CORPUS_N;
				double trans_kj = transition[k][j] * (CORPUS_N * initial[k]);
				double trans_ji = transition[j][i] * visit_j;
				double fq = trigram[k][j][i] * trans_kj;

				double p1 = trans_kj > 1 && fq >= 1 ? (fq-1)/(trans_kj-1) : 0;
				double p2 = visit_j >= 1 && trans_ji>=1 ? (trans_ji - 1)/(visit_j-1) : 0; 
				double p3 = visit_i >= 1 ? (visit_i-1)/(CORPUS_N-1) : 0;
				if(p1 >= p2 && p1 >= p3){
					lbd[2] += fq;
				}else if(p2 >= p1 && p2 >= p3){
					lbd[1] += fq;
				}else{
					lbd[0] += fq;
				}
			}
		}
	}

	assert(lbd[0]>=0 && lbd[1] >= 0 && lbd[2] >= 0 && lbd[0] + lbd[1] + lbd[2] > 0);

	for(int k = 0; k < TAGNUM; k++){
		for(int j = 0; j < TAGNUM; j++){
			for(int i = 0; i < TAGNUM; i++){
				trigram[k][j][i] = (lbd[0]*initial[i] + lbd[1]*transition[j][i] + lbd[2]*trigram[k][j][i]) / (lbd[0] + lbd[1] + lbd[2]);
			}
		}
	}
}
void smooth_emission(){
// TODO: implement Kneser–Ney smoothing
	std::vector<double> P_observ(EQCLASS, 0);	
	for(int o = 0; o < EQCLASS; o++){
		for(int i = 0; i < TAGNUM; i++){
			P_observ[o] += emission[i][o] * initial[i];
		}
	}

	double lbd[] = {0, 0};
	for(int i = 0; i < TAGNUM; i++){
		for(int o = 0; o < EQCLASS; o++){			
			if(emission[i][o]==0) continue;			
			double N_obs = P_observ[o] * CORPUS_N;
			double N_tag = initial[i] * CORPUS_N;
			double fq = emission[i][o] * N_tag;

			double p1 = N_tag > 1 && fq >= 1 ? (fq-1)/(N_tag-1) : 0;
			double p2 = N_obs >= 1 ? (N_obs-1)/(CORPUS_N-1) : 0;

			if(p1 >= p2){
				lbd[1] += fq;
			}else{
				lbd[0] += fq;
			}
		}
	}

	assert(lbd[0]>=0 && lbd[1] >= 0 && lbd[0] + lbd[1] > 0);

	for(int i = 0; i < TAGNUM; i++){
		for(int o = 0; o < EQCLASS; o++){
			emission[i][o] = (lbd[0]*P_observ[o] + lbd[1]*emission[i][o]) / (lbd[0] + lbd[1]);
		}
	}
}



void load_lexicon(std::string const& name){
	std::ifstream source(name, std::ifstream::in);
	std::string buffer;
	
	int freq;
	CORPUS_N = 0;

	std::vector<std::vector<int> > f_S_Eqv;
	std::vector<int> sum_of_Suffix(EQCLASS, 0);
	int d, wcount=0;
	int max_suf = -1;
	// while(!source.eof()){
	// while(source >> buffer >> freq >> d){
		// source >> buffer >> freq >> d;
	std::string sentence;
	while(std::getline(source, sentence)){
		if(!(std::stringstream(sentence) >> buffer >> freq >> d)) continue;
		std::string vocab(buffer);
		// dict_map[vocab] = d;
		dict_map.insert(vocab, d);
		assert(wordIndex(vocab) == d);
		wcount++;
		CORPUS_N += freq;

		if(freq < RARETHRES){
			for(int i = 0; i < MAXSUFFIX && i < vocab.size(); i++){
				std::string suf = vocab.substr(vocab.size()-1-i,i+1);
				// int suf_index = suffIndex(suf);
				int suf_index = suffix_map.wordIndex(suf);
				if(suf_index==-1){					
					f_S_Eqv.push_back(std::vector<int>(EQCLASS, 0));
					suf_index = f_S_Eqv.size() - 1;
					// suffix_map[suf] = suf_index;
					suffix_map.insert(suf, suf_index);
					if(suf_index > max_suf)
						max_suf = suf_index;
				}
				f_S_Eqv[suf_index][d] += freq;
				sum_of_Suffix[d] += freq;
			}
		}
		// std::cerr << "(" << wcount << ") words loaded.\r";
		// std::cerr << "word: " << vocab << ", word freq " << freq << ", Equivalent Class " << d << "\r";
	}
	std::cerr << "[Info] Dictionary Size: " << wcount << std::endl;
	std::cerr << "[Info] Suffix Map Size: " << max_suf + 1 << std::endl;

	// smooth_transition();
	// smooth_emission();

	// TODO: remove suffixes that appeared only once?
	// suffix smoothing	
	suf_emission.reshape({TAGNUM, max_suf + 1});

	for(int tg = 0; tg < TAGNUM; tg++){
		for(int s = 0; s < max_suf + 1; s++){
			double prob = 0.;
			for(int q = 0; q < EQCLASS; q++){
				double P_S_Eqv = sum_of_Suffix[q] ? ((double)f_S_Eqv[s][q])/(double)sum_of_Suffix[q] : 0.;
				prob += emission[tg][q] * P_S_Eqv;
			}
			suf_emission[{tg, s}] = prob;
		}		
	}
}

void load_model(std::string const& name, std::vector<std::string>& tagname){
	std::ifstream source(name, std::ifstream::in);
	std::string buffer;

	source >> buffer;
	assert(buffer == "#initial");
	source >> TAGNUM;
	std::cerr << "[Info] Tagset Size: " << TAGNUM << std::endl;
	tagname.resize(TAGNUM);
	
	for(int i = 0; i < TAGNUM; i++){
		source >> tagname[i];
	}
	for(int i = 0; i < TAGNUM; i++){
		source >> initial[i];
	}

	source >> buffer;
	assert(buffer == "#transition");

	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			source >> transition[i][j];
		}
	}

	source >> buffer;
	assert(buffer == "#trigram");

	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			for(int k = 0; k < TAGNUM; k++){		
				source >> trigram[i][j][k];
			}
		}
	}

	source >> buffer;
	assert(buffer == "#emission");
	
	source >> EQCLASS;
	for(int i = 0; i < TAGNUM; i++){
		emission[i] = new double[EQCLASS]();		
		for(int j = 0; j < EQCLASS; j++){
			source >> emission[i][j];
		}
	}
	std::cerr << "[Info] Equivalent Classes: " << EQCLASS << std::endl;
}

void save_model(std::string const& name, std::vector<std::string> const& tagname){
	std::ofstream destin(name, std::ofstream::out);
	std::string buffer;

	destin << std::fixed << std::setprecision(9);// << std::endl;

	destin << "#initial " << std::endl;
	destin << TAGNUM << " " << std::endl;

	for(int i = 0; i < TAGNUM; i++){
		destin << tagname[i] << " ";
	}
	destin << std::endl;
	for(int i = 0; i < TAGNUM; i++){
		destin << initial[i] << " ";
	}

	destin << std::endl << "#transition " << std::endl;

	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			destin << transition[i][j] << " ";
		}
		destin << std::endl;
	}

	destin << std::endl << "#trigram " << std::endl;

	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < TAGNUM; j++){		
			for(int k = 0; k < TAGNUM; k++){		
				destin << trigram[i][j][k] << " ";
			}
			destin << std::endl;
		}
	}

	destin << std::endl << "#emission " << std::endl;
	
	destin << EQCLASS << " " << std::endl;
	for(int i = 0; i < TAGNUM; i++){		
		for(int j = 0; j < EQCLASS; j++){
			destin << emission[i][j] << " ";
		}
		destin << std::endl;
	}
}