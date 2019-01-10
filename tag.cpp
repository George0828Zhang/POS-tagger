#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <cassert>
#include <array>
#include <vector>
#include <limits>// for inf
#include <cmath>// for log
#include <iomanip>// for set precision
#include "Array.hpp"
#include "tag.hpp"

constexpr int MAXTAGNUM = 62;
constexpr int RARETHRES = 10;
constexpr int MAXSUFFIX = 10;
constexpr int BEAMSIZE = 100;
constexpr double inf = std::numeric_limits<double>::infinity();

using Pair = std::array<int, 2>;

std::unordered_map<std::string, int> dict_map;
std::unordered_map<std::string, int> suffix_map;
double initial[MAXTAGNUM];
double transition[MAXTAGNUM][MAXTAGNUM];
double trigram[MAXTAGNUM][MAXTAGNUM][MAXTAGNUM];
double* emission[MAXTAGNUM];
DyArray<double> suf_emission;

int TAGNUM;
int EQCLASS;

int wordIndex(std::string const& word){
	auto iter = dict_map.find(word);
	return iter==dict_map.end() ? -1 : iter->second;
}
int suffIndex(std::string const& word){
	auto iter = suffix_map.find(word);
	return iter==suffix_map.end() ? -1 : iter->second;
}

double oovEmission(std::string const& word, int tg){
	int wlen = word.size();
	for(int l = std::min(MAXSUFFIX, wlen); l > 0; l--){
		int index = suffIndex(word.substr(wlen-l, l));
		if(index!=-1){
			return suf_emission[{tg, index}];
		}
	}
	return (1.0/TAGNUM);
}


DyArray<double> delta;
DyArray<int> phi;
Kmax<double, Pair> beamspace(BEAMSIZE);
std::vector<Pair> beam;
void POStag3(std::vector<std::string> const& sentence, std::vector<int>& tag){
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
				beamspace.insert(delta[{i,j,1}], {i,j});
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
						double logprob = delta[{j,k,t-1}] + log(trigram[k][j][i]) + log(emi_i);
						if(logprob > delta[{i,j,t}]){
							delta[{i,j,t}] = logprob;
							phi[{i,j,t}] = k;							
						}
					}
					beamspace.insert(delta[{i,j,t}], {i,j});
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






DyArray<float> alpha;
DyArray<float> beta;
DyArray<float> Gamma;
DyArray<float> epsilon;
DyArray<float> epsilon_3;

void initEmission(DyArray<float>& ehat, std::vector<std::string> const& sentence);
void makeAlpha(std::vector<std::string> const& sentence, DyArray<float> & ehat);
void makeBeta(std::vector<std::string> const& sentence, DyArray<float>& ehat);
void makeGamma(int T);
void makeEpsilons(int T, DyArray<float>& ehat);

void refineHMM(std::vector<std::vector<std::string>> const& sentences){
	// enhance speed by calculating oovemission proprocessing
	int samples = sentences.size();
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

	for(auto const& sent : sentences ){
		int T = sent.size();
		initEmission(emit_hat, sent);
		makeAlpha(sent, emit_hat);
		makeBeta(sent, emit_hat);
		makeGamma(T);
		makeEpsilons(T, emit_hat);

		for(int t = 0; t < T - 2; t++){
			// cummulate for initial
			for(int i = 0; i < TAGNUM; i++){
				initial_cummu[{i}] += Gamma[{0, i}];
			}

			// cummulate for transition
			for(int i = 0; i < TAGNUM; i++){
				for(int j = 0; j < TAGNUM; j++){
					transition_cummu[{i, j}] += epsilon[{t, i, j}];
				}
				visit_cummu[{i}] += Gamma[{t, i}];
			}

			// cummulate for trigram
			for(int i = 0; i < TAGNUM; i++){
				for(int j = 0; j < TAGNUM; j++){
					for(int k = 0; k < TAGNUM; k++){
						trigram_cummu[{i, j, k}] += epsilon_3[{t, i, j, k}];
					}
				}
			}

			// cummulate for emission
			int index = wordIndex(sent[t]);
			if(index != -1){
				for(int i = 0; i < TAGNUM; i++){
					observ_cummu[{i, index}] += Gamma[{t, i}];
				}
			}
		}
	}

	// reestimate model
	for(int i = 0; i < TAGNUM; i++){
		// initial
		initial[i] = initial_cummu[{i}] / samples;
	}

	for(int i = 0; i < TAGNUM; i++){
		// transition
		for(int j = 0; visit_cummu[{i}] > 0 && j < TAGNUM; j++){
			transition[i][j] = transition_cummu[{i, j}] / visit_cummu[{i}];
			for(int k = 0; transition_cummu[{i, j}] > 0 && k < TAGNUM; k++){
				trigram[i][j][k] = trigram_cummu[{i, j, k}] / transition_cummu[{i, j}];
			}
		}

		for(int j = 0; visit_cummu[{i}] > 0 && j < EQCLASS; j++){
			emission[i][j] = observ_cummu[{i, j}] / visit_cummu[{i}];
		}
	}


}
void initEmission(DyArray<float>& ehat, std::vector<std::string> const& sentence){
	int T = sentence.size();
	ehat.reshape({T, TAGNUM});
	for(int t = 0; t < T; t++){
		int index = wordIndex(sentence[t]);
		for(int i = 0; i < TAGNUM; i++){
			ehat[{t, i}] = index==-1 ? (float)oovEmission(sentence[t], i) : (float)emission[i][index];
		}
	}
}
void makeAlpha(std::vector<std::string> const& sentence, DyArray<float>& ehat){
	int T = sentence.size();

	alpha.reshape({T, TAGNUM});
	alpha.clear(0);

	// t = 0
	for(int i = 0; i < TAGNUM; i++){
		float emi_i = ehat[{0, i}];
		if(emi_i>0){
			alpha[{0, i}] = initial[i] * emi_i;
		}
	}

	if(T > 1){
		// t = 1
		// index = indices[1];
		for(int i = 0; i < TAGNUM; i++){
			float emi_i = ehat[{1, i}];
			for(int j = 0; j < TAGNUM; j++){				
				if(emi_i>0 && transition[j][i]>0){
					alpha[{1, i}] += alpha[{0, j}] * transition[j][i] * emi_i;	
				}
			}
		}

		// t > 1
		for(int t = 2; t < T;t++){
			// index = indices[t];
			for(int i = 0; i < TAGNUM; i++){
				float emi_i = ehat[{t, i}];
				if(emi_i>0){
					// int index2 = indices[t-1];
					for(int j = 0; j < TAGNUM; j++){
						float emi_j = ehat[{t-1, j}];
						if(emi_j>0){
							float prior = 0.;
							for(int k = 0; k < TAGNUM; k++){
								if(trigram[k][j][i] > 0){
									prior += alpha[{t-2, j, k}] * trigram[k][j][i];
								}
							}
							alpha[{t, i}] = prior * emi_i * emi_j;
						}
					}
				}
			}
		}
	}	
}
void makeBeta(std::vector<std::string> const& sentence, DyArray<float>& ehat){
	int T = sentence.size();

	beta.reshape({T, TAGNUM});
	beta.clear(0);

	// t = T - 1
	for(int i = 0; i < TAGNUM; i++){
		beta[{T-1, i}] = 1.;
	}

	if(T > 1){
		// t = T - 2
		for(int j = 0; j < TAGNUM; j++){			
			for(int i = 0; i < TAGNUM; i++){	
				float emi_i = ehat[{T-1, i}];			
				if(emi_i>0 && transition[j][i]>0){
					beta[{T-2, j}] += beta[{T-1, i}] * transition[j][i] * emi_i;	
				}
			}
		}

		// t > 1
		for(int t = T-3; t >= 0;t--){
			for(int k = 0; k < TAGNUM; k++){
				for(int j = 0; j < TAGNUM; j++){	
					float emi_j = ehat[{t+1, j}];	
					if(emi_j>0){	
						for(int i = 0; i < TAGNUM; i++){	
							float emi_i = ehat[{t+2, i}];			
							if(emi_i>0){
								beta[{t, k}] += beta[{t+2, i}] * trigram[k][j][i] * emi_i * emi_j;	
							}
						}
					}
				}
			}
		}
	}	
}
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

	epsilon_3.reshape({T, TAGNUM, TAGNUM, TAGNUM});
	
	for(int t = 0; t < T - 2; t++){
		float sum = 0;
		for(int k = 0; k < TAGNUM; k++){
			for(int j = 0; j < TAGNUM; j++){
				for(int i = 0; i < TAGNUM; i++){
					epsilon_3[{t, k, j, i}] = alpha[{t, k}] * trigram[k][j][i] * ehat[{t+1, j}] * ehat[{t+2, i}] * beta[{t+2, i}];
					sum += epsilon_3[{t, k, j, i}];
				}
			}
		}
		if(sum > 0){
			for(int k = 0; k < TAGNUM; k++){
				for(int j = 0; j < TAGNUM; j++){
					for(int i = 0; i < TAGNUM; i++){
						epsilon_3[{t, k, j, i}] /= sum;
					}
				}
			}
		}
	}
}
// void makeAlpha(std::vector<std::string> const& sentence, std::vector<int> const& indices){
// 	int T = sentence.size();

// 	alpha.reshape({T, TAGNUM, TAGNUM});
// 	alpha.clear(0);

// 	// t = 0
// 	int index = indices[0];
// 	for(int i = 0; i < TAGNUM; i++){
// 		double emi_i = index==-1 ? oovEmission(sentence[0], i) : emission[i][index];
// 		if(emi_i>0){
// 			for(int j = 0; j < TAGNUM; j++){
// 				alpha[{0, i, j}] = initial[i] * emi_i;
// 			}
// 		}
// 	}

// 	if(T > 1){
// 		// t = 1
// 		index = indices[1];
// 		for(int i = 0; i < TAGNUM; i++){
// 			double emi_i = index==-1 ? oovEmission(sentence[1], i) : emission[i][index];
// 			for(int j = 0; j < TAGNUM; j++){				
// 				if(emi_i>0 && transition[j][i]>0){
// 					alpha[{1, i,j}] = alpha[{0, j,0}] * transition[j][i] * emi_i;	
// 				}
// 			}
// 		}

// 		// t > 1
// 		for(int t = 2; t < T;t++){
// 			index = indices[t];
// 			for(int i = 0; i < TAGNUM; i++){
// 				double emi_i = index==-1 ? oovEmission(sentence[t], i) : emission[i][index];
// 				if(emi_i>0){					
// 					for(int j = 0; j < TAGNUM; j++){
// 						double prior = 0.;
// 						for(int k = 0; k < TAGNUM; k++){
// 							if(trigram[k][j][i] > 0){
// 								prior += alpha[{t-1, j,k}] * trigram[k][j][i];
// 							}
// 						}
// 						alpha[{t, i, j}] = prior * emi_i; 
// 					}
// 				}
// 			}
// 		}
// 	}	
// }
// void makeBeta(std::vector<std::string> const& sentence, std::vector<int> const& indices){
// 	int T = sentence.size();

// 	beta.reshape({TAGNUM, T});
// 	beta.clear(0);

// 	// t = T-1
// 	int index = indices[T-1];
// 	for(int i = 0; i < TAGNUM; i++){
// 		beta[{i,T-1}] = 1.;
// 	}

// 	// t < T-1
// 	for(int t = T-2; t >= 0 ; t--){
// 		index = indices[t+1];
// 		for(int j = 0; j < TAGNUM; j++){
// 			for(int i = 0; i < TAGNUM; i++){
// 				double emi_i = index==-1 ? oovEmission(sentence[t+1], i) : emission[i][index];
// 				beta[{j,t}] += transition[j][i] * emi_i * beta[{i,t+1}];
// 			}
// 		}
// 	}
// }









void load_lexicon(std::string const& name){
	std::ifstream source(name, std::ifstream::in);
	std::string buffer;
	
	int freq;

	std::vector<std::vector<int> > f_S_Eqv;
	std::vector<int> sum_of_Suffix(EQCLASS, 0);
	int d;
	while(!source.eof()){
		source >> buffer >> freq >> d;
		std::string vocab(buffer);
		dict_map[vocab] = d;

		if(freq < RARETHRES){
			for(int i = 0; i < MAXSUFFIX && i < vocab.size(); i++){
				std::string suf = vocab.substr(vocab.size()-1-i,i+1);
				int suf_index = suffIndex(suf);
				if(suf_index==-1){					
					f_S_Eqv.push_back(std::vector<int>(EQCLASS, 0));
					suf_index = f_S_Eqv.size() - 1;
					suffix_map[suf] = suf_index;
				}
				f_S_Eqv[suf_index][d] += freq;
				sum_of_Suffix[d] += freq;
			}
		}
	}
	std::cerr << "[Info] Dictionary Size: " << dict_map.size() << std::endl;
	std::cerr << "[Info] Suffix Map Size: " << suffix_map.size() << std::endl;

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