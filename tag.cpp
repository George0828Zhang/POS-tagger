#include <iostream>
#include <fstream>
#include <sstream>//parsing line
#include <string>
#include <cassert>
#include <array>
#include <vector>
#include <numeric>// for accumulate
#include <limits>// for inf
#include <cmath>// for log
#include <iomanip>// for set precision

#include "Tree.hpp"
#include "Array.hpp"
#include "tag.hpp"

constexpr int MAXTAGNUM = 62;
constexpr int RARETHRES = 20;
constexpr int MAXSUFFIX = 10;
constexpr int BEAMSIZE = 50;
constexpr double inf = std::numeric_limits<double>::infinity();

using Pair = std::array<int, 2>;

lexiTree dict_map; 
suffixTree suffix_map;
double initial[MAXTAGNUM];
double begin[MAXTAGNUM];
double end[MAXTAGNUM];
double transition[MAXTAGNUM][MAXTAGNUM];
double trigram[MAXTAGNUM][MAXTAGNUM][MAXTAGNUM];
double* emission[MAXTAGNUM];
DyArray<double> suf_emission;

int TAGNUM;
int EQCLASS;
int CORPUS_N;
int t_DET;

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
				// delta[{0, i,j}] = ( log(initial[i]) + log(emi_i) );
				delta[{0, i,j}] = ( log(begin[i]) + log(emi_i) );
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
		double prob = delta[{slen-1,i,j}] * end[j];	
		if(prob>best_prob){
			best_tail_i = i;
			best_tail_j = j;
			best_prob = prob;
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
		// delta[{0,i}] = initial[i] * emi_i;
		delta[{0,i}] = begin[i] * emi_i;
	}

	// t >= 1
	for(int t = 1; t < slen;t++){
		index = dict_map.wordIndex(sentence[t]);

		for(int i = 0; i < TAGNUM; i++){
			double emi_i = index==-1 ? oovEmission(sentence[t], i) : emission[i][index];
			for(int j = 0; j < TAGNUM; j++){

				double prob = delta[{t-1,j}] * transition[j][i] * emi_i;

				// if(t>1 && phi[{t,j}]==t_DET){
				// 	prob = delta[{t-1,j}] * trigram[t_DET][j][i] * emi_i;
				// }

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
		double prob = delta[{slen-1,i}] * end[i];
		if(prob>best_prob){
			best_tail_i = i;
			best_prob = prob;
		}
	}


	// back tracking	
	for(int t = slen-1; t >= 0; t--){
		tag[t] = best_tail_i;
		if (t==0) break;
		best_tail_i = phi[{t,best_tail_i}];
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
	for(int i = 0; i < TAGNUM; i++){
		begin[i] = (lbd[0]*initial[i] + lbd[1]*begin[i]) / (lbd[0] + lbd[1]);
		end[i] = (lbd[0]*initial[i] + lbd[1]*end[i]) / (lbd[0] + lbd[1]);
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
void smooth_probabilities(){
	smooth_transition();
	smooth_trigram();
	smooth_emission();
}

struct compare {
    bool operator()(const std::string& first, const std::string& second) {
        return first.size() < second.size();
    }
} slenComp;

void load_lexicon(std::string const& name){
	std::ifstream source(name, std::ifstream::in);
	std::string buffer;
	
	int freq;
	CORPUS_N = 0;

	std::vector<std::vector<int> > f_S_Eqv;
	// std::vector<int> suffFreq;
	// std::vector<std::string> Unique_Suffix;
	std::vector<int> classFreq(EQCLASS, 0);
	int d, wcount=0;
	int max_suf = -1;

	std::string sentence;
	while(std::getline(source, sentence)){
		if(!(std::stringstream(sentence) >> buffer >> freq >> d)) continue;
		std::string vocab(buffer);
		dict_map.insert(vocab, d);
		assert(wordIndex(vocab) == d);
		wcount++;
		CORPUS_N += freq;
		classFreq[d] += freq;
		if(freq < RARETHRES){
			for(int i = 0; i < MAXSUFFIX && i < vocab.size(); i++){
				std::string suf = vocab.substr(vocab.size()-1-i,i+1);
				int suf_index = suffix_map.wordIndex(suf);
				if(suf_index==-1){
					// Unique_Suffix.push_back(suf);
					// suffFreq.push_back(0);
					f_S_Eqv.push_back(std::vector<int>(EQCLASS, 0));
					suf_index = f_S_Eqv.size() - 1;
					suffix_map.insert(suf, suf_index);
					max_suf = suf_index;
				}
				f_S_Eqv[suf_index][d] += freq;
				// suffFreq[suf_index] += freq;
			}
		}
	}
	std::cerr << "[Info] Dictionary Size: " << wcount << std::endl;
	std::cerr << "[Info] Suffix Map Size: " << max_suf + 1 << std::endl;

	// suffix smoothing	
	suf_emission.reshape({TAGNUM, max_suf + 1});

	for(int tg = 0; tg < TAGNUM; tg++){
		for(int s = 0; s < max_suf + 1; s++){
			double prob = 0.;
			for(int q = 0; q < EQCLASS; q++){
				double P_S_Eqv = classFreq[q] ? ((double)f_S_Eqv[s][q])/(double)classFreq[q] : 0.;
				prob += emission[tg][q] * P_S_Eqv;
			}
			suf_emission[{tg, s}] = prob;
		}		
	}

	// smoothing
	// double meanP = std::accumulate(initial, initial + TAGNUM, 0) / TAGNUM;
	// double varP = 0;
	// for(int i = 0; i < TAGNUM; i++){
	// 	varP += (initial[i]-meanP)*(initial[i]-meanP);
	// }
	// varP /= (TAGNUM - 1);


	// std::sort(Unique_Suffix.begin(), Unique_Suffix.end(), slenComp);
	// assert(Unique_Suffix[0].size()<Unique_Suffix.back().size());
	// for(auto& suffix : Unique_Suffix){
	// 	int l = suffix.size();
	// 	int index = suffix_map.suffixIndex(suffix);
	// 	int next = l>1 ? suffix_map.suffixIndex(suffix.substr(1,l-1)) : -1;

	// 	for(int t = 0; t < TAGNUM; t++){
	// 		double tag_freq = CORPUS_N * initial[t];
	// 		double next_emit_inv = 0;
	// 		if(l==1){
	// 			next_emit_inv = initial[t];
	// 		}
	// 		else{
	// 			next_emit_inv = suf_emission[{t, next}] * tag_freq / suffFreq[next];
	// 		}

	// 		double emit_inv = (1 + varP * next_emit_inv)/(1 + varP);
	// 		suf_emission[{t, index}] = emit_inv * suffFreq[index] / tag_freq;
	// 	}
	// }
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
		if(tagname[i]=="DET")
			t_DET = i;
	}
	for(int i = 0; i < TAGNUM; i++){
		source >> initial[i];
	}

	source >> buffer;
	assert(buffer == "#begin");
	for(int i = 0; i < TAGNUM; i++){
		source >> begin[i];
	}

	source >> buffer;
	assert(buffer == "#end");
	for(int i = 0; i < TAGNUM; i++){
		source >> end[i];
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