#include <vector>
#include <string>
#include <iostream>// for cin,cout,cerr,endl
#include <fstream>// for input file
#include <sstream>// for parsing ord
#include <cassert>
#include <chrono>// for time evaluation
#include "tag.hpp"

using uchar = unsigned char;

void parseInput(
	int argc, 
	char** argv, 
	std::string& fmodel, 
	std::string& flexicon,
	std::string& finput,
	std::string& foutput,
	int& order);

int main(int argc, char** argv){
	auto start = std::chrono::system_clock::now();
	std::string modelname;
	std::string lexiname;
	std::string input;
	std::string output;
	int order;

	parseInput(argc,argv,modelname,lexiname,input,output,order);

	std::ifstream readfromf(input, std::fstream::in);
	std::ofstream writetof(output, std::fstream::out);

	std::istream& in = input.size() ? readfromf : std::cin;
	std::ostream& out = output.size() ? writetof : std::cout;

	std::cerr << "[Info] Mode: " << (input.size() ? "File" : "Input") << \
	" to " << (output.size() ? "file" : "output")  << ", order = " << order << std::endl;

	std::vector<std::string> tagname;
	load_model(modelname, tagname);
	load_lexicon(lexiname);
	smooth_probabilities();

	std::string sentence;
	while(getline(in, sentence)){
		// std::cerr << "[Debug] sent: " << sentence << "\n";
		std::vector<std::string> tokenized;
		if(input.size())
			tokenizer(sentence, tokenized);
		else
			tokenizer2(sentence, tokenized);
		int T = tokenized.size();

		// Capitalization
		// std::string Cap(tokenized[0]);
		// if(tokenized[0][0]>='A' && tokenized[0][0]<='Z') 
		// tokenized[0][0] = tolower(tokenized[0][0]);

		// POS tagging
		std::vector<int> tags(T);
		if(order==3)
			POStag3(tokenized, tags);
		else
			POStag(tokenized, tags);
		// tokenized[0] = Cap;
		for(int i = 0; i < T; i++){
			out << tokenized[i] << "/" << tagname[tags[i]] << " ";
		}
		out << std::endl;
	}

	auto end = std::chrono::system_clock::now(); 
    std::chrono::duration<double> elapsed_seconds = end-start; 
    std::cerr << "[Info] Elapsed time: " << elapsed_seconds.count() << "s\n";
}

void parseInput(
	int argc, 
	char** argv, 
	std::string& fmodel, 
	std::string& flexicon,
	std::string& finput,
	std::string& foutput,
	int& order)
{
	order = 2;
	for(int i = 1; i < argc - 1; i++){
		std::string input(argv[i]);
		if(input=="-m"){
			fmodel = std::string(argv[i+1]);
		}else if(input=="-l"){
			flexicon = std::string(argv[i+1]);
		}else if(input=="-i"){
			finput = std::string(argv[i+1]);
		}else if(input=="-o"){
			foutput = std::string(argv[i+1]);
		}else if(input=="-ord"){
			std::stringstream(argv[i+1]) >> order;
			order = std::min(3, std::max(2, order));
		}
	}
}

void tokenizer2(std::string const& sentence, std::vector<std::string>& tokenized){
	const std::string common_ch = "~!@#$%%^&*()_+=`\\/.,<>?\":;";
	// "-" might be hyphenized text, and ' and . might be abreviation
	int slen = sentence.size();
	int tlen = 0;
	int c_at;
	for(int i = 0; i < slen; i++){
		uchar c = sentence[i];
		if(common_ch.find(c) != std::string::npos){
			if(tlen){
				// process sub-word tokenizing TODO
				std::string word = sentence.substr(i-tlen, tlen);
				c_at = word.find('\'');
				if(c_at != std::string::npos && c_at != 0){
					tokenized.push_back(word.substr(0, c_at));
					word.erase(0, c_at);
				}				
				tokenized.push_back(word);
			}
			tokenized.push_back(std::string(1, c));
			tlen = 0;
		}else if(tlen && (c == ' ' || c == '\n' || c == '\t' || c == '\r')){
			// process sub-word tokenizing TODO
			std::string word = sentence.substr(i-tlen, tlen);
			c_at = word.find('\'');
			if(c_at != std::string::npos && c_at != 0){
				tokenized.push_back(word.substr(0, c_at));
				word.erase(0, c_at);
			}				
			tokenized.push_back(word);
			tlen = 0;
		}else if(i == slen - 1){
			// process sub-word tokenizing TODO
			std::string word = sentence.substr(i-tlen, 1+tlen);
			c_at = word.find('\'');
			if(c_at != std::string::npos && c_at != 0){
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

void tokenizer(std::string const& sentence, std::vector<std::string>& tokenized){

	int slen = sentence.size();
	int tlen = 0;
	int c_at;
	for(int i = 0; i < slen; i++){
		uchar c = sentence[i];
		if(tlen && (c == ' ' || c == '\n' || c == '\t' || c == '\r')){
			std::string word = sentence.substr(i-tlen, tlen);							
			tokenized.push_back(word);
			tlen = 0;
		}else if(i == slen - 1){
			std::string word = sentence.substr(i-tlen, 1+tlen);							
			tokenized.push_back(word);
			tlen = 0;
		}
		else
			tlen++;
	}
}