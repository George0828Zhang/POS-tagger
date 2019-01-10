#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include "tag.hpp"

int main(int argc, char** argv){
	std::vector<std::string> tagname;
	assert(argc==3);
	load_model(argv[1], tagname);
	load_lexicon(argv[2]);
	std::string sentence;
	while(getline(std::cin, sentence)){
		// cerr << "input received." << endl;
		// if(sentence.size()<2) continue;
		std::vector<std::string> tokenized;
		tokenizer(sentence, tokenized);
		// Capitalization
		std::string Cap(tokenized[0]);
		if(tokenized[0][0]>='A' && tokenized[0][0]<='Z') 
			tokenized[0][0] = tolower(tokenized[0][0]);
		// POS tagging
		std::vector<int> tags(tokenized.size());
		POStag3(tokenized, tags);
		tokenized[0] = Cap;
		for(int i = 0; i < tags.size(); i++){
			std::cout << tokenized[i] << "/" << tagname[tags[i]] << " ";
		}
		std::cout << std::endl;
	}
}



void tokenizer2(std::string const& sentence, std::vector<std::string>& tokenized){
	const std::string common_ch = "~!@#$%%^&*()_+=`\\/.,<>?\":;";
	// "-" might be hyphenized text, and ' and . might be abreviation
	int slen = sentence.size();
	int tlen = 0;
	int c_at;
	for(int i = 0; i < slen; i++){
		char c = sentence[i];
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
		char c = sentence[i];
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