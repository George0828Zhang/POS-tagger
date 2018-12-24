#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#define MAXTAGCNT 200
#define MAXWDLEN 50
#define MAXTAGLEN 10

using namespace std;

vector<string> tag_name;
unordered_map<string, int> tag_map;
int tag_uni[MAXTAGCNT]={};
int tag_bi[MAXTAGCNT][MAXTAGCNT]={};
int tag_tri[MAXTAGCNT][MAXTAGCNT][MAXTAGCNT]={};

vector<string> dictionary;
vector<int> wordFreq;
vector<string> canTags;
unordered_map<string, int> word_map;

void parse(char* filename);

int main(int argc, char** argv){
	for(int i = 1; i < argc; i++){
		parse(argv[i]);
	}
}


int tagIndex(string const& t){
	auto iter = tag_map.find(t);
	return iter==tag_map.end() ? -1 : iter->second;
}
int wordIndex(string const& w){
	auto iter = word_map.find(w);
	return iter==word_map.end() ? -1 : iter->second;
}
bool isSpace(char c){
	return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}
void trim(string& s){
    s.erase(remove_if(s.begin(), s.end(), isSpace), s.end());
}
void modifyBrownTag(string& t){
	//'FW-''-HL''-TL''-NC''-T''-N''*'
	// auto plus = {"FW-", "-HL", "-TL", "-NC", "-T", "-N"};
	size_t x = t.find("fw-");
	if(x != string::npos){
		t.erase(x, 3);
	}
	x = t.find("-hl");
	if(x != string::npos){
		t.erase(x, 3);
	}
	x = t.find("-tl");
	if(x != string::npos){
		t.erase(x, 3);
	}
	x = t.find("-nc");
	if(x != string::npos){
		t.erase(x, 3);
	}
	x = t.find("-t");
	if(x != string::npos){
		t.erase(x, 2);
	}
	x = t.find("-n");
	if(x != string::npos){
		t.erase(x, 2);
	}
	if(t.size()<=1) return;
	x = t.find('*');
	if(x != string::npos){
		t.erase(x, 1);
	}
	return;
}


void parse(char* filename){
	fstream file(filename, fstream::in);
	string line;
	while(!file.eof()){		
		getline(file, line);
		vector<string> tokens;
		string word_tag;
		// cout << line << endl;
		istringstream lstream(line);
		int cx1_index = -1;
		int cx2_index = -1;
		// cx2 -> cx1 -> cur
		while(!lstream.eof()){
			getline(lstream, word_tag, ' ');
			trim(word_tag);
			if(word_tag.size()>0){
				size_t slh = word_tag.find('/');
				if(slh == string::npos) continue;
				string word = word_tag.substr(0, slh);
				string tag = word_tag.substr(slh + 1, word_tag.size() - slh);				
				modifyBrownTag(tag);

				// handling tag
				int t_index = tagIndex(tag);
				if(t_index==-1){
					tag_name.push_back(tag);
					t_index = tag_name.size() - 1;
					tag_map[tag] = t_index;
				}

				tag_uni[t_index] += 1;
				if(cx1_index!=-1){
					tag_bi[cx1_index][t_index] += 1;
					if(cx2_index!=-1)
						tag_tri[cx2_index][cx1_index][t_index] += 1;
				}

				cx2_index = cx1_index;
				cx1_index = t_index;

				// handling word
				
			}
		}
	}	
}
