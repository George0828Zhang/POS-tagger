#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
#include <limits>
#include <cmath>

#define MAXWDLEN 50
#define MAXTAGLEN 10

using namespace std;



void parse(char* filename);

int main(int argc, char** argv){
	for(int i = 1; i < argc; i++){
		parse(argv[i]);
	}
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
		while(!lstream.eof()){
			getline(lstream, word_tag, ' ');
			cout << word_tag << endl;
		}
	}	
}
