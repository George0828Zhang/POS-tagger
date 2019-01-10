#ifndef _TAG_HPP_
#define _TAG_HPP_
#include <string>
#include <vector>
void load_model(std::string const& name, std::vector<std::string>& tagname);
void save_model(std::string const& name, std::vector<std::string> const& tagname);
void load_lexicon(std::string const& name);
int wordIndex(std::string const& word);
void POStag3(std::vector<std::string> const& sentence, std::vector<int>& tag);

void tokenizer2(std::string const& sentence, std::vector<std::string>& tokenized);
void tokenizer(std::string const& sentence, std::vector<std::string>& tokenized);
#endif