#!/usr/bin/env python3
tagged_sentences = []

from nltk.corpus import masc_tagged as corpus
tagged_sentences += corpus.tagged_sents();
# untag_sentences += corpus.sents();

data = open("testdata.txt", "wt")
res = open("answer.txt", "wt")

amount = 100

for sentences in tagged_sentences:
	for (word, tag) in sentences:
		data.write("{} ".format(word))
		res.write("{}/{} ".format(word, tag))
	data.write("\n")
	res.write("\n")
	amount -= 1
	if amount == 0:
		break