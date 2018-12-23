#!/usr/bin/env python3
# from nltk.corpus import treebank
from nltk.corpus import brown as corpus

tagged_sentences = corpus.tagged_sents();

initial = {}
# words = set()
words = {}
wordcount = {}
# emission = {}
transition = {}


# preparation
for sentence in tagged_sentences:
	for i, (word, tag) in enumerate(sentence):
		tag = tag.replace('FW-', '').replace('-HL', '').replace('-TL', '').replace('-NC', '').replace('-T', '').replace('-N', '')
		if len(tag) > 1:
			tag = tag.replace('*', '')

		# words.add(word)
		if word not in words:
			words[word] = set()
			wordcount[word] = 0
		words[word].add(tag)
		wordcount[word] += 1

		if tag not in initial:
			# emission[tag] = {}
			initial[tag] = 0
		# if word not in emission[tag]:
		# 	emission[tag][word] = 0
		# emission[tag][word] += 1
		initial[tag] += 1
		if i > 0:
			ptag = sentence[i-1][1].replace('FW-', '').replace('-HL', '').replace('-TL', '').replace('-NC', '').replace('-T', '').replace('-N', '')
			if len(ptag) > 1:
				ptag = ptag.replace('*', '')
			if ptag not in transition:
				transition[ptag] = {}
			if tag not in transition[ptag]:
				transition[ptag][tag] = 0
			transition[ptag][tag] += 1

# create equivalence classes
Eqv = []
best100 = sorted(wordcount, key=lambda x: wordcount[x])[-100:]

# output probability Model
file = open("model.txt", "wt")

# initial probs
file.write("#initial \n")
file.write("{} \n".format(len(initial)))
denom = sum(initial.values())

for tag in initial:
	file.write("{} ".format(tag))
file.write("\n")
for tag in initial:
	logprob = initial[tag] / denom
	file.write("{:.9f} ".format(logprob))

file.write("\n#transition \n")
# transition probs
for ptag in initial:
	denom = sum(transition[ptag].values())
	for tag in initial: 
		domi = 0 if tag not in transition[ptag] else transition[ptag][tag]
		logprob = domi / denom
		file.write("{:.9f} ".format(logprob))
	file.write("\n")

file.write("\n#vocab \n")
for word in words:
	cl = words[word]
	if word in best100:
		Eqv.append(word)
		file.write("{} {} ".format(word,  len(Eqv)-1))
		continue
	elif cl not in Eqv:
		Eqv.append(cl)
	file.write("{} {} ".format(word,  Eqv.index(cl)))
file.write("\n")
# print(len(Eqv)) 571


emission = {}
for sentence in tagged_sentences:
	for i, (word, tag) in enumerate(sentence):
		tag = tag.replace('FW-', '').replace('-HL', '').replace('-TL', '').replace('-NC', '').replace('-T', '').replace('-N', '')
		if len(tag) > 1:
			tag = tag.replace('*', '')

		if tag not in emission:
			emission[tag] = [0]*len(Eqv)
		
		if word in best100:
			eqclass = Eqv.index(word)
		else:
			eqclass = Eqv.index(words[word])
		
		emission[tag][eqclass] += 1


file.write("\n#emission \n")
file.write("{}\n".format(len(Eqv)))
# emission probs
for tag in initial:
	denom = sum(emission[tag])
	# file.write("{} ".format(len(emission[tag])))
	for i in range(len(Eqv)): 
		domi = emission[tag][i]
		logprob = domi / denom
		file.write("{:.9f} ".format(logprob))
	file.write("\n")



