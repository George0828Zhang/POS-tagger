#!/usr/bin/env python3
tagged_sentences = []

# from nltk.corpus import masc_tagged as corpus
# tagged_sentences += corpus.tagged_sents(tagset='universal');


from nltk.corpus import treebank as corpus
tagged_sentences += corpus.tagged_sents(tagset='universal');


initial = {}
# words = set()
words = {}
wordcount = {}
# emission = {}
transition = {}

trigram = {}
corpus_N = 0

mod_brown = False
longest = 0

def tprocess(rawtag):	
	if mod_brown:
		rawtag = rawtag.replace('FW-', '').replace('-HL', '').replace('-TL', '').replace('-NC', '').replace('-T', '').replace('-N', '')
		if len(rawtag) > 1:
			rawtag = rawtag.replace('*', '')
	if rawtag == None:
		rawtag = "-NONE-"
	return rawtag

def process_token(word, tag, ptag=None, pptag=None):
	global corpus_N

	corpus_N += 1
		
	tag = tprocess(tag)

	if word not in words:
		words[word] = set()
		wordcount[word] = 0
	words[word].add(tag)
	wordcount[word] += 1

	if tag not in initial:
		initial[tag] = 0
	if tag not in transition:
			transition[tag] = {}

	# if "+" in tag and "\'" in word:
	# 	print (word, tag)
	# 	print (sentence)

	initial[tag] += 1
	if ptag != None:
		ptag = tprocess(ptag)	
		if tag not in transition[ptag]:
			transition[ptag][tag] = 0
		transition[ptag][tag] += 1
		if pptag != None:
			pptag = tprocess(pptag)	
			if (pptag, ptag) not in trigram:
				trigram[(pptag, ptag)] = {}
			if tag not in trigram[(pptag, ptag)]:
				trigram[(pptag, ptag)][tag] = 0
			trigram[(pptag, ptag)][tag] += 1

# preparation
for sentence in tagged_sentences:
	for i, (word, tag) in enumerate(sentence):
		ptag, pptag = None, None
		if i > 0:
			ptag = sentence[i-1][1]
			if i > 1:
				pptag = sentence[i-2][1]
		process_token(word, tag, ptag, pptag);
# def cut(w_t):
# 	cut = w_t.rfind('_')
# 	return w_t[:cut], w_t[cut+1:]
# with open("browntag_nolines.txt", "r") as fobj:
# 	for sentence in fobj:		
# 		tokens = sentence.strip().split()
# 		tagged_sentences.append(tokens)
# 		for i, word_tag in enumerate(tokens):
# 			word, tag = cut(word_tag)

# 			ptag, pptag = None, None
# 			if i > 0:
# 				trash, ptag = cut(tokens[i-1])
# 				if i > 1:
# 					trash, pptag = cut(tokens[i-2])
# 			process_token(word, tag, ptag, pptag);

best100 = sorted(wordcount, key=lambda x: wordcount[x])[-100:]

# output probability Model
file = open("model.txt", "wt")
file_lex = open("lexicon.txt", "wt")

# initial probs
file.write("#initial \n")
file.write("{} \n".format(len(initial)))
denom = sum(initial.values())
assert(denom==corpus_N)
for tag in initial:
	file.write("{} ".format(tag))
file.write("\n")
for tag in initial:
	prob = initial[tag] / denom
	file.write("{:.9f} ".format(prob))

file.write("\n#transition \n")
# transition probs
# lbd = [0, 0]
# for t1 in transition:
# 	for t2 in transition[t1]:		
# 		fq = transition[t1][t2]
# 		p1 = (fq - 1)/(initial[t1]-1) if initial[t1]>1 else 0
# 		p2 = (initial[t2]-1)/(corpus_N-1) if corpus_N>1 else 0
# 		mx = p1 if p1 > p2 else p2
# 		if mx == p1:
# 			lbd[1] += fq
# 		elif mx == p2:
# 			lbd[0] += fq
# sumlbd = sum(lbd)
sumlbd = 1
lbd = [0, 1]

for ptag in initial:
	denom = sum(transition[ptag].values())
	for tag in initial:
		domi = 0 if tag not in transition[ptag] else transition[ptag][tag]
		transition[ptag][tag] = domi
		bi = 0 if denom==0 else (domi / denom)
		term1 = lbd[0]*initial[tag]/corpus_N
		term2 = lbd[1]*bi
		file.write("{:.9f} ".format((term1+term2)/sumlbd))
	file.write("\n")

file.write("\n#trigram \n")
# traigram transition probs
# lbd = [0, 0, 0]
# for key in trigram:
# 	for t3 in trigram[key]:
# 		(t1, t2) = key
# 		fq = trigram[key][t3]
# 		p1 = (fq - 1)/(transition[t1][t2]-1) if transition[t1][t2]>1 else 0
# 		p2 = (transition[t2][t3]-1)/(initial[t2]-1) if initial[t2]>1 else 0
# 		p3 = (initial[t3] - 1)/(corpus_N-1) if corpus_N>1 else 0
# 		mx = p3 if p3 > p2 else p2
# 		mx = mx if mx > p1 else p1
# 		if mx == p1:
# 			lbd[2] += fq
# 		elif mx == p2:
# 			lbd[1] += fq
# 		elif mx == p3:
# 			lbd[0] += fq
# sumlbd = sum(lbd)
lbd = [0, 0, 1]

for pptag in initial:
	for ptag in initial:
		hastri = (pptag, ptag) in trigram
		if hastri:
			denom = sum(trigram[(pptag, ptag)].values())
		for tag in initial: 
			term1 = lbd[0]*initial[tag]/corpus_N
			term2 = lbd[1]*transition[ptag][tag]/initial[ptag]
			if hastri and (tag in trigram[(pptag, ptag)]):
				term3 = trigram[(pptag, ptag)][tag] / denom
			else:
				term3 = 0
			file.write("{:.9f} ".format((term1+term2+term3)/sumlbd))
			# if((term1+term2+term3)/sumlbd <= 0):
			# 	print(lbd, initial[tag], transition[ptag][tag])
		file.write("\n")

Eqv = []
for word in words:	
	if word in best100:
		words[word] = word
	cl = words[word]
	if cl not in Eqv:
		Eqv.append(cl)


emission = {}
for sentence in tagged_sentences:
	for i, (word, tag) in enumerate(sentence):
		tag = tprocess(tag)
		if tag not in emission:
			emission[tag] = [0]*len(Eqv)
		eqclass = Eqv.index(words[word])
		
		emission[tag][eqclass] += 1
# for tokens in tagged_sentences:
# 	for word_tag in tokens:
# 		word, tag = cut(word_tag)
# 		tag = tprocess(tag)
# 		if tag not in emission:
# 			emission[tag] = [0]*len(Eqv)
# 		eqclass = Eqv.index(words[word])
		
# 		emission[tag][eqclass] += 1
		


file.write("\n#emission \n")
file.write("{} \n".format(len(Eqv)))
# emission probs
for tag in initial:
	denom = sum(emission[tag])
	# file.write("{} ".format(len(emission[tag])))
	for i in range(len(Eqv)): 
		domi = emission[tag][i]
		logprob = domi / denom
		file.write("{:.9f} ".format(logprob))
	file.write("\n")


# file_lex.write("\n#vocab_freq_Eqv\n")
for word in words:
	file_lex.write("{} {} {}\n".format(word, wordcount[word], Eqv.index(words[word])))	
file_lex.write("\n")



