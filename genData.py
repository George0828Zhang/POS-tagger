#!/usr/bin/env python3
untagged_sentences = []

# import nltk
# from nltk.corpus import brown as corpus
# untagged_sentences += corpus.sents()

# data = open("testdata.txt", "wt")
# data2 = open("traindata.txt", "wt")
# res = open("answer.txt", "wt")

# tru_amount = len(untagged_sentences)
# amount = 5000
# amount2 = 5000
# print("[Data] Extracted {} out of {} ({:.2f}%)".format(amount2, tru_amount, amount2/tru_amount*100))

# for sentence in untagged_sentences:
# 	# print(sentence)
# 	tagged = nltk.pos_tag(sentence)
# 	for word,tag in tagged:
# 		if(amount>0):
# 			data.write("{} ".format(word))
# 		else:
# 			data2.write("{} ".format(word))
# 		res.write("{}/{} ".format(word, tag))
# 	if(amount>0):
# 		data.write("\n")
# 	else:
# 		data2.write("\n")
# 	res.write("\n")
	
# 	amount -= 1
# 	amount2 -= 1
# 	if amount2 == 0:
# 		break


#!/usr/bin/env python3
tagged_sentences = []

from nltk.corpus import conll2000 as corpus
tagged_sentences += corpus.tagged_sents(tagset='universal')

import nltk
untagged_sentences += corpus.sents()

data = open("testdata.txt", "wt")
res = open("answer.txt", "wt")
control = open("control.txt", "wt")

tru_amount = len(tagged_sentences)
amount = tru_amount
print("[Data] Extracted {} out of {} ({:.2f}%)".format(amount, tru_amount, amount/tru_amount*100))


for sentences, untagged in zip(tagged_sentences, untagged_sentences):
	# if "" in [a[0] for a in sentences]:
	# 	continue
	tagged = nltk.pos_tag(untagged, tagset='universal')
	for (word, tag), (word2, tag2) in zip(sentences, tagged):
		data.write("{} ".format(word))
		res.write("{}/{} ".format(word, tag))
		control.write("{}/{} ".format(word2, tag2))
	data.write("\n")
	res.write("\n")
	control.write("\n")
	amount -= 1
	if amount == 0:
		break