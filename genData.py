#!/usr/bin/env python3
untagged_sentences = []

import nltk
from nltk.corpus import brown as corpus
untagged_sentences += corpus.sents()

data = open("testdata.txt", "wt")
res = open("answer.txt", "wt")

amount = 500
tru_amount = len(untagged_sentences)
print("[Data] Extracted {} out of {} ({:.2f}%)".format(amount, tru_amount, amount/tru_amount*100))

for sentences in untagged_sentences:
	tagged = nltk.pos_tag(sentences)
	for word,tag in tagged:
		data.write("{} ".format(word))
		res.write("{}/{} ".format(word, tag))
	data.write("\n")
	res.write("\n")
	amount -= 1
	if amount == 0:
		break


# #!/usr/bin/env python3
# tagged_sentences = []

# from nltk.corpus import brown as corpus
# tagged_sentences += corpus.tagged_sents()

# data = open("testdata.txt", "wt")
# res = open("answer.txt", "wt")

# amount = 3000

# for sentences in tagged_sentences:
# 	for (word, tag) in sentences:
# 		data.write("{} ".format(word))
# 		res.write("{}/{} ".format(word, tag))
# 	data.write("\n")
# 	res.write("\n")
# 	amount -= 1
# 	if amount == 0:
# 		break