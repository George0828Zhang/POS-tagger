#!/usr/bin/env python3
import sys

resname = sys.argv[1]
ansname = sys.argv[2] 

res = open(resname, "r")
ans = open(ansname, "r")
# res = open("result.txt", "r")
# ans = open("answer.txt", "r")

blkcorrect = 0
blktotal = 0

sent_res = res.readline()

while sent_res != '':
	sent_res = sent_res.strip().split()
	sent_ans = ans.readline().strip().split()
	if len(sent_res) != len(sent_ans):
		print(type(sent_res))
		print(type(sent_ans))
	for i, blk_res in enumerate(sent_res):
		blk_ans = sent_ans[i]
		if blk_res == blk_ans:
			blkcorrect += 1
		blktotal += 1
	sent_res = res.readline()
	

print("{} words out of {} are correct. Accuracy: {:.2f}%".format(blkcorrect, blktotal, blkcorrect/blktotal*100))