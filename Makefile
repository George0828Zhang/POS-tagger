all: tag
run: model.txt tag
	./tag model.txt
model.txt:
	python3 makemodel.py
	# ./genmodel
tag: tag.cpp Array.h
	g++ -std=c++17 tag.cpp -o tag
genmodel: genmodel.cpp
	g++ genmodel.cpp -o genmodel