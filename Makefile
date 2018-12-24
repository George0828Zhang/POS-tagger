all: tag
run: model.txt tag
	./tag model.txt
model.txt: genmodel
	python3 makemodel.py
	# ./genmodel
tag: tag.cpp
	g++ tag.cpp -o tag
genmodel: genmodel.cpp
	g++ genmodel.cpp -o genmodel