all: tag
run: model.txt tag
	./tag model.txt
model.txt:
	python3 makemodel.py
tag: tag.cpp
	g++ tag.cpp -o tag