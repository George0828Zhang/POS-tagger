all: tag
eval: model.txt tag
	python3 genData.py
	@./tag model.txt < testdata.txt > result.txt
	python3 evaluation.py
run: model.txt tag
	@./tag model.txt
model:
	@python3 makemodel.py
tag: tag.cpp Array.h
	@g++ -std=c++17 tag.cpp -o tag
genmodel: genmodel.cpp
	@g++ genmodel.cpp -o genmodel
clean:
	rm result.txt answer.txt