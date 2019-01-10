all: tag
eval: model.txt tag
	python3 genData.py
	@./tag model.txt < testdata.txt > result.txt
	python3 evaluation.py
run: model.txt lexicon.txt tag
	@./tag model.txt lexicon.txt
model:
	@python3 makemodel.py
tag: main.cpp Array.hpp tag.hpp tag.cpp
	@g++ -std=c++17 main.cpp tag.hpp tag.cpp -o tag
genmodel: genmodel.cpp
	@g++ genmodel.cpp -o genmodel
clean:
	rm result.txt answer.txt