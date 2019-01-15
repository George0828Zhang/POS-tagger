MOD=model.txt
LEX=lexicon.txt
ORD=2

all: tag
eval: testdata.txt answer.txt ${MOD} ${LEX} tag	
	@./tag -m ${MOD} -l ${LEX} -ord ${ORD} -i testdata.txt -o result.txt
	python3 evaluation.py result.txt answer.txt
data: 
	@python3 genData.py
run: ${MOD} ${LEX} tag
	@./tag -m ${MOD} -l ${LEX} -ord ${ORD}
model:
	@python3 makemodel.py
tag: main.cpp Array.hpp Tree.hpp tag.hpp tag.cpp
	@g++ -std=c++17 main.cpp tag.cpp -o tag
clean:
	rm result.txt answer.txt