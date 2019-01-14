MOD=model.txt
LEX=lexicon.txt
ORD=2

all: tag
eval: testdata.txt answer.txt ${MOD} ${LEX} tag
	@# python3 genData.py
	@# ./tag model.txt lexicon.txt < testdata.txt > result.txt
	@#./tag -m model.txt -l lexicon.txt < testdata.txt > result.txt
	@./tag -m ${MOD} -l ${LEX} -ord ${ORD} -i testdata.txt -o result.txt
	python3 evaluation.py result.txt answer.txt
run: ${MOD} ${LEX} tag
	@./tag -m ${MOD} -l ${LEX} -ord ${ORD}
model:
	@python3 makemodel.py
tag: main.cpp Array.hpp Tree.hpp tag.hpp tag.cpp
	@g++ -std=c++17 main.cpp tag.hpp tag.cpp -o tag
# train: model.txt lexicon.txt bwestimate
# 	cp model.txt model_og.txt
# 	@./bwestimate model.txt lexicon.txt traindata.txt 30
# bwestimate: bwestimate.cpp Array.hpp tag.hpp tag.cpp
# 	@g++ -std=c++17 bwestimate.cpp tag.hpp tag.cpp -o bwestimate
# genmodel: genmodel.cpp
# 	@g++ genmodel.cpp -o genmodel
clean:
	rm result.txt answer.txt