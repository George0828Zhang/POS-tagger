# Statistical Part-of-Speech Tagger
DSP2018 Final Project, B05902064 張致強

### Features
- can correctly tag most sentences with simple structure.
- can handle out-of-volcabulary words as well as non-sense sentences.
- supports order 2 and order 3 models.
- acceptible speed & accuracy

### Usage
##### Compile
```bash=
make
```
##### Generate Model
```bash=
make model
```
##### Run Interactive Mode
```bash=
make run 
```
##### Run Evaluation
```bash=
make eval
```

### Implementation Details
- Underlying model: HMM
    * The program adopts Hidden-Markov-Model, with tags as states and word equivalence classes as observations. Transitions and emission probabilities are estimated from statistics of the corpus. 
    * Both first and second order models are implemented.
- Equivalence classes
    * In order to increase speed as well as generalization, equivalence classes are used to represent words instead of its text. 
    * For any 2 words, we group them as equivalent if they share the same possible tags in the corpus. 
    * For example, the corpus used here has 12408 different words, but are generalized into 152 classes.

- Probability smoothing
    * Since the model is statistics based, smoothing is used to combat data sparsity
    * Bi-grams are smoothed using unigrams, and trigrams with bigrams and unigrams. Emission is smooth by $P(Eqv)$.
    * All smoothing used deleted interpolation.

- Handling unknown words
    * Another way to combat data sparsity is handling of OOV. The technique used here is suffix analysis. 
    * Suffix of words with less than 20 appearences in corpus are analyzed. Suffixes of length up to 10 are considered.
    * For an OOV word with longest matched suffix $S_i$, instead of estimating with $P(Eqv|tag)$, $P(S_i|tag)$ is used.
    * $P(S_i|tag)=\sum_j P(S_i|Eqv_j)P(Eqv_j|tag)$

- Other small details
    * Use of lexicon tree and suffix tree to significantly increase query speed
    * Beam search is used on second order model to increase speed
    * Self implemented tokenizer can handle basic word-word and word-punctuation separations. 

### Evaluation
The corpus used for estimation is `treebank` corpus, available from nltk. The tagset used is `Universal` tagset. Brown corpus is used as evaluation. The entire corpus (`57,340` sentences, `1,161,192` word tokens) is fed into the program and the result is compared with its universal tags. 
The result is: `1,029,429` words out of `1,161,192` are correct. Accuracy: `88.65%`. The program took `71.7372s` to tag the whole corpus. (The nltk pos_tagger acheives `91.84%` using averaged perceptron tagger.)

### Problems And Future Work
As seen in last section, the accuracy is below 90%, which is acceptible, but not satisfying. The reference paper acheived 96.7% accuracy on Penn Treebank corpus (Wall Street Journal). There are features and improvements to be adopted from this paper in order to improve accuracy. On the other hand, the speed is satisfiable, with near instantaneous response for interactive mode. Another problem is generalization. I tried to implement Baum-Welch algorithm for re-estimation and adaption, but the trainning speed is unacceptibly slow and the outcome is not as good as expected. In the future, I shall look into more papers on unsupervised approaches to further improve the portability. 

### References
- Thorsten Brants. 2000. TnT --- A Statistical Part-ofSpeech Tagger. In *Proceedings of the Sixth Applied Natural Language Processing Conference ANLP-2000, April 29-May 3, 2000, Seattle, WA*.

- Julian Kupiec. 1992. Robust part-of-speech tagging using a hidden Markov model.

- Natural Language Toolkit https://www.nltk.org/
- Brown corpus
- Penn Treebank Corpus (Sample)
