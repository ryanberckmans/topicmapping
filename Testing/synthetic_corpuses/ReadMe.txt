
Using this package you can create two synthetic datasets. The language datasets and the fully synthetic datasets.


*************************************
1. Language datasets

From a terminal, type:

python synthetic_datasets/syn_from_languages.py 1000 0.2

*** Input:
 1000 is the number of documents
 0.2 is the total probability of the two biggest languages
     in the paper we used 0.2 (all languages are equiprobable) and 0.6

*** Output:
 The program creates 4 files:
 1) syn_corpus.txt: each line is a document, each number is a word
 2) syn_corpus.corpus: same information as before, but it is written in the same corpus format as the one used in LDA: http://www.cs.princeton.edu/~blei/lda-c/index.html
 3) syn_thetas.txt each line is a document: each number is the probability of using a topic: p(topic|doc). The first number is for the first topic, and so on. 
 4) syn_betas.txt each line is a topic: each number is the log-probability of using a word: log p(word|topic). First number is for the  word '0', second for word '1' and so on. 

You can also change the words per document from synthetic_datasets/syn_from_languages.py (at the bottom of the file).

*************************************
2. Synthetic datasets


From a terminal, type:

python synthetic_datasets/syn_corpus_alpha.py 1e-3 0.1 0.2


*** Arguments:
 1e-3 is alpha
 0.1 is the fraction of generic words
 0.2 is the probability of the top 20% topic

*** Output:
 Same as above

The other parameters (number of documents, document length etc.) can be tuned from synthetic_datasets/syn_corpus_alpha.py (at the bottom of the file).


