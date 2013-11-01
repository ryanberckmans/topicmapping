

#python ../testing-topic-modeling/synthetic_corpuses/synthetic_datasets/syn_corpus_alpha.py 0.001 0.5

rm a.out

g++ -O3 -funroll-loops -Wall  ./Sources/TopicMapping/docmap_tmp.cpp
time ./a.out -f syn_corpus.txt

#time ./a.out -f syn_corpus.txt -part infomap.part -minf 0.0 -maxf 0.0001 -nos
#time ./a.out -fullout -f syn_corpus.txt -part infomap.part -minf 0.0 -maxf 0.0001 -nos







