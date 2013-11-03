

#python ../testing-topic-modeling/synthetic_corpuses/synthetic_datasets/syn_corpus_alpha.py 0.001 0.5

rm a.out
g++ -O3 -funroll-loops -Wall  ./Sources/TopicMapping/docmap_tmp.cpp
time ./a.out -f run_exp/syn_corpus.txt -part infomap.part -minf 0.15 -maxf 0.17


#time ./a.out -f syn_corpus.txt
#time ./a.out -f syn_corpus.txt -part infomap.part -minf 0.0 -maxf 0.0001 -nos
#time ./a.out -fullout -f syn_corpus.txt -part infomap.part -minf 0.0 -maxf 0.0001 -nos



compare_models lda_gammas.txt lda_gammas.bak2 run_exp/syn_corpus.txt 