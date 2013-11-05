
# include "../standard_package/standard_include.cpp"

// global variables (really minor things)
bool global_verbose=false;
const string INFOMAP_PATH="/Users/andrealancichinetti/Desktop/work/2013/LDA_eval/topicmapping/bin/Infomap";
//const string INFOMAP_PATH="/home/staff/andrea/Desktop/ESA/esa_bangkok/TreatCorpusPhoenix/topicmapping/Sources/Infomap-0.11.5/Infomap";


#include "../SingleSliceInfomap/infomap_single_function.cpp"
#include "word_corpus.h"
#include "set_parameters.h"


int main(int argc, char * argv[]) {
        
    parameters P;
	set_parameters_for_docmap(P, argc, argv);
    
    
    // setting corpus from file
    word_corpus C(P.double_ps.at("-minf"),
                  P.double_ps.at("-maxf"),
                  P.double_ps.at("-p"), 
                  P.int_ps.at("-t"),  
                  P.string_ps.at("-part"), 
                  P.string_ps.at("-f") );
    

    // doc_topic[doc] is p(t|doc)
    deque<mapid> doc_topic_best;
    // topic_word_best[topic][word] is p(w|t)
    map<int, mapid> topic_word_best;
    // doc_assignments[doc][wn] is the topic to which the word has been assigned
    deque<mapii> doc_assignments;
    // pt[topic] is p(t)
    mapid pt;
    
    
    double eff_ntopics=C.dimap(P.int_ps.at("-r"), \
                               P.double_ps.at("-step"),
                               pt,doc_topic_best, \
                               topic_word_best,
                               doc_assignments);
    
    cout<<"Effective number topics: "<<eff_ntopics<<endl;
    
    // writing p(t|doc) and p(w|t) in files thetas.txt and betas.txt
    C.write_short_beta_and_theta_files(doc_topic_best, topic_word_best, \
                                       "doc_topics.txt", "topic_words.txt", "topic_summary.txt", pt);
    C.write_beta_and_theta_files(doc_topic_best, topic_word_best, "thetas.txt", "betas.txt");
    
    // optimizing LDA
    C.lda_model(doc_topic_best, topic_word_best, P.double_ps.at("-alpha"));
    
    
    

    return 0;
}


