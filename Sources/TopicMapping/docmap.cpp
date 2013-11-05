
# include "../standard_package/standard_include.cpp"

// global variables (really minor things)
bool global_verbose=false;
const string INFOMAP_PATH="___INSERT_INFOMAP_PATH_HERE___";

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
    
    

    // topic_word[topic][word] is p(w|t)
    map<int, mapid> topic_word;

    if(P.string_ps.at("-model").size()==0) {
        
        /*
         getting initial guess from Infomap
         this is the **novel** part of the algorithm
         the rest is LDA (variational inference)
         with asymmetric alpha priors
        */
        
        // doc_topic[doc] is p(t|doc)
        deque<mapid> doc_topic_best;
        // pt[topic] is p(t)
        mapid pt;
        
        
        double eff_ntopics=C.dimap(P.int_ps.at("-r"), \
                                   P.double_ps.at("-step"),
                                   pt,doc_topic_best, \
                                   topic_word);
        
        cout<<"Effective number topics: "<<eff_ntopics<<endl;
        
        // TODO: I should probably remove this (intermediate step)
        // writing p(t|doc) and p(w|t) in files thetas.txt and betas.txt
        C.write_short_beta_and_theta_files(doc_topic_best, topic_word, \
                                           "doc_topics.txt", "topic_words.txt", "topic_summary.txt", pt);    
        C.write_beta_and_theta_files(doc_topic_best, topic_word, "thetas.txt", "betas.txt");
        //
    } else {
        
        read_topic_model_from_file(topic_word, P.string_ps.at("-model"));
    }
    
    
    // optimizing LDA
    C.lda_model(topic_word, P.double_ps.at("-alpha"), P.bool_ps.at("-skip_opt_al"));
    
    
    

    return 0;
}


