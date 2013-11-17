
# include "../standard_package/standard_include.cpp"

// global variables (really minor things)
bool global_verbose=false;
const string INFOMAP_PATH="/Users/andrealancichinetti/Desktop/work/2013/LDA_eval/topicmapping/bin/Infomap";

#include "../SingleSliceInfomap/infomap_single_function.cpp"
#include "word_corpus.h"
#include "set_parameters.h"


int main(int argc, char * argv[]) {
        
    parameters P;
	set_parameters_for_docmap(P, argc, argv);
    
    
    /*
     setting corpus from file
     lots of parameter are thrown here because
     they might be used in multiple functions
     (mostly likelihood filtering though)
    */
    word_corpus C(P.double_ps.at("-minf"),
                  P.double_ps.at("-maxf"),
                  P.double_ps.at("-p"), 
                  P.int_ps.at("-t"),  
                  P.string_ps.at("-part"));
    
    // setting docs_ and basic structures
    if(P.string_ps.at("-word_wn").size()==0) {
        cout<<"*** corpus file: "<<P.string_ps.at("-f")<<endl;
        C.set_from_file(P.string_ps.at("-f"));
    } else {
        C.set_from_file(P.string_ps.at("-f"), P.string_ps.at("-word_wn"));        
    }
    cout<<"corpus was set"<<endl;
    
    if(P.bool_ps.at("-corpus")) {
        cout<<"writing corpus file CORPUS.corpus"<<endl;
        C.write_corpus_file();
    }
    
    if(P.string_ps.at("-parall").size()!=0) {
        cout<<"running null model only"<<endl;
        C.null_model(P.string_ps.at("-parall"));
        return 0;
    }

    // topic_word[topic][word] is p(w|t)
    map<int, mapid> topic_word;

    if(P.string_ps.at("-model").size()==0) {
        
        /*
         getting initial guess from Infomap
         this is the **novel** part of the algorithm
         the rest is LDA (variational inference)
         with asymmetric alpha priors
         and sparse data structures
        */
        
        //doc_topic[doc] is p(t|doc)
        deque<mapid> doc_topic_best;
        // pt[topic] is p(t)
        mapid pt;
        cout<<"getting a model clustering words"<<endl;
        double eff_ntopics=C.dimap(P.int_ps.at("-r"), \
                                   P.double_ps.at("-step"),
                                   P.bool_ps.at("-write_net"),
                                   pt, doc_topic_best, \
                                   topic_word);

        cout<<"effective number topics: "<<eff_ntopics<<endl;
        
        // writing p(t|doc) and p(w|t) in files thetas.txt and betas.txt
        C.write_short_beta_and_theta_files(doc_topic_best, topic_word, \
                                           "plsa_thetas_sparse.txt", \
                                           "plsa_betas_sparse.txt", \
                                           "plsa_summary.txt", pt);    
        C.write_theta_file(doc_topic_best, topic_word, "plsa_thetas.txt");
    
    } else {
        // skipping all the previous part because we are loading a model from file
        read_topic_model_from_file(P.string_ps.at("-model"), topic_word);
    }
    
    // optimizing LDA
    C.lda_model(topic_word,
                P.double_ps.at("-alpha"), \
                P.bool_ps.at("-skip_opt_al"), \
                P.bool_ps.at("-infer"),\
                P.int_ps.at("-lag"), \
                P.string_ps.at("-alpha_file"));
    
    cout<<"--done!--"<<endl;

    return 0;
}

