


# include "../standard_package/standard_include.cpp"

// global variables (really minor things)
bool global_verbose=false;
const string INFOMAP_PATH="___INSERT_INFOMAP_PATH_HERE___";


#include "../SingleSliceInfomap/infomap_single_function.cpp"
#include "sig_topics.h"
#include "set_parameters.h"


int main(int argc, char * argv[]) {
    
    
    parameters P;
	set_parameters_for_docmap(P, argc, argv);
    
    string corpus_file=P.string_ps.at("-f");
    double p_value=P.double_ps.at("-p");
    double min_filter=P.double_ps.at("-minf");
    double max_filter=P.double_ps.at("-maxf");
    int min_docs=P.int_ps.at("-t");
    int Nruns=P.int_ps.at("-r");
    string partition_file=P.string_ps.at("-part");
    
    // checking parameters
    max_filter=min(0.51, max_filter);
    min_filter=max(0., min_filter);
    if(p_value<=0) {
        cerr<<"p-value cannot be smaller than 0: "<<p_value<<endl;
        exit(-1);
    }
    if(min_filter>max_filter) {
        cerr<<"ERROR: max_filter is smaller than min_filter"<<endl;
        exit(-1);
    }
    
    cout<<"*** corpus file: "<<corpus_file<<endl;
    cout<<"*** p-value: "<<p_value<<endl;
    cout<<"*** min filter: "<<min_filter<<endl;
    cout<<"*** max filter: "<<max_filter<<endl;
    cout<<"*** min docs per topic: "<<min_docs<<endl;
    cout<<"*** convergence parameter: "<<P.double_ps.at("-conv")<<endl;
    if (partition_file.size()>0)
        cout<<"***  partition file: "<<partition_file<<endl;
    
    // setting corpus from file
    word_corpus C;
	C.set_from_file(corpus_file);

    // doc_topic[doc] is p(t|doc)
    deque<mapid> doc_topic_best;
    // topic_word_best[topic][word] is p(w|t)
    map<int, mapid> topic_word_best;
    // doc_assignments[doc][wn] is the topic to which the word has been assigned
    deque<mapii> doc_assignments;
    // pt[topic] is p(t)
    mapid pt;
    
    double eff_ntopics=C.dimap(min_filter, max_filter, min_docs, \
            Nruns, p_value, partition_file, \
            pt,
            doc_topic_best, \
            topic_word_best, \
            doc_assignments, 0, P.double_ps.at("-conv"));
    
    cout<<"Effective number topics: "<<eff_ntopics<<endl;
    
    // writing p(t|doc) and p(w|t) in files thetas.txt and betas.txt
    C.write_short_beta_and_theta_files(doc_topic_best, topic_word_best, \
                                       "doc_topics.txt", "topic_words.txt", "topic_summary.txt", pt);
        C.write_beta_and_theta_files(doc_topic_best, topic_word_best, "thetas.txt", "betas.txt");
    
    
    return 0;
}


