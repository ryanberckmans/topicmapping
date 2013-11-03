

typedef deque< pair<int, double> > did;
typedef map<int, deque< pair<int, double> > > mapiid;
typedef map<int, mapid> mapimapid;

#include "degree_block.h"
#include "utilities.h"
#include "log_table.h"


// TODO
// in class doc, wn_occurences_ is not efficient


class doc {

public:
	
	doc(){num_words_=0;};
	~doc(){};
	
	void set_from_string(string s, mapsi & word_number_all, mapis & number_word_all);
    void compute_thetas(int_matrix & word_partition, mapid & theta_doc, int black_module);

	UI num_words_;			// num_words_ in the document
	mapii wn_occurences_;	// word-number - occurences
    
};





class word_corpus {

public:

	word_corpus(double min_filter_p, double max_filter_p, 
                double p_value_p, int min_docs_p, 
                string partition_file_p, string corpus_file_p);
	~word_corpus(){};


    // setting the corpus
	void set_from_file(string filename);
	void write_corpus_file();

    // computing the network 
	void null_model(DI & links1, DI & links2, DD & weights);
    
    
    void write_beta_and_theta_files(deque<mapid> & doc_topic, map<int, mapid> & topic_word, \
                                    string theta_file, string beta_file);
    void write_short_beta_and_theta_files(deque<mapid> & doc_topic, map<int, mapid> & topic_word, \
                                          string theta_file, string beta_file, string beta_file_short, \
                                          mapid & pt);
    double dimap(int Nruns, \
                 mapid & pt, \
                 deque<mapid> & doc_topic_best, \
                 map<int, mapid> & topic_word_best, \
                 deque<mapii> & doc_assignments);
    
    double lda_model(deque<mapid> & doc_topic, 
                     map<int, mapid> & topic_word);

    
private:
    
    void write_partition(mapii & hard_mems);

    void update_doc_num_thetas(DI & doc_numbers,
                                bool just_one_noise_topic, 
                                map<int, mapid> & doc_number_thetas, 
                                int_matrix & word_partition, 
                                int & doc_number_multipart);

    void initial_ptopic(deque<mapid> & doc_topic, map<int, mapii> & word_topic, const mapii & hard_mems, const DI & doc_prevalent_topics);
    void get_betas(map<int, mapii> & word_topic, map<int, mapid> & topic_word, mapid & pt);
    double compute_likelihood(map<int, mapid> & topic_word, map<int, mapii> & word_topic, deque<mapid> & doc_topic);
    double likelihood_filter(map<int, mapii> & word_topic, deque<mapid> & doc_topic, \
                             map<int, mapid> & topic_word, mapid & pt, \
                             const mapii & hard_mems, double filtering_par, const DI & doc_prevalent_topics, \
                             deque<mapii> & doc_assignments);
    void get_rid_of_non_prevalent_topics(mapii & hard_mems, DI & doc_prevalent_topics);
    void get_prevalent_topics(DI & doc_prevalent_topics, mapii & hard_mems);

    double optimal_filtering(mapii & hard_mems, \
                             mapid & pt_best, \
                             deque<mapid> & doc_topic_best, \
                             map<int, mapid> & topic_word_best, \
                             deque<mapii> & doc_assignments_best, bool verbose);
    
    
    
    // basic data structures
    deque<string> word_strings_;		// word_strings_[wn] = "word_in_text"	
	DI word_occurrences_;            // word_occurrences_[wn] = occurences
	deque<doc> docs_;                // documents
    
    
    // parameters
    double max_filter_;
    double min_filter_;;
    double p_value_;
    int min_docs_;;
    string partition_file_;

    
    // lda functions
    double lda_inference(int doc_number);
    double compute_likelihood(int doc_number, DD & var_gamma);
    double run_em();
    double run_em_sparse();
    void set_class_words_to_zeros_map();
    void set_class_words_to_zeros();
    void initialize_lda_data(deque<mapid> & doc_topic, map<int, mapid> & topic_word);
    void optimize_alpha();
    void optimize_alpha_sparse();
    
    double compute_likelihood_sparse(int doc_number);
    double lda_inference_sparse(int doc_number);

    
    // lda data structures
    // (add explanation here)
    
    deque<DD> phis_ldav_;
    deque<DD> betas_ldav_;
    DD alphas_ldav_;
    deque<DD> class_word_ldav_;
    DD class_total_ldav_;
    deque<DD> gammas_ldav_;
    int num_topics_ldav_;
    
    // same thing as before but this is done
    // with sparse data structures
    deque<mapid> phis_ldav_map_;
    deque<mapid> betas_ldav_map_;
    deque<mapid> class_word_ldav_map_;
    mapid class_total_ldav_map_;
    deque<mapid> gammas_ldav_map_;



    
    
};




word_corpus::word_corpus(double min_filter_p, double max_filter_p, 
                         double p_value_p, int min_docs_p, 
                         string partition_file_p, string corpus_file_p) {
    
    
    max_filter_=min(0.51, max_filter_p);
    min_filter_=max(0., min_filter_p);
    p_value_=p_value_p;
    min_docs_=min_docs_p;
    partition_file_=partition_file_p;
    
    
    // checking parameters
    if(p_value_<=0) {
        cerr<<"p-value cannot be smaller than 0: "<<p_value_<<endl;
        exit(-1);
    }
    
    if(min_filter_>max_filter_) {
        cerr<<"ERROR: max_filter is smaller than min_filter"<<endl;
        exit(-1);
    }
    
    
    cout<<"*** corpus file: "<<corpus_file_p<<endl;
    cout<<"*** p-value: "<<p_value_<<endl;
    cout<<"*** min filter: "<<min_filter_<<endl;
    cout<<"*** max filter: "<<max_filter_<<endl;
    cout<<"*** min docs per topic: "<<min_docs_<<endl;
    if(partition_file_.size()>0)
        cout<<"***  partition file: "<<partition_file_<<endl;
    
    
    // setting docs_ and basic structures
    set_from_file(corpus_file_p);
}


#include "doc.cpp"
#include "word_corpus.cpp"
#include "sig_topics.cpp"
#include "filter_topics.cpp"
#include "lda_model.h"

