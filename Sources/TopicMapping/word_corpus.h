

typedef deque< pair<int, double> > did;
typedef map<int, deque< pair<int, double> > > mapiid;
typedef map<int, mapid> mapimapid;

// values below this threshold are 
// not included in betas
#define SPARSE_limit 0

#include "degree_block.cpp"
#include "utilities.cpp"
#include "log_table.cpp"


class doc {

public:
	
	doc(){num_words_=0;};
	~doc(){};
	
	void set_from_string(string s, mapsi & word_number_all, mapis & number_word_all);
    void set_from_string_given_wn(string buffer, mapsi & word_wn_all);


	UI num_words_;			// num_words_ (not unique) in the document
	deqii wn_occs_;         // word-number - occurences
    
};


class word_corpus {

public:

	word_corpus(double min_filter_p, double max_filter_p, 
                double p_value_p, int min_docs_p, 
                string partition_file_p);
	~word_corpus(){};


    // setting the corpus
	void set_from_file_first_time(string filename, string out_dir);
    void set_from_file(string filename, string wn_file);
	void write_corpus_file(string out_dir);


    // getters
    int get_num_terms() { return int(word_occurrences_.size()); };

    // computing the network
    void null_model(string parall_str, string out_dir);
    void write_theta_file(deque<mapid> & doc_topic,\
                          map<int, mapid> & topic_word, \
                          string theta_file);
    void write_short_beta_and_theta_files(deque<mapid> & doc_topic,\
                                          map<int, mapid> & topic_word, \
                                          string theta_file, \
                                          string beta_file, \
                                          string beta_file_short, \
                                          mapid & pt);
    double dimap(int Nruns, \
                 double step,
                 bool print_sig_words,
                 mapid & pt, \
                 deque<mapid> & doc_topic_best, \
                 map<int, mapid> & topic_word_best,\
                 string out_dir);
    
    double lda_model(map<int, mapid> & topic_word, \
                     double alpha_init,\
                     bool skip_alpha_opt, \
                     bool infer_flag, int print_lag,\
                     string alpha_file, string out_dir);
    
private:
    
    
    // ======================= member private methods ======================

    
    // null model
    void dotpr_similarity_of_connected_words(map<pair<int, int> , int> & cooc);
    void dotpr_similarity_of_connected_words_parallel(map<pair<int, int> , int> & cooc, \
                                                      int par_a, int par_b, int max_ab);
    void null_model(DI & links1, DI & links2, DD & weights, \
                    int par_a, int par_b, int max_ab, bool print_sig_words, string out_dir);
        
    // likelihood filtering
    void write_partition(mapii & hard_mems, string out_dir);
    void initial_ptopic(deque<mapid> & doc_topic, \
                        map<int, mapii> & word_topic, \
                        const mapii & hard_mems, \
                        const DI & doc_prevalent_topics);
    void get_betas(map<int, mapii> & word_topic, \
                   map<int, mapid> & topic_word, mapid & pt);
    double compute_likelihood(map<int, mapid> & topic_word, \
                              map<int, mapii> & word_topic, deque<mapid> & doc_topic);
    double likelihood_filter(map<int, mapii> & word_topic, deque<mapid> & doc_topic, \
                             map<int, mapid> & topic_word, mapid & pt, \
                             const mapii & hard_mems, \
                             double filtering_par, const DI & doc_prevalent_topics,\
                             deque<mapii> & doc_assignments);
    void get_rid_of_non_prevalent_topics(mapii & hard_mems, DI & doc_prevalent_topics);
    void get_prevalent_topics(DI & doc_prevalent_topics, mapii & hard_mems);
    double optimal_filtering(mapii & hard_mems, \
                             mapid & pt_best, \
                             deque<mapid> & doc_topic_best, \
                             map<int, mapid> & topic_word_best, \
                             bool verbose, double step, string out_dir);
    
    
    // lda functions
    double lda_inference(int doc_number);
    double compute_likelihood(int doc_number, DD & var_gamma);
    double E_step(bool verbose, string out_dir, int iter);
    double run_em();
    double run_em_sparse(bool skip_alpha_opt, bool infer_flag, \
                         int print_lag, string out_dir);
    void set_class_words_to_zeros_map();
    void set_class_words_to_zeros();
    void initialize_lda_data(map<int, mapid> & topic_word,\
                             double alphas_init, string alpha_file);
    void compute_non_sparse_gammas(deque<DD> & gammas_ldav);
    void optimize_alpha_sparse();
    double compute_likelihood_sparse(int doc_number, \
                                     mapid & var_gamma, mapid & digamma_gam, 
                                     const double & digsum, 
                                     const double & likelihood_alpha,
                                     const double & likelihood_const);
                                     
    double compute_likelihood_sparse_constants(int doc_number, \
                                               mapid & var_gamma, \
                                               double & digsum, \
                                               const double & sum_alphas);
    double compute_likelihood_alpha_terms(double & sum_alphas);
    double lda_inference_sparse(int doc_number, const double & sum_alphas,\
                                const double & likelihood_alpha,\
                                bool print_assignments, ostream & asgout);
    
    void print_lda_results(string out_dir, int iter);
    void print_betas(string outfile);
    
    
    // ======================= member variables ======================
    
    // basic data structures
    deque<string> word_strings_;        // word_strings_[wn] = "word_in_text"	
	DI word_occurrences_;               // word_occurrences_[wn] = occurences
	deque<doc> docs_;                   // documents
    
    
    // parameters
    double max_filter_;
    double min_filter_;;
    double p_value_;
    int min_docs_;;
    string partition_file_;

    
    // lda data structures
    // (add explanation here)
    
    /* this are for the non-sparse implementation 
        and are not in use anymore
    deque<DD> phis_ldav_;
    deque<DD> betas_ldav_;
    deque<DD> class_word_ldav_;
    DD class_total_ldav_;
    deque<DD> gammas_ldav_;
    */
    
    
    int num_topics_ldav_;
    DD alphas_ldav_;
    
    // same thing as before but this is done
    // with sparse data structures
    deque<deqid> phis_ldav_map_;
    deque<mapid> betas_ldav_map_;
    deque<mapid> class_word_ldav_map_;
    mapid class_total_ldav_map_;
    deque<mapid> gammas_ldav_map_;
    
    
};




word_corpus::word_corpus(double min_filter_p, double max_filter_p, 
                         double p_value_p, int min_docs_p, 
                         string partition_file_p) {
    
    
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
    
    
    cout<<"*** p-value: "<<p_value_<<endl;
    cout<<"*** min filter: "<<min_filter_<<endl;
    cout<<"*** max filter: "<<max_filter_<<endl;
    cout<<"*** min docs per topic: "<<min_docs_<<endl;
    if(partition_file_.size()>0)
        cout<<"***  partition file: "<<partition_file_<<endl;
    
    
}


#include "doc.cpp"
#include "word_corpus.cpp"
#include "sig_topics.cpp"
#include "filter_topics.cpp"
#include "lda_model.h"

