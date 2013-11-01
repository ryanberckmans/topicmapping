

typedef deque< pair<int, double> > did;
typedef map<int, deque< pair<int, double> > > mapiid;
typedef map<int, mapid> mapimapid;

#include "degree_block.h"
#include "utilities.h"
#include "log_table.h"

// TODO
// in class doc, wn_occurences is not efficient




class doc {

public:
	
	doc(){num_words=0;};
	~doc(){};
	
	void set_from_string(string s, mapsi & word_number_global, mapis & number_word_global);
    void compute_thetas(int_matrix & word_partition, mapid & theta_doc, int black_module);

	UI num_words;			// num_words in the document
	mapii wn_occurences;	// word-number - occurences
    
};





class word_corpus {

public:
	
	word_corpus(){};
	~word_corpus(){};


    // setting the corpus
	void set_from_file(string filename);
	void write_corpus_file();

    // computing the network 
	void null_model(double p_value, \
                    DI & links1, DI & links2, DD & weights, \
                    bool use_dotproduct, bool verbose);
    
    
    void write_beta_and_theta_files(deque<mapid> & doc_topic, map<int, mapid> & topic_word, \
                                    string theta_file, string beta_file);
    void write_short_beta_and_theta_files(deque<mapid> & doc_topic, map<int, mapid> & topic_word, \
                                          string theta_file, string beta_file, string beta_file_short, \
                                          mapid & pt);


    double dimap(double min_filter, double max_filter, int min_docs, \
               int Nruns, double p_value, string partition_file, \
               mapid & pt, \
               deque<mapid> & doc_topic_best, \
               map<int, mapid> & topic_word_best, \
               deque<mapii> & doc_assignments,\
               int level, const double & convergence_precision);
    
    string get_topic_title(const mapid & topic_distr);
    
    
    
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
    void get_prevalent_topics(DI & doc_prevalent_topics, mapii & hard_mems, int min_docs);

    double optimal_filtering(mapii & hard_mems, double min_filter, double max_filter, int min_docs, \
                             mapid & pt_best, \
                             deque<mapid> & doc_topic_best, \
                             map<int, mapid> & topic_word_best, \
                             deque<mapii> & doc_assignments_best, bool verbose);
    
    
    
    // basic data structures
    deque<string> word_strings;		// word_strings[wn] = "word_in_text"	
	DI word_occurrences;            // word_occurrences[wn] = occurences
	deque<doc> docs;                // documents
    
    
    // lda functions
    double lda_inference(int doc_number);
    double compute_likelihood(int doc_number, DD & var_gamma);
    void gibbs_sampling(deque<mapii> & doc_assignments);    
    double run_em();
    void set_class_words_to_zeros();
    void initialize_lda_data(deque<mapid> & doc_topic, map<int, mapid> & topic_word);
    void optimize_alpha();
    
    // lda data structures
    // (add explanation here)
    deque<DD> phis_ldav;
    deque<DD> betas_ldav;
    DD alphas_ldav;
    deque<DD> class_word_ldav;
    DD class_total_ldav;
    deque<DD> gammas_ldav;
    int num_topics_ldav;
    
};

#include "doc.cpp"
#include "sig_topics.cpp"
#include "lda_model.h"
