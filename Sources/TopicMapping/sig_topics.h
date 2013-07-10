

typedef deque< pair<int, double> > did;
typedef map<int, deque< pair<int, double> > > mapiid;
typedef map<int, mapid> mapimapid;

#include "degree_block.h"
#include "utilities.h"
#include "log_table.h"

void update_map_wordnum(mapsi & word_number, mapis & number_word, string s) {
	
	if(word_number.find(s)==word_number.end()) {
		int num=word_number.size();
		word_number[s]=num;
		number_word[num]=s;
	}
}

class doc {

public:
	
	doc(){num_words=0;};
	~doc(){};
	
	void set_from_string(string s, mapsi & word_number_global, mapis & number_word_global);
    void compute_thetas(int_matrix & word_partition, mapid & theta_doc, int black_module);

	UI num_words;			// num_words in the document
	mapii wn_occurences;	// word-number - occurences
    
};


void doc::set_from_string(string s, mapsi & word_number_global, mapis & number_word_global) {
	
	stringstream ss(s);
    deque<string> tokens;
	string buf;
    while (ss >> buf)
        tokens.push_back(buf);
	wn_occurences.clear();
	RANGE_loop(i, tokens) {
		update_map_wordnum(word_number_global, number_word_global, tokens[i]);
		int_histogram(word_number_global[tokens[i]], wn_occurences);		
	}
	num_words=tokens.size();
	
}


void doc::compute_thetas(int_matrix & word_partition, mapid & theta_doc, int black_module) {
    
    
    theta_doc.clear();
    SI covered_words;
    
    RANGE_loop(i, word_partition) {
        int module=word_partition[i][0];
        word_partition[i].pop_front();

        RANGE_loop(j, word_partition[i]) {
        
            if(wn_occurences.count(word_partition[i][j])==0) {
                cerr<<"error in compute_thetas"<<endl;
                exit(-1);
            }
            int_histogram(module, theta_doc, double(wn_occurences[word_partition[i][j]])/num_words);
            covered_words.insert(word_partition[i][j]);
        }
    }

    // words which are not covered 
    IT_loop(mapii, itm, wn_occurences) if(covered_words.find(itm->first)==covered_words.end()) {
        int_histogram(black_module, theta_doc, double(itm->second)/num_words);
    }
    
    // check that the norm is one
    double norm=0.;
    IT_loop(mapid, itm, theta_doc) norm+=itm->second;
    if(fabs(norm-1)>1e-4) {
        cerr<<"error in norm"<<endl;
        exit(-1);
    }

}


class word_corpus {

public:
	
	word_corpus(){total_words=0;};
	~word_corpus(){};

	void set_from_file(string filename);
	void write_corpus_file();
    void write_theta_file_from_multipart(string multipart_file, string doc_number_file);
	void null_model(double p_value, bool compute_only_first_order,\
                    DI & links1, DI & links2, DD & weights, bool use_dotproduct, bool verbose);
    
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
    
    void fix_data_structs();
    string get_topic_title(const mapid & topic_distr);
    
    // this is only used from files (it could be cleaner)
    // therefore it makes sense to use it
    // only from oC
	mapsi word_number_global;		// word - wn
	mapis number_word_global;		// wn - word
	
	mapii wn_occurences_global;                 // wn - occurences
	deque<doc> docs;                            // documents
	int total_words;
    
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
};


#include "sig_topics.cpp"


void check_distr_equal(mapid & distr1, mapii & distr2) {
    
    mapid distr2norm;
    IT_loop(mapii, itm, distr2) distr2norm[itm->first]=double(itm->second);
    normalize_mapid(distr2norm);
    
    cout<<"check_distr_equal"<<endl;
    IT_loop(mapid, itm, distr1) {
        //cout<<itm->first<<" "<<itm->second<<" "<<distr2norm.at(itm->first)<<endl;
        if(fabs(itm->second-distr2norm.at(itm->first))>1e-5) {
            cerr<<"error in check_distr_equal"<<endl;
            exit(-1);
        }
    }
    cout<<"done"<<endl;
    
}


void second_level(double min_filter, double max_filter, int min_docs, \
                  int Nruns, double p_value, word_corpus & oC, \
                  deque<mapii> & doc_assignments_old, \
                  const deque<mapid> & doc_topic_old, \
                  const map<int, mapid> & topic_word_old, const double & convergence_precision, 
                  int subtopic_docs, bool fullout, mapid & pt_old) {
    
    /* for each old topic, we are going to have a list of new topics
    for each doc, there must be a list of ids of that doc in the new corpuses */
    
    
    // for each doc, doc_new_topics[doc] is p(t|doc) with new labels
    deque<mapid> doc_new_topics;
    RANGE_loop(i, doc_topic_old) {
        mapid mapid_;
        doc_new_topics.push_back(mapid_);
    }
    // for each old topic, how this relates to the new topics
    deque<deque<int> > old_topic_new_topic_labels;
    // for each new topic, betas
    map<int, mapid> new_topics;
    
    
    
    //=============== One new corpus per topic ==================
    // for each topic, what are the original labels of the docs
    map<int, mapii> topic_old_docs;
    // defining a new corpus for each topic
    map<int, word_corpus> old_topic_corpus;
    // setting corpus from topics
    set_corpuses_from_topics(doc_assignments_old, old_topic_corpus, topic_old_docs, oC);
    //==================================================

    /*
     each sub topic has a label starting from 0
     topic 0 has subtopics 0, 1, 2, x
     topic 1 has subtopics x+1, x+2 ...
    */
    int topic_label=0;
    
    ofstream pout("subtopic_summary.txt");
    pout<<"#docs "<<oC.docs.size()<<endl;

    
    for(map<int, word_corpus>::iterator top_c_itm = old_topic_corpus.begin(); top_c_itm!=old_topic_corpus.end(); top_c_itm++) {
        
        cout<<"================ finding subtopics for topic =========== "<<top_c_itm->first<<" pt: "<<pt_old.at(top_c_itm->first)<<endl;
        cout<<oC.get_topic_title(topic_word_old.at(top_c_itm->first))<<endl;
        // subtopics summary 
        pout<<"****topic: "<<top_c_itm->first<<" pt: "<<pt_old.at(top_c_itm->first)<<endl;
        pout<<oC.get_topic_title(topic_word_old.at(top_c_itm->first))<<endl;
        word_corpus & C =top_c_itm->second;
        C.fix_data_structs();
        
        // ======================= running dimap ==================
        deque<mapid> doc_topic_sub;
        map<int, mapid> topic_word_sub;
        deque<mapii> doc_assignments_subtopics;
        mapid pt_sub;
        double eff_ntopics=C.dimap(min_filter, max_filter, min_docs, \
                                   Nruns, p_value, "", pt_sub, \
                                   doc_topic_sub, \
                                   topic_word_sub, \
                                   doc_assignments_subtopics, 1, convergence_precision);
        
        // ======================= running dimap ==================

        mapid doc_pr;
        RANGE_loop(i, C.docs) doc_pr[i]=double(C.docs[i].num_words);
        normalize_mapid(doc_pr);
        double eff_docs= compute_eff_num_topics(doc_pr);
        double eff_words = compute_eff_num_topics(topic_word_old.at(top_c_itm->first));
        
        //cout<<"Number of sub-topics: "<<topic_word_sub.size()<<" vs "<<doc_pr.size()<<endl;
        cout<<"Effective #topics: "<<eff_ntopics<<" #docs "<<eff_docs<<" #words "<<eff_words<<endl;
        cout<<"#docs per sub-topic: "<<eff_docs/eff_ntopics<<" #words per sub-topic: "<<eff_words/eff_ntopics<<endl;
        /* if there are actually no sub-topics ====================
         to say that there are subtopics, 
         the average subtopic needs to 
         cover at least subtopic_docs docs.*/
        //if(eff_ntopics<2. or eff_docs/eff_ntopics<subtopic_docs or eff_words/eff_ntopics<subtopic_words) {
        if(eff_ntopics<1.1 or eff_docs/eff_ntopics<subtopic_docs) {
            mapid mapid_;
            mapid_[0]=1.;
            RANGE_loop(i, doc_topic_sub) doc_topic_sub[i]=mapid_;
            topic_word_sub.clear();
            // the old topic
            topic_word_sub[0]=topic_word_old.at(top_c_itm->first);
            pt_sub=mapid_;
            //check_distr_equal(topic_word_sub[0], C.wn_occurences_global);
        } else {
            cout<<"Subtopics found!"<<endl;
        }
        
        // setting a new branch
        DI di_;
        old_topic_new_topic_labels.push_back(di_);
        
        // mapping the topic labels here to the global ones
        mapii label_loc_to_glob;

        for(map<int, mapid>::iterator itm=topic_word_sub.begin(); itm!=topic_word_sub.end(); itm++) {
            // topic relations
            old_topic_new_topic_labels[old_topic_new_topic_labels.size()-1].push_back(topic_label);
            // labels
            label_loc_to_glob[itm->first]=topic_label;
            // betas
            new_topics[topic_label]=itm->second;
            if(topic_word_sub.size()>1) {
                pout<<"  -> subtopic: "<<topic_label<<" #words: "<<itm->second.size()<<" pt: "<<pt_old.at(top_c_itm->first)*pt_sub.at(itm->first)<<endl;
                pout<<"  "<<oC.get_topic_title(itm->second)<<endl;
            }
            ++topic_label;
        }
        if(topic_word_sub.size()<=1)
            pout<<"  no subtopics"<<endl;
        pout<<"-------------------------"<<endl;
        RANGE_loop(i, doc_topic_sub) {
            // which doc is this?
            int original_doc = topic_old_docs.at(top_c_itm->first).at(i);
            //cout<<"local doc: "<<i<<" original_doc: "<<original_doc<<endl;
            // how much the doc is using the old topic            
            double old_topic_pr= doc_topic_old[original_doc].at(top_c_itm->first);
            // the topic distribution of this doc
            mapid & topic_distr = doc_new_topics[original_doc];
            IT_loop(mapid, itm, doc_topic_sub[i]) {
                topic_distr[label_loc_to_glob.at(itm->first)]=old_topic_pr*itm->second;
            }
        } 
    }

    // pt is used only in topic summary, which is not written here
    // (it has been written already)
    mapid empty_mapid;
    oC.write_short_beta_and_theta_files(doc_new_topics, new_topics, \
                                        "doc_topics2.txt", "topic_words2.txt", "", empty_mapid);
    if(fullout)
        oC.write_beta_and_theta_files(doc_new_topics, new_topics, "thetas2.txt", "betas2.txt");

}






