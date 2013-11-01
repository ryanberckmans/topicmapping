

void word_corpus::set_from_file(string filename) {
    
	ofstream pout("word_wn_count.txt");
	ifstream gin(filename.c_str());
	
    multimap<int, int> occurences_wn;			// occurences - wn
    mapsi word_number_all;
    mapis number_word_all;
    mapii wn_occurences__all;
    
    
	string buffer;
	while(getline(gin, buffer)) {
        
		doc newdoc;
		newdoc.set_from_string(buffer, word_number_all, number_word_all);
		IT_loop(mapii, itm, newdoc.wn_occurences_) {			
			int_histogram(itm->first, wn_occurences__all, itm->second);
		}
		if(newdoc.num_words_>0)
			docs_.push_back(newdoc);
        else {
            cerr<<"Some line is empty in "<<filename<<endl;
            cerr<<"Please avoid that, because it is messing with the document indexing"<<endl;
            cerr<<"Exiting..."<<endl;
            exit(-1);
        }
	}
	
	int total_words=0;
	IT_loop(mapii, itm, wn_occurences__all) {			
		occurences_wn.insert(make_pair(itm->second,itm->first));
		total_words+=itm->second;
	}	
    
	cout<<"total_words: "<<total_words<<" total unique words: "<<wn_occurences__all.size()<<endl;
	cout<<"#docs_ "<<docs_.size()<<endl;
	for(multimap<int, int>::iterator itm=occurences_wn.begin(); itm!=occurences_wn.end(); itm++) {			
		pout<<number_word_all[itm->second]<<" "<<itm->second<<" "<<itm->first<<endl;
	}
    
    
    
    // word_occurrences_ is copied from wn_occurences__all
	word_occurrences_.clear();
    IT_loop(mapii, itm, wn_occurences__all) {
        word_occurrences_.push_back(itm->second);
        assert_ints(itm->first, word_occurrences_.size()-1);
    }

    // word_strings_ is copied from number_word_all
    word_strings_.clear();	
    IT_loop(mapis, itm, number_word_all) {
        word_strings_.push_back(itm->second);
        assert_ints(itm->first, word_strings_.size()-1);
    }
    

	
    
}


void word_corpus::write_corpus_file() {
	
	ofstream pout("CORPUS.corpus");
	
	RANGE_loop(i, docs_) {
		
		pout<<docs_[i].wn_occurences_.size()<<" ";
		IT_loop(mapii, itm, docs_[i].wn_occurences_) {
			pout<<itm->first<<":"<<itm->second<<" ";
		}
		pout<<endl;
	}
}


void word_corpus::null_model(double p_value, \
                             DI & links1, DI & links2, DD & weights, \
                             bool use_dotproduct, bool verbose) {
	
	if(verbose) cout<<"null model. p-value: "<<p_value<<endl;
    // this is not efficient
	// EFFICIENCY-ALERT
    
    // word_norm2[wn] is the sum over docs_ of wn_occurences_**2
    mapid word_norm2;
    int_matrix original_corpus;
	RANGE_loop(i, docs_) {
		DI new_word_list;
		IT_loop(mapii, itm, docs_[i].wn_occurences_) {
			for(int occ=0; occ<itm->second; occ++) {
				new_word_list.push_back(itm->first);
            }
            int_histogram(itm->first, word_norm2, double(itm->second)*double(itm->second));
		}
		original_corpus.push_back(new_word_list);
	}
    IT_loop(mapid, itm, word_norm2) itm->second=sqrt(itm->second);
    
	degree_block DB;
	DB.random_similarities_poissonian_set_data(original_corpus);
	    
	map<pair<int, int> , int> cooc;
	map<pair<int, int> , int> cooc_nullterm;
	dotpr_similarity_of_connected_words(original_corpus, cooc);
    
    
    // included_words[word_id] is true if the word is not isolated
    // this requires that words are labelled from 0 and it is only
    // true is corpus was set from file
    //deque<bool> included_words;
    //included_words.assign(word_occurrences_.size(), false);
    
	int qfive_links=0;
	int total_possible_links=0;
    
    // words which have at least a significant connection
    SI included_words_si;
    //map<int, pair<int, int> > word_most_similar_word_weights;
    
    // if you want to see the network, uncomment sigfile
    // and check that there is no second level (or change file names for second level)
    //ofstream sigfile("similar_words.txt");
	
    links1.clear(); links2.clear(); weights.clear();
    // computing the significant links and writing the first order network
	for(map<pair<int, int> , int>::iterator it=cooc.begin(); it!=cooc.end(); it++) {
		
		++total_possible_links;
		
		int k1=min(word_occurrences_[it->first.first], word_occurrences_[it->first.second]);
		int k2=max(word_occurrences_[it->first.first], word_occurrences_[it->first.second]);
		int qfive=DB.qfive(k1,k2,p_value);
		
        
        //update_most_similar_words(word_most_similar_word_weights, it->first.first, it->first.second, it->second- qfive);
        //update_most_similar_words(word_most_similar_word_weights, it->first.second, it->first.first, it->second- qfive);
                
        // significant links
        if(it->second- qfive>0) {
            
            /*if(compute_only_first_order==false) {
                // you should remove this
                included_words[it->first.first]=true;
                included_words[it->first.second]=true;
            }*/

            included_words_si.insert(it->first.first);
            included_words_si.insert(it->first.second);

            links1.push_back(it->first.first);
            links2.push_back(it->first.second);
            if(use_dotproduct) {
                weights.push_back(double(it->second- qfive));
            } else {
                //double w_= double(it->second- qfive) / (word_norm2[it->first.first] * word_norm2[it->first.second]);
                double w_= double(it->second- qfive) / (qfive+1.);
                weights.push_back(w_);
                //cout<<"words "<<it->first.first<<" "<<it->first.second<<" "<<double(it->second- qfive);
                //cout<<" "<< (word_norm2[it->first.first] * word_norm2[it->first.second])<<" "<<w_<<endl;
            }
            
            //sigfile<<it->first.first<<" "<<it->first.second<<" "<<double(it->second- qfive)<<endl;
            ++qfive_links;
            cooc_nullterm.insert(make_pair(it->first, qfive));
		}
	}
    
    
    
	if(verbose) {
        cout<<"p-value: "<<p_value<<" significant link fraction: "<<double(qfive_links)/total_possible_links<<endl;
        cout<<"total possible links "<<total_possible_links<<endl;
    }
}



void print_theta(mapid & theta_doc, mapii & new_names_for_topics, ostream & pout) {
    
    DD thetas;
    thetas.assign(new_names_for_topics.size(), 0.);
    IT_loop(mapid, itm, theta_doc) {
        thetas[new_names_for_topics[itm->first]]=itm->second;
    }
    prints(thetas, pout);
}



void word_corpus::update_doc_num_thetas(DI & doc_numbers,
                                        bool just_one_noise_topic, 
                                        map<int, mapid> & doc_number_thetas, 
                                        int_matrix & word_partition, 
                                        int & doc_number_multipart) {
    if(word_partition.size()>0) {
        
        // getting original label
        int doc_number_original=doc_numbers[doc_number_multipart];
        mapid theta_doc;
        // getting theta
        int missing_topic=-doc_number_original-1;
        if(just_one_noise_topic)
            missing_topic=-1;
        docs_[doc_number_original].compute_thetas(word_partition, theta_doc, missing_topic);
        doc_number_thetas[doc_number_original]=theta_doc;
        word_partition.clear();
        ++doc_number_multipart;
    }
}



void word_corpus::get_rid_of_non_prevalent_topics(mapii & hard_mems, DI & doc_prevalent_topics) {
    
    // this function is deleting from hard_mems words which appear in topics
    // which are not prevalent in any document
    
    
    //cout<<"total words in partition: "<<hard_mems.size()<<endl;
    
    SI prevalent_topics;
    RANGE_loop(i, docs_) {
        prevalent_topics.insert(doc_prevalent_topics[i]);
    }
    
    
    //cout<<"prevalent topics: "<<prevalent_topics.size()<<endl;
    mapii new_hard_mems;
    IT_loop(mapii, itm, hard_mems) {
        // if the topic is prevalent in at least one doc
        if(prevalent_topics.find(itm->second)!=prevalent_topics.end()) {
            new_hard_mems.insert(*itm);
        }
    }
    
    hard_mems=new_hard_mems;
    //cout<<"words in partition after removing small topics: "<<hard_mems.size()<<endl;
    
    
    
}


void word_corpus::initial_ptopic(deque<mapid> & doc_topic, map<int, mapii> & word_topic, \
                                 const mapii & hard_mems, const DI & doc_prevalent_topics) {
    
    // we set doc_topic[doc] as p(t|doc)
    // and word_topic[word] as n(w,t)
    // pt[topic] is n(t)
    
    doc_topic.clear();
    word_topic.clear();
    
    RANGE_loop(i, docs_) {
        
        mapid topic_distr;
        
        IT_loop(mapii, itm, docs_[i].wn_occurences_) {
            
            if(hard_mems.count(itm->first)>0) {
                // if the word is in the partition
                int topic_num=hard_mems.at(itm->first);
                int_histogram(topic_num, topic_distr, double(itm->second)/docs_[i].num_words_);
                if(word_topic.count(itm->first)==0) {
                    mapii newmapii;
                    word_topic[itm->first]=newmapii;
                }
                int_histogram(topic_num, word_topic[itm->first], itm->second);
            }            
        }
        
        int prevalent_topic = doc_prevalent_topics[i];
        
        IT_loop(mapii, itm, docs_[i].wn_occurences_) {
            
            if(hard_mems.count(itm->first)==0) {
                
                // if the word is not here, I use the most significant topic in this doc
                int_histogram(prevalent_topic, topic_distr, double(itm->second)/docs_[i].num_words_);
                if(word_topic.count(itm->first)==0) {
                    mapii newmapii;
                    word_topic[itm->first]=newmapii;
                }
                int_histogram(prevalent_topic, word_topic[itm->first], itm->second);  
            }            
        }
        
        doc_topic.push_back(topic_distr);
        
    }
    
}


void word_corpus::get_betas(map<int, mapii> & word_topic, map<int, mapid> & topic_word, mapid & pt) {
    
    pt.clear();
    topic_word.clear();
    
    for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
        
        IT_loop(mapii, itm2, itm->second) if(itm2->second>1e-10) {
            if(topic_word.count(itm2->first)==0) {
                mapid newmapid;
                topic_word[itm2->first]=newmapid;
            }
            int_histogram(itm->first, topic_word[itm2->first], double(itm2->second));
        }
    }
    
    // 
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        pt[itm->first]=normalize_mapid(itm->second);
    }
    normalize_mapid(pt);
}


double word_corpus::compute_likelihood(map<int, mapid> & topic_word, map<int, mapii> & word_topic, deque<mapid> & doc_topic) {
    
    double loglik=0.;
    RANGE_loop(i, docs_) {
        
        
        IT_loop(mapii, itm, docs_[i].wn_occurences_) {
            double pr=0;
            mapii & possible_topics = word_topic[itm->first];
            
            // intersection of possible_topics & doc_topic[i]
            mapid intermap;
            mapid::iterator out_itr( intermap.begin() );
            set_intersection(doc_topic[i].begin(), doc_topic[i].end(), 
                             possible_topics.begin(), possible_topics.end(), 
                             inserter( intermap, out_itr ), doc_topic[i].value_comp());

            IT_loop(mapid, itm2, intermap) {
                if(topic_word.count(itm2->first)>0) {
                    pr+= itm2->second * get_from_mapid(topic_word.at(itm2->first), itm->first);
                }
            }
            
            if(pr==0) {
                
                cout<<"doc "<<i<<endl;
                cout<<"word "<<itm->first<<endl;
                prints(doc_topic[i]);
                prints(word_topic[itm->first]);                 
                cerr<<"error in pr"<<endl;
                exit(-1);
            }
            loglik+= itm->second * log(pr);
        }
    }
    return loglik;
}



double word_corpus::likelihood_filter(map<int, mapii> & word_topic, deque<mapid> & doc_topic, \
                                      map<int, mapid> & topic_word, mapid & pt, \
                                      const mapii & hard_mems, double filtering_par, const DI & doc_prevalent_topics, \
                                      deque<mapii> & doc_assignments) {
    
    topic_word.clear();
    pt.clear();
    
    // we set doc_topic[doc] as p(t|doc)
    // and word_topic[word] as n(w,t)
    
    // doc_assignments[doc][wn] is the topic to which the word has been assigned
    doc_assignments.clear();
    RANGE_loop(i, docs_) {
        mapii void_mapii_;
        doc_assignments.push_back(void_mapii_);
    }
    
    
    RANGE_loop(i, docs_) {
        
        int prevalent_topic = doc_prevalent_topics[i];
        
        IT_loop(mapii, itm, docs_[i].wn_occurences_) {
            
            // topic num is either the prevalent_topic or the topic
            // the word belongs to originally 
            int topic_num = prevalent_topic;
            if(hard_mems.count(itm->first)>0) {
                topic_num=hard_mems.at(itm->first);
            }
            // if the topic is small, we move words to prevalent_topic
            if(doc_topic[i].at(topic_num)<filtering_par and topic_num!=prevalent_topic) {
                word_topic[itm->first][topic_num]-=itm->second;
                int_histogram(prevalent_topic, word_topic[itm->first], itm->second);
            }
            // now topic_num is prevalent_topic, if it is small
            if(doc_topic[i].at(topic_num)<filtering_par) {
                topic_num=prevalent_topic;
            }
            // we insert the word assignement
            doc_assignments[i][itm->first]=topic_num;
        }
        
        IT_loop(mapid, itm, doc_topic[i]) if(itm->second<filtering_par and itm->first!=prevalent_topic) {
            doc_topic[i][prevalent_topic]+=itm->second;
            itm->second=0;
        }
    }
    
    get_betas(word_topic, topic_word, pt);
    
    /**
    // checking word_topic is correct
    for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
        int occ=0; IT_loop(mapii, itm2, itm->second) occ+=itm2->second;
        if(occ!=word_occurrences_[itm->first]) {
            cerr<<"occ is wrong for word  "<<itm->first<<endl;
            exit(-1);
        }
    }*/
    
    if(false) {
        
        for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
            cout<<"word:: "<<itm->first<<endl;
            prints(itm->second);
        }
        RANGE_loop(i, docs_) {
            cout<<"doc:: "<<i<<endl;
            prints(doc_topic[i]);
        }
        for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
            cout<<"topic:: "<<itm->first<<endl;
            prints(itm->second);
        }
        cout<<"pt"<<endl;
        prints(pt);
    }
    
    return compute_likelihood(topic_word, word_topic, doc_topic);
    
}


string word_corpus::get_topic_title(const mapid & topic_distr) {
    
    // returning the most probable 20 words in this topic
    
    deque<pair<double, int> > pr_word;
    IT_const_loop(mapid, itm2, topic_distr) {
        pr_word.push_back(make_pair(-itm2->second, itm2->first));
    }
    sort(pr_word.begin(), pr_word.end());
    
    string title="";
    for(int i=0; i<min(20,int(pr_word.size())); i++){
        title+=word_strings_[pr_word[i].second]+" ";
    }
    return title;

}


void word_corpus::write_short_beta_and_theta_files(deque<mapid> & doc_topic, map<int, mapid> & topic_word, \
                                             string theta_file, string beta_file, string beta_file_short, \
                                             mapid & pt) {
    
    
    // writing handy files:
    /*
     one for the topics of a document
     one for the word distribution of a topic
     and one for the topic names
     */
    
    // for each topic, I use consecutive names
    mapii topic_names;
    
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        int tp=topic_names.size();
        topic_names[tp]=itm->first;
    }
    
    // asserting labels are correct
    IT_loop(mapii, itm, topic_names) {
        if (itm->first!=itm->second) {
            cerr<<"error in topic names"<<endl;
            exit(-1);
        }
    }
    
    ofstream pout1(theta_file.c_str());
    RANGE_loop(i, doc_topic) {
        IT_loop(mapid, itm, doc_topic[i])
            pout1<<itm->first<<":"<<itm->second<<" ";
        pout1<<endl;
    }
    pout1.close();
    pout1.open(beta_file.c_str());
    ofstream pout2;
    // the file is opened only if beta_file_short is passed
    if(beta_file_short.size()>0) {
        pout2.open(beta_file_short.c_str());
        pout2<<"#docs_ "<<docs_.size()<<endl;
    }
    
    IT_loop(mapii, itm, topic_names) {
        // printing topic distribution over words
        mapid & topic_distr = topic_word[itm->second];
        deque<pair<double, int> > pr_word;
        IT_loop(mapid, itm2, topic_distr) {
            pr_word.push_back(make_pair(-itm2->second, itm2->first));
        }
        sort(pr_word.begin(), pr_word.end());
        if(pout2.is_open()) pout2<<"topic: "<<itm->first<<" #words: "<<pr_word.size()<<" pt: "<<pt.at(itm->first)<<endl;
        RANGE_loop(i, pr_word){
            pout1<<word_strings_[pr_word[i].second]<<":"<<-pr_word[i].first<<" ";
            if (i<20 and pout2.is_open()) {
                pout2<<word_strings_[pr_word[i].second]<<" ";
            } 
        }
        pout1<<endl;
        if(pout2.is_open()) pout2<<endl;
    }
    pout1.close();
    pout2.close();
}




void word_corpus::write_beta_and_theta_files(deque<mapid> & doc_topic, map<int, mapid> & topic_word, \
                                             string theta_file, string beta_file) {
    
    // for each topic, I use consecutive names
    mapii topic_names;
    
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        int tp=topic_names.size();
        topic_names[tp]=itm->first;
    }
    
    // asserting labels are correct
    for(UI i=0; i<topic_names.size(); i++) {
        if(topic_names.count(i)==0) {
            cerr<<"error in label names"<<endl;
            exit(-1);
        }
    } 

    ofstream pout1(theta_file.c_str());
    RANGE_loop(i, docs_) {
        IT_loop(mapii, itm, topic_names) {
            if(doc_topic[i].count(itm->second)>0) {
                pout1<<doc_topic[i][itm->second]<<" ";
            }
            else {
                pout1<<"0 ";
            }
        }
        pout1<<endl;
    }
    
    double smoothing_par=1e-4;
    
    deque<DD> betas;
    set_matrix_to_zero(topic_names.size(), word_occurrences_.size(), betas);
    
    double smoothing_par_per_word=smoothing_par/(double(word_occurrences_.size()));
    
    IT_loop(mapii, itm, topic_names) {
        
        DD & row= betas[itm->first];
        mapid & topic_distr = topic_word[itm->second];
        IT_loop(mapid, itm2, topic_distr) {
            double smoothed_prob= (itm2->second+smoothing_par_per_word) / (1+smoothing_par) ;
            row[itm2->first]=log(smoothed_prob);
        }
        RANGE_loop(w, row) if(row[w]==0) row[w]=log((smoothing_par_per_word) / (1+smoothing_par));
    }
    
    print_matrix_to_file_or_screen(betas, beta_file);
}

double get_z_score(int N, double p, int real) {
    
    /* this is the z-score of a binomial distribution
     we could actually compute the right cumularive.
     but the approximation should hold pretty well */
    
    
    double expected_usage =  p * N;
    double std= sqrt(max(0., p * N * (1.-p) ));
    double z_score=0;
    if(std<1e-10) {
        z_score= 1e300;
    } else {
        z_score = (double(real) - expected_usage)/std;
    }
    return z_score;
    
}


void word_corpus::get_prevalent_topics(DI & doc_prevalent_topics, mapii & hard_mems, int min_docs) {
    
    /* first I need to compute the probability of using a topic
     in the null model, where we draw words at random according to 
     simply their probability */
    
    
    log_fact_table log_fact_table_;
    int max_doc_size=0;
    RANGE_loop(i, docs_) max_doc_size=max(max_doc_size, int(docs_[i].num_words_));
    log_fact_table_._set_(max_doc_size);
    // for each topic, the probability we use it in the null model
    mapid topic_pr;
    
    RANGE_loop(wn, word_occurrences_) {
        if (hard_mems.count(wn)>0) {
            int topic_num=hard_mems.at(wn);
            int_histogram(topic_num, topic_pr, double(word_occurrences_[wn]));
        }
    }
    
    // so far, for each topic I am counting the occurrences of 
    // the words in that topic
    normalize_mapid(topic_pr);
    
    double min_number_of_docs=double(min_docs);
    // topics which are used in a very few documents
    SI banned_topics;
    IT_loop(mapid, itm, topic_pr) if( itm->second < min_number_of_docs / docs_.size() ) {
        banned_topics.insert(itm->first);
    }
    
    
    RANGE_loop(i, docs_) {
        
        mapii topic_usage;
        // getting topic_usage for this doc
        IT_loop(mapii, itm, docs_[i].wn_occurences_) {
            if(hard_mems.count(itm->first)>0) {
                int topic_num=hard_mems.at(itm->first);
                int_histogram(topic_num, topic_usage, itm->second);
            }
        }
        
        double smallest_binpvalue=2;
        double smallest_binpvalue_non_banned=2;
        int prevalent_topic=-1;
        int prevalent_topic_non_banned=-1;
        
        IT_loop(mapii, itm, topic_usage) {            
            // topic_id is itm->first
            // usage is itm->second            
            
            double bin_pvalue= log_fact_table_.cum_binomial_right(itm->second, docs_[i].num_words_, topic_pr.at(itm->first));
            if(bin_pvalue<smallest_binpvalue_non_banned and banned_topics.find(itm->first) == banned_topics.end()) {
                prevalent_topic_non_banned=itm->first;
                smallest_binpvalue_non_banned=bin_pvalue;
            }
            if(bin_pvalue<smallest_binpvalue) {
                prevalent_topic=itm->first;
                smallest_binpvalue=bin_pvalue;
            }
        }
        // prevalent_topic==-1 means that the document had just filtered out words
        if(prevalent_topic_non_banned!=-1) {
            doc_prevalent_topics.push_back(prevalent_topic_non_banned);
        } else {
            doc_prevalent_topics.push_back(prevalent_topic);
        }
    }
    
    
    get_rid_of_non_prevalent_topics(hard_mems, doc_prevalent_topics);
    int number_of_topics=0; IT_loop(mapii, itm, hard_mems) number_of_topics=max(number_of_topics, itm->second);
    
    RANGE_loop(i, doc_prevalent_topics) {
        int prevalent_topic=doc_prevalent_topics[i];
        if(prevalent_topic==-1) {
            ++number_of_topics;
            prevalent_topic=number_of_topics;
            doc_prevalent_topics[i]=prevalent_topic;
        }
    }
}


double compute_eff_num_topics(const mapid & pt) {
    double h=0.;
    IT_const_loop(mapid, itm, pt) h+=itm->second*log2(itm->second);
    return pow(2,-1*h);
}


double word_corpus::optimal_filtering(mapii & hard_mems, double min_filter, double max_filter, int min_docs, \
                                      mapid & pt_best, \
                                      deque<mapid> & doc_topic_best, \
                                      map<int, mapid> & topic_word_best, \
                                      deque<mapii> & doc_assignments_best, bool verbose) {
    
    
    // for each document, the topic which we believe
    // the document mostly belongs to
    DI doc_prevalent_topics;
    get_prevalent_topics(doc_prevalent_topics, hard_mems, min_docs);    
    
    // this is the starting point (the same for all values of filtering)
    deque<mapid> doc_topic_initial;
    map<int, mapii> word_topic_initial;
    if(verbose) cout<<"initial step"<<endl;
    initial_ptopic(doc_topic_initial, word_topic_initial, hard_mems, doc_prevalent_topics);
    if(verbose) cout<<"initial step done"<<endl;
    
    // these are the parameters of the best model
    double max_ll=-1e300;
    double optimal_par=0.;
    
    double eff_ntopics=0.;
    
    for(double filtering_par=min_filter; filtering_par<max_filter+0.01; filtering_par+=0.01) {
        
        map<int, mapii> word_topic; // for each word, {topic:usage}
        deque<mapid> doc_topic;     // for each doc, {topic:probability}
        map<int, mapid> topic_word;
        mapid pt;
        deque<mapii> doc_assignments;
        
        word_topic= word_topic_initial;
        doc_topic= doc_topic_initial;
        
        // filtering 
        double loglikelihood=likelihood_filter(word_topic, doc_topic, topic_word, pt, hard_mems, filtering_par, doc_prevalent_topics, doc_assignments);
        
        //double aic= (docs_.size()+word_occurrences_.size())*(topic_word.size()-1) - loglikelihood;
        if(verbose) {
            cout<<"filtering: "<<filtering_par<<" loglikelihood: "<<loglikelihood<<" #topics: "<<topic_word.size()<<endl;
        }
        
        if(loglikelihood>max_ll) {
            max_ll=loglikelihood;
            optimal_par=filtering_par;
            doc_topic_best=doc_topic;
            topic_word_best=topic_word;
            doc_assignments_best=doc_assignments;
            pt_best=pt;
            eff_ntopics=compute_eff_num_topics(pt);
            if(verbose) cout<<"best filtering so far: "<<optimal_par<<endl;
        }
        
    }

    // sorting topic names so that they start from zero and there are no gaps
    make_topic_names_consecutive(doc_topic_best, topic_word_best, doc_assignments_best, pt_best);
    
    if(verbose) cout<<"optimal filtering: "<<optimal_par<<endl;

    return eff_ntopics;
    
}



void word_corpus::write_partition(mapii & hard_mems) {
    
    map<int, DI> partition;
    IT_loop(mapii, itm, hard_mems) {
        if(partition.count(itm->first)==0) {
            DI voiddi;
            partition.insert(make_pair(itm->second, voiddi));
        }
        partition[itm->second].push_back(itm->first);
    }
    ofstream pout("infomap.part");
    for(map<int, DI>::iterator itm= partition.begin(); itm!=partition.end(); itm++) {
        prints(itm->second, pout);
    }
    pout.close();
    pout.open("infomap-words.part");
    for(map<int, DI>::iterator itm= partition.begin(); itm!=partition.end(); itm++) {
        deque<pair<int, int> > occ_wn;
        RANGE_loop(i, itm->second) {
            occ_wn.push_back(make_pair(-word_occurrences_.at(itm->second[i]), itm->second[i]));
        }
        sort(occ_wn.begin(), occ_wn.end());
        RANGE_loop(i, occ_wn) {
            pout<<word_strings_.at(occ_wn[i].second)<<" ";
        }
        pout<<endl;
    }
    pout.close();
    
    
}



double word_corpus::dimap(int Nruns, \
                          mapid & pt, \
                          deque<mapid> & doc_topic_best, \
                          map<int, mapid> & topic_word_best, \
                          deque<mapii> & doc_assignments) {
    
    doc_topic_best.clear();
    topic_word_best.clear();
    doc_assignments.clear();
    pt.clear();
    
    mapii hard_memberships;
    if (partition_file_.size()==0) {
        
        DI links1;
        DI links2;
        DD weights;
        // running null model
        null_model(p_value_, links1, links2, weights, true);
        
        // collecting infomap initial partition
        if(links1.size()==0 and level==0) {
            cerr<<"empty graph"<<endl;
            cerr<<"Please try again with a bigger p-value"<<endl;
            exit(-1);
        }
        get_infomap_partition_from_edge_list(Nruns, irand(100000000), \
                                             links1, links2, weights, \
                                             hard_memberships);
        links1.clear();
        links2.clear();
        weights.clear();
        write_partition(hard_memberships);
    } else {
        cout<<"reading partition from file: "<<partition_file_<<endl;
        int_matrix ten;
        get_partition_from_file(partition_file, ten);
        RANGE_loop(i, ten) RANGE_loop(j, ten[i]) hard_memberships[ten[i][j]]=i;
    }
    
    // max-likelihood filter
    double eff_ntopics=optimal_filtering(hard_memberships, 
                                         pt, doc_topic_best,
                                         topic_word_best, doc_assignments);
                        
    return eff_ntopics;

}



void invert_doc_topic_newid(deque<mapii> & doc_topic_newid, map<int, mapii> & topic_old_docs) {
    
    topic_old_docs.clear();
    RANGE_loop(i, doc_topic_newid) {
        IT_loop(mapii, itm, doc_topic_newid[i]) {

            int topic=itm->first;
            int new_id=itm->second;
            
            if(topic_old_docs.count(topic)==0) {
                mapii mapii_;
                topic_old_docs[topic]=mapii_;
            }
            // in this topic, here is how the new id maps to the original doc
            topic_old_docs[topic][new_id]=i;
        }
    }
}







int sample_topic(DD & probs, DI & topics, double & sum) {
    
    // check this
    //prints(topics);
    //prints(probs);
    int nn=lower_bound(probs.begin(), probs.end(), ran4()*sum)-probs.begin();
    //cout<<"nn:: "<<nn<<" "<<topics.size()<<endl;
    return topics[nn];
    
}

void word_corpus::gibbs_sampling(deque<mapii> & doc_assignments) {
    
    
    // you probably need arrays
    
    
    double alpha_hyper=0.01;
    double beta_hyper=0.01;
    
    deque<mapii> n_d_topic;
    map<int, mapii> n_w_topic;
    mapii n_topic;
    
    if(doc_assignments.size()!=docs_.size()) {  cerr<<"doc_assignments size does not match"<<endl; exit(-1);  }
    
    
    // initialization
    RANGE_loop(doc_number, doc_assignments) {
        
        mapii d_topic;
        IT_loop(mapii, itm, doc_assignments[doc_number]) {            
            
            //cout<<"-------------- "<<doc_number<<endl;
            //prints(doc_assignments[i]);
            int_histogram(itm->second, d_topic);
            if(n_w_topic.count(itm->first)==0) {
                mapii new_mapii;
                n_w_topic[itm->first]=new_mapii;
            }
            int_histogram(itm->second, n_w_topic[itm->first]);
            int_histogram(itm->second, n_topic);
        }
        n_d_topic.push_back(d_topic);
    }
    
    // sampling
    for(int iter=0; iter<1000; iter++) RANGE_loop(doc_number, doc_assignments) {
        
        //cout<<"-------------- "<<doc_number<<" "<<n_d_topic[doc_number].size()<<endl;
        
        IT_loop(mapii, itm, doc_assignments[doc_number]) {            
            
            int topic= itm->second;
            int word= itm->first;
            // make sure they are positive
            n_d_topic[doc_number].at(topic)-=1;
            n_w_topic[word].at(topic)-=1;
            n_topic.at(topic)-=1;
            
            if(n_d_topic[doc_number][topic]==0) n_d_topic[doc_number].erase(topic);
            if(n_w_topic[word][topic]==0) n_w_topic[word].erase(topic);
            if(n_topic[topic]==0) n_topic.erase(topic);
            
            // looping over doc_topics
            DI topics;
            DD probs;
            double sum=0.;
            IT_loop(mapii, itm2, n_d_topic[doc_number]) {
                double pr= (itm2->second + alpha_hyper) * (n_w_topic[word][itm2->first]+ beta_hyper);
                probs.push_back(sum+pr);
                topics.push_back(itm2->first);
                sum+=pr;
            }
            // check  norm
            topic=sample_topic(probs, topics, sum);
            itm->second=topic;
            int_histogram(topic, n_d_topic[doc_number]);
            int_histogram(topic, n_w_topic[word]);
            int_histogram(topic, n_topic);            
        }
    }
    
    
    // topics can be empty... check this
    
    deque<mapid> n_d_topic_d;
    map<int, mapid> n_topic_w_d;
    
    
    
    RANGE_loop(doc_number, n_d_topic) {
        mapid new_mapid;
        IT_loop(mapii, itm, n_d_topic[doc_number]) new_mapid[itm->first]=itm->second;
        cout<<doc_number<<" ---- "<<n_d_topic[doc_number].size()<<endl;
        normalize_mapid(new_mapid);
        n_d_topic_d.push_back(new_mapid);
    }
    
    
    for (map<int, mapii>::iterator word_itm= n_w_topic.begin(); 
         word_itm!=n_w_topic.end(); word_itm++) {
        
        IT_loop(mapii, itm2, word_itm->second) {
            
            int word= word_itm->first;
            int topic= itm2->first;
            int count= itm2->second;
            if(n_topic_w_d.count(topic)==0) {
                mapid new_mapid;
                n_topic_w_d[topic]=new_mapid;
            }
            n_topic_w_d[topic][word]=count;
        }
    }
    
    for(map<int, mapid>::iterator topic_itm= n_topic_w_d.begin(); 
        topic_itm!=n_topic_w_d.end(); topic_itm++) {
        
        normalize_mapid(topic_itm->second);
        
    }
    
    
    write_beta_and_theta_files(n_d_topic_d, n_topic_w_d, "thetas_gs.txt", "betas_gs.txt");
    
}


