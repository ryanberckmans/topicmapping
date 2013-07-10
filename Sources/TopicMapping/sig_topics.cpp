

void word_corpus::set_from_file(string filename) {
    
	ofstream pout("word_wn_count.txt");
	ifstream gin(filename.c_str());
	
    multimap<int, int> occurences_wn;			// occurences - wn
    
	string buffer;
	while(getline(gin, buffer)) {
        
		doc newdoc;
		newdoc.set_from_string(buffer, word_number_global, number_word_global);
		IT_loop(mapii, itm, newdoc.wn_occurences) {			
			int_histogram(itm->first, wn_occurences_global, itm->second);
		}
		if(newdoc.num_words>0)
			docs.push_back(newdoc);
        
	}
	
	total_words=0;
	IT_loop(mapii, itm, wn_occurences_global) {			
		occurences_wn.insert(make_pair(itm->second,itm->first));
		total_words+=itm->second;
	}	
    
	cout<<"total_words: "<<total_words<<" total unique words: "<<wn_occurences_global.size()<<endl;
	cout<<"#docs "<<docs.size()<<endl;
	for(multimap<int, int>::iterator itm=occurences_wn.begin(); itm!=occurences_wn.end(); itm++) {			
		pout<<number_word_global[itm->second]<<" "<<itm->second<<" "<<itm->first<<endl;
	}	
    
	
    
}


void word_corpus::write_corpus_file() {
	
	ofstream pout("CORPUS.corpus");
	
	RANGE_loop(i, docs) {
		
		pout<<docs[i].wn_occurences.size()<<" ";
		IT_loop(mapii, itm, docs[i].wn_occurences) {
			pout<<itm->first<<":"<<itm->second<<" ";
		}
		pout<<endl;
	}
}

/*
void update_most_similar_words(map<int, pair<int, int> > & word_most_similar_word_weights, int w1, int w2, int weight) {

  
    if(word_most_similar_word_weights.count(w1)==0) {
        word_most_similar_word_weights[w1]=make_pair(-1,0);
    }
    
    if(word_most_similar_word_weights[w1].first==-1 or word_most_similar_word_weights[w1].second<weight) {
        word_most_similar_word_weights[w1]=make_pair(w2, weight);
    }
}*/

void word_corpus::null_model(double p_value, bool compute_only_first_order,\
                             DI & links1, DI & links2, DD & weights, bool use_dotproduct, bool verbose) {
	
	if(verbose) cout<<"null model. p-value: "<<p_value<<endl;
    // this is not efficient
	// EFFICIENCY-ALERT
    
    // word_norm2[wn] is the sum over docs of wn_occurences**2
    mapid word_norm2;
    int_matrix original_corpus;
	RANGE_loop(i, docs) {
		DI new_word_list;
		IT_loop(mapii, itm, docs[i].wn_occurences) {
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
    //included_words.assign(wn_occurences_global.size(), false);
    
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
		
		int k1=min(wn_occurences_global[it->first.first], wn_occurences_global[it->first.second]);
		int k2=max(wn_occurences_global[it->first.first], wn_occurences_global[it->first.second]);
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
    
    
    //cout<<"qfive_links: "<<qfive_links<<endl; 
    /*
    for(map<int, pair<int, int> >:: iterator itm= word_most_similar_word_weights.begin(); \
        itm!=word_most_similar_word_weights.end(); itm++) {
        
        // if the word was never included
        if(included_words_si.find(itm->first)==included_words_si.end()) {
            
            if(itm->second.second>0) {
                cerr<<"this should not happen "<<itm->second.second<<" "<<itm->first<<" "<<itm->second.first<<endl;
            }
            
            links1.push_back(itm->first);
            links2.push_back(itm->second.first);
            weights.push_back(1.);
            sigfile<<itm->first<<" "<<itm->second.first<<" "<<1.<<endl;
            // cooc_nullterm is not updated_here
            // also included_words should be removed
        }
    }
    */
    
	if(verbose) {
        cout<<"p-value: "<<p_value<<" significant link fraction: "<<double(qfive_links)/total_possible_links<<endl;
        cout<<"total possible links "<<total_possible_links<<endl;
    }
    
    if(compute_only_first_order==false) {
        
        if(number_word_global.size()==0) {
            cerr<<"empty number_word_global. bug!"<<endl;
            exit(-1);
        }
        
        //ofstream filtered_out("filtered_out.txt");
        //RANGE_loop(i, included_words) if(included_words[i]==false) {
        //    filtered_out<<i<<" "<<number_word_global[i]<<endl;
        //}
        
        // writing the second order network
        ofstream out_second("second_order.txt");
        //  only docs which have some significant words will be numbered
        ofstream doc_numbers("doc_numbers.txt");
        
        RANGE_loop(i, docs) {
            
            bool significant_doc=false;
            // a deque of just the word_numbers 
            DI words_in_docs;        
            IT_loop(mapii, itm, docs[i].wn_occurences) {
                words_in_docs.push_back(itm->first);
            }
            
            sort(words_in_docs.begin(), words_in_docs.end());
            
            RANGE_loop(j, words_in_docs) for(UI k=0; k<j; k++) {
                
                pair<int, int> word_pair(words_in_docs[k], words_in_docs[j]);
                
                if (word_pair.first>=word_pair.second) {
                    cerr<<"error in sorting words"<<endl;
                    exit(-1);
                }
                if(cooc.count(word_pair)==0 or cooc[word_pair]<=0) {
                    cerr<<"error in word_pair"<<endl;
                    exit(-1);
                }
                if(cooc_nullterm.count(word_pair)>0) {
                    if(significant_doc==false) {
                        significant_doc=true;
                        out_second<<"Slice:"<<endl;
                        doc_numbers<<i<<endl;
                    }
                    double K= 1. - double(cooc_nullterm[word_pair])/cooc[word_pair]; 
                    out_second<<word_pair.first<<" "<<word_pair.second<<" ";
                    out_second<<K*docs[i].wn_occurences[word_pair.first]*docs[i].wn_occurences[word_pair.second]<<endl;
                }
            }
        }
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
        docs[doc_number_original].compute_thetas(word_partition, theta_doc, missing_topic);
        doc_number_thetas[doc_number_original]=theta_doc;
        word_partition.clear();
        ++doc_number_multipart;
    }
}

void word_corpus::write_theta_file_from_multipart(string multipart_file, string doc_number_file) {
    
    
    
    // if this is true, all missing words
    // are assigned to one single topic (-1)
    bool just_one_noise_topic=true;
    
    
    /*
     the problem with the labels is that some 
     docs might be missing because they are just noise
     */
    
    //getting original document labels
    
    DI doc_numbers;
    ifstream gin(doc_number_file.c_str());
    int a;
    while(gin>>a) {
        doc_numbers.push_back(a);
    }
    gin.close();
    
    // for each doc (original label) a map of the topics used
    map<int, mapid> doc_number_thetas;
    
    gin.open(multipart_file.c_str());
    
    string s;
    int doc_number_multipart=0;
    int_matrix word_partition;
    
    while(getline(gin, s)) {
        
        if(s!="Slice:") {
            DI v;
            cast_string_to_doubles(s, v);
            // adding module to doc partition
            word_partition.push_back(v);
        } else {
            // doc partition is done
            update_doc_num_thetas(doc_numbers, just_one_noise_topic, doc_number_thetas, \
                                  word_partition, doc_number_multipart);
        }
    }
    // repeating for last doc
    update_doc_num_thetas(doc_numbers, just_one_noise_topic, doc_number_thetas, \
                          word_partition, doc_number_multipart);
    
    // setting theta doc for missing documents    
    RANGE_loop(i, docs) if(doc_number_thetas.find(i)==doc_number_thetas.end()) {
        mapid theta_doc;
        if(just_one_noise_topic)
            theta_doc[-1]=1;
        else
            theta_doc[-i-1]=1;
        doc_number_thetas[i]=theta_doc;
    }
    
    
    // getting new names for topics
    mapii new_names_for_topics;
    for(map<int, mapid>::iterator itm= doc_number_thetas.begin(); itm!=doc_number_thetas.end(); itm++) {
        IT_loop(mapid, itm2, itm->second) {
            if(new_names_for_topics.count(itm2->first)==0) {
                new_names_for_topics.insert(make_pair(itm2->first, new_names_for_topics.size()));
            }
        }
    }
    
    cout<<"new names"<<endl;
    prints(new_names_for_topics);
    
    // printing
    ofstream pout("multi_part_thetas.txt");
    for(map<int, mapid>::iterator itm= doc_number_thetas.begin(); itm!=doc_number_thetas.end(); itm++) {
        //cout<<"printing "<<itm->first<<endl;
        print_theta(itm->second, new_names_for_topics, pout);
    }
    
}


void word_corpus::get_rid_of_non_prevalent_topics(mapii & hard_mems, DI & doc_prevalent_topics) {
    
    // this function is deleting from hard_mems words which appear in topics
    // which are not prevalent in any document
    
    
    //cout<<"total words in partition: "<<hard_mems.size()<<endl;
    
    SI prevalent_topics;
    RANGE_loop(i, docs) {
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
    
    RANGE_loop(i, docs) {
        
        mapid topic_distr;
        
        IT_loop(mapii, itm, docs[i].wn_occurences) {
            
            if(hard_mems.count(itm->first)>0) {
                // if the word is in the partition
                int topic_num=hard_mems.at(itm->first);
                int_histogram(topic_num, topic_distr, double(itm->second)/docs[i].num_words);
                if(word_topic.count(itm->first)==0) {
                    mapii newmapii;
                    word_topic[itm->first]=newmapii;
                }
                int_histogram(topic_num, word_topic[itm->first], itm->second);
            }            
        }
        
        int prevalent_topic = doc_prevalent_topics[i];
        
        IT_loop(mapii, itm, docs[i].wn_occurences) {
            
            if(hard_mems.count(itm->first)==0) {
                
                // if the word is not here, I use the most significant topic in this doc
                int_histogram(prevalent_topic, topic_distr, double(itm->second)/docs[i].num_words);
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
    RANGE_loop(i, docs) {
        
        
        IT_loop(mapii, itm, docs[i].wn_occurences) {
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
    RANGE_loop(i, docs) {
        mapii void_mapii_;
        doc_assignments.push_back(void_mapii_);
    }
    
    
    RANGE_loop(i, docs) {
        
        int prevalent_topic = doc_prevalent_topics[i];
        
        IT_loop(mapii, itm, docs[i].wn_occurences) {
            
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
        if(occ!=wn_occurences_global[itm->first]) {
            cerr<<"occ is wrong for word  "<<itm->first<<endl;
            exit(-1);
        }
    }*/
    
    if(false) {
        
        for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
            cout<<"word:: "<<itm->first<<endl;
            prints(itm->second);
        }
        RANGE_loop(i, docs) {
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
        title+=number_word_global[pr_word[i].second]+" ";
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
        pout2<<"#docs "<<docs.size()<<endl;
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
            pout1<<number_word_global[pr_word[i].second]<<":"<<-pr_word[i].first<<" ";
            if (i<20 and pout2.is_open()) {
                pout2<<number_word_global[pr_word[i].second]<<" ";
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
    RANGE_loop(i, docs) {
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
    set_matrix_to_zero(topic_names.size(), wn_occurences_global.size(), betas);
    
    double smoothing_par_per_word=smoothing_par/(double(wn_occurences_global.size()));
    
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
    RANGE_loop(i, docs) max_doc_size=max(max_doc_size, int(docs[i].num_words));
    log_fact_table_._set_(max_doc_size);
    // for each topic, the probability we use it in the null model
    mapid topic_pr;
    
    IT_loop(mapii, itm, wn_occurences_global) {
        if (hard_mems.count(itm->first)>0) {
            int topic_num=hard_mems.at(itm->first);
            int_histogram(topic_num, topic_pr, double(itm->second));
        }
    }
    
    // so far, for each topic I am counting the occurrences of 
    // the words in that topic
    normalize_mapid(topic_pr);
    
    double min_number_of_docs=double(min_docs);
    // topics which are used in a very few documents
    SI banned_topics;
    IT_loop(mapid, itm, topic_pr) if( itm->second < min_number_of_docs / docs.size() ) {
        banned_topics.insert(itm->first);
    }
    
    
    RANGE_loop(i, docs) {
        
        mapii topic_usage;
        // getting topic_usage for this doc
        IT_loop(mapii, itm, docs[i].wn_occurences) {
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
            
            double bin_pvalue= log_fact_table_.cum_binomial_right(itm->second, docs[i].num_words, topic_pr.at(itm->first));
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
        
        //double aic= (docs.size()+wn_occurences_global.size())*(topic_word.size()-1) - loglikelihood;
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
            occ_wn.push_back(make_pair(-wn_occurences_global.at(itm->second[i]), itm->second[i]));
        }
        sort(occ_wn.begin(), occ_wn.end());
        RANGE_loop(i, occ_wn) {
            pout<<number_word_global.at(occ_wn[i].second)<<" ";
        }
        pout<<endl;
    }
    pout.close();
    
    
}



double word_corpus::dimap(double min_filter, double max_filter, int min_docs, \
                        int Nruns, double p_value, string partition_file, \
                        mapid & pt, \
                        deque<mapid> & doc_topic_best, \
                        map<int, mapid> & topic_word_best, \
                        deque<mapii> & doc_assignments, int level, const double & convergence_precision) {
    
    doc_topic_best.clear();
    topic_word_best.clear();
    doc_assignments.clear();
    pt.clear();
    
    mapii hard_memberships;
    if (partition_file.size()==0) {
        
        DI links1;
        DI links2;
        DD weights;
        // running null model
        null_model(p_value, true,\
                     links1, links2, weights, true, level==0);
        
        // collecting infomap initial partition
        if(links1.size()==0 and level==0) {
            cerr<<"empty graph"<<endl;
            cerr<<"Please try again with a bigger p-value"<<endl;
            exit(-1);
        }
        get_infomap_partition_from_edge_list(Nruns, irand(100000000), \
                                             links1, links2, weights, \
                                             hard_memberships, level==0, convergence_precision);
        links1.clear();
        links2.clear();
        weights.clear();
        if(level==0)
            write_partition(hard_memberships);
    }
    
    else {
        cout<<"reading partition from file: "<<partition_file<<endl;
        int_matrix ten;
        get_partition_from_file(partition_file, ten);
        RANGE_loop(i, ten) RANGE_loop(j, ten[i]) hard_memberships[ten[i][j]]=i;
    }
    
    // max-likelihood filter
    double eff_ntopics=optimal_filtering(hard_memberships, min_filter, max_filter, min_docs, \
                                         pt, doc_topic_best, topic_word_best, doc_assignments, \
                                         level==0);
    
                        
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


void set_corpuses_from_topics(deque<mapii> & doc_assignments, \
                              map<int, word_corpus> & old_topic_corpus, \
                              map<int, mapii> & topic_old_docs, word_corpus & oC) {
    
    
    
    // doc_topic_newid[doc][topic] is the new id of the doc in that corpus
    deque<mapii> doc_topic_newid;
    
    RANGE_loop(i, doc_assignments) {
        
        mapii topic_newid;
        IT_loop(mapii, itm, doc_assignments[i]) {
            
            // itm->second is the topic
            // itm->first is the word
            // i is the doc
            
            // adding this topic to the corpuses
            if(old_topic_corpus.count(itm->second)==0) {
                word_corpus word_corpus_;
                old_topic_corpus[itm->second]=word_corpus_;
            }
            
            // adding this doc to the corpus
            // it is the first time, we see this topic for this doc
            if(topic_newid.count(itm->second)==0) {
                topic_newid[itm->second]=old_topic_corpus[itm->second].docs.size();
                doc doc_;
                old_topic_corpus[itm->second].docs.push_back(doc_);
            }
            
            old_topic_corpus[itm->second].docs[topic_newid[itm->second]].wn_occurences.insert(make_pair(itm->first, oC.docs[i].wn_occurences[itm->first]));
        }
        // storing the doc ids
        doc_topic_newid.push_back(topic_newid);
    }
    
    invert_doc_topic_newid(doc_topic_newid, topic_old_docs);

    // asserting old_topic_corpus have labels in order
    for(UI i=0; i<old_topic_corpus.size(); i++) {
        if(old_topic_corpus.count(i)==0) {
            cerr<<"error in second level"<<endl;
            exit(-1);
        }
    } 
    
}



void word_corpus::fix_data_structs() {
    
    total_words=0;
    RANGE_loop(i, docs) {
        
        //cout<<"doc "<<i<<endl;
        //prints(docs[i].wn_occurences);
        
        docs[i].num_words=0;
        IT_loop(mapii, itm, docs[i].wn_occurences) {
            docs[i].num_words+=itm->second;
            int_histogram(itm->first, wn_occurences_global, itm->second);
        }
        total_words+=docs[i].num_words;
    }
}

