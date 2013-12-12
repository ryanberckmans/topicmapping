
void word_corpus::set_from_file_first_time(string filename, string out_dir) {
    
	ofstream pout((out_dir+"/word_wn_count.txt").c_str());
	ifstream gin(filename.c_str());
	
    multimap<int, int> occurences_wn;			// occurences - wn
    mapsi word_number_all;
    mapis number_word_all;
    mapii wn_occurences__all;
    
    
	string buffer;
	while(getline(gin, buffer)) {
        
		doc newdoc;
		newdoc.set_from_string(buffer, word_number_all, number_word_all);
		IT_loop(deqii, itm, newdoc.wn_occs_) {
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
    gin.close();
	
	int total_words=0;
	IT_loop(mapii, itm, wn_occurences__all) {			
		occurences_wn.insert(make_pair(itm->second,itm->first));
		total_words+=itm->second;
	}	
    
	cout<<"total_words: "<<total_words<<" total unique words: "<<wn_occurences__all.size()<<endl;
	cout<<"#docs "<<docs_.size()<<endl;
	for(multimap<int, int>::iterator itm=occurences_wn.begin(); itm!=occurences_wn.end(); itm++) {			
		pout<<number_word_all[itm->second]<<" "<<itm->second<<" "<<itm->first<<endl;
	}
    pout.close();
    
    
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


void word_corpus::set_from_file(string filename, string wn_file) {
    
    /*
     only words which are found in this string will be considered
     word_strings_ and word_occurrences_ are set from file
     docs gets their word ids from word_occurrences_
     words which are not found in word_strings_ are discarded
     there might be words in word_strings_ which never appear in the datasets
     that's fine unless we conpute p(topic|word)
     */
    

    // ==== getting word_strings_ and word_occurrences_ =============
    
    word_strings_.clear();          // word_strings_[wn] = "word_in_text"	
	word_occurrences_.clear();      // word_occurrences_[wn] = occurences
    map<int, pair<string, int> > wn_word_occ_map;
    string word_str;
    int wn;
    int occ;
    ifstream gin_wn(wn_file.c_str());
    while(gin_wn>>word_str) {
        gin_wn>>wn; gin_wn>>occ;
        wn_word_occ_map[wn]=make_pair(word_str, occ);
    }
    int counter=0;
    mapsi word_wn_all;
    for(map<int, pair<string, int> >::iterator itm= wn_word_occ_map.begin(); \
        itm!=wn_word_occ_map.end(); itm++) {
        assert_ints(itm->first, counter, "error in word_wn_count.txt. word_ids are not consecutive!");
        // itm->second.first is a string
        word_strings_.push_back(itm->second.first);
        // itm->second.second is the occs
        word_occurrences_.push_back(itm->second.second);
        ++counter;
        word_wn_all[itm->second.first]=itm->first;
    }
    gin_wn.close();
    
    // ============= setting documents ===============================

    ifstream gin(filename.c_str());
	string buffer;
	while(getline(gin, buffer)) {
        
		doc newdoc;
		newdoc.set_from_string_given_wn(buffer, word_wn_all);
		if(newdoc.num_words_>0)
			docs_.push_back(newdoc);
        else {
            cerr<<"Some line is empty in "<<filename<<endl;
            cerr<<"Please avoid that, because it is messing with the document indexing"<<endl;
            cerr<<"Exiting..."<<endl;
            exit(-1);
        }
	}
    gin.close();
    cout<<"#docs "<<docs_.size()<<endl;

	
}



void word_corpus::write_corpus_file(string out_dir) {
	
	ofstream pout((out_dir+"/CORPUS.corpus").c_str());
	RANGE_loop(i, docs_) {
		
		pout<<docs_[i].wn_occs_.size()<<" ";
		IT_loop(deqii, itm, docs_[i].wn_occs_) {
			pout<<itm->first<<":"<<itm->second<<" ";
		}
		pout<<endl;
	}
    pout.close();
}




void word_corpus::write_partition(mapii & hard_mems, string out_dir) {
    
    map<int, DI> partition;
    IT_loop(mapii, itm, hard_mems) {
        if(partition.count(itm->first)==0) {
            DI voiddi;
            partition.insert(make_pair(itm->second, voiddi));
        }
        partition[itm->second].push_back(itm->first);
    }
    ofstream pout((out_dir+"/infomap.part").c_str());
    for(map<int, DI>::iterator itm= partition.begin(); itm!=partition.end(); itm++) {
        prints(itm->second, pout);
    }
    pout.close();
    pout.open((out_dir+"/infomap-words.part").c_str());
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


void word_corpus::write_theta_file(deque<mapid> & doc_topic,\
                                   map<int, mapid> & topic_word, \
                                   string theta_file) {
    
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
                pout1<<doc_topic[i][itm->second] * docs_[i].num_words_<<" ";
            }
            else {
                pout1<<"0 ";
            }
        }
        pout1<<endl;
    }
    pout1.close();
}


void word_corpus::write_short_beta_and_theta_files(deque<mapid> & doc_topic,
                                                   map<int, mapid> & topic_word, \
                                                   string theta_file, 
                                                   string beta_file, 
                                                   string beta_file_short, \
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
    
    // theta file in short format
    ofstream pout1(theta_file.c_str());
    RANGE_loop(i, doc_topic) {
        IT_loop(mapid, itm, doc_topic[i]) {
            pout1<<itm->first<<" "<<itm->second<<" ";
        }
        pout1<<endl;
    }
    pout1.close();
    
    // gammas[doc][topic] is how many words in doc come from topic
    deque<DD> gammas;
    RANGE_loop(i, doc_topic) {
        DD vs;
        vs.assign(topic_names.size(), 0.);
        IT_loop(mapid, itm, doc_topic[i]) {
            vs[itm->first]=itm->second * docs_[i].num_words_;
        }
        gammas.push_back(vs);
    }
    
    DD ptopic;
    get_ptopic_distr(ptopic, gammas);
    // writing topics in sparse format
    print_topic_sparse_format_complete(topic_word, beta_file, beta_file_short, \
                                       word_strings_, ptopic, 100);

}



