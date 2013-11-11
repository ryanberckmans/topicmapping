
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
		
		pout<<docs_[i].wn_occs_.size()<<" ";
		IT_loop(deqii, itm, docs_[i].wn_occs_) {
			pout<<itm->first<<":"<<itm->second<<" ";
		}
		pout<<endl;
	}
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
                pout1<<doc_topic[i][itm->second]<<" ";
            }
            else {
                pout1<<"0 ";
            }
        }
        pout1<<endl;
    }
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
    
    // writing topics in sparse format
    print_topic_sparse_format(topic_word, beta_file, beta_file_short, word_strings_);

}



