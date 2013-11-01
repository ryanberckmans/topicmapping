

void doc::set_from_string(string s, mapsi & word_number_all, mapis & number_word_all) {
	
	stringstream ss(s);
    deque<string> tokens;
	string buf;
    while (ss >> buf)
        tokens.push_back(buf);
	wn_occurences_.clear();
	RANGE_loop(i, tokens) {
		update_map_wordnum(word_number_all, number_word_all, tokens[i]);
		int_histogram(word_number_all[tokens[i]], wn_occurences_);		
	}
	num_words_=tokens.size();
	
}


void doc::compute_thetas(int_matrix & word_partition, mapid & theta_doc, int black_module) {
    
    
    theta_doc.clear();
    SI covered_words;
    
    RANGE_loop(i, word_partition) {
        int module=word_partition[i][0];
        word_partition[i].pop_front();
        
        RANGE_loop(j, word_partition[i]) {
            
            if(wn_occurences_.count(word_partition[i][j])==0) {
                cerr<<"error in compute_thetas"<<endl;
                exit(-1);
            }
            int_histogram(module, theta_doc, double(wn_occurences_[word_partition[i][j]])/num_words_);
            covered_words.insert(word_partition[i][j]);
        }
    }
    
    // words which are not covered 
    IT_loop(mapii, itm, wn_occurences_) if(covered_words.find(itm->first)==covered_words.end()) {
        int_histogram(black_module, theta_doc, double(itm->second)/num_words_);
    }
    
    // check that the norm is one
    double norm=0.;
    IT_loop(mapid, itm, theta_doc) norm+=itm->second;
    if(fabs(norm-1)>1e-4) {
        cerr<<"error in norm"<<endl;
        exit(-1);
    }
    
}
