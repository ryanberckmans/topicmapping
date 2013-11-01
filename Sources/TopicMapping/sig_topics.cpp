



void word_corpus::null_model(DI & links1, DI & links2, DD & weights) {
	
    
    bool verbose=true;
    
	if(verbose) cout<<"null model. p-value: "<<p_value_<<endl;
    // this is not efficient
	// EFFICIENCY-ALERT
    
    int_matrix original_corpus;
	RANGE_loop(i, docs_) {
		DI new_word_list;
		IT_loop(mapii, itm, docs_[i].wn_occurences_) {
			for(int occ=0; occ<itm->second; occ++) {
				new_word_list.push_back(itm->first);
            }
		}
		original_corpus.push_back(new_word_list);
	}
    
	degree_block DB;
	DB.random_similarities_poissonian_set_data(original_corpus);
	    
	map<pair<int, int> , int> cooc;
	map<pair<int, int> , int> cooc_nullterm;
	dotpr_similarity_of_connected_words(original_corpus, cooc);
    
    
	int qfive_links=0;
	int total_possible_links=0;
	
    links1.clear(); links2.clear(); weights.clear();
    // computing the significant links and writing the first order network
	for(map<pair<int, int> , int>::iterator it=cooc.begin(); it!=cooc.end(); it++) {
		
		++total_possible_links;
		
		int k1=min(word_occurrences_[it->first.first], word_occurrences_[it->first.second]);
		int k2=max(word_occurrences_[it->first.first], word_occurrences_[it->first.second]);
		int qfive=DB.qfive(k1,k2,p_value_);        
                
        // significant links
        if(it->second- qfive>0) {
        
            links1.push_back(it->first.first);
            links2.push_back(it->first.second);
            weights.push_back(double(it->second- qfive));
            
            ++qfive_links;
            cooc_nullterm.insert(make_pair(it->first, qfive));
		}
	}
    
    
	if(verbose) {
        cout<<"p-value: "<<p_value_<<" significant link fraction: "<<double(qfive_links)/total_possible_links<<endl;
        cout<<"total possible links "<<total_possible_links<<endl;
    }
}








