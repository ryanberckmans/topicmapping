



void word_corpus::null_model(DI & links1, DI & links2, DD & weights) {
	
    
    bool use_dotproduct=true;
    bool verbose=true;
    
	if(verbose) cout<<"null model. p-value: "<<p_value_<<endl;
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
		int qfive=DB.qfive(k1,k2,p_value_);
		
        
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
        cout<<"p-value: "<<p_value_<<" significant link fraction: "<<double(qfive_links)/total_possible_links<<endl;
        cout<<"total possible links "<<total_possible_links<<endl;
    }
}








