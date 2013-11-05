
void word_corpus::dotpr_similarity_of_connected_words(map<pair<int, int> , int> & cooc) {

    cooc.clear();
	RANGE_loop(doc_number, docs_) {
        
        if(doc_number%1000==0 and doc_number!=0){
            cout<<"computing connected words:: "<<doc_number<<endl;
        }
        deqii & wn_occs_doc= docs_[doc_number].wn_occs_;
        
        // loop over all pairs of documents
        RANGE_loop(j, wn_occs_doc) for(UI k=j+1; k<wn_occs_doc.size(); k++) {

            //general_assert(wn_occs_doc[j].first < wn_occs_doc[k].first, "error in word order");
            pair<int, int> pp(wn_occs_doc[j].first, wn_occs_doc[k].first);
            int value= wn_occs_doc[j].second * wn_occs_doc[k].second;
            pair<map<pair<int, int> , int>::iterator, bool> ret = cooc.insert(make_pair(pp, value));
            if(ret.second==false)
                ret.first->second+=value;
        }
    }
}


void word_corpus::null_model(DI & links1, DI & links2, DD & weights) {
	
    
    bool verbose=true;
    
	if(verbose) cout<<"null model. p-value: "<<p_value_<<endl;
    
    DI doc_sizes;
	RANGE_loop(i, docs_) {
        doc_sizes.push_back(docs_[i].num_words_);
	}
    
    // initialize with doc_sizes
	degree_block DB;
	DB.random_similarities_poissonian_set_data(doc_sizes);

	map<pair<int, int> , int> cooc;
	map<pair<int, int> , int> cooc_nullterm;
    // co-occurrences matrix
	dotpr_similarity_of_connected_words(cooc);
    
	int qfive_links=0;
	int total_possible_links=0;	
    links1.clear(); links2.clear(); weights.clear();
    
    // computing the significant links
    cout<<"pair to evaluate for significance:: "<<cooc.size()<<endl;
	for(map<pair<int, int> , int>::iterator it=cooc.begin(); it!=cooc.end(); it++) {
		
		++total_possible_links;
        if(total_possible_links%100000==0) {
            cout<<"pairs done: "<<total_possible_links<<endl;
        }
		
		int k1=min(word_occurrences_[it->first.first], word_occurrences_[it->first.second]);
		int k2=max(word_occurrences_[it->first.first], word_occurrences_[it->first.second]);
		int qfive=DB.qfive(k1,k2,p_value_);
        
        // significant links
        if(it->second - qfive>0) {
            
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



/*



void word_corpus::dotpr_similarity_of_connected_words_parallel(map<pair<int, int> , int> & cooc, \
                                                               int par_a, \
                                                               int par_b, \
                                                               int max_a, int max_b) {
    
    // 
    cooc.clear();
	RANGE_loop(doc_number, docs_) {
        
        if(doc_number%1000==0 and doc_number!=0){
            cout<<"computing connected words:: "<<doc_number<<endl;
        }
        deqii & wn_occs_doc_orig= docs_[doc_number].wn_occs_;
        
        // only considering words in par_a and par_b
        DI first_part;
        IT_loop(deqii, itm, wn_occs_doc_orig) if(itm->first%max_a==par_a)
            first_part.push_back(itm->first);
        DI second_part;
        IT_loop(deqii, itm, wn_occs_doc_orig) if(itm->first%max_b==par_b)
            second_part.push_back(itm->first);
        
        // sort out what to do for the order!
        // loop over all pairs of documents
        ////////////////////////////////
        RANGE_loop(j, partial_wn_occ) for(UI k=j+1; k<partial_wn_occ.size(); k++) {
            
            //general_assert(partial_wn_occ[j].first < partial_wn_occ[k].first, "error in word order");
            pair<int, int> pp(partial_wn_occ[j].first, partial_wn_occ[k].first);
            int value= partial_wn_occ[j].second * partial_wn_occ[k].second;
            pair<map<pair<int, int> , int>::iterator, bool> ret = cooc.insert(make_pair(pp, value));
            if(ret.second==false)
                ret.first->second+=value;
        }
    }
}


void word_corpus::null_model_parallel(DI & links1, DI & links2, DD & weights,\
                                      int par_a, int par_b, 
                                      int max_a, int max_b) {
	
    // same function as null_model
    // but it works only for 
    // wn%par_a and wn%par_b
    
    bool verbose=true;
    
	if(verbose) cout<<"null model. p-value: "<<p_value_<<endl;
    
    DI doc_sizes;
	RANGE_loop(i, docs_) {
        doc_sizes.push_back(docs_[i].num_words_);
	}
    
    // initialize with doc_sizes
	degree_block DB;
	DB.random_similarities_poissonian_set_data(doc_sizes);
    
	map<pair<int, int> , int> cooc;
	map<pair<int, int> , int> cooc_nullterm;
    // co-occurrences matrix
	dotpr_similarity_of_connected_words(cooc, par_a, par_b, max_a, max_b);
    
	int qfive_links=0;
	int total_possible_links=0;	
    links1.clear(); links2.clear(); weights.clear();
    
    // computing the significant links
    cout<<"pair to evaluate for significance:: "<<cooc.size()<<endl;
	for(map<pair<int, int> , int>::iterator it=cooc.begin(); it!=cooc.end(); it++) {
		
		++total_possible_links;
        if(total_possible_links%100000==0) {
            cout<<"pairs done: "<<total_possible_links<<endl;
        }
		
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
    
    
    char outfile[100];
    sprintf(outfile, "sig_links_%d_%d.txt", par_a, par_b);
    
    
    
}

*/



