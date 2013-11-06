


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





void word_corpus::dotpr_similarity_of_connected_words_parallel(map<pair<int, int> , int> & cooc, \
                                                               int par_a, \
                                                               int par_b, int max_ab) {
    
    // this is to compute the network using multiple cpus.
    /*
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
        // this is to be finished
        RANGE_loop(j, partial_wn_occ) for(UI k=j+1; k<partial_wn_occ.size(); k++) {
            
            //general_assert(partial_wn_occ[j].first < partial_wn_occ[k].first, "error in word order");
            pair<int, int> pp(partial_wn_occ[j].first, partial_wn_occ[k].first);
            int value= partial_wn_occ[j].second * partial_wn_occ[k].second;
            pair<map<pair<int, int> , int>::iterator, bool> ret = cooc.insert(make_pair(pp, value));
            if(ret.second==false)
                ret.first->second+=value;
        }
    }
     */
}



void word_corpus::null_model(DI & links1, DI & links2, DD & weights,\
                             int par_a, int par_b, int max_ab, bool print_sig_words) {
	
    
    //  this is the general function
    //  for computing the word graph
     
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
    // co-occurrences matrix
    if(max_ab==1) {
        dotpr_similarity_of_connected_words(cooc);
    } else {
        dotpr_similarity_of_connected_words_parallel(cooc, par_a, par_b, max_ab);        
    }
    
	int qfive_links=0;
	int total_possible_links=0;	
    links1.clear(); links2.clear(); weights.clear();
    
    // computing the significant links
    if(verbose) cout<<"pair to evaluate for significance:: "<<cooc.size()<<endl;
	for(map<pair<int, int> , int>::iterator it=cooc.begin(); it!=cooc.end(); it++) {
		
		++total_possible_links;
        if(verbose and total_possible_links%1000000==0) {
            cout<<"pairs done: "<<double(total_possible_links)/cooc.size()<<endl;
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
		}
	}
    
	if(verbose) {
        cout<<"p-value: "<<p_value_<<" significant link fraction: ";
        cout<<double(qfive_links)/total_possible_links<<endl;
        cout<<"total possible links "<<total_possible_links<<endl;
    }
    
    // if max_ab>1, this is always true
    if(print_sig_words or max_ab>1) {
        char outfile[200];
        if(max_ab>1) {
            sprintf(outfile, "sig_words_%d_%d_%d.edges.txt", par_a, par_b, max_ab);
        } else {  
            sprintf(outfile, "sig_words.edges.txt");
        }
        ofstream sigout(outfile);
        RANGE_loop(i, links1)
            sigout<<links1[i]<<" "<<links2[i]<<" "<<weights[i]<<endl;
    }

}





void word_corpus::null_model(string parall_str) {
    
    // public interface for parallelization

    // parsing string
    DI param_parall;
    cast_string_to_doubles(parall_str, param_parall);
    assert_ints(param_parall.size(), 3, "-parall should be followed by i:j:n");
    int parall_i=param_parall[0];
    int parall_j=param_parall[1];
    int parall_n=param_parall[2];
    general_assert(parall_i>0 and parall_i<parall_n, "i out of range respect to n");
    general_assert(parall_j>0 and parall_j<parall_n, "j out of range respect to n");

    cout<<"running null model. Parallel parameters: "<<parall_i<<" "<<parall_j<<" "<<parall_n<<endl;
    DI links1;
    DI links2;
    DD weights;
    
    // default is no parallelization 
    null_model(links1, links2, weights, parall_i, parall_j, parall_n, true);
    
    
}

