


void word_corpus::dotpr_similarity_of_connected_words(map<pair<int, int> , int> & cooc) {

    cooc.clear();
	RANGE_loop(doc_number, docs_) {
        
        if(doc_number%1000==0 and doc_number!=0){
            cout<<"computing word co-occurrences. doc_number::"<<doc_number<<endl;
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
                                                               int par_a, int par_b, int max_ab) {
    
    // this is to compute the network using multiple cpus.
    
    cooc.clear();
	RANGE_loop(doc_number, docs_) {
        
        if(doc_number%1000==0 and doc_number!=0){
            cout<<"computing connected words... docs processed::  "<<doc_number<<endl;
        }
        deqii & wn_occs_doc_orig= docs_[doc_number].wn_occs_;
        
        // only considering words in par_a and par_b
        deqii first_part;
        IT_loop(deqii, itm, wn_occs_doc_orig) if(itm->first%max_ab==par_a)
            first_part.push_back(*itm);
        deqii second_part;
        IT_loop(deqii, itm, wn_occs_doc_orig) if(itm->first%max_ab==par_b)
            second_part.push_back(*itm);
        
        ///cout<<"---------------------"<<endl;
        ///prints_deqii(first_part);
        ///prints_deqii(second_part);
        //cout<<"---------------------"<<endl;
        
        // both first_part and second_part are sorted
        // if you are in first_part[j]
        // you only need to move from k 
        // such that second_part[k] > first_part[i]
        
        // this is the index in second_part from which we left
        UI k_zero=0;
        RANGE_loop(j, first_part) {
            int & first_word=first_part[j].first;
            // getting the right k_zero
            while(k_zero<second_part.size() and first_word >= second_part[k_zero].first) {
                ++k_zero;
            }
            ///cout<<"k_zero "<<k_zero<<endl;
            ///cout<<"first_word "<<first_word<<endl;
            // now we are sure that second_part[k_zero]>first_part[j]
            // then also second_part[k_zero+something] is > first_part[j]
            for(UI k=k_zero; k<second_part.size(); k++) {
                ///cout<<"second_word "<< second_part[k].first<<" K: "<<k<<endl;
                general_assert(first_word < second_part[k].first, "error in word order");
                pair<int, int> pp(first_word, second_part[k].first);
                int value= first_part[j].second * second_part[k].second;
                pair<map<pair<int, int> , int>::iterator, bool> ret = cooc.insert(make_pair(pp, value));
                if(ret.second==false)
                    ret.first->second+=value;
            }            
        }
    }
}



void word_corpus::null_model(DI & links1, DI & links2, DD & weights,\
                             int par_a, int par_b, int max_ab,\
                             bool print_sig_words, string out_dir) {
	
    
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
            sprintf(outfile, "sig_words_%d_%d_%d.edges", par_a, par_b, max_ab);
            cout<<"writing to file: "<<outfile<<endl;
        } else {  
            sprintf(outfile, "sig_words.edges");
        }
        string outfile_s(outfile);
        outfile_s = out_dir+"/"+outfile_s;
        ofstream sigout(outfile_s.c_str());
        RANGE_loop(i, links1)
            sigout<<links1[i]<<" "<<links2[i]<<" "<<weights[i]<<endl;
        sigout.close();
    }

}





void word_corpus::null_model(string parall_str, string out_dir) {
    
    // public interface for parallelization

    // parsing string
    DI param_parall;
    cast_string_to_doubles(parall_str, param_parall);
    assert_ints(param_parall.size(), 3, "-parall should be followed by i:j:n");
    int parall_i=param_parall[0];
    int parall_j=param_parall[1];
    int parall_n=param_parall[2];
    
    cout<<"running null model. Parallel parameters: "<<parall_i<<" "<<parall_j<<" "<<parall_n<<endl;
    general_assert(parall_i>=0 and parall_i<parall_n, "i out of range respect to n");
    general_assert(parall_j>=0 and parall_j<parall_n, "j out of range respect to n");

    DI links1;
    DI links2;
    DD weights;
    
    // default is no parallelization 
    null_model(links1, links2, weights, parall_i, parall_j, parall_n, true, out_dir);
    
    
}

