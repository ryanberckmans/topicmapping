
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
