

void doc::set_from_string(string s, mapsi & word_number_all, mapis & number_word_all) {
	
	stringstream ss(s);
    deque<string> tokens;
	string buf;
    while (ss >> buf)
        tokens.push_back(buf);
	mapii wn_occurences;
	RANGE_loop(i, tokens) {
		update_map_wordnum(word_number_all, number_word_all, tokens[i]);
		int_histogram(word_number_all[tokens[i]], wn_occurences);		
	}
    
    IT_loop(mapii, itm, wn_occurences) {
        // storing data in deqii
        wn_occs_.push_back(*itm);
    }
	num_words_=tokens.size();
	
}


void doc::set_from_string_given_wn(string s, mapsi & word_wn_all) {
    
    stringstream ss(s);
    deque<string> tokens;
	string buf;
    while(ss >> buf)
        tokens.push_back(buf);
	mapii wn_occurences;
    int num_words=0;
	RANGE_loop(i, tokens) {
        // skips words which are not in word_wn_all (not in training set)
        if(word_wn_all.count(tokens[i])>0) {
            int_histogram(word_wn_all[tokens[i]], wn_occurences);
            ++num_words;	
        }
	}
    IT_loop(mapii, itm, wn_occurences) {
        // storing data in deqii
        wn_occs_.push_back(*itm);
    }
	num_words_=num_words;


}







