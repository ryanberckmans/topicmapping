


#include <stdlib.h>     /* atoi */


bool separate_strings_tree_file(string &b, deque<string> & v) {		
    	
	v.clear();
	string s1;
	
	for(int i=0; i<int(b.size()); i++) {
		
		if(b[i]==' ' || b[i]=='\t' || b[i]=='\n' || b[i]==',' || b[i]==':' || b[i]=='"') {
			if(s1.size()>0)
				v.push_back(s1);
			s1.clear();
            
		} else
			s1.push_back(b[i]);
		
		if(i==int(b.size())-1) {
        
			if(s1.size()>0)
				v.push_back(s1);			
			s1.clear();
            
		}
	}
    
	return true;
}




void get_the_partition(string tree_file, mapii & hard_memberships) {
    
    
    hard_memberships.clear();
    
    // getting the partition
    ifstream tree_in(tree_file.c_str());
    
    string gins;
    while(getline(tree_in, gins)) if(gins.size()>0 and gins[0]!='#') {
        deque<string> vss;
        separate_strings_tree_file(gins,  vss);
        if(vss.size()!=4) {
            cerr<<"error in tree file"<<endl;
            exit(-1);
        }
        int node_from_tree_file= atoi(vss[3].c_str());
        int cluster= atoi(vss[0].c_str());
        hard_memberships[node_from_tree_file] = cluster;
    }
    
    tree_in.close();


}

double get_infomap_partition_from_edge_list(int Ntrials, int random_seed, \
                                            const deque<int> & links1, const deque<int> & links2,
                                            const deque<double> & weights,
                                            mapii & hard_memberships, bool verbose, string out_dir) {

    //
    // same function as below
    // but the input is the edge list
    // links1, links2, weights
    // links1 and link2 can start from 
    //
    
    hard_memberships.clear();
    
    // old id -> new id
    mapii labels;
    //new id -> old id
    DI new_labels;
    // {(node1, node2):weight}
    map<pair<int, int>, double> edges;
    
    RANGE_loop(li, links1) {
                
        int node1=links1[li];
        int node2=links2[li];
        double weight=weights[li];
        
        if (labels.count(node1)==0) {
            labels.insert(make_pair(node1, new_labels.size()));
            new_labels.push_back(node1);
        }
        
        if (labels.count(node2)==0) {
            labels.insert(make_pair(node2, new_labels.size()));
            new_labels.push_back(node2);
        }
        
        node1=labels[node1];
        node2=labels[node2];
        
        pair<int, int> link(min(node1,node2), max(node1, node2));
        if (edges.count(link)==0) {
            edges[link]=weight;
        } else {
            edges[link]+=weight;
        }
        
        if(edges.size()%1000000==0) {
            cout<<"read "<<edges.size()/1000000<<" million links"<<endl;
        }
        
    }
    
    assert_ints(labels.size(), new_labels.size(), "labels old and new have diff. sizes");
    
    
    cout<<"writing word graph to sig_words.net"<<endl;
    ofstream sigout((out_dir+"/sig_words.net").c_str());
    sigout<<"*Vertices "<<labels.size()<<endl;
    RANGE_loop(i, new_labels) {
        sigout<<i+1<<" \""<<new_labels[i]<<"\""<<endl;
    }
    sigout<<"*Edges"<<endl;
    for(map<pair<int, int>, double>::iterator itm=edges.begin(); itm!=edges.end(); itm++) {
        sigout<<itm->first.first+1<<" "<<itm->first.second+1<<" "<<itm->second<<endl;
    }
    sigout.close();
    
    
    char num_trials_ch[200];
    sprintf(num_trials_ch, "  --num-trials %d ", Ntrials);
    char seed_ch[200];
    sprintf(seed_ch, " --seed %d ", random_seed);
    string option_file= " "+out_dir+"/sig_words.net "+out_dir+" --two-level --undirected  ";
    string option_log= " > "+out_dir+"/infomap.log";
    string command_line = INFOMAP_PATH + option_file + string(num_trials_ch) + string(seed_ch) + option_log;
    
    // running the code
    cout<<"Running Infomap"<<endl;
    int sy=system(command_line.c_str());
    cout<<"Infomap's call returned:: "<<sy<<" "<<endl;
    
    get_the_partition(out_dir+"/sig_words.tree", hard_memberships);
    
    return 0.;
}




