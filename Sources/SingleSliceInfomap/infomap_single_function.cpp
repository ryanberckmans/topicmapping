


#include "infomap.h"
#include <stdlib.h>     /* atoi */

unsigned stou(char *s){
    return strtoul(s,(char **)NULL,10);
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials, const double & convergence_precision);
void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip);
void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent, const double & convergence_precision);


#include "infomap_utilities.cpp"


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

double get_infomap_partition_from_edge_list(int Ntrials, int random_seed, \
                                            deque<int> links1, deque<int> links2,
                                            const deque<double> & weights,
                                            mapii & hard_memberships, bool verbose) {
    /*
     same function as below
     but the input is the edge list
     links1, links2, weights
     links1 and link2 can start from 
     whatever number
     links1 and link2 are reset, 
     that's why we pass them by value
     */
    
    hard_memberships.clear();

    if (links1.size()!=links2.size() or weights.size()!=links1.size()) {
        cerr<<"sizes do not match in get_infomap_partition"<<endl;
        exit(-1);
    }
    
    if(links1.size()==0) {
        return 0.;
    }
    
    
    // relabelling 
    mapii old_labels_new_labels;
    
    RANGE_loop(i, links1) {
        int n1= links1[i];
        int n2= links2[i];
        if(old_labels_new_labels.count(n1)==0) {
            old_labels_new_labels.insert(make_pair(n1, old_labels_new_labels.size()));
        }
        if(old_labels_new_labels.count(n2)==0) {
            old_labels_new_labels.insert(make_pair(n2, old_labels_new_labels.size()));
        }
        links1[i]=old_labels_new_labels[n1];
        links2[i]=old_labels_new_labels[n2];
    }
    
    //cout<<"----------------------------"<<endl;
    //prints(old_labels_new_labels);
    // getting Links in map format
    map<int, map<int, double> > Links;
    for(unsigned int i=0; i<links1.size(); i++) {
        if(Links.count(links1[i])==0) {
            map<int, double > new_map;
            Links.insert(make_pair(links1[i], new_map));
        }
        if(Links.count(links2[i])==0) {
            map<int, double > new_map;
            Links.insert(make_pair(links2[i], new_map));
        }
        if(Links[links1[i]].count(links2[i])==0) {
            Links[links1[i]][links2[i]]=0;
        }
        if(Links[links2[i]].count(links1[i])==0) {
            Links[links2[i]][links1[i]]=0;
        }
        Links[links1[i]][links2[i]]+=weights[i];
        Links[links2[i]][links1[i]]+=weights[i];

    }

    // back dictionary
    mapii new_labels_old_labels;
    IT_loop(mapii, itm, old_labels_new_labels) new_labels_old_labels[itm->second]=itm->first;
    
    
    // printing network
    ofstream tempout("tmp_net.net");
    tempout<<"*Vertices "<<new_labels_old_labels.size()<<endl;
    IT_loop(mapii, itm, new_labels_old_labels) {
        tempout<<itm->first+1<<" \""<<itm->first+1<<"\""<<endl;
    }
    tempout<<"*Edges"<<endl;
    for(map<int, map<int, double> >::iterator itm_l =  Links.begin(); itm_l!=Links.end(); itm_l++) {
        IT_loop(mapid, itm_l2, itm_l->second) {
            if(itm_l->first < itm_l2->first) {
                tempout<<itm_l->first+1<<" "<<itm_l2->first+1<<" "<<itm_l2->second<<endl;
            }
        }
    }
    tempout.close();
    
    char num_trials_ch[200];
    sprintf(num_trials_ch, "  --num-trials %d ", Ntrials);
    char seed_ch[200];
    sprintf(seed_ch, " --seed %d ", random_seed);
    string option_file(" tmp_net.net ./ --two-level --undirected  ");
    string option_log(" > infomap.log");
    string command_line = INFOMAP_PATH + option_file + string(num_trials_ch) + string(seed_ch) + option_log;
    
    // running the code
    system(command_line.c_str());
    
    // getting the partition    
    ifstream tree_in("tmp_net.tree");
    string gins;
    while(getline(tree_in, gins)) if(gins.size()>0 and gins[0]!='#') {
        deque<string> vss;
        separate_strings_tree_file(gins,  vss);
        //prints(vss);
        if(vss.size()!=4) {
            cerr<<"error in tree file"<<endl;
            exit(-1);
        }
        int node_new_label= atoi(vss[3].c_str())-1;
        int cluster= atoi(vss[0].c_str());
        hard_memberships[new_labels_old_labels.at(node_new_label)] = cluster;
    }
    tree_in.close();
    return 0.;
}

double get_infomap_partition_from_file(int Ntrials, int random_seed, \
                                       string filename, mapii & hard_memberships, 
                                       bool verbose) {
    
    ifstream gin(filename.c_str());
    string s;
    DI links1;
    DI links2;
    DD weights;
    
    while(getline(gin, s)) {
        DD link_s;
        cast_string_to_doubles(s, link_s);
        int n1= cast_int(link_s[0]);
        int n2= cast_int(link_s[1]);
        links1.push_back(n1);
        links2.push_back(n2);
        weights.push_back(link_s[2]);
    }
    
    
    return get_infomap_partition_from_edge_list(Ntrials, \
                            random_seed, links1, links2, weights, \
                            hard_memberships, verbose);
}




