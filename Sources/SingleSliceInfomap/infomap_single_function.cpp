


#include "infomap.h"
#include <stdlib.h>     /* atoi */

unsigned stou(char *s){
    return strtoul(s,(char **)NULL,10);
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials, const double & convergence_precision);
void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip);
void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent, const double & convergence_precision);


#include "infomap_utilities.cpp"

/*
double get_infomap_partition(int Ntrials, int random_seed, \
                          map<int,map<int,double> > & Links, \
                          deque<int> & memberships, bool verbose, const double & convergence_precision ) {
    
    ////
     this function takes 
     Links
     the graph in the 
     list of neighbors format
     *** Links shoud be symmetrical
     *** and nodes should be labelled starting from zero
     and sets a partition
     which optimizes the code length
     the function returns the 
     relative gain in code length
     * /
    
    memberships.clear();
    
    char random_seed_ch[100];
    sprintf(random_seed_ch, "%d", random_seed);
    MTRand *R = new MTRand(stou(random_seed_ch));
    
    int Nnode = Links.size();    
    string *nodeNames = new string[Nnode];    
    char node_name[100];
    for(int i=0;i<Nnode;i++) {
        sprintf(node_name, "%d", i);
        string node_name_str(node_name);
        nodeNames[i]= node_name_str;
    }
    
    /////////// Partition network /////////////////////
    double totalDegree = 0.0;
    vector<double> degree(Nnode);
    Node **node = new Node*[Nnode];
    for(int i=0;i<Nnode;i++){
        node[i] = new Node(i);
        degree[i] = 0.0;
    }
    
    int NselfLinks = 0;
    for(map<int,map<int,double> >::iterator fromLink_it = Links.begin(); fromLink_it != Links.end(); fromLink_it++){
        for(map<int,double>::iterator toLink_it = fromLink_it->second.begin(); toLink_it != fromLink_it->second.end(); toLink_it++){
            
            int from = fromLink_it->first;
            int to = toLink_it->first;
            double weight = toLink_it->second;
            if(weight > 0.0){
                if(from == to){
                    NselfLinks++;
                }
                else{
                    node[from]->links.push_back(make_pair(to,weight));
                    node[to]->links.push_back(make_pair(from,weight));
                    totalDegree += 2*weight;
                    degree[from] += weight;
                    degree[to] += weight;
                }
            }
        }
    }
   
    
    //Swap maps to free memory
    for(map<int,map<int,double> >::iterator it = Links.begin(); it != Links.end(); it++)
        map<int,double>().swap(it->second);
    map<int,map<int,double> >().swap(Links);
    
    // Initiation
    GreedyBase* greedy;
    greedy = new Greedy(R,Nnode,totalDegree,node);
    greedy->initiate();
    
    double uncompressedCodeLength = -greedy->nodeDegree_log_nodeDegree;
    repeated_partition(R, &node, greedy, !verbose, Ntrials, convergence_precision);
    //int Nmod = greedy->Nnode;
    //cout << "Done! Code length " << greedy->codeLength << " in " << Nmod << " modules." << endl;
    //cout << "Compressed by " << 100.0*(1.0-greedy->codeLength/uncompressedCodeLength) << " percent." << endl;
    
    double relative_gain=1.0-greedy->codeLength/uncompressedCodeLength;
    
    // Order modules by size
    multimap<double,treeNode,greater<double> > treeMap;
    multimap<double,treeNode,greater<double> >::iterator it_tM;
    for(int i=0;i<greedy->Nnode;i++){
        
        int Nmembers = node[i]->members.size();
        treeNode tmp_tN;
        it_tM = treeMap.insert(make_pair(node[i]->degree/totalDegree,tmp_tN));
        for(int j=0;j<Nmembers;j++)
            it_tM->second.members.insert(make_pair(degree[node[i]->members[j]]/totalDegree,make_pair(node[i]->members[j],nodeNames[node[i]->members[j]]))); 
        
    }
    
    memberships.assign(Nnode, 0);
    int clusterNr = 0;  
    for(multimap<double,treeNode,greater<double> >::iterator mod = treeMap.begin(); mod != treeMap.end(); mod++){
        for(multimap<double,pair<int,string>,greater<double> >::iterator mem = mod->second.members.begin(); mem != mod->second.members.end(); mem++){
            memberships[mem->second.first] = clusterNr;
        }
        clusterNr++;
    }
    
    // freeing memory ------------------------
    delete [] nodeNames;
    for(int i=0;i<greedy->Nnode;i++){
        delete node[i];
    }
    delete [] node;
    delete greedy;
    delete R;
    
    return relative_gain;
    
}


*/

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
                                            mapii & hard_memberships, bool verbose,
                                            const double & convergence_precision) {
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
        hard_memberships[new_labels_old_labels[node_new_label]] = cluster;
    }
    tree_in.close();
    return 0.;
}

double get_infomap_partition_from_file(int Ntrials, int random_seed, \
                                       string filename, mapii & hard_memberships, 
                                       bool verbose, const double & convergence_precision) {
    
    
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
                            hard_memberships, verbose, convergence_precision);
}




