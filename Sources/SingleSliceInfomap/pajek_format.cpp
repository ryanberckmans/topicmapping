

# include "../standard_package/standard_include.cpp"

int main(int argc, char * argv[]) {
	
	
    if(argc<2) {
        
        cerr<<argv[0]<<" [edge_list_file]"<<endl;
        cerr<<"output is pajek.net"<<endl;
        return -1;
    }
        
    // old id -> new id
    mapii labels;
    //new id -> old id
    DI new_labels;
    // {(node1, node2):weight}
    map<pair<int, int>, double> edges;
    
    string gins;
    ifstream gin(argv[1]);
    
    while (getline(gin, gins)) {
        
        DD ds;
        cast_string_to_doubles(gins, ds);

        if(ds.size()>=2) {
            
            int node1=cast_int(ds[0]);
            int node2=cast_int(ds[1]);
            
            double weight=1.;
            if (ds.size()>2)
                weight=ds[2];
                
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
        }
        
        if(edges.size()%1000000==0) {
            cout<<"read "<<edges.size()/1000000<<" million links"<<endl;
        }
                
    }
    
    assert_ints(labels.size(), new_labels.size(), "labels old and new have diff. sizes");
     
    ofstream pout("pajek.net"); 
    
    pout<<"*Vertices "<<labels.size()<<endl;
    RANGE_loop(i, new_labels) {
        pout<<i+1<<" \""<<new_labels[i]<<"\""<<endl;
    }
    pout<<"*Edges"<<endl;
    for(map<pair<int, int>, double>::iterator itm=edges.begin(); itm!=edges.end(); itm++) {
        pout<<itm->first.first+1<<" "<<itm->first.second+1<<" "<<itm->second<<endl;
    }
    pout.close();
    
    return 0;
}

