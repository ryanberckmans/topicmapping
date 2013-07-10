
void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent, const double & convergence_precision){
    
    int Nnode = greedy->Nnode;
    Node **cpy_node = new Node*[Nnode];
    for(int i=0;i<Nnode;i++){
        cpy_node[i] = new Node();
        cpyNode(cpy_node[i],(*node)[i]);
    }
    
    int iteration = 0;
    double outer_oldCodeLength;
    do{
        outer_oldCodeLength = greedy->codeLength;
        
        if((iteration > 0) && (iteration % 2 == 0) && (greedy->Nnode > 1)){  // Partition the partition
            
            
            if(!silent)
                cout << "Iteration " << iteration+1 << ", moving " << flush;
            
            Node **rpt_node = new Node*[Nnode];
            for(int i=0;i<Nnode;i++){
                rpt_node[i] = new Node();
                cpyNode(rpt_node[i],cpy_node[i]);
            }
            vector<int> subMoveTo(Nnode);
            vector<int> moveTo(Nnode);
            int subModIndex = 0;
            
            for(int i=0;i<greedy->Nnode;i++){
                
                int sub_Nnode = (*node)[i]->members.size();
                
                if(sub_Nnode > 1){
                    
                    Node **sub_node = new Node*[sub_Nnode]; 
                    set<int> sub_mem;
                    for(int j=0;j<sub_Nnode;j++)
                        sub_mem.insert((*node)[i]->members[j]);
                    set<int>::iterator it_mem = sub_mem.begin();
                    int *sub_renumber = new int[Nnode];
                    int *sub_rev_renumber = new int[sub_Nnode];
                    double totalDegree = 0.0;
                    for(int j=0;j<sub_Nnode;j++){
                        
                        //    fprintf(stderr,"%d %d\n",j,(*it_mem));
                        int orig_nr = (*it_mem);
                        int orig_Nlinks = cpy_node[orig_nr]->links.size(); // ERROR HERE
                        sub_renumber[orig_nr] = j;
                        sub_rev_renumber[j] = orig_nr;
                        sub_node[j] = new Node(j);
                        for(int k=0;k<orig_Nlinks;k++){
                            int orig_link = cpy_node[orig_nr]->links[k].first;
                            int orig_link_newnr = sub_renumber[orig_link];
                            double orig_weight = cpy_node[orig_nr]->links[k].second;
                            if(orig_link < orig_nr){
                                if(sub_mem.find(orig_link) != sub_mem.end()){
                                    sub_node[j]->links.push_back(make_pair(orig_link_newnr,orig_weight));
                                    sub_node[orig_link_newnr]->links.push_back(make_pair(j,orig_weight));
                                    totalDegree += 2.0*orig_weight;
                                }
                            }
                        }
                        it_mem++;
                    }
                    
                    GreedyBase* sub_greedy;
                    sub_greedy = new Greedy(R,sub_Nnode,totalDegree,sub_node);
                    sub_greedy->initiate();
                    partition(R,&sub_node,sub_greedy,true, convergence_precision);
                    for(int j=0;j<sub_greedy->Nnode;j++){
                        int Nmembers = sub_node[j]->members.size();
                        for(int k=0;k<Nmembers;k++){
                            subMoveTo[sub_rev_renumber[sub_node[j]->members[k]]] = subModIndex;
                        }
                        moveTo[subModIndex] = i;
                        subModIndex++;
                        delete sub_node[j];
                    }
                    
                    delete [] sub_node;
                    delete sub_greedy;
                    delete [] sub_renumber;
                    delete [] sub_rev_renumber;
                    
                }
                else{
                    
                    subMoveTo[(*node)[i]->members[0]] = subModIndex;
                    moveTo[subModIndex] = i;
                    
                    subModIndex++;
                    
                }
            }
            
            for(int i=0;i<greedy->Nnode;i++)
                delete (*node)[i];
            delete [] (*node);
            
            greedy->Nnode = Nnode;
            greedy->Nmod = Nnode;
            greedy->node = rpt_node;
            greedy->initiate();
            greedy->determMove(subMoveTo);
            greedy->level(node,false); 
            greedy->determMove(moveTo);
            (*node) = rpt_node;
            
            outer_oldCodeLength = greedy->codeLength;
            
            if(!silent)
                cout << greedy->Nnode << " modules, looping " << flush;
            
        }
        else if(iteration > 0){
            
            if(!silent)
                cout << "Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping " << flush;
            
            
            Node **rpt_node = new Node*[Nnode];
            for(int i=0;i<Nnode;i++){
                rpt_node[i] = new Node();
                cpyNode(rpt_node[i],cpy_node[i]);
            }
            
            vector<int> moveTo(Nnode);
            for(int i=0;i<greedy->Nnode;i++){
                int Nmembers = (*node)[i]->members.size();
                for(int j=0;j<Nmembers;j++){
                    moveTo[(*node)[i]->members[j]] = i;
                }
            }
            
            for(int i=0;i<greedy->Nnode;i++)
                delete (*node)[i];
            delete [] (*node);
            
            greedy->Nnode = Nnode;
            greedy->Nmod = Nnode;
            greedy->node = rpt_node;
            greedy->initiate();
            greedy->determMove(moveTo);
            
            (*node) = rpt_node;
        }
        else{
            
            if(!silent)
                cout << "Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping " << flush;
            
        }
        
        double oldCodeLength;
        do{
            oldCodeLength = greedy->codeLength;
            bool moved = true;
            int Nloops = 0;
            int count = 0;
            while(moved){
                moved = false;
                double inner_oldCodeLength = greedy->codeLength;
                greedy->move(moved);
                Nloops++;
                count++;
                if(fabs(inner_oldCodeLength-greedy->codeLength) < convergence_precision)
                    moved = false;
                
                if(count == 10){	  
                    greedy->tune();
                    count = 0;
                }
                // 	if(!silent){
                // 	  cerr << Nloops;
                // 	  int loopsize = to_string(Nloops).length();
                // 	  for(int i=0;i<loopsize;i++)
                // 	    cerr << "\b";
                // 	}
            }
            
            greedy->level(node,true);
            
            if(!silent)
                cout << Nloops << " " << flush;
            
        } while(oldCodeLength - greedy->codeLength >  convergence_precision);
        
        iteration++;
        if(!silent)
            cout << "times between mergings to code length " <<  greedy->codeLength << " in " << greedy->Nmod << " modules." << endl;
        
    } while(outer_oldCodeLength - greedy->codeLength > convergence_precision);
    
    for(int i=0;i<Nnode;i++)
        delete cpy_node[i];
    delete [] cpy_node;
    
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent, 
                        int Ntrials, const double & convergence_precision) {
    
    double shortestCodeLength = 1000.0;
    int Nnode = greedy->Nnode;
    vector<int> cluster(Nnode);
    
    for(int trial = 0; trial<Ntrials;trial++) {
        
        if(!silent)
            cout << "Attempt " << trial+1 << "/" << Ntrials << endl;
        
        Node **cpy_node = new Node*[Nnode];
        for(int i=0;i<Nnode;i++){
            cpy_node[i] = new Node();
            cpyNode(cpy_node[i],(*node)[i]);
        }
        
        greedy->Nnode = Nnode;
        greedy->Nmod = Nnode;
        greedy->node = cpy_node;
        greedy->initiate();
        
        partition(R,&cpy_node,greedy,silent, convergence_precision);
        
        if(greedy->codeLength < shortestCodeLength){
            
            shortestCodeLength = greedy->codeLength;
            
            // Store best partition
            for(int i=0;i<greedy->Nnode;i++){
                for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
                    cluster[(*mem)] = i;
                }
            }
        }
        
        for(int i=0;i<greedy->Nnode;i++){
            delete cpy_node[i];
        }
        delete [] cpy_node;
        
    }
    
    // Commit best partition
    greedy->Nnode = Nnode;
    greedy->Nmod = Nnode;
    greedy->node = (*node);
    greedy->initiate();
    greedy->determMove(cluster);
    greedy->level(node,true);
    
}

void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip){
    
    multimap<double,treeNode,greater<double> >::iterator it;
    if(it_tM->second.nextLevel.size() > 0){
        int i=1;
        for(it = it_tM->second.nextLevel.begin(); it != it_tM->second.nextLevel.end(); it++){
            string cpy_s(s + to_string(i) + ":");
            printTree(cpy_s,it,outfile,flip);
            i++;
        }
    }
    else{
        int i = 1;
        for(multimap<double,pair<int,string>,greater<double> >::iterator mem = it_tM->second.members.begin(); mem != it_tM->second.members.end(); mem++){
            if(flip){
                string cpy_s(s + to_string(i) + " \"" + mem->second.second + "\" " + to_string(mem->first));
                (*outfile) << cpy_s << endl;
            }
            else{
                string cpy_s(s + to_string(i) + " " + to_string(mem->first) + " \"" + mem->second.second + "\"");
                (*outfile) << cpy_s << endl;
            }
            i++;
        } 
    }  
}
