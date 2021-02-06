#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include "head/queue.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/bnet.h"

using namespace std;
extern int sample_num;

/*
functions in this file:
1. Graph::Graph(){}
2. Graph::~Graph()
3. void Graph::init(ifstream &fin, int &N, int &PO, int ite)
4. void Graph::addelement(int ni1, int ni2, int no, string ipattern1, string ipattern2, int &indEdge)
5. void Graph::traverse(int flag)
6. void Graph::topSort(vector<int> &sort_list)
7. double Graph::COP(int target, vector<int> &sort_list)
8. void Graph::WAA(vector<int> &sort_list)
9. double Graph::findmin(double p)
10. void Graph::sort_cutnodes()
11. void Graph::MFFC(vector<int> &sort_list)
*/

//Global variables and external variables


/*topSort()*/
void topSort(BnetNetwork *net, vector<char*> &sort_list)
{
	BnetNode *nd, *auxnd, *tmp;
    queue sort_queue;
    map<char*, int> indegree_list; 
    map<char*, int>::iterator itrm_si;
    
    nd = net->nodes;
    while (nd != NULL) 
    {
//    	cout << "cnode: " << nd->name << endl;
    	indegree_list.insert(pair<char*, int>(nd->name, nd->ninp));
    	if(nd->ninp == 0)
    		sort_queue.push(nd->name);
		nd = nd->next;
    }
    
    while(!sort_queue.empty())
    {
    	sort_queue.traverse();
        char *pnode = sort_queue.pop(); 
//        cout << endl << "pop: " << pnode << endl;
        st_lookup(net->hash, pnode, &tmp);
        sort_list.push_back(tmp->name);
        for(int i = 0; i < tmp->nfo; i++)
        {        	
            char *outnode = tmp->fanouts[i];
            st_lookup(net->hash, outnode, &auxnd);
//            cout << "output: " << auxnd->name;
            itrm_si = indegree_list.find(auxnd->name);
            if(itrm_si == indegree_list.end())
            {
            	cout << " not exist!";
            	exit(1);
            }
            else 
            {
            	itrm_si->second = itrm_si->second - 1;
//            	cout << ", num: " << itrm_si->second << endl;
            	if(itrm_si->second == 0)
            	{
                	sort_queue.push(auxnd->name);
//                	cout << "node " << auxnd->name << " is enqueued. " << endl;
            	}
            }
        }         
    }     
}



/*get_MFFC()*/
void get_MFFC(BnetNetwork *net, vector<char*> &sort_list, map<char*, set<char*> > &TFI_set_char, map<char*, set<char*> > &MFFC_set)
{
	set<char*>::iterator itrs, itrs1, itrs2;
	BnetNode *nd, *tmp, *auxnd;	
	map<char*, set<char*> >::iterator itrm_cs, itrm_cs1;	
	map<char*, map<int, char*> >::iterator itrm_cm, itrm_cm1;
	map<char*, map<int, char*> > TFI_set;
	map<int, char*> ind_sort;
	map<char*, int> sort_ind; 
	map<int, char*>::iterator itrm_ic, itrm_ic1, itrm_ic2, itrm_ic3;
	map<char*, int>::iterator itrm_ci;
	for(int i = 0; i < sort_list.size(); i++)
	{
		ind_sort.insert(pair<int, char*>(i, sort_list[i]));
		sort_ind.insert(pair<char*, int>(sort_list[i], i));
	}
/*	cout << "ind_sort" << endl;
	for(itrm_ic = ind_sort.begin(); itrm_ic != ind_sort.end(); itrm_ic++)
		cout << itrm_ic->first << ", " << itrm_ic->second << endl;
*/	
	map<int, char*> this_TFI;
	map<int, char*> pred_TFI;
	set<char*> this_TFI_set;
	for(itrm_ic = ind_sort.begin(); itrm_ic != ind_sort.end(); itrm_ic++)
    {
        char *cnode = itrm_ic->second;
        st_lookup(net->hash, cnode, &nd);
		int cind = itrm_ic->first;
//		cout << "cnode = " << cnode << endl;
		
        this_TFI.clear();
        this_TFI.insert(pair<int, char*>(cind, nd->name));
        for(int j = 0; j < nd->ninp; j++)
        {
        	char *innode = nd->inputs[j];
        	st_lookup(net->hash, innode, &tmp);  	
        	itrm_cm = TFI_set.find(tmp->name);
        	if(itrm_cm == TFI_set.end())
        	{
        		cout << "this innode is not in TFI_set" << endl;
        		exit(1);
        	}
        	pred_TFI = itrm_cm->second;
			for(itrm_ic1 = pred_TFI.begin(); itrm_ic1 != pred_TFI.end(); itrm_ic1++)
        	{
        		st_lookup(net->hash, itrm_ic1->second, &auxnd);  
        		this_TFI.insert(pair<int, char*>(itrm_ic1->first, auxnd->name));
        	}
        }
        TFI_set.insert(pair<char*, map<int, char*> >(cnode, this_TFI));
        this_TFI_set.clear();
        for(itrm_ic1 = this_TFI.begin(); itrm_ic1 != this_TFI.end(); itrm_ic1++)
        {
        	char *node = itrm_ic1->second;
        	st_lookup(net->hash, node, &tmp);  	
        	this_TFI_set.insert(tmp->name);
        }
        MFFC_set.insert(pair<char*, set<char*> >(nd->name, this_TFI_set));
        TFI_set_char.insert(pair<char*, set<char*> >(nd->name, this_TFI_set));
    }
    
    
    for(int i = 0; i < sort_list.size(); i++)
    {
        char *cnode = sort_list[i];
//        cout << endl << "cnode: " << cnode << endl;
        st_lookup(net->hash, cnode, &nd);
        itrm_cm = TFI_set.find(nd->name);
        map<int, char*> current_MFFC = itrm_cm->second;
        map<int, char*> delete_list;
        itrm_ic = current_MFFC.end();
        itrm_ic--;
        if(current_MFFC.size() > 1)
	    {
	        while(1)
	        {
	        	char *mnode = itrm_ic->second;        	
	        	if(!strcmp(cnode, mnode))
	        	{
	        		itrm_ic--;
	        		if(itrm_ic == current_MFFC.begin())
	        			break;
	        		else
	        			continue;
	        	}
//	        	cout << "checking node " << mnode << endl;
	        	st_lookup(net->hash, mnode, &tmp);
	        	if(tmp->type == BNET_OUTPUT_NODE)
	        	{
//	        		cout << "node " << mnode << " is to be deleted." << endl;
					itrm_cm1 = TFI_set.find(tmp->name);
					if(itrm_cm1 == TFI_set.end())
					{
						cout << "current node " << mnode << " is not in TFI_set" << endl;
						exit(1);
					}
	        		map<int, char*> delete_tfi = itrm_cm1->second;
	        		for(itrm_ic3 = delete_tfi.begin(); itrm_ic3 != delete_tfi.end(); itrm_ic3++)
	        			delete_list.insert(pair<int, char*>(itrm_ic3->first, itrm_ic3->second));	        		        			
	        		itrm_ic--;
	        		if(itrm_ic == current_MFFC.begin())
	        			break;		
					else continue;        
	        	}
	        	for(int j = 0; j < tmp->nfo; j++)
	        	{
	        		char *outnode = tmp->fanouts[j];
	        		st_lookup(net->hash, outnode, &auxnd);
	        		itrm_ci = sort_ind.find(auxnd->name);
	        		int ind = itrm_ci->second;
	        		itrm_ic1 = current_MFFC.find(ind);
	        		itrm_ic2 = delete_list.find(ind);
	        		if((itrm_ic1 == current_MFFC.end()) || (itrm_ic2 != delete_list.end()))
	        		{
//	        			cout << "node " << mnode << " is to be deleted." << endl;
						itrm_cm1 = TFI_set.find(tmp->name);
						if(itrm_cm1 == TFI_set.end())
						{
							cout << "current node " << outnode << " is not in TFI_set" << endl;
							exit(1);
						}
	        			map<int, char*> delete_tfi = itrm_cm1->second;
	        			for(itrm_ic3 = delete_tfi.begin(); itrm_ic3 != delete_tfi.end(); itrm_ic3++)
	        				delete_list.insert(pair<int, char*>(itrm_ic3->first, itrm_ic3->second));	        			
						break;
	        		}	        		
	        	}
	        	itrm_ic--;
	        	if(itrm_ic == current_MFFC.begin())
	        		break;
	        }//while(1)
	    }//if(current_MFFC.size() > 1)
	    
        //deal with the first element
        itrm_ic = current_MFFC.begin();
        char *mnode = itrm_ic->second;        	
        if(strcmp(cnode, mnode))
		{
//			cout << "checking node " << mnode << endl;
	        st_lookup(net->hash, mnode, &tmp);
	        for(int j = 0; j < tmp->nfo; j++)
	        {
	        	char *outnode = tmp->fanouts[j];
//	        	cout << "outnode " << outnode << endl;
	        	st_lookup(net->hash, outnode, &auxnd);
	        	itrm_ci = sort_ind.find(auxnd->name);
	        	int ind = itrm_ci->second;
	        	itrm_ic1 = current_MFFC.find(ind);
	        	itrm_ic2 = delete_list.find(ind);
	        	if((itrm_ic1 == current_MFFC.end()) || (itrm_ic2 != delete_list.end()))
	        	{
//	        		cout << "node " << mnode << " is to be deleted." << endl;
					itrm_cm1 = TFI_set.find(tmp->name);
					if(itrm_cm1 == TFI_set.end())
					{
						cout << "current fanout node " << outnode << " is not in TFI_set" << endl;
						exit(1);
					}
	        		map<int, char*> delete_tfi = itrm_cm1->second;
	        		for(itrm_ic3 = delete_tfi.begin(); itrm_ic3 != delete_tfi.end(); itrm_ic3++)
	        			delete_list.insert(pair<int, char*>(itrm_ic3->first, itrm_ic3->second));
/*	        		cout << "delete_list: " << endl;
					for(itrm_ic3 = delete_list.begin(); itrm_ic3 != delete_list.end(); itrm_ic3++)
						cout << itrm_ic3->second << " ";
					cout << endl;	
*/
					break;      
	        	}	        	
	        }
	    }

        itrm_cs1 = MFFC_set.find(cnode);
        for(itrm_ic = delete_list.begin(); itrm_ic != delete_list.end(); itrm_ic++)
        {
        	char *dnode = itrm_ic->second;
            itrs = itrm_cs1->second.find(dnode);
            if(itrs != itrm_cs1->second.end())
                itrm_cs1->second.erase(itrs);
        }
        
/*        cout << "final MFFC "<< endl;
        set<char*> this_MFFC = itrm_cs1->second;
        for(itrs = this_MFFC.begin(); itrs != this_MFFC.end(); itrs++)
        	cout << *itrs << " ";
        cout << endl;
*/
  	}//for(int i = 0; i < sort_list.size(); i++)      
            
}


/*findmin()*/
double findmin(double p)
{
	if (p <= 0.5) return p;
	else return (1-p);
}

/*sort()*/
void asc_sort(BnetNetwork *net, vector<char*> &cutnodes)
{

	BnetNode *tmp1, *tmp2, *tmp3;
	char *temp1, *temp2;
	int i, j;
	for(i = 1; i < cutnodes.size(); i++) // Insertion Sort
	{
	    temp1 = cutnodes[i];
	    st_lookup(net->hash, temp1, &tmp1);
	    j = i - 1;
	    temp2 = cutnodes[j];
	    st_lookup(net->hash, temp2, &tmp2);
	    double sp_v0 = tmp1->p;
	    double sp_v1 = tmp2->p;
	    while(j >= 0 && (sp_v0 < sp_v1))
	    {	    
		    strcpy(cutnodes[j+1], cutnodes[j]);
		    j = j-1;		 
		    if(j >= 0) 
		    {
		    	st_lookup(net->hash, cutnodes[j], &tmp3);
		    	sp_v1 = tmp3->p;  
		    }
		    else
		    	break;
	    }
	    strcpy(cutnodes[j+1], tmp1->name);
	}

}


//permute()
int permute(vector<string> &dont_care, vector<string> &insig_string, vector<char*> &cutnodes)
{
//	cout << endl << "Coming into permute!" << endl;
	map<string, int> insig_map;
	map<string, int>::iterator itrm_si;
	map<char*, int> insig_org_map;
	map<char*, int>::iterator itrm_ci;
	
	
/*	cout << "insig_string: " << endl;
	for(int i = 0; i < insig_string.size(); i++)
		cout << insig_string[i] << " ";
	cout << endl;
	cout << "unsort_cutnodes: " << endl;
	for(int i = 0; i < cutnodes.size(); i++)
		cout << cutnodes[i] << " ";
	cout << endl;
*/	
	int index;
	for(int i = 0; i < insig_string.size(); i++)
	{
		string str = insig_string[i];
		if(str == "(null)")
			return 1;
		insig_map.insert(pair<string, int>(str, i));
//		cout << str << ": " << i << endl;		
	}
	
	vector<int> sort_order;
	for(int i = 0; i < cutnodes.size(); i++)
	{
		string str(cutnodes[i]);
		itrm_si = insig_map.find(str);
		if(itrm_si == insig_map.end())
			return 1;
		int index = itrm_si->second;
		sort_order.push_back(index);
	}
	
	vector<string> dont_care_p1;
	for(int i = 0; i < dont_care.size(); i++)
	{
		string dc = dont_care[i];
		string dc_p1;
		for(int j = 0; j < dc.size(); j++)
		{
			int order = sort_order[j];
//			cout << "order = " << order << endl;
			char ch = dc[order];
			dc_p1.append(1, ch);
		}
		dont_care_p1.push_back(dc_p1);
	}

	dont_care = dont_care_p1;
	
	return 0;
}


//permute()
void permute_v2(vector<string> &dont_care, vector<char*> &insig_string, vector<char*> &cutnodes)
{
//	cout << endl << "Coming into permute!" << endl;
	map<string, int> insig_map;
	map<string, int>::iterator itrm_si;
	map<char*, int> insig_org_map;
	map<char*, int>::iterator itrm_ci;
	
	
	int index;
	for(int i = 0; i < insig_string.size(); i++)
	{
		char *name = insig_string[i];
		string str(name);
		if(str == "(null)")
			return;
		insig_map.insert(pair<string, int>(str, i));
//		cout << str << ": " << i << endl;		
	}
	
	vector<int> sort_order;
	for(int i = 0; i < cutnodes.size(); i++)
	{
		string str(cutnodes[i]);
		itrm_si = insig_map.find(str);
		int index = itrm_si->second;
		sort_order.push_back(index);
	}

	
	vector<string> dont_care_p1;
	for(int i = 0; i < dont_care.size(); i++)
	{
		string dc = dont_care[i];
		string dc_p1;
		for(int j = 0; j < dc.size(); j++)
		{
			int order = sort_order[j];
//			cout << "order = " << order << endl;
			char ch = dc[order];
			dc_p1.append(1, ch);
		}
		dont_care_p1.push_back(dc_p1);
	}

	dont_care = dont_care_p1;
}




/*topSort()*/
void topSort(BnetNetwork **net, vector<char*> &sort_list)
{
	BnetNode *nd, *auxnd, *tmp;
    queue sort_queue;
    map<char*, int> indegree_list; 
    map<char*, int>::iterator itrm_si;
    
    nd = (*net)->nodes;
    int i = 0;
    string name_str;
    while (nd != NULL) 
    {
    	indegree_list.insert(pair<char*, int>(nd->name, nd->ninp));
    	if(nd->ninp == 0)
    		sort_queue.push(nd->name);
		nd = nd->next;
    }
    
    while(!sort_queue.empty())
    {
    	sort_queue.traverse();
        char *pnode = sort_queue.pop(); 
//        cout << endl << "pop: " << pnode << endl;
        st_lookup((*net)->hash, pnode, &tmp);
        sort_list.push_back(tmp->name);
        for(int i = 0; i < tmp->nfo; i++)
        {        	
            char *outnode = tmp->fanouts[i];
            st_lookup((*net)->hash, outnode, &auxnd);
//            cout << "output: " << auxnd->name;
            itrm_si = indegree_list.find(auxnd->name);
            if(itrm_si == indegree_list.end())
            {
            	cout << " not exist!";
            	exit(1);
            }
            else 
            {
            	itrm_si->second = itrm_si->second - 1;
//            	cout << ", num: " << itrm_si->second << endl;
            	if(itrm_si->second == 0)
            	{
                	sort_queue.push(auxnd->name);
//                	cout << "node " << auxnd->name << " is enqueued. " << endl;
            	}
            }
        }         
    }     
}





/*6. update_ckt()*/
void update_ckt(BnetNetwork *net, DdManager *dd)
{
	cout << "update_ckt : " << endl;
	
	ifstream fin;
	ofstream fout;
	FILE *fp;
    string str;
	BnetNode *nd;
	BnetTabline *f;
	
	//0. delete "_out" in ckt_sim.blif
	char com[100];	
    sprintf(com, "sed -i 's/_out/ /g' ./blif_files/ckt_sim.blif");
    system(com);
    

    //1. read ckt_sim.blif and store it in net   
    fp = fopen("./blif_files/ckt_sim.blif", "r");
    net = Bnet_ReadNetwork(fp);
    fclose(fp);

    
    //2. write net to ckt_sim_buffer.blif: add buffers to output nodes    
    fout.open("./blif_files/ckt_sim_buffer.blif", iostream::out);
    fout << ".model ckt_sim" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    for(int i = 0; i < net->npos; i++)
    {
        fout << net->outputs[i] <<  "_out ";
    }
    fout  << endl;        
    nd = net->nodes;
    while(nd != NULL)
    {
        if(nd->type == BNET_INPUT_NODE)
        {
        	nd = nd->next;
            continue;
        }
        fout << ".names ";
        for(int i = 0; i < nd->ninp; i++)
            fout << nd->inputs[i] << " ";
        fout << nd->name << endl;
        f = nd->f;
        while(f != NULL)
        {
        	if (f->values != NULL) 
		    	fout << f->values << " " << 1 - nd->polarity << endl;
			else 
			    fout <<  1 - nd->polarity << endl;
        	f = f->next;
        }
        nd = nd->next;
    }    
    fout  << endl;
    for(int i = 0; i < net->npos; i++)
    {
    	fout << ".names " << net->outputs[i] <<  " " << net->outputs[i] << "_out" << endl;
    	fout << "1 1" << endl;
    }        
    fout << ".end" << endl;
    fout.close();
    
    Bnet_FreeNetwork_Bdd(net, dd);
           
    sprintf(com, "rm -rf ckt_sim.blif");
    system(com);
	sprintf(com, "cp ./blif_files/ckt_sim_buffer.blif ./blif_files/ckt_sim.blif");
    system(com);
        
    //read ckt_sim_buffer and store it in net    
    fp = fopen("./blif_files/ckt_sim.blif", "r");
    net = Bnet_ReadNetwork(fp);
    fclose(fp);
    
/*    cout << "in update: " << endl;
    Bnet_PrintNetwork(net);
*/    
//    net = net_new;
  
}


double comp_pdiff_simu(map<string, string> &node_signatures, string &snd1, string &snd2)
{
	map<string, string>::iterator itrm_ss; 
	
	itrm_ss = node_signatures.find(snd1);
	string nd1_sig = itrm_ss->second;
	
	itrm_ss = node_signatures.find(snd2);
	string nd2_sig = itrm_ss->second;
	
	int num_diff = 0;
	for(int i = 0; i < nd1_sig.size(); i++)
	{
		char c1 = nd1_sig[i];
		char c2 = nd2_sig[i];
		if(c1 != c2)
			num_diff++;
	}
	double pdiff = num_diff/(double)sample_num;
	return pdiff;
}
