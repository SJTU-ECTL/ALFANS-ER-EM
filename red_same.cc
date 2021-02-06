#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include "head/queue.h"
#include "head/stack.h"
#include "head/basics.h"
#include "head/helper.h"
#include "head/simu_ckt_ss.h"
#include "head/sub_sim_ckt.h"
#include "head/write_func.h"


using namespace std;

extern int numPI_ini;

/*
functions in this file:

*/

//Global variables and external variables 


int red_same_helper(BnetNetwork *net, vector<char*> &sort_list, map<char*, set<char*> > &MFFC_set, map<char*, set<char*> > &TFI_set, map<string, string> &node_signatures, string &ts_best, string &ss_best, int &pol_best)
{
	BnetNode *nd, *tmp, *auxnd;
	set<char*>::iterator itrs, itrs1;
	set<string>::iterator itrs_str, itrs_str1;
	map<char*, set<char*> >::iterator itrm_cs, itrm_cs1, itrm_cs2;
	map<string, set<char*> > node_supp;
	map<string, set<char*> >::iterator itrm_ss;
	multimap<set<string>, string> supp_node;	
	typedef multimap<set<string>, string>::iterator ITE;
	ITE itrmm_ss;		
	
	multimap<int, char*> ts_cand_top;
	multimap<int, char*>::iterator itrmm_ic, itrmm_ic1;	
	for(int i = 0; i < sort_list.size(); i++)
	{
		char *csig = sort_list[i];
		st_lookup(net->hash, csig, &nd);
    	if(nd->type == BNET_INPUT_NODE)
    	{
    		nd = nd->next;
    		continue;
    	}
    	itrm_cs = MFFC_set.find(nd->name);
		set<char*> this_MFFC = itrm_cs->second;
		int indlog_size = this_MFFC.size(); 
    	ts_cand_top.insert(pair<int, char*>(indlog_size, nd->name));
    }
    
    int max_score = -1;
	int flag = 0;
	string ts, ss;
    for(itrmm_ic1 = ts_cand_top.begin(); itrmm_ic1 != ts_cand_top.end(); itrmm_ic1++)
	{	
        char *cnode = itrmm_ic1->second;
        string nd1_str(cnode);
        st_lookup(net->hash, cnode, &nd);
        cout << endl << "*current node = " << nd1_str << ", indlog_size = " << itrmm_ic1->first << endl;
        
        //obtain current node's pi_supp
        set<char*> tfi;
        set<string> pi_supp;
        itrm_cs = TFI_set.find(nd->name); 
        tfi = itrm_cs->second;
        for(itrs = tfi.begin(); itrs != tfi.end(); itrs++)
        {
        	char *tfi_node = *itrs;
        	st_lookup(net->hash, tfi_node, &tmp);
        	if(tmp->type == BNET_INPUT_NODE)
        	{
        		string snode(tfi_node);
        		pi_supp.insert(snode);
        	}
        }    
        itrmm_ss = supp_node.find(pi_supp);
        if(itrmm_ss == supp_node.end())
        {
        	supp_node.insert(make_pair(pi_supp, nd1_str));
        	cout << "#node " << nd1_str << " is inserted " << endl;
        	nd = nd->next;
        	cout << "--no previous nodes with same pi_supp!" << endl;
        	continue;
        }
        else
        {
        	cout << "--nodes with same pi_supp exist!" << endl;        	
        	itrm_cs1 = MFFC_set.find(nd->name);
        	set<char*> nd1_MFFC = itrm_cs1->second;
        	itrm_cs1 = TFI_set.find(nd->name);
			set<char*> nd1_tfi = itrm_cs1->second;

        	pair<ITE, ITE> ret;
        	ret = supp_node.equal_range(pi_supp);
        	int num = supp_node.count(pi_supp);
        	cout << "num_cand = " << num << endl;
        	for(itrmm_ss = ret.first; itrmm_ss != ret.second; itrmm_ss++)
        	{
        		string nd2_str = itrmm_ss->second;
        		cout << "cand_node " << nd2_str << endl;
        		char *same_node = new char[nd2_str.size()+1];
        		strcpy(same_node, nd2_str.c_str());
        		st_lookup(net->hash, same_node, &auxnd);
        		delete []same_node;
        		
        		//skip case: nd1_str = nd2_str
        		if(nd2_str == nd1_str)
        			continue;
        		//skip case: nd1_str or nd2_str contains "_ss"
        		size_t found1 = nd1_str.find("_ss");		
				if(found1 != string::npos)
					continue;
        		size_t found2 = nd2_str.find("_ss");		
				if(found2 != string::npos)
					continue;
				//skip case: nd1 has only one input which is nd2
				if(nd->ninp == 1)
				{
					char *onlyin = nd->inputs[0];
					string onlyin_str(onlyin);
					if(onlyin_str == nd2_str)
						continue;
				}
				//skip case: nd2 has only one input which is nd1
				if(auxnd->ninp == 1)
				{
					char *onlyin = auxnd->inputs[0];
					string onlyin_str(onlyin);
					if(onlyin_str == nd1_str)
						continue;
				}

        		double pdiff = comp_pdiff_simu(node_signatures, nd1_str, nd2_str);
        		cout << "pdiff = " << pdiff << endl;
        		if(pdiff == 0 || pdiff == 1)
        		{        		        			   			
        			itrm_cs2 = MFFC_set.find(auxnd->name);
        			set<char*> nd2_MFFC = itrm_cs2->second;
        			if(nd->type == BNET_OUTPUT_NODE && nd1_MFFC.size() == 1)
        				if(auxnd->type == BNET_OUTPUT_NODE && nd2_MFFC.size() == 1)
							continue;
        			
        			//Judge whether nd1 is in tfi of nd2 or nd2 is in tfi of nd1        
					itrm_cs1 = TFI_set.find(auxnd->name);
					set<char*> nd2_tfi = itrm_cs1->second;
					itrs = nd1_tfi.find(auxnd->name);
					int flag_type1= 0, flag_type2 = 0;
					if(itrs != nd1_tfi.end())
					{
						cout << "nd2 is in tfi of nd1" << endl;
						flag_type1 = 1;
					}
					else
						flag_type1 = 0;
					itrs = nd2_tfi.find(nd->name);
					if(itrs != nd2_tfi.end())
					{
						cout << "nd1 is in tfi of nd2" << endl;
						flag_type2 = 1;
					}
					else
						flag_type2 = 0;	
						
        			int num_ts, num_diff = 0, num_common = 0;
					if(flag_type1 == 0 && flag_type2 == 0)
					{
						if(nd->type == BNET_OUTPUT_NODE && nd1_MFFC.size() == 1)
						{
        					ts = nd2_str;
        					ss = nd1_str;
        					num_ts = nd2_MFFC.size();
        				}
        				else if(auxnd->type == BNET_OUTPUT_NODE && nd2_MFFC.size() == 1)
        				{        					
        					ts = nd1_str;
        					ss = nd2_str;
        					num_ts = nd1_MFFC.size();
        				}
						else
						{
							if(nd1_MFFC.size() >= nd2_MFFC.size())
        					{        					
        						ts = nd1_str;
        						ss = nd2_str;
        						num_ts = nd1_MFFC.size();
        					}
        					else
        					{
        						ts = nd2_str;
        						ss = nd1_str;
        						num_ts = nd2_MFFC.size();
        					}
        				}
					}
					else if(flag_type1 == 0 && flag_type2 == 1) //nd1 is in tfi of nd2
					{
						if(auxnd->type == BNET_OUTPUT_NODE && nd2_MFFC.size() == 1)
							continue;
						ts = nd2_str;
        				ss = nd1_str;
        				num_ts = nd2_MFFC.size();
					}
					else if(flag_type1 == 1 && flag_type2 == 0) //nd2 is in tfi of nd1
        			{
        				if(nd->type == BNET_OUTPUT_NODE && nd1_MFFC.size() == 1)
							continue;
						ts = nd1_str;
        				ss = nd2_str;
        				num_ts = nd1_MFFC.size();
					}        			
        					       			
					for(itrs = nd1_MFFC.begin(); itrs != nd1_MFFC.end(); itrs++)
					{
						char *node = *itrs;
						st_lookup(net->hash, node, &auxnd);
						itrs1 = nd2_MFFC.find(auxnd->name);
						if(itrs1 != nd2_MFFC.end())
							num_common++;
					}
					flag = 1;     
					num_diff = num_ts - num_common;
					int score = num_diff;
					cout << "score = " << score << endl;
					if(score > max_score)
					{
						ts_best = ts;
						ss_best = ss;
						cout << "ts_best = " << ts_best << ", ss_best = " << ss_best << endl;
						if(pdiff == 0)
        					pol_best = 1;
        				else
        					pol_best = 0;
						max_score = score;
					}
        		}//if(pdiff == 0 || pdiff == 1)
        		
        	}//end of for(itrmm_ss = ret.first; itrmm_ss != ret.second; itrmm_ss++)
        	   
        	supp_node.insert(make_pair(pi_supp, nd1_str));
        	cout << "#node " << nd1_str << " is inserted " << endl;  	
        }//else (itrmm_ss != supp_node.end())
        	
    }//end of for loop
    
    if(!flag)
    	return 1;
    else 
    	return 0;
}



void red_same(BnetNetwork *net)
{
	struct timeb startTime, endTime; 
	ftime(&startTime);
	BnetNode *nd, *tmp;
	map<char*, set<char*> >::iterator itrm_cs, itrm_cs1, itrm_cs2;	
	
	map<string, string> node_signatures;
	map<string, double> node_sp;
	double real_er = 0;
	string ts_best, ss_best;
	int pol_best;
	int iIndex = 0;
	while(real_er <= 0)
	{		
		//step0. run simulation
		if(iIndex == 0)
		{
			cout << "running simulation for the first iteration: " << endl;
			simu_ckt_ss(net, node_signatures, node_sp, real_er, iIndex, 0);
		}
	
		//step1. topological sort	
		cout << endl << "step1. topological sort: " << endl;
		vector<char*> sort_list;
		sort_list.clear();
		topSort(&net, sort_list);
		
		//step2. get MFFC for each signal
		cout << endl << "step2. get MFFC: " << endl;
		map<char*, set<char*> > MFFC_set;
		map<char*, set<char*> > TFI_set;
		map<char*, set<char*> >::iterator itrm_cs, itrm_cs1, itrm_cs2;	
		set<char*>::iterator itrs, itrs1;
		TFI_set.clear();
		MFFC_set.clear();		
		get_MFFC(net, sort_list, TFI_set, MFFC_set);
			
		//step3. red_same_helper
		cout << endl << "step3. red_same_helper: " << endl;
		int res = red_same_helper(net, sort_list, MFFC_set, TFI_set, node_signatures, ts_best, ss_best, pol_best);
		cout << "summary: " << endl;
		if(res)
		{
			cout << "no more same node pairs!" << endl;
			break;
		}
		else
		{
			cout << "make substitution: " << endl;
			cout << "ts_best = " << ts_best << ", ss_best = " << ss_best << endl;
			char *ts = new char[ts_best.size()+1];
			strcpy(ts, ts_best.c_str());
			st_lookup(net->hash, ts, &nd);
			itrm_cs = MFFC_set.find(nd->name);			
			set<char*> ts_MFFC = itrm_cs->second;
			
			char *ss = new char[ss_best.size()+1];
			strcpy(ss, ss_best.c_str());
			st_lookup(net->hash, ss, &nd);
			itrm_cs = TFI_set.find(nd->name);
			set<char*> ss_TFI = itrm_cs->second;
			sub_sim_ckt(net, ts_MFFC, ss_TFI, ts, ss, pol_best, iIndex);
			delete []ts;
			delete []ss;
		}
		
		cout << "compute real error rate after this substitution: " << endl;
		//sweep ckt_sim.blif
		char com[200];
		sprintf(com, "sis -t none -f ./script/sim_ckt_sweep_whole.rug > ./script/sis.txt"); 
	    	system(com);
		//report area
		cout << "mapped area: " << endl;
	    	sprintf(com, "sis -t none -f ./script/map.rug");	    
	    	system(com);	
		
	    	//write ckt_org_sim.blif	
		write_ckt_comb(net);	 //ckt_org_sim.blif

		cout << "running simulation for checking real error rate: " << endl;
		node_signatures.clear();
		node_sp.clear();
		FILE *fp = fopen("./blif_files/ckt_org_sim.blif", "r");
	    	BnetNetwork *net_comb = Bnet_ReadNetwork(fp); 
	    	fclose(fp); 
		simu_ckt_ss(net_comb, node_signatures, node_sp, real_er, iIndex, 1);
		cout << "report: real_er = " << real_er << endl; 
		
		if(real_er > 0)	
	    		break;
	    	
	    	Bnet_FreeNetwork(net);
		fp = fopen("./blif_files/ckt_sim.blif", "r");
	    	net = Bnet_ReadNetwork(fp); 
	    	fclose(fp);  

		iIndex++;
	//	if(iIndex == 1)
	//		break;
	}
	
	ftime(&endTime);
    	double rt_rem_red = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    	cout << "total runtime for rem_red: " << rt_rem_red << endl;
}
