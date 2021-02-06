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
#include "head/basics.h"
#include "head/exdc_factor_new_v2.h"
#include "head/helper.h"
#include "head/read_file.h"
#include "head/write_func.h"
#include "head/loc_sim_main.h"
#include "head/simu_ckt.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/cudd_build.h"
#include "cudd/cudd_comp.h"
#include "cudd/cudd_dst.h"

using namespace std;

extern int numPI_ini;
extern double ini_threshold_er, ini_threshold_em;

/*
functions in this file:

*/

//Global variables and external variables 


/*1. find_ave_sp()*/
void find_ave_sp(BnetNetwork *net, multimap<double, char*> &ave_sp)
{
	BnetNode *nd, *tmp;

	nd = net->nodes;
	while(nd != NULL)
	{
		if(nd->ninp < 2 || nd->ninp > 15)
//		if(nd->ninp < 2)
		{
			nd = nd->next;
			continue;
		}
		double prod = 1;
//		cout << "current node: " << nd->name << endl;
		for(int i = 0; i < nd->ninp; i++)
		{
			char *innode = nd->inputs[i];
			if(!st_lookup(net->hash, innode, &tmp))
			{
				cout << "this node doesn't exist in hash!" << endl;
				exit(1);
			}
			prod = prod * findmin(tmp->rp);
		}
//		cout << "prod = " << prod << endl;
		prod = pow(prod, 1.0/nd->ninp);
		ave_sp.insert(pair<double, char*>(prod, nd->name));
//		cout << "node " << nd->name << ", ave_sp: " << prod << endl;
		nd = nd->next;
	}
}


//2. find_exdc_sim()
char *find_exdc_sim(BnetNetwork *net, DdManager **dd, BnetNetwork *net_comb, DdManager **dd_comb,  vector<string> &final_pla, map<string, struct score_pla> &sim_record, map<string, int> &internal_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int &num_output, map<string, map<string, double> > &node_pattern_rate, double &max_score, int &ave_em, double threshold_er, double threshold_em, double real_em, int &min_modified_po, int &max_modified_po, int iIndex)
{
	//variables
    struct timeb st_es, et_es;
    double total_exdc = 0;
    BnetNode *nd, *tmp, *last_nd, *auxnd, *nd1, *nd2;
    ifstream fin;
    FILE *fp;
    string str;
    
    //iterators
    multimap<double, char*>::iterator itrm_dc; 
    multimap<double, string>::iterator itrmm_ds;
    map<string, struct score_pla>::iterator itrm_ss;
    map<string, multimap<double, string> >:: iterator itrm_cm;
    set<char*>::iterator itrs;
    map<string, int>::iterator itrm_si;
   
    
    char com[100];
    multimap<double, char*> ave_sp;
    find_ave_sp(net, ave_sp);
    itrm_dc = ave_sp.begin();
    char *max_score_node = itrm_dc->second;  
    int max_score_status;          
    int index = 0;  
	        	
    //start the for loop
    multimap<double, string> score_record;
    int cur_min_modified_po, cur_max_modified_po;
    for(itrm_dc = ave_sp.begin(); itrm_dc != ave_sp.end(); itrm_dc++, index++)
    {
        char *cnode = itrm_dc->second;
        st_lookup(net->hash, cnode, &nd);

        cout << endl <<"--------------------------------------------" << endl;  
        cout << "%%check for " << cnode << ", index = " << index << ", rp = " << nd->rp << endl; 
    //    cout << "# of fanouts: " << nd->nfo << endl;
        
		//get unsort_cutnodes: input nodes of cnode
        vector<char*> unsort_cutnodes;  
        vector<struct index_flag> input_index;
        struct index_flag indf;
        for(int i = 0; i < nd->ninp; i++)
        {
        	cout << nd->inputs[i] << " " << endl;
			char *cp1 = new char[50];			
			strcpy(cp1, nd->inputs[i]);
			unsort_cutnodes.push_back(cp1);
			string str(nd->inputs[i]);
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if(tmp->type != BNET_INPUT_NODE && iIndex > 0) str.append("sim");
			itrm_si = internal_index.find(str);
			if(itrm_si == internal_index.end()) 
			{
				cout << "no such input node " << " in internal_index!" << endl;
				exit(1);
			}
			indf.index = itrm_si->second;
			if(tmp->type == BNET_INPUT_NODE) indf.flag = 1;
			else indf.flag = 2;
			input_index.push_back(indf);
        }
       	cout << endl;

  	 	//call ABC to obtain the don't cares for the current node
  	 	cout << "**********************************" << endl;
    	cout << endl << "a1. Run ABC to get don't cares!" << endl;
    	ftime(&st_es);
    	char com_abc[100]; 
    	write_abc_rug(iIndex, cnode);                                 
    	sprintf(com_abc, "abc -f ./script/abc_dc.rug > abc.txt");
    	system(com_abc);    
   	ftime(&et_es);
    	double rt_abc = ((et_es.time - st_es.time)*1000 + (et_es.millitm - st_es.millitm))/1000.0;
    	cout << "@runtime for a1. abc for don't cares: " << rt_abc << endl;    
    	vector<string> dont_care;
	string filename = "abc.txt";
    	read_abc_result(dont_care, filename); 
    	cout << "dont_care size: " << dont_care.size() << endl;  
        for (int i = 0; i < dont_care.size(); i++)
	    cout << dont_care[i] << " ";
        cout << endl;
       
        vector<string> org_pla; 
    	BnetTabline *t = nd->f;
		while(t != NULL) 
		{
			string str(t->values);
			org_pla.push_back(str);
			t = t->next;
		}
				
		//find all affected primary outputs of current node
       	nd1 = net->nodes;
    	while(nd1 != NULL)
    	{
    		nd1->visited = 0;
    		nd1 = nd1->next;
    	}
    	set<char*> po_set;
    	find_tranfanout_po(net, nd->name, po_set);
    	cout << endl << "current node " << nd->name << "'s affected pos: " << po_set.size() << endl;
    	for(itrs = po_set.begin(); itrs != po_set.end(); itrs++)
    		cout << *itrs << " ";
    	cout << endl;
    	
    	//call ABC to obtain the don't cares for the current node in terms of each affected PO
    	map<string, vector<string> > po_dont_care_map;
/*    	for(itrs = po_set.begin(); itrs != po_set.end(); itrs++)
    	{
    		char *po = *itrs;
			write_abc_rug_po(po, cnode);
			char com[100];
			sprintf(com, "echo /dev/null > abc_po.txt; abc -f ./script/abc_dc_po.rug > abc_po.txt");
    		system(com);  
    		vector<string> dont_care_po;
    		string filename = "abc_po.txt";
    		read_abc_result(dont_care_po, filename); 	
			vector<string> dont_care_po_filter;
			for(int i = 0; i < dont_care_po.size(); i++)
				if(!isIncludeVec(dont_care, dont_care_po[i]))
					dont_care_po_filter.push_back(dont_care_po[i]);
			string str_po(po);
			po_dont_care_map.insert(make_pair(str_po, dont_care_po_filter)); 
		}
*/
		
		//Find exdcs from previous iterations
		vector<string> local_dc_set;	 //empty   	
	    //Find all qualified exdcs
	    cout << "**********************************" << endl;
		cout << endl << "a2. exdc_factor_new_v2: " << endl;
		ftime(&st_es);	
		struct score_pla max_sp;
		cout << "threshold_em = " << threshold_em << endl;
		int org_min_modified_po = min_modified_po;
		int org_max_modified_po = max_modified_po;
		exdc_factor_new_v2(net, dd, cnode, unsort_cutnodes, sim_record, org_pla, dont_care, local_dc_set, po_set, po_dont_care_map, internal_index, input_index, sim_output_wi, pis, rand, simu_res, num_output, node_pattern_rate, max_sp, threshold_er, threshold_em, real_em, min_modified_po, max_modified_po, iIndex);
		double each_score = max_sp.score;
		cout << endl << "$cnode = " << cnode << ", each_score = " << each_score << endl;
        ftime(&et_es);
        double rt_exdc = ((et_es.time - st_es.time)*1000 + (et_es.millitm - st_es.millitm))/1000.0;
        total_exdc += rt_exdc;
		cout << "@runtime for a2. exdc_new_v2: " << rt_exdc << endl;
	
		//Insert the current simplification result into sim_record
		cout << "max_sp: " << endl;
		cout << "score : " << max_sp.score << endl;
		cout << "real_er: " << max_sp.real_er << ", ave_em: " << max_sp.ave_em << endl;
		cout << "pla: " << endl;
		for(int j = 0; j < max_sp.pla.size(); j++)
			cout << max_sp.pla[j] << endl;
		string cnode_str(cnode);
		itrm_ss = sim_record.find(cnode_str);
		if(itrm_ss == sim_record.end() && max_sp.score > 0)
  		{
  			cout << "inserted!" << endl;
  			sim_record.insert(pair<string, score_pla>(cnode_str, max_sp));
  			score_record.insert(pair<double, string>(max_sp.score, cnode_str));
  		}
  		//Update max_score and max_score_node
		if(each_score > max_score)
		{
			max_score = each_score;
            max_score_node = cnode;   
            max_score_status = max_sp.status;    
            final_pla = max_sp.pla;
            ave_em = max_sp.ave_em;
            cout << "min_modified_po = " << min_modified_po << ", max_modified_po = " << max_modified_po << endl;
            cur_min_modified_po = min_modified_po;
            cur_max_modified_po = max_modified_po;
		}		
		min_modified_po = org_min_modified_po;
        max_modified_po = org_max_modified_po;
		for(int i = 0; i < unsort_cutnodes.size(); i++)
			delete []unsort_cutnodes[i]; 
		
    }//for loop
    cout << "0. max_score_node = " <<  max_score_node << endl;  
    min_modified_po = cur_min_modified_po;
	max_modified_po = cur_max_modified_po;
    
/*    cout << "all simplifications: " << sim_record.size() << endl;
    index = 0;
    for(itrm_ss = sim_record.begin(); itrm_ss != sim_record.end(); itrm_ss++)
    {	
    	string node = itrm_ss->first;
    	struct score_pla sp = itrm_ss->second;
    	cout << index << ". cnode: " << node << endl;
    	cout << "score = " << sp.score << ", lit_save = " << sp.lit_save << ", real_er = " << sp.real_er << endl;
    	cout << "sim_pla: " << endl;
    	for(int j = 0; j < sp.pla.size(); j++)
    		cout << sp.pla[j] << endl;
    	index++;
    }
*/
    
    cout << "score_record: " << endl;
    for(itrmm_ds = score_record.begin(); itrmm_ds != score_record.end(); itrmm_ds++)
    {
    	double score = itrmm_ds->first;
    	string node = itrmm_ds->second;
    	cout << "node: " << node << ", score = " << score << endl;
    }
    
/*    if(!sim_record.empty())
    {
    	index = 0;
    	max_score = -1;
    	string max_node;
    	cout << "top 3 nodes: " << endl;
    	string filename_org = "./blif_files/ckt_org.blif";
    	fp = fopen(filename_org.c_str(), "r");
    	BnetNetwork *net_org = Bnet_ReadNetwork(fp); 
    	fclose(fp); 
    	itrmm_ds = score_record.end();
    	itrmm_ds--;
    	for(; index < score_record.size(); itrmm_ds--)
    	{
    		if(index > 2)	
    			break;    	
    		string snode = itrmm_ds->second;    		
    		itrm_ss = sim_record.find(snode);
    		struct score_pla sp = itrm_ss->second;
    		int lit_save = sp.lit_save;
    		vector<string> sim_pla = sp.pla;
    		double er_assure;
    		simu_assure(net_org, net, snode, sim_pla, er_assure);
    		cout << endl << "assure node: " << snode << endl;
    		cout << "score = " << sp.score << ", lit_save = " << sp.lit_save << ", real_er = " << sp.real_er << endl;
    		cout << "er_assure = " << er_assure << endl;
    		double er_incre = sp.real_er;
//    		double score_assure = (lit_save/(er_assure*0.8+er_incre*0.2))*(ini_threshold-er_assure);
			double score_assure = lit_save/(er_assure*0.8+er_incre*0.2);
    		cout << "max_score = " << max_score << ", score_assure = " << score_assure << endl;
    		if(score_assure > max_score)
    		{
    			max_score = score_assure;
    			max_node = snode;
    			final_pla = sim_pla;
    		}	    	
    		index++;
    	}
    	char *mn = new char[max_node.size()+1];
    	strcpy(mn, max_node.c_str());
   		st_lookup(net->hash, mn, &tmp);
    	max_score_node = tmp->name;
    	delete []mn;
    	cout << "1. max_score_node = " <<  max_score_node << endl;   
	}
*/	
    sim_record.clear();

#ifdef use_record	    
    //Remove the tfo nodes of current node from sim_record
    tmp = net->nodes;
	while(tmp != NULL)
	{
		tmp->visited = 0;
		tmp = tmp->next;
	}
	set<char*> tfo;
	set<char*>::iterator itrs;
    find_tranfanout(net, max_score_node, tfo);
    for(itrs = tfo.begin(); itrs != tfo.end(); itrs++)
    {
    	char *node = *itrs;
    	string node_str(node);
		itrm_ss = sim_record.find(node_str);
		if(itrm_ss != sim_record.end())
		{
			cout << "tfo node " << node << " is removed!" << endl;
			sim_record.erase(itrm_ss);
		}
    }
#endif
    
    
    cout << endl << "summary: " << endl;
    cout << "max_score_node = " << max_score_node << endl; 
    cout << "max_score = " << max_score << endl;
    cout << "final_pla: " << endl;
    for(int i = 0; i < final_pla.size(); i++)
    	cout << final_pla[i] << endl;
    cout << "total_exdc = " << total_exdc << endl; 
    
    return max_score_node;

}


void extract_sim_nodes(vector<string> &sim_node_lines, vector<string> &sim_node_nodes)
{    

	string str, s, s1;
	str = sim_node_lines[0];
	istringstream ss(str);
	ss >> s;
	if(s == ".names")
	{
		sim_node_nodes = sim_node_lines;
		return;
	}
		
	vector<string> sim_node_nodes_new;	
	string input_names(".names ");
	for(int i = 0; i < sim_node_lines.size(); i++)
	{
		string str = sim_node_lines[i];
		istringstream ss(str);
		ss >> s;		
		if(s[0] == '.')
		{
			if(s == ".ilb")
			{					
				while(ss >> s)
				{
					input_names.append(s);
					input_names.append(" ");
				}
			}
			else if(s == ".ob")
			{
				ss >> s;
				input_names.append(s);
				sim_node_nodes.push_back(input_names);
			}
			continue;
		}					
		ss >> s1;
		if(s1 == "2")
			break;
		sim_node_nodes.push_back(str);				
	}
	
	if(sim_node_nodes.size() == 1)
    {
    	stringstream ss(sim_node_nodes[0]);
    	string s;
    	while(ss >> s);
    	string str(".names ");
    	str.append(s);    	
    	sim_node_nodes_new.push_back(str);     //.names node
    	str.clear();
    	str.append(1, '0');
    	sim_node_nodes_new.push_back(str);     //0
    	sim_node_nodes = sim_node_nodes_new;
    }   
    
}
