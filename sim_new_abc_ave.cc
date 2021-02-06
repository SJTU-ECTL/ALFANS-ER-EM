#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include "head/queue.h"
#include "head/basics.h"
#include "head/exdc_factor_new_v2_ave.h"
#include "head/exdc_factor_new_v3_ave.h"
#include "head/helper.h"
#include "head/write_func.h"
#include "head/loc_sim_main_v2.h"
#include "head/simu_ckt.h"
#include "head/read_file.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"

using namespace std;

//Global variables and external variables 
static int consecu_num_large = 0, consecu_num_small = 0;
extern int em_class;
extern double T_em0;
extern int numPI_ini;
extern double ini_threshold_er, ini_threshold_em;
extern double er_weight, em_weight;
extern double ave_error_max, ave_error_max_cur, max_weight_max, max_weight_max_cur, real_er_max, real_er_max_cur;
double real_er_min, ave_error_min;
#define e 2.718

int score_mode = 6;
/*
functions in this file:

*/



//1. find_exdc_sim_ave(
char *find_exdc_sim_ave(BnetNetwork *net, DdManager **dd, BnetNetwork *net_comb, DdManager **dd_comb,  vector<string> &final_pla, map<string, struct score_pla> &sim_record, multimap<double, struct score_pla> &sim_record_top, map<string, int> &internal_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int &num_output, map<string, map<string, double> > &node_pattern_rate, double &max_score, double threshold_er, double threshold_em, double &real_er, double &real_em, int &min_modified_po, int &max_modified_po, double &T_em, double T_em1, double &T_er, double T_er1, int iIndex)
{
	//variables
    struct timeb st_es, et_es;
    double total_exdc = 0;
    BnetNode *nd, *tmp, *last_nd, *auxnd, *nd1, *nd2;
    ifstream fin;
    FILE *fp;
    string str;
    double weight_limit = log(threshold_em)/log(2);
    cout << "weight_limit = " << weight_limit << endl;    
    double real_em_pre = real_em;
    double real_er_pre = real_er;

    //iterators
    multimap<double, char*>::iterator itrm_dc; 
    multimap<double, string>::iterator itrmm_ds;
    map<string, struct score_pla>::iterator itrm_ss;
    map<string, multimap<double, string> >:: iterator itrm_cm;
    multimap<double, struct score_pla>::iterator itrmm_dss;
    set<char*>::iterator itrs;
    map<string, int>::iterator itrm_si;
   
    
    char com[100];
    multimap<double, char*> ave_sp;
    find_ave_sp(net, ave_sp);
    itrm_dc = ave_sp.begin();
    char *max_score_node = itrm_dc->second;  
    int max_score_status, max_score_min_weight = -1, max_score_max_weight = -1, max_score_lit_save = -1;
    double max_score_real_er, max_score_ave_em;
    max_score = -1;          
    int index = 0;  
    int flag_add_const = 0;
    ave_error_max = 0;
    real_er_max = 0;
	        	
    //start the for loop
    multimap<double, struct score_pla> score_record;
    int cur_min_modified_po, cur_max_modified_po;
    real_er_min = 1;
    ave_error_min = ini_threshold_em;
    for(itrm_dc = ave_sp.begin(); itrm_dc != ave_sp.end(); itrm_dc++, index++)
    {
        char *cnode = itrm_dc->second;
        st_lookup(net->hash, cnode, &nd);

        cout << endl <<"--------------------------------------------" << endl;  
        cout << "%%check for " << cnode << ", index = " << index << ", rp = " << nd->rp << endl; 
	cout << "currently, real_er = " << real_er << endl;
        
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
       
        //call abc to get the level of cnode
/*      if (iIndex == 0)
        	sprintf (com_abc, "abc -c 'read_blif ./blif_files/ckt_org.blif; print_level %s' > abc_level.out", cnode);
        else
        	sprintf (com_abc, "abc -c 'read_blif ./blif_files/ckt_sim.blif; print_level %s' > abc_level.out", cnode);
        system(com_abc);
	string abc_level_file = "abc_level.out";
        int level = read_abc_level(abc_level_file);
	cout << "level of " << cnode << ": " << level << endl;
*/
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
	set<char*> tfo;
	int affected_po_size = 0;
        find_tranfanout(net, nd->name, tfo, po_set, affected_po_size);
    	cout << endl << "current node " << nd->name << "'s affected pos: " << affected_po_size << endl;
	cout << "po_set: " << endl;
	for (itrs = po_set.begin(); itrs != po_set.end(); itrs++)
		cout << *itrs << " ";
	cout << endl;
//      double score_factor = tfo.size()/(double)affected_po_size;
//	cout << "tfo.size: " << tfo.size() << ", affected_po_size = " << affected_po_size << ", score_factor = " << score_factor << endl;
//	double score_factor = tfo.size();
//	cout << "tfo.size(): " << tfo.size() << endl;

//	find_tranfanout_po(net, nd->name, po_set);
	double score_factor = 1;

	//Find all qualified exdcs
	cout << "**********************************" << endl;
	cout << endl << "a2. exdc_factor_new_v2: " << endl;
	ftime(&st_es);	
	struct score_pla max_sp;
	cout << "threshold_em = " << threshold_em << endl;
       	int org_min_modified_po = min_modified_po;
	int org_max_modified_po = max_modified_po;
	exdc_factor_new_v3_ave(net, dd, cnode, score_factor, unsort_cutnodes, score_record, org_pla, dont_care, po_set, internal_index, input_index, sim_output_wi, pis, rand, simu_res, num_output, node_pattern_rate, max_sp, max_score, threshold_er, threshold_em, real_em, min_modified_po, max_modified_po, T_em, T_er, flag_add_const, iIndex);
        ftime(&et_es);
        double rt_exdc = ((et_es.time - st_es.time)*1000 + (et_es.millitm - st_es.millitm))/1000.0;
        total_exdc += rt_exdc;
	cout << "@runtime for a2. exdc_new_v2: " << rt_exdc << endl;

	if (score_mode == 3)
	{
		double each_score = max_sp.score;
		if (each_score > max_score)
			max_score = each_score;
	}
	if (score_mode == 6)
	{
		double each_score = max_sp.score;
		cout << "each_score = " << each_score << ", max_score = " << max_score << endl;
		if (each_score > max_score)
		{
			max_score = each_score;
			max_score_node = cnode;   
			max_score_lit_save = max_sp.lit_save;
			max_score_real_er = max_sp.real_er;
			max_score_ave_em = max_sp.ave_em;
			max_score_min_weight = max_sp.min_weight;
			max_score_max_weight = max_sp.max_weight;
			final_pla = max_sp.pla;
			cur_min_modified_po = min_modified_po;
			cur_max_modified_po = max_modified_po;
		//	real_er = max_score_real_er + real_er_pre;
			real_er = max_score_real_er;
			real_em = max_score_ave_em;
		}
	}

	min_modified_po = org_min_modified_po;
       	max_modified_po = org_max_modified_po;

	for(int i = 0; i < unsort_cutnodes.size(); i++)
		delete []unsort_cutnodes[i]; 

    }//for loop

    //************************************************************************************
    //update the scores and find the ASE with highest score
/*    cout << "after this iteration, real_er_max = " << real_er_max << ", ave_error_max = " << ave_error_max << endl;
    cout << "score_record of this iteration:" << score_record.size() << endl;
    double ls, er, em, er_part, em_part;
    max_score = -1;
    double em_min = ave_error_max;
    for (itrmm_dss = score_record.begin(); itrmm_dss != score_record.end(); itrmm_dss++)
    {
	ls = itrmm_dss->second.lit_save;
	er = itrmm_dss->second.real_er; 
	em = itrmm_dss->second.ave_em;
	if (em < em_min) em_min = em;
    }
    double add_em_part = -log(em_min/ave_error_max)/log(2);
    cout << "em_min = " << em_min << ", add_em_part = " << add_em_part << endl;
*/
    //************************************************************************************
    //update the scores and find the ASE with highest score
    cout << "after this iteration, real_er_min = " << real_er_min << ", ave_error_min = " << ave_error_min << endl;
    cout << "score_record of this iteration:" << score_record.size() << endl;
    double ls, er, em, er_part, em_part, min_w, max_w;
    double decoup_em, decoup_em_min = ini_threshold_em * 2;
    double max_ave = 0;
    double em_p;
    for (itrmm_dss = score_record.begin(); itrmm_dss != score_record.end(); itrmm_dss++)
    {
	ls = itrmm_dss->second.lit_save;
	er = itrmm_dss->second.er_part; 
	em = itrmm_dss->second.ave_em;
	min_w = itrmm_dss->second.min_weight;
	max_w = itrmm_dss->second.max_weight;
	em_p = pow(e, (log(em)/log(2))/T_em);	
	cout << setw(12) << itrmm_dss->second.score << ": " << setw(10) << itrmm_dss->second.node << ", l: " << setw(3) << ls << ", er: " << setw(12) << er;
	cout << ", ave_error: " << setw(10) << em << ", em_p: " << em_p << ", min_weight: " << setw(3) << min_w << ", max_weight: " << setw(3) << max_w << endl;	
	if (em > max_ave) max_ave = em;
    }
    cout << "max_ave = " << max_ave << endl;

    //summary of this iteration
    cout << endl << "summary:" << endl;
    cout << "max_score_node = " << max_score_node << endl; 
    cout << "final_pla: " << endl;
    for(int i = 0; i < final_pla.size(); i++)
    	cout << final_pla[i] << endl;

    min_modified_po = cur_min_modified_po;
    max_modified_po = cur_max_modified_po;

    	double ratio_1 = 0.5, ratio_2 = 0.1;
	// current T_em adjusting mechanism
//	 if (max_score_weight <= weight_limit)
/*	if (iIndex > 0)
	{
		if (max_score_ave_em > real_em_pre)
		{
			double tmp = log(max_score_ave_em/real_em_pre)/log(2);
			T_em = T_em * pow(ratio_1, tmp);
			if (T_em < T_em1) T_em = T_em1;
		}
	 	cout << "after adjusting, T_em = " << T_em << endl;
	}
*/
//	double tmp_pow = log(ini_threshold_em/max_score_ave_em)/log(2);
/*	if (max_score_ave_em > real_em_pre)
	{
		double tmp_pow = 1;
	//	if (iIndex > 0)
	//	{
	//		tmp_pow = (max_score_ave_em/real_em_pre);
	//		tmp_pow = log(tmp_pow)/log(2);
	//	}
	//	tmp_pow *= pow(2, max_score_ave_em/ini_threshold_em);
	//	cout << "  tmp_pow = " << tmp_pow << endl;
	//	T_em = T_em * pow(ratio_1, tmp_pow);
		T_em = T_em * pow((1 - max_score_ave_em/ini_threshold_em), 1);
		if (T_em < T_em1) T_em = T_em1;
	 	cout << "after adjusting, T_em = " << T_em << endl;
	}
*/
/*	if (iIndex > 0)
	{
		if (max_score_real_er > real_er_pre)
		{
			double tmp = max_score_real_er/real_er_pre;
			T_er = T_er * pow(ratio_1, tmp);
			if (T_er < T_er1) T_er = T_er1;
		}
	 	cout << "after adjusting, T_er = " << T_er << endl;
	}
*/
        //2017.6.20
/*	 if (max_score_max_weight <= weight_limit)
	 {
		consecu_num_small++;
		double pow_val = log(1 + consecu_num_small)/log(2);
		T_em = T_em * pow(0.5, pow_val);
		if (T_em < T_em1) T_em = T_em1;
	 	cout << "after adjusting, T_em = " << T_em << endl;
	 }
	else
		consecu_num_small = 0;
*/

	//2017.6.21
/*	if (iIndex > 0)
	{
		if (max_score_ave_em > real_em_pre)
		{
		//	double tmp = log(max_score_ave_em/real_em_pre)/log(2);
			T_em = T_em * pow(ratio_1, 1);
			if (T_em < T_em1) T_em = T_em1;
		}
	 	cout << "after adjusting, T_em = " << T_em << endl;
	}
//end of 2017.6.21
*/
	//2017.6.22
	if (iIndex > 0)
	{
		ratio_1 = ini_threshold_er * 0.8;
		if (max_score_ave_em > real_em_pre)
		{
			consecu_num_small++;
			T_em = T_em * pow(ratio_1, log(1+consecu_num_small)/log(2));
			if (T_em < T_em1) T_em = T_em1;
		}
		else
			consecu_num_small = 0;
	 	cout << "after adjusting, T_em = " << T_em << endl;
	}
//end of 2017.6.22

    cout << "total_exdc = " << total_exdc << endl; 
    return max_score_node;

}


