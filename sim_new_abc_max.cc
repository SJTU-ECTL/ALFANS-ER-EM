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
#include "head/exdc_factor_new_v2_max.h"
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
extern int numPI_ini, numPO_ini;
extern double ini_threshold_er, ini_threshold_em;
extern vector<string> sub_abs_ckt, sub_abs_pi, sub_abs_po;
extern vector<string> comparator_ckt, comparator_pi, comparator_po;
extern vector<int> comp_number;
extern double er_weight, em_weight;
extern double ave_error_max, ave_error_max_cur, max_weight_max, max_weight_max_cur, real_er_max, real_er_max_cur;

int min_weight_pre = 0;
int num_candidate_pre = 0;
int num_within_em_over = 0;
/*
functions in this file:
char *find_exdc_sim_max();

*/


//1. find_exdc_sim()
char *find_exdc_sim_max(BnetNetwork *net, DdManager **dd, BnetNetwork *net_comb, DdManager **dd_comb,  vector<string> &final_pla, map<string, struct score_pla> &sim_record, multimap<double, struct score_pla> &sim_record_top, map<string, int> &internal_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int &num_output, map<string, map<string, double> > &node_pattern_rate, double &max_score, double threshold_er, double threshold_em, double &real_er, double &real_em, int &min_modified_po, int &max_modified_po, double &T_em, double T_em1, int iIndex)
{
	//variables
    struct timeb st_es, et_es;
    double total_exdc = 0;
    BnetNode *nd, *tmp, *last_nd, *auxnd, *nd1, *nd2;
    ifstream fin;
    FILE *fp;
    string str;
    int weight_limit = log(threshold_em)/log(2);
    cout << "weight_limit = " << weight_limit << endl;    
    double real_er_pre = real_er;

    //iterators
    multimap<double, char*>::iterator itrm_dc; 
    multimap<double, string>::iterator itrmm_ds;
    map<string, struct score_pla>::iterator itrm_ss;
    map<string, multimap<double, string> >:: iterator itrm_cm;
    multimap<double, struct score_pla>::iterator itrmm_dss;
    set<char*>::iterator itrs;
    map<string, int>::iterator itrm_si;
    map<string, struct wi_pair>::iterator itrm_sw;
   
    
    char com[100];
    multimap<double, char*> ave_sp;
    find_ave_sp(net, ave_sp);
    itrm_dc = ave_sp.begin();
    char *max_score_node = itrm_dc->second;  
    int max_score_status, max_score_min_weight = -1, max_score_max_weight = -1, max_score_lit_save = -1;
    double  max_score_real_er = 1;
    max_score = -1;          
    double max_score_second = -1;
    int index = 0;  
    ave_error_max_cur = 0;
    max_weight_max_cur = 0;
    real_er_max_cur = 0;
	        	
    //start the for loop
    multimap<double, struct score_pla> score_record;
    multimap<double, struct score_pla>::reverse_iterator re_itrmm_dss;
    int cur_min_modified_po, cur_max_modified_po;
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
    	set<char*> po_set, tfo;
	map<char*, int> po_level;
	map<char*, int>::iterator itrm_ci;
	int affected_po_size = 0;
//      find_tranfanout(net, nd->name, tfo, po_set, affected_po_size);
	ftime(&st_es);
      	find_tranfanout_level(net, nd->name, 0, tfo, po_level, affected_po_size);
	ftime(&et_es);
        double rt_fanout = ((et_es.time - st_es.time)*1000 + (et_es.millitm - st_es.millitm))/1000.0;
	cout << "runtime for find_tranfanout: " << rt_fanout << endl;
    	cout << endl << "current node " << nd->name << "'s affected pos: " << affected_po_size << endl;
//	if (affected_po_size > numPO_ini * 0.8) continue;
/*	cout << "po_set: " << endl;
	for (itrs = po_set.begin(); itrs != po_set.end(); itrs++)
		cout << *itrs << " ";
	cout << endl;
*/
	double weight_factor = 0, score_factor = 0;
	for (itrm_ci = po_level.begin(); itrm_ci != po_level.end(); itrm_ci++)
	{
		cout << itrm_ci->first << ": " << itrm_ci->second << endl;
		po_set.insert(itrm_ci->first);
	/*	string po(itrm_ci->first);
		if(iIndex > 0) po.append("sim");
		itrm_sw = sim_output_wi.find(po);
		int weight = itrm_sw->second.weight;
		int dist = itrm_ci->second;
		weight_factor += (double)weight/dist;
	*/
	}
//	cout << "weight_factor = " << weight_factor << endl; 
//	score_factor = log(score_factor)/log(2);
//	if (score_factor == 0) score_factor = 0.5;
//	score_factor = 1/score_factor;
//	score_factor = pow(2, -weight_factor/T_em);
//	cout << "score_factor = " << score_factor << endl; 
	
	
//      double score_factor = tfo.size()/(double)affected_po_size;
//	cout << "tfo.size: " << tfo.size() << ", affected_po_size = " << affected_po_size << ", score_factor = " << score_factor << endl;
//	double score_factor = tfo.size();
//	cout << "tfo.size(): " << tfo.size() << endl;

//	find_tranfanout_po(net, nd->name, po_set);

	//Find all qualified exdcs
	cout << "**********************************" << endl;
	cout << endl << "a2. exdc_factor_new_v2: " << endl;
	ftime(&st_es);	
	struct score_pla max_sp;
	exdc_factor_new_v2_max(net, dd, cnode, score_factor, unsort_cutnodes, score_record,  org_pla, dont_care, po_set, internal_index, input_index, sim_output_wi, pis, rand, simu_res, num_output, node_pattern_rate, max_sp, max_score, threshold_er, threshold_em, real_em, T_em, iIndex);
	cout << endl << "$cnode = " << cnode << ", score = " << max_sp.score << endl;
        ftime(&et_es);
        double rt_exdc = ((et_es.time - st_es.time)*1000 + (et_es.millitm - st_es.millitm))/1000.0;
        total_exdc += rt_exdc;
	cout << "@runtime for a2. exdc_new_v2: " << rt_exdc << endl;

	cout << "max_sp: " << endl;
	cout << "score : " << max_sp.score << endl;
	cout << "real_er: " << max_sp.real_er << endl;
	cout << "max_weight: " << max_sp.max_weight << endl;
	if (max_sp.score > max_score)
		max_score = max_sp.score;

 	//Update max_score and max_score_node
/*	if (em_class == 1)
	{
		cout << "start pick_max_sp!" << endl;
		pick_max_sp(sim_record_top, T_em, max_sp, max_score_second);
		cout << "new max_sp: " << max_sp.node << endl;
		cout << "score : " << max_sp.score << " ";
		cout << "real_er: " << max_sp.real_er << " ";
		cout << "max_score_second: " << max_score_second << endl;
		max_score_node = max_sp.node;
		max_score = max_sp.score;
		final_pla = max_sp.pla;
		cout << "max_weight: " << max_sp.max_weight << endl;
		max_score_weight = max_sp.max_weight;
	}
	else if (em_class == 2)
*/
/*	double each_score = max_sp.score;
	if(each_score > max_score)
	{
		max_score = each_score;
		max_score_node = cnode;   
		max_score_weight = max_sp.max_weight;
		max_score_lit_save = max_sp.lit_save;
		max_score_real_er = max_sp.real_er;
		final_pla = max_sp.pla;
		real_er = max_sp.real_er;
	}		
*/

	for(int i = 0; i < unsort_cutnodes.size(); i++)
		delete []unsort_cutnodes[i]; 

    }//for loop

    //************************************************************************************
    //summary of this iteration
/*    cout << endl << "summary: T_em = " << T_em << endl;
    cout << "max_score_node = " << max_score_node << endl; 
    cout << "max_score = " << max_score << ", max_score_lit_save = " << max_score_lit_save << ", max_score_real_er = " << max_score_real_er << endl;
    cout << "max_score_weight = " << max_score_weight << endl;
    real_er_max = max_score_real_er;
    max_weight_max = max_score_weight;
    cout << "final_pla: " << endl;
    for(int i = 0; i < final_pla.size(); i++)
    	cout << final_pla[i] << endl;
*/
    if (score_record.empty())
    {
	cout << "score_record is empty!" << endl;
	max_score = -1;
	return max_score_node;
    }
    cout << "num_within_em_over = " << num_within_em_over << endl;
    int num_total = score_record.size();
    int num_check = 0;
    for (re_itrmm_dss = score_record.rbegin(); ; re_itrmm_dss++)
    {
	char *cnode = re_itrmm_dss->second.node;
	cout << "# " << re_itrmm_dss->second.score << ": " << cnode << ", sf: " << re_itrmm_dss->second.sf <<  ", l: " << re_itrmm_dss->second.lit_save << ", er_part: " << re_itrmm_dss->second.er_part << ", min_w = " << re_itrmm_dss->second.min_weight << endl;
	num_check++;
	if (num_check >= num_total) break;
    }

    double real_er_acc;
    num_check = 0;
    struct score_pla max_sp;
    int flag_find = 0;
    cout << endl << "score_record of this iteration:" << score_record.size() << endl;
    for (re_itrmm_dss = score_record.rbegin(); ; re_itrmm_dss++)
    {
	char *cnode = re_itrmm_dss->second.node;
	double er_single = re_itrmm_dss->second.real_er;
	cout << endl << "# " << re_itrmm_dss->second.score << ": " << cnode << ", l: " << re_itrmm_dss->second.lit_save << ", er_part: " << re_itrmm_dss->second.er_part << "er_single: " << er_single;
	cout <<  ", min_w: " << re_itrmm_dss->second.min_weight << ", max_w: " << re_itrmm_dss->second.max_weight << endl;	
	vector<string> cur_final_pla = re_itrmm_dss->second.pla;
	get_real_er_acc(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number, cur_final_pla, real_er_acc);
//	real_er_acc = real_er + er_single - real_er * er_single;
	if (real_er_acc <= ini_threshold_er)
	{
		max_sp = re_itrmm_dss->second;
		max_sp.real_er = real_er_acc;
		flag_find = 1;
		break;
	}
	num_check++;
	if (num_check >= num_total) break;
    }
    if (flag_find == 0)
    {
	cout << "all beyond ini_threshold_er!" << endl;
	max_score = -1;
	return max_score_node;
    }

    max_score = max_sp.score;
    max_score_node = max_sp.node;   
    max_score_lit_save = max_sp.lit_save;
    max_score_real_er = max_sp.real_er;
    max_score_min_weight = max_sp.min_weight;
    max_score_max_weight = max_sp.max_weight;
    final_pla = max_sp.pla;
    real_er = max_sp.real_er;
    cout << endl << "summary: T_em = " << T_em << endl;
    cout << "max_score_node = " << max_score_node << endl; 
    cout << "max_score = " << max_score << ", max_score_lit_save = " << max_score_lit_save << ", max_score_real_er = " << max_score_real_er << endl;
//    cout << "max_score_min_weight = " << max_score_min_weight << endl;
    cout << "final_pla: " << endl;
    for(int i = 0; i < final_pla.size(); i++)
    	cout << final_pla[i] << endl;

    //adjusting T_em	
    	double ratio_1 = 1.0, ratio_2 = 0.5, ratio_3 = 0.9;
	// current T_em adjusting mechanism
	cout << "weight_limit = "<< weight_limit << endl;
	cout << "before adjusting, T_em = " << T_em << endl;
	double pow_val; 

	//2017.6.19
/*	 if (max_score_max_weight <= weight_limit)
	 {
		consecu_num_small++;
		pow_val = log(1 + consecu_num_small)/log(2);
		T_em = T_em * pow(ratio_2, pow_val);
		if (T_em < T_em1) T_em = T_em1;
	 }
	else
		consecu_num_small = 0;
//end of 2017.6.19
*/
	//2017.6.22
/*	 if (max_score_max_weight <= weight_limit)
	 {
		consecu_num_small++;
		pow_val = log(1 + consecu_num_small)/log(2);
		ratio_2 = ini_threshold_er/pow_val;
		cout << "ratio_2 = " << ratio_2 << endl;
		T_em = T_em * pow(ratio_2, 1);
		if (T_em < T_em1) T_em = T_em1;
	 }
*/
	ratio_2 = ini_threshold_er * 0.8;
	 if (max_score_max_weight <= weight_limit)
	 {
		consecu_num_small++;
		pow_val = log(1 + consecu_num_small)/log(2);
		T_em = T_em * pow(ratio_2, pow_val);
		if (T_em < T_em1) T_em = T_em1;
	 }
	else
		consecu_num_small = 0;


//end of 2017.6.22
/*	else
	{
		consecu_num_small = 0;
	//	pow_val = max_score_min_weight/weight_limit;
	//	T_em = T_em * pow(ratio_2, pow_val);
		if (max_score_min_weight <= min_weight_pre)
		{
		//	consecu_num_large++;
		//	T_em = T_em * pow(ratio_3, log(1+consecu_num_large)/log(2));
		//	min_weight_pre = max_score_min_weight;
			T_em = T_em * ratio_2;
		}
	}
*/

	//2017.6.21
/*	if (max_score_min_weight > min_weight_pre)
	{
		consecu_num_large++;
		pow_val = log(1 + consecu_num_large)/log(2);
		T_em = T_em * pow(0.8, pow_val);
		if (T_em < T_em1) T_em = T_em1;
	}
	else
		consecu_num_large = 0;
	min_weight_pre = max_score_min_weight;
*/
//end of 2017.6.21

//	double er_margin = (ini_threshold_er - max_score_real_er)/(1 - max_score_real_er);
//	double ratio_er = 2 - er_margin/ini_threshold_er;
//	cout << "ratio_er = " << ratio_er << endl;
/*	double er_margin_ratio = (ini_threshold_er - max_score_real_er)/(1-max_score_real_er);
	er_margin_ratio = er_margin_ratio/ini_threshold_er;
	cout << "er_margin_ratio = " << er_margin_ratio << endl;
	int num_candidate_now = score_record.size();
*/
/*	if (max_score_min_weight > min_weight_pre)
	{
	//	T_em = T_em * pow(ratio_2, max_score_min_weight - min_weight_pre);
	//	T_em = (weight_limit - max_score_min_weight)/log(1.5 * ratio_er);	
	//	T_em = weight_limit/log(1.5 * max_score_min_weight);	
	//	T_em = T_em * pow(ratio_2, 1);
		double ratio_em = (double)max_score_min_weight/weight_limit;
		T_em = T_em * (1 - ratio_em) * er_margin_ratio;
		if (T_em < T_em1) T_em = T_em1;
		min_weight_pre = max_score_min_weight;
	}
*/
/*	cout << "num_candidate_pre = " << num_candidate_pre << ", num_candidate_now = " << num_candidate_now << endl;
	if (num_candidate_now < num_candidate_pre)
	{
		consecu_num_small++;
		T_em = T_em * pow(ratio_2, log(1+consecu_num_small)/log(2));
		if (T_em < T_em1) T_em = T_em1;
	}
//	else
//		consecu_num_small = 0;
	num_candidate_pre = num_candidate_now;
*/
//	double ratio_em = 1.5 / (1 - min_weight_pre/weight_limit); 
//	T_em = weight_limit / log(ratio_em);


	 cout << "after adjusting, T_em = " << T_em << endl;

    //adjusting er_weight
  //  er_weight = 1 /pow(1.2, max_score_real_er/ini_threshold_er);
    // er_weight = pow(2, 1+max_score_real_er/ini_threshold_er);
//	 cout << "after adjusting, er_weight = " << er_weight << endl;

    cout << "total_exdc = " << total_exdc << endl; 

    return max_score_node;

}
