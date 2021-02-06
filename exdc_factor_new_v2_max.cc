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
#include "head/read_file.h"
#include "head/write_func.h"
#include "head/exdc_helper.h"
#include "head/btree.h"
#include "head/simu_ckt.h"
#include "head/loc_sim_main_v2.h"
#include "cudd/bnet.h"
#include "cudd/cudd_build_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"

using namespace std;

extern int numPI_ini;
extern double total_deviation;
extern int sample_num;
extern double ini_threshold_em, ini_threshold_er;
extern int num_one_po_equal, num_one_po_unequal;
extern map<string, set<char*> > po_tfi;
extern map<string, set<string> > po_inputs;
extern map<string, vector<string> > po_cone_string;
extern vector<string> sub_abs_ckt, sub_abs_pi, sub_abs_po;
extern vector<string> comparator_ckt, comparator_pi, comparator_po;
extern vector<int> comp_number;
extern double max_weight_max, max_weight_max_cur, real_er_max, real_er_max_cur;
extern int lit_save_max, lit_save_max_cur;
extern double er_weight, em_weight;
extern int num_within_em_over;

static int error_rate_unconstraint;
#define mode 2    ////mode = 2 represents maximum error magnitude; mode = 1 represents average error magnitude
#define BDD_SAT_mode 2 // 1: BDD, 2: SAT
#define exp_ratio 0.005
#define score_mode 3 
#define e 2.718
//#define e 2

#define flag_er_part 1 
#define flag_score_factor 1 
#define flag_score_factor_er 0 
#define flag_beyond 0 
#define flag_round 0
#define flag_pick_max 1 

#define beyond_ratio_1 0.3
#define beyond_ratio_2 0.3 
#define ignore_po_size_one 1
#define mul_ratio 1 


/*
functions in this file:

*/

//Global variables and external variables 


//exdc_factor
void exdc_factor_new_v2_max(BnetNetwork *net, DdManager **dd, char *cnode, double score_factor, vector<char*> &unsort_cutnodes, multimap<double, struct score_pla> &score_record, vector<string> &org_pla, vector<string> &dont_care, set<char*> &po_set, map<string, int> &internal_index, vector<struct index_flag> &input_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int num_output, map<string, map<string, double> > &node_pattern_rate, struct score_pla &max_sp, double max_score,  double threshold_er, double threshold_em, double real_em, double &T_em, int iIndex)
{
    //iterators 
    struct timeb st, et, st1, et1, st2, et2;   
    map<string, double>::iterator itrm_sd, itrm_sd1;  
    map<string, struct score_pla>::iterator itrm_ss; 
    multimap<double, struct score_pla>::iterator itrmm_ds;
    set<int>::iterator itrs0, itrs1, itrs2;
    set<string>::iterator itrss, itrss0, itrss1;
    set<char*>::iterator itrs_char;
    map<string, set<char*> >::iterator itrm_sc;
    map<string, set<string> >::iterator itrm_sss;
    map<string, vector<string> >::iterator itrm_sv;
    map<string, struct wi_pair>::iterator itrm_sw;
    map<int, int>::iterator itrmi;
    
    	//variables
    	//initialize max_sp
    	max_sp.node = cnode;
    	max_sp.score = -(1e+9);
	max_sp.lit_save = 0;
	max_sp.real_er = 1;
	max_sp.max_em = 1e+9;
	max_sp.status = 0;
	int lit_save, status = 0;
  	double score, this_error_mag = 0;
  	double min_error_mag = 0;
        cout << "threshold_er = " << threshold_er << endl;
	double er_margin = threshold_er/(1 - (ini_threshold_er-threshold_er));
        cout << "er_margin = " << er_margin << endl;
	double cur_real_er = ini_threshold_er - threshold_er;
  	vector<string> final_pla;	
  	char com[100];
    	string str, s;
    	BnetNode *nd, *tmp, *auxnd;	
	double est_er_acc, est_er_margin;
	double lit_part, er_part, em_part;
	double rt_one_check_total = 0;
    
    	if (iIndex == 0)
	{
		if (threshold_er == 1)
			error_rate_unconstraint = 1;
		else
			error_rate_unconstraint = 0;
	}
	
	cout << endl << "$current node: " << cnode << endl;	
	if (error_rate_unconstraint)
		threshold_er = 1;	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* step1. obtain pattern_rate for cnode and see if there are patterns with probabilities within threshold_er */
	ftime(&st1);
	set<string> dont_care_set;
	for (int i = 0; i < dont_care.size(); i++) 
		dont_care_set.insert(dont_care[i]);
	string sname(cnode);
	map<string, map<string, double> >::iterator itrm_ssd;
	itrm_ssd = node_pattern_rate.find(sname);
	map<string, double> pattern_rate = itrm_ssd->second;
	if (error_rate_unconstraint == 0)
        {
		int num_useful = 0;
		if (dont_care.empty())
		{
			for(itrm_sd = pattern_rate.begin(); itrm_sd != pattern_rate.end(); itrm_sd++)
			{
				double er = itrm_sd->second;
				if ( er <= mul_ratio * er_margin )
				{
					num_useful++;
					break;
				}
			}
			if(!num_useful)
			{
				cout << "no useful patterns!" << endl;
				return; 
			}
		}
        }	
	ftime(&et1);
	double rt_step1 = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
	cout << "exdc: runtime for step1: " << rt_step1 << endl;
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* step2. check the range of weights of the affected POs to exclude some impossible cases */
	//if current node has only one affected po, then first get its weight and compute its maximum error magnitude
	double weight_limit = log(threshold_em)/log(2);
	int min_weight = 1000, max_weight = -1;
	if (po_set.size() == 1) 
	{
		cout << "Current node has one affected po!" << endl;
		itrs_char = po_set.begin();
		string po(*itrs_char);
		if(iIndex > 0) po.append("sim");
		itrm_sw = sim_output_wi.find(po);
		int weight = itrm_sw->second.weight;
		min_weight = max_weight = weight;
		if (iIndex == 0 && min_weight > weight_limit) return;

	/*	st_lookup(net->hash, cnode, &nd);
		if (nd->nfo == 1 && iIndex == 0)
		{
			char *fanout = nd->fanouts[0];
			st_lookup(net->hash, fanout, &tmp);
			if (tmp->type == BNET_OUTPUT_NODE)
			{
				cout << "weight = " << weight << ", weight_limit = " << weight_limit << endl;
				min_error_mag = pow(2.0, weight);
				if (min_error_mag > threshold_em)
				{
					cout << "cnode is PO itself and its error mag is larger than threshold_em!" << endl;
					return;
				}
			}
		}			
	*/
	}
	else
	{
		cout << "Current node has multiple affected pos!" << endl;
		int weight;
		for(itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
		{
			string po(*itrs_char);
			if(iIndex > 0) po.append("sim");
			itrm_sw = sim_output_wi.find(po);
			weight = itrm_sw->second.weight;
			if (weight < min_weight) min_weight = weight;
			if (weight > max_weight) max_weight = weight;
		}
	//	if (min_weight > weight_limit) return;
	}
/*
//	int pre_min_modified_po = min_modified_po;
//	int pre_max_modified_po = max_modified_po;
//	cout << "min_modified_po = " << min_modified_po << ", max_modified_po = " << max_modified_po << endl;
	cout << "min_weight = " << min_weight << ", max_weight = " << max_weight << endl;
//	int this_min_modified_po = (min_weight < min_modified_po)? min_weight: min_modified_po;
//	int this_max_modified_po = (max_weight > max_modified_po)? max_weight: max_modified_po;
//	cout << "this_min_modified_po = " << this_min_modified_po << ", this_max_modified_po = " << this_max_modified_po << endl;
//	min_modified_po = this_min_modified_po;
//	max_modified_po = this_max_modified_po;	
//	if (this_min_modified_po > weight_limit) return;
*/
	if (flag_score_factor)	
	{
	//	double ratio = -min_weight/(double)T_em;
		double ratio = -max_weight/(double)T_em;
		double control= pow(e, ratio);
		score_factor = control;
		cout << "T_em = " << T_em << ", min_weight = " << min_weight << ", score_factor = " << score_factor << endl;
	}
	else
		score_factor = 1;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* step4. get the factored form of cnode and build its binary tree */
  	//get num_input and name_pos
  	ftime(&st1);
  	int num_input = unsort_cutnodes.size();
	cout << "unsort_cutnodes: " << endl;
	set<string> insig;
    	map<string, int> name_pos;
    	map<string, int>::iterator itrm_si;
	for(int i = 0; i < unsort_cutnodes.size(); i++)
	{
	    cout << unsort_cutnodes[i] << " ";
	    string name(unsort_cutnodes[i]);
	    insig.insert(name);
	    name_pos.insert(pair<string, int>(name, i));
	}
	cout << endl; 	
	
	//print and read the factored form of cnode
	cout << endl << "$factor form of bignode: " << endl;		
	write_bignode_pla(net, cnode);   		                        
	sprintf(com, "sis -t none -f ./script/print_factor_org.rug > factor.txt");
	system(com);
	string ffe;
	read_factor_v2(ffe);
	
	//build a binary tree for the factored-form-expression
	btNode *ini_root;
	string ffe_space;
	add_space_star(ffe, ffe_space);
	build_tree_from_exp(ffe_space, &ini_root);
	cout << "ffe_space = " << ffe_space << endl;
	ftime(&et1);
	double rt_step4 = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
	cout << "exdc: runtime for step4: " << rt_step4 << endl;
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* step5. prepare for the simplification process */  	
  	//pick exdcs for ffe(factor-form-expression) 
	cout << endl << "**************************************" << endl;
	cout << "$pick exdcs for ffe(factor-form-expression): " << endl;  	
  	cout << "current factor: " << ffe << endl; 
  	btNode *root;  		
  	//find all leaf nodes and set their indexes: leaf_set
  	map<int, btNode*> ini_leaf_set, leaf_set;
  	map<int, btNode*>::iterator itrm_ib;
  	int start_ind = 0;
	visitleaf(ini_root, ini_leaf_set, start_ind);	
	leaf_set = ini_leaf_set;	
	int num_lit = ini_leaf_set.size();
	cout << "num_lit = " << num_lit << endl;	
	comp_exp(&ini_root);  //compute expression at each node  		
		
	//obtain ini_node_inv_pla & ini_current_inv_pla
	map<int, vector<string> > ini_index_inv_pla, index_inv_pla;
	map<int, vector<string> >::iterator itrm_iv;
	multimap<string, int> ini_node_index, node_index;
	multimap<string, int>::iterator itrmm_si;
	set<string> ini_current_inv_pla, current_inv_pla;
	update_index_inv_pla(ini_root, ini_leaf_set, name_pos, num_input, ini_index_inv_pla, ini_node_index, ini_current_inv_pla);
	index_inv_pla = ini_index_inv_pla;
	current_inv_pla = ini_current_inv_pla;
	node_index = ini_node_index;
	 
	 
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* step6. start check every possible simplification	 */  
	set<string> exdc_cubes, all_exdc_minterms;
	vector<string> sim_org_pla, sim_org_pla_tmp, final_org_pla;
	int num_lit_org = num_lit;
	if(num_lit > 5) num_lit = 5;           //set the limit for literal save
        int flag_continue = 0; 
	int num_loop = 0;
	//ASPDAC modification: ignore the case where lit_save = 1
	for(int i = 1; i <= num_lit; i++)
	{
		if(i == num_lit)
		{
			cout << endl << "##lit_save = " << i << endl;
			num_loop++;
			ftime(&st);
			cout << "case 1. removing this whole factor: " << endl;
			ftime(&st2);
			exdc_cubes.clear();
			for(itrss = ini_current_inv_pla.begin(); itrss != ini_current_inv_pla.end(); itrss++)
				exdc_cubes.insert(*itrss);
			all_exdc_minterms.clear();
			for(itrss = exdc_cubes.begin(); itrss != exdc_cubes.end(); itrss++)
			{
				string cube = *itrss;
				exp_cube_set(cube, all_exdc_minterms);
			}

			double this_real_er_case1, this_real_er_case2;	
			status = simu_real_er(all_exdc_minterms, dont_care_set, pattern_rate, this_real_er_case1);
			cout << "this_real_er_case1 = " << this_real_er_case1 << ", er_margin = " << er_margin << endl;	
			ftime(&et2);
			double rt_part1 = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
			cout << "runtime for part1: " << rt_part1 << endl;

			if (this_real_er_case1 > mul_ratio * er_margin) cout << "!error rate is beyond er_margin" << endl;	
			else
			{
				cout << "!error rate is within er_margin" << endl;	
			    	// compute lit_save
			    	ftime(&st2);
			    	sim_org_pla.clear();
				for(int k = 0; k < org_pla.size(); k++)
				{
					itrss = ini_current_inv_pla.find(org_pla[k]);
					if(itrss == ini_current_inv_pla.end())
						sim_org_pla.push_back(org_pla[k]);
				}	
				final_org_pla.clear();			
				lit_save = get_save_new(net, cnode, 0, sim_org_pla, final_org_pla, iIndex);

				//compute score using lit_save and this_real_er
				flag_continue = 0;
				if (score_mode == 3)
				{
				/*	if(max_weight > weight_limit)
						score_factor = pow(e, -weight_limit/T_em);
					else
						score_factor = pow(e, -max_weight/T_em);
				*/
					if (this_real_er_case1 == 0 && status == 0)
					{
						est_er_acc = cur_real_er;
						score = lit_save * score_factor * (1e+5);
					}
					else if (this_real_er_case1 == 0 && status == 1)
					{	
						est_er_acc = cur_real_er + exp_ratio - cur_real_er * exp_ratio;
					//	est_er_margin = (ini_threshold_er - est_er_acc)/(1 - est_er_acc);
					//	er_part = est_er_acc/est_er_margin;
					//	er_part = exp_ratio;
						if (flag_er_part == 1)
							er_part = est_er_acc/(ini_threshold_er-est_er_acc);
						else if (flag_er_part == 2)
							er_part = exp_ratio;
						else if (flag_er_part == 3)
							er_part = exp_ratio * est_er_acc;
					//	er_part = est_er_acc - cur_real_er;
					//	er_part = er_weight * exp_ratio + (1-er_weight) * est_er_acc; 
					//	er_part = est_er_acc * exp_ratio;
					//	er_part = pow(est_er_acc, er_weight);
						score = lit_save * score_factor / er_part;
					}
					else 
					{
						est_er_acc = cur_real_er + this_real_er_case1 - cur_real_er * this_real_er_case1;
					//	est_er_margin = (ini_threshold_er - est_er_acc)/(1 - est_er_acc);
					//	er_part = est_er_acc/est_er_margin;
					//	er_part = this_real_er_case1;
						if (flag_er_part == 1)
							er_part = est_er_acc/(ini_threshold_er-est_er_acc);
						else if (flag_er_part == 2)
							er_part = this_real_er_case1;
						else if (flag_er_part == 3)
							er_part = this_real_er_case1 * est_er_acc;
					//	er_part = est_er_acc - cur_real_er;
					//	er_part = er_weight * this_real_er_case1 + (1-er_weight) * est_er_acc; 
					//	er_part = est_er_acc * this_real_er_case1;
					//	er_part = pow(est_er_acc, er_weight);
						score = lit_save * score_factor / er_part;
					}
					cout << "er_part = " << er_part << ", this_real_er_case1 = " << this_real_er_case1 << endl;
					cout << "predict: score = " << score << endl;
				}
                                if (score > 0.5 * max_score) flag_continue = 1;
				ftime(&et2);
				double rt_part2 = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
				cout << "runtime for part2: " << rt_part2 << endl;

				ftime(&st2);
                                if (flag_continue)	
				{
                                        int res_case1;
					if (iIndex == 0 && min_error_mag != 0)
					{
						if (min_error_mag <= threshold_em)
							res_case1 = 1;
						else
							res_case1 = 0;
					}
					else
					{
						ftime(&st1);
						res_case1 = check_valid_error_mag(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number, final_pla, BDD_SAT_mode);
						ftime(&et1);
						double rt_part31 = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
						cout << "runtime for part31: " << rt_part31 << endl;
					}

					ftime(&st1);
					if (res_case1)
					{		
						cout << "within threshold_em!" << endl;	
						cout << "min_weight = " << min_weight << endl;
						cout << "inserted score: " << score << endl;
						struct score_pla can_sp;
						can_sp.node = cnode;
						can_sp.score = score;
						can_sp.lit_save = lit_save;
						can_sp.real_er = this_real_er_case1;
					//	can_sp.real_er = est_er_acc;
						can_sp.er_part = er_part;
						can_sp.max_weight = max_weight;	
						can_sp.min_weight = min_weight;	
						can_sp.pla = final_org_pla;
						can_sp.sf = score_factor;
						score_record.insert(make_pair(score, can_sp));
				
						if (flag_pick_max && score > max_score)
						{
							cout << "find new max_score: " << score << endl;
							cout << "lit_save = " << lit_save << ", real_er = " << this_real_er_case1 << ", min_weight = " << min_weight << endl;
							cout << "est_er_acc = " << est_er_acc << endl;
							max_score = score;
							max_sp = can_sp;
						/*	
							max_sp.score = score;
							max_sp.lit_save = lit_save;
							max_sp.real_er = this_real_er_case1;
							max_sp.max_weight = max_weight;
							max_sp.pla = final_org_pla;
						*/
						}
					}
					ftime(&et1);
					double rt_part32 = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
					cout << "runtime for part32: " << rt_part32 << endl;
				}//if (flag_continue == 1)
				else
				{
					cout << "max_score = " << max_score << endl;
					cout << "less than max_score!" << endl;
				}
				ftime(&et2);
				double rt_part3 = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
				cout << "runtime for part3: " << rt_part3 << endl;
			}
			ftime(&et);
			double rt_one_check_case1 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
			rt_one_check_total += rt_one_check_case1;
			cout << "exdc: runtime for one check = " << rt_one_check_case1 << endl;
			
			cout << "case 2. make this whole factor become 1: " << endl;
			num_loop++;
			ftime(&st);
			string cdc(num_input, '-');
			vector<string> org_cubes;
			org_cubes.push_back(cdc);
			all_exdc_minterms.clear();
			minus_cubes(org_cubes, org_pla, all_exdc_minterms);

			ftime(&st2);
			status = simu_real_er(all_exdc_minterms, dont_care_set, pattern_rate, this_real_er_case2);
			ftime(&et2);
			double rt_real_er_case2 = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
			cout << "runtime for simu_real_er: " << rt_real_er_case2 << endl;
			cout << "status = " << status << endl;
			cout << "this_real_er_case2 = " << this_real_er_case2 << ", er_margin = " << er_margin << endl; 
			if (flag_round)
			{
				round_error_rate(this_real_er_case2, 2);
				cout << "after round, this_real_er_case2 = " << this_real_er_case2 << endl;
			}

			if(this_real_er_case2 > mul_ratio * er_margin)	cout << "!error rate is beyond er_margin" << endl;	
			else
			{
				cout << "!error rate is within er_margin" << endl;	
			    	//compute lit_save
			    	sim_org_pla.clear();
				final_org_pla.clear();
				lit_save = get_save_new(net, cnode, 1, sim_org_pla, final_org_pla, iIndex);
		
				//compute score using lit_save and this_real_er
				flag_continue = 0;
				if (score_mode == 3)
				{
				/*	if(max_weight > weight_limit)
						score_factor = pow(e, -weight_limit/T_em);
					else
						score_factor = pow(e, -max_weight/T_em);
				*/
					if (this_real_er_case2 == 0 && status == 0)
					{
						est_er_acc = cur_real_er;
						score = lit_save * score_factor * (1e+5);
					}
					else if (this_real_er_case2 == 0 && status == 1)
					{
						est_er_acc = cur_real_er + exp_ratio - cur_real_er * exp_ratio;
					//	est_er_margin = (ini_threshold_er - est_er_acc)/(1 - est_er_acc);
					//	er_part = est_er_acc/est_er_margin;
					//	er_part = exp_ratio;
						if (flag_er_part == 1)
							er_part = est_er_acc/(ini_threshold_er-est_er_acc);
						else if (flag_er_part == 2)
							er_part = exp_ratio;
						else if (flag_er_part == 3)
							er_part = exp_ratio * est_er_acc;
					//	er_part = est_er_acc - cur_real_er;
					//	er_part = er_weight * exp_ratio + (1-er_weight) * est_er_acc; 
					//	er_part = est_er_acc * exp_ratio;
					//	er_part = pow(est_er_acc, er_weight);
						score = lit_save * score_factor / er_part;
					}
					else 
					{
						est_er_acc = cur_real_er + this_real_er_case2 - cur_real_er * this_real_er_case2;
					//	est_er_margin = (ini_threshold_er - est_er_acc)/(1 - est_er_acc);
					//	er_part = est_er_acc/est_er_margin;
					//	er_part = this_real_er_case2;
						if (flag_er_part == 1)
							er_part = est_er_acc/(ini_threshold_er-est_er_acc);
						else if (flag_er_part == 2)
							er_part = this_real_er_case2;
						else if (flag_er_part == 3)
							er_part = this_real_er_case2 * est_er_acc;
					//	er_part = est_er_acc - cur_real_er;
					//	er_part = er_weight * this_real_er_case2 + (1-er_weight) * est_er_acc; 
					//	er_part = est_er_acc * this_real_er_case2;
					//	er_part = pow(est_er_acc, er_weight);
						score = lit_save * score_factor / er_part;
					}
					cout << "er_part = " << er_part << ", this_real_er_case2 = " << this_real_er_case2 << endl;
				}
				else if (score_mode == 4)
				{
					if (this_real_er_case2 == 0 && status == 0) er_part = 1e-8;
					else if (this_real_er_case2 == 0 && status == 1) 
						er_part = er_weight * exp_ratio/real_er_max;
					else
						er_part = er_weight * this_real_er_case2/real_er_max;
					score = lit_save/er_part;
				}
                                if (score >  0.5 * max_score) flag_continue = 1;
                                
				if (flag_continue)
				{
					int res_case2;
					if (iIndex == 0 && min_error_mag != 0)
					{
						if (min_error_mag <= threshold_em)
							res_case2 = 1;
						else
							res_case2 = 0;
					}
					else
						res_case2 = check_valid_error_mag(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number, final_pla, BDD_SAT_mode);
					if (res_case2)
					{
						cout << "within threshold_em!" << endl;
						cout << "min_weight = " << min_weight << endl;
						if (score_mode == 4)
						{
							if (this_real_er_case2 == 0 && status == 0)	
								score = 1e+8;
							else
								score = lit_save/(er_part + em_part);
							cout << "lit_save = " << lit_save << ", er_part = " << er_part << ", em_part = " << em_part << endl;
						}
						struct score_pla can_sp;
						can_sp.node = cnode;
						can_sp.score = score;
						can_sp.lit_save = lit_save;
						can_sp.real_er = this_real_er_case2;
					//	can_sp.real_er = est_er_acc;
						can_sp.er_part = er_part;
						can_sp.max_weight = max_weight;	
						can_sp.min_weight = min_weight;	
						can_sp.pla = final_org_pla;
						can_sp.sf = score_factor;
						score_record.insert(make_pair(score, can_sp));
						cout << "inserted score: " << score << endl;
				
						if (flag_pick_max && score > max_score)
						{
							cout << "find new max_score: " << score << endl;
							cout << "lit_save = " << lit_save << ", real_er = " << this_real_er_case2 << ", min_weight = " << min_weight << endl;
							cout << "est_er_acc = " << est_er_acc << endl;
							max_score = score;
							max_sp = can_sp;
						/*
							max_sp.score = score;
							max_sp.lit_save = lit_save;
							max_sp.real_er = this_real_er_case2;
							max_sp.max_weight = max_weight;
							max_sp.pla = final_org_pla;
						*/
						}
					}
				}
				else
					cout << "max_score = " << max_score << ", less than max_score!" << endl;
			}			
			ftime(&et);
			double rt_one_check_case2 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
			rt_one_check_total += rt_one_check_case2;
			cout << "exdc: runtime for one check = " << rt_one_check_case2 << endl;
			break;
		}//if(i == num_lit)
		
		
		//case 3: lit_save != num_lit
		int lit_save = i;
		cout << endl << "##lit_save = " << i << endl;
		double num_choice = (factorial(num_lit_org) / factorial(i) ) / factorial(num_lit_org - i);
		int index_vec[num_lit];
		int index_limit[num_lit];
		set<int> over_er_index;	
		//initialize index_vec
		int q = 0, p = 0;
		for(int j = 0; j < i; j++)
		{	
			index_vec[q++] = j;
		    index_limit[p++] = num_lit_org - i + j;
		}
		    
		double num_total = 0;
		int flag_ignore = 0;			
		while(num_total < num_choice)
		{	    	
			num_loop++;
			ftime(&st);
			//step1. obtain index_str  	
		    	cout << endl << "**index_vec: ";
		    	vector<string> node_vec;
			for(int k = 0; k < lit_save; k++)
			{
				cout << index_vec[k];
			    itrm_ib = ini_leaf_set.find(index_vec[k]);
			    string node = itrm_ib->second->data;
			    node_vec.push_back(node);
			}
			cout << endl;
			    
			//step2. restore varilabes and containers
			sim_org_pla.clear();
			for(itrss = ini_current_inv_pla.begin(); itrss != ini_current_inv_pla.end(); itrss++)
				sim_org_pla.push_back(*itrss);
			leaf_set = ini_leaf_set;
			index_inv_pla = ini_index_inv_pla;
			node_index = ini_node_index;
			btNode *root = new btNode;
			root->parent = NULL;
			copyTree(ini_root, root);
			
			//step3. check this simplification
			vector<string> add_pla;
			set<string> delete_pla;
			int lit_save = i;
			int flag_reduce = 0;
			int flag_same_node = 0;
			ftime(&st1);
			for(int k = 0; k < lit_save; k++)
			{
				int index;
				string leaf_node;
				if(k == 0)
				{
					index = index_vec[k];
					leaf_node = node_vec[k];
				}
				else if(k > 0)
				{
					leaf_node = node_vec[k];
				//	cout << "leaf_node: " << leaf_node << endl;
					int num_same = node_index.count(leaf_node);
					if(num_same > 1)
					{
						cout << "same node appears more than once!" << endl;
						flag_same_node = 1;
						break;
					}
					else if(num_same == 0)
						break;
					else
					{
						itrmm_si = node_index.find(leaf_node);
						index = itrmm_si->second;
					}
				}
				itrm_ib = leaf_set.find(index);
				btNode *leaf = itrm_ib->second;    	
				itrm_iv = index_inv_pla.find(index);		
				vector<string> inv_pla = itrm_iv->second;
				add_pla.clear();
				delete_pla.clear();					
				btNode *parent = leaf->parent;
				if(parent == NULL || parent->data == "+")
				{
				//	cout << "remove node " << leaf->data << " by reducing: " << endl;
					flag_reduce = 1;
				}
				else if(parent->data == "*")
				{
				//	cout << "remove node " << leaf->data << " by expanding: " << endl;
					//update add_pla
					vector<string> exp_cubes;
					string true_name;
					int sign;
					get_true_name(leaf->data, true_name, sign);
					itrm_si = name_pos.find(true_name);
					int pos = itrm_si->second;
				    for(int q = 0; q < inv_pla.size(); q++)
				    {
				    	string cstr = inv_pla[q];
				    	cstr[pos] = '-';
				    	exp_cubes.push_back(cstr);
				    }
				    add_pla.insert(add_pla.end(), exp_cubes.begin(), exp_cubes.end());
				}
						
				//set delete_pla = inv_pla
				for(int q = 0; q < inv_pla.size(); q++)
					delete_pla.insert(inv_pla[q]);		
					
				//update sim_org_pla by removing delete_pla and add add_pla
				sim_org_pla_tmp.clear();
				for(int q = 0; q < sim_org_pla.size(); q++)
				{
					itrss = delete_pla.find(sim_org_pla[q]);
					if(itrss == delete_pla.end())
						sim_org_pla_tmp.push_back(sim_org_pla[q]);
				}
				sim_org_pla_tmp.insert(sim_org_pla_tmp.begin(), add_pla.begin(), add_pla.end());
				sim_org_pla = sim_org_pla_tmp;
					
				//update tree for this simplified factor	
				ftime(&st2);		
				int res = removeLeafNode(&root, leaf_node);
				if(res) break;						
				if(k == lit_save - 1)
					break;	
						
				map<int, btNode*> new_leaf_set;
			  	int start_ind = 0;
				visitleaf(root, new_leaf_set, start_ind);	
				comp_exp(&root);
				leaf_set = new_leaf_set;				    
				ftime(&et2);
				double rt_tree = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
    				cout << "@runtime for update_tree: " << rt_tree << endl;    
    				
				index_inv_pla.clear();
				node_index.clear();
				update_index_inv_pla(root, new_leaf_set, name_pos, num_input, index_inv_pla, node_index, current_inv_pla); 
				
			}//for(int k = 0; k < i; k++)
			
			freeTree(root);
			ftime(&et1);
			double rt_step2 = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
    			cout << "@runtime for check_this_sim: " << rt_step2 << endl;      	
			if (flag_same_node)
			{
				//update index_vec
				num_total++;
				int cindex = i - 1;	    	  	
				if(cindex < 0)
					break;	
			    	while(index_vec[cindex] == index_limit[cindex])
				{
					cindex--;
					if(cindex < 0)
						break;
				} 
				index_vec[cindex] += 1;		
				for(int k = cindex+1; k < i; k++)
					index_vec[k] = index_vec[k-1] + 1;	
				continue;
			} 
    			
    			//step4.If all removed nodes cause expansion, then check if this case can be ignored.
    			int flag_ignore = 0;				
			if(!flag_reduce)
				for(int q = 0; q < lit_save; q++)
				{
					itrs0 = over_er_index.find(index_vec[q]);
					if(itrs0 != over_er_index.end())
					{
						cout << "covers some node in over_er_index!" << endl;
						flag_ignore = 1;
						break;
					}
				}				
			if(flag_ignore)
			{
				//update index_vec
				num_total++;
				int cindex = i - 1;	    	  	
				if(cindex < 0)
					break;	
			    while(index_vec[cindex] == index_limit[cindex])
			    {
			    	cindex--;
			    	if(cindex < 0)
			    		break;
			    }
				index_vec[cindex] += 1;		
				for(int k = cindex+1; k < i; k++)
					index_vec[k] = index_vec[k-1] + 1;	
				continue;
			}
								
			//step5. obtain exdc_cubes & sim_org_pla
			ftime(&st1);				
			for(int q = 0; q < org_pla.size(); q++)
			{
				itrss = ini_current_inv_pla.find(org_pla[q]);
				if(itrss == ini_current_inv_pla.end())
					sim_org_pla.push_back(org_pla[q]);
			}
			exdc_cubes.clear();
			find_diff_cubes(org_pla, sim_org_pla, exdc_cubes);
			ftime(&et1);
			double rt_step3 = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
    			cout << "@runtime for obtain_sim_pla: " << rt_step3 << endl;    	
				
			//step6. compute this_real_er & obtain lit_save, score
			double this_real_er_case3;
			if(exdc_cubes.empty())
				this_real_er_case3 = 0;
			else
			{
				all_exdc_minterms.clear();
				for(itrss = exdc_cubes.begin(); itrss != exdc_cubes.end(); itrss++)
				{
					string cube = *itrss;
					exp_cube_set(cube, all_exdc_minterms);
				}
			}	
			
			ftime(&st2);
			status = simu_real_er(all_exdc_minterms, dont_care_set, pattern_rate, this_real_er_case3);
			ftime(&et2);
			double rt_real_er_case3 = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
			cout << "runtime for simu_real_er: " << rt_real_er_case3 << endl;
			cout << "this_real_er_case3 = " << this_real_er_case3 << ", er_margin = " << er_margin << endl; 

			if(this_real_er_case3 > mul_ratio * er_margin)
			{
				cout << "!error rate is beyond er_margin" << endl;
				if(lit_save == 1) over_er_index.insert(index_vec[0]);    		
			}
			else
			{
				cout << "!error rate is within er_margin" << endl;					
				if(sim_org_pla == org_pla)
				{
					cout << "sim_org_pla == org_pla" << endl;
					score = -1;
				}
				else
				{
			    		//compute lit_save
			    		final_org_pla.clear();
					lit_save = get_save_new(net, cnode, this_real_er_case3, sim_org_pla, final_org_pla, iIndex);
						
					//compute score using lit_save and this_real_er
					flag_continue = 0;
					if (score_mode == 3)
					{
					/*	if(max_weight > weight_limit)
							score_factor = pow(e, -weight_limit/T_em);
						else
							score_factor = pow(e, -max_weight/T_em);
					*/
						if (this_real_er_case3 == 0 && status == 0)
						{
							est_er_acc = cur_real_er;
							score = lit_save * score_factor * (1e+5);
						}
						else if (this_real_er_case3 == 0 && status == 1)
						{
							est_er_acc = cur_real_er + exp_ratio - cur_real_er * exp_ratio;
						//	est_er_margin = (ini_threshold_er - est_er_acc)/(1 - est_er_acc);
						//	er_part = est_er_acc/est_er_margin;
						//	er_part = exp_ratio;
							if (flag_er_part == 1)
								er_part = est_er_acc/(ini_threshold_er-est_er_acc);
							else if (flag_er_part == 2)
								er_part = exp_ratio;
							else if (flag_er_part == 3)
								er_part = exp_ratio * est_er_acc;
						//	er_part = est_er_acc - cur_real_er;
						//	er_part = er_weight * exp_ratio + (1-er_weight) * est_er_acc; 
						//	er_part = est_er_acc * exp_ratio;
						//	er_part = pow(est_er_acc, er_weight);
							score = lit_save * score_factor / er_part;
						}
						else 
						{
							est_er_acc = cur_real_er + this_real_er_case3 - cur_real_er * this_real_er_case3;
						//	est_er_margin = (ini_threshold_er - est_er_acc)/(1 - est_er_acc);
						//	er_part = est_er_acc/est_er_margin;
						//	er_part = this_real_er_case3;
							if (flag_er_part == 1)
								er_part = est_er_acc/(ini_threshold_er-est_er_acc);
							else if (flag_er_part == 2)
								er_part = this_real_er_case3;
							else if (flag_er_part == 3)
								er_part = this_real_er_case3 * est_er_acc;
						//	er_part = est_er_acc - cur_real_er;
						//	er_part = er_weight * this_real_er_case3 + (1-er_weight) * est_er_acc; 
						//	er_part = est_er_acc * this_real_er_case3;
						//	er_part = pow(est_er_acc, er_weight);
							score = lit_save * score_factor / er_part;
						}
						cout << "er_part = " << er_part << ", this_real_er_case3 = " << this_real_er_case3 << endl;
					}
					else if (score_mode == 4)
					{
						if (this_real_er_case3 == 0 && status == 0) er_part = 1e-8;
						else if (this_real_er_case3 == 0 && status == 1) 
							er_part = er_weight * exp_ratio/real_er_max;
						else
							er_part = er_weight * this_real_er_case3/real_er_max;
						score = lit_save/er_part;
					}
                                	if (score > 0.5 * max_score) flag_continue = 1;
				
					if (flag_continue)	
				        {
						int res_case3;
						if (iIndex == 0 && min_error_mag != 0)
						{
							if (min_error_mag <= threshold_em)
								res_case3 = 1;
							else
								res_case3 = 0;
						}
						else
			 				res_case3 = check_valid_error_mag(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number, final_pla, BDD_SAT_mode);
						if (res_case3)
						{		
							cout << "within threshold_em!" << endl;			
							cout << "min_weight = " << min_weight << endl;
							if (score_mode == 4)
							{
								if (this_real_er_case3 == 0 && status == 0)	
									score = 1e+8;
								else
									score = lit_save/(er_part + em_part);
								cout << "lit_save = " << lit_save << ", er_part = " << er_part << ", em_part = " << em_part << endl;
							}
							struct score_pla can_sp;
							can_sp.node = cnode;
							can_sp.score = score;
							can_sp.lit_save = lit_save;
							can_sp.real_er = this_real_er_case3;
						//	can_sp.real_er = est_er_acc;
							can_sp.er_part = er_part;
							can_sp.max_weight = max_weight;	
							can_sp.min_weight = min_weight;	
							can_sp.pla = final_org_pla;
							can_sp.sf = score_factor;
							score_record.insert(make_pair(score, can_sp));
							cout << "inserted score: " << score << endl;
						
							if (flag_pick_max && score > max_score)
							{
								cout << "find new max_score: " << score << endl;
								cout << "lit_save = " << lit_save << ", real_er = " << this_real_er_case3 << ", min_weight = " << min_weight << endl; 				
								cout << "est_er_acc = " << est_er_acc << endl;
								max_score = score;
								max_sp = can_sp;
							/*
								max_sp.score = score;
								max_sp.lit_save = lit_save;
								max_sp.real_er = this_real_er_case3;
								max_sp.max_weight = max_weight;
								max_sp.pla = final_org_pla;
							*/
							}
						}
					}//if (flag_continue == 1)
					else
						cout << "max_score = " << max_score << ", less than max_score!" << endl;
				}
			}
				
			//update index_vec
			num_total++;
			int cindex = i - 1;	    	  	
			if(cindex < 0)
				break;	
			while(index_vec[cindex] == index_limit[cindex])
			{
				cindex--;
				if(cindex < 0)
					break;
			}
			index_vec[cindex] += 1;		
			for(int k = cindex+1; k < i; k++)
				index_vec[k] = index_vec[k-1] + 1;								   			      

			ftime(&et);
			double rt_one_check_case3 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
			rt_one_check_total += rt_one_check_case3;
			cout << "exdc: runtime for one check = " << rt_one_check_case3 << endl;
		}//while(num_total < num_choice)		    
	}//for(int i = 1; i <= num_lit; i++)	
	
//	max_modified_po = log(max_sp.max_em)/log(2);    

	cout << "num_loop = " << num_loop << endl;
	cout << "rt_one_check_total = " << rt_one_check_total << endl;
	cout << endl << "max_sp: " << cnode << endl;
	cout << "score : " << max_sp.score << ", real_er: " << max_sp.real_er << ", max_em = " << max_sp.max_em << ", status = " << max_sp.status << endl;
	
}

