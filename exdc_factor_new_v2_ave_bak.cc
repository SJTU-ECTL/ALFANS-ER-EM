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
#include "head/sim_new_abc.h"
#include "cudd/bnet.h"
#include "head/loc_sim_main.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/cudd_build_v2.h"
//#include "cudd/ntr_cudd.h"

using namespace std;

extern int numPI_ini;
extern double total_deviation;
extern int num_deviation;
extern int sample_num;
extern double ini_threshold_em;
extern int num_one_po_equal, num_one_po_unequal;
extern map<string, set<char*> > po_tfi;
extern map<string, set<string> > po_inputs;
extern map<string, vector<string> > po_cone_string;

int sample_size = 1;
int conf_prob = 0.95;  //confidence level

#define mode 2    ////mode = 2 represents maximum error magnitude; mode = 1 represents average error magnitude

#define inacc_ratio 0
#define exp_ratio 0.005
#define em_ratio 1
#define em_ratio_single 2
#define score_mode 2
#define ignore_flag 1
#define ignore_po_size_one 1

#define lit_weight 0.1
#define er_weight 0.6
#define em_weight 0.3

/*
functions in this file:

*/

//Global variables and external variables 



int exdc_cubes_real_er(BnetNetwork *net, DdManager **dd, vector<string> &dont_care, vector<string> &local_dc_set, vector<char*> &cutnodes, set<string> &all_exdc_minterms, map<set<string>, double> &pattern_er, double &this_real_er, double threshold)
{
	struct timeb st, et;
	set<string>::iterator itrs, itrs1;
	map<string, double>::iterator itrm_sd;
	map<set<string>, double>::iterator itrm_vd;
	BnetNode *auxnd;
	
	cout << "In exdc_cubes_real_er, all_exdc_minterms: " << endl;
	for(itrs = all_exdc_minterms.begin(); itrs != all_exdc_minterms.end(); itrs++)
		cout << *itrs << endl;
	
	itrm_vd = pattern_er.find(all_exdc_minterms);
	if(itrm_vd != pattern_er.end())
	{
		this_real_er = itrm_vd->second;
		cout << "found match in pattern_er!" << endl;
		return 0;
	}

	ftime(&st);
	ofstream fout;
	fout.open("./pla_files/exdc_cubes.pla", ios::out);
	itrs = all_exdc_minterms.begin();
	string str = *itrs;
	int numin = str.size();
	int numout = 1;
	fout << ".i " << numin << endl;
	fout << ".o " << numout << endl;
	fout << ".ilb ";
	for(int i = 0; i < numin; i++)
		fout << "n" << i << " ";
	fout << endl;
	fout << ".ob out"  << endl;
//	cout << "non-dc exdc_cubes: " << endl;
	for(itrs = all_exdc_minterms.begin(); itrs != all_exdc_minterms.end(); itrs++)
	{
		string minterm = *itrs;
		if(isIncludeVec(dont_care, minterm) || isIncludeVec(local_dc_set, minterm))
		{
			cout << "minterm " << minterm  << " is dc." << endl;
			continue;
		}
		fout << *itrs << " 1" << endl;
	}
	fout << ".e" << endl;
	fout.close();
	char com[100];
	sprintf(com, "sis -t none -f ./script/sim_pla_cube.rug > ./pla_files/pla.txt");
	system(com);
	
	vector<string> exdc_cubes_sim;
	ifstream fin;
	fin.open("./pla_files/exdc_cubes_sim.pla", ios::in);
	string s;
//	cout << "exdc_cubes_sim: " << endl;
	while(getline(fin, str))
	{
		istringstream ss(str);
		ss >> s;
		if(s[0] != '.')
		{
			exdc_cubes_sim.push_back(s);
//			cout << s << endl;
		}
	}
	fin.close();
	ftime(&et);
    double rt_sim_pla = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
//    cout << "@runtime for sim_pla: " << rt_sim_pla << endl;   

    /*Build bdds for don't cares in dc_include*/ 
    int num_exdc_cubes = exdc_cubes_sim.size(); 	
	ftime(&st);
	DdNode *tmp;
    DdNode *func, *prod, *var, *var_tmp;
    func = Cudd_ReadLogicZero(*dd);
	Cudd_Ref(func);
	cout << "exdc_cubes_sim: " << endl;
	for(int i = 0; i < exdc_cubes_sim.size(); i++)
	{		
		prod = Cudd_ReadOne(*dd);
	    Cudd_Ref(prod);
		string dc = exdc_cubes_sim[i];
		cout << dc << endl;	
		for(int i = 0; i < dc.size(); i++)
		{
			if(dc[i] == '-')
				continue;
			char *cnode = cutnodes[i];
			if(!st_lookup((net)->hash, cnode, &auxnd))
			{
				cout << "current node doesn't exixt in st_table!" << endl;
				exit(1);
			}
			if (dc[i] == '1')
				var = auxnd->dd;
		    else 
				var = Cudd_Not(auxnd->dd);

			tmp = Cudd_bddAnd(*dd, prod, var);
			if (tmp == NULL) 
			{
				cout << "tmp is NULL!" << endl;
				exit(1);
			}
			Cudd_Ref(tmp);
			Cudd_IterDerefBdd(*dd, prod);
			prod = tmp;
		}
		tmp = Cudd_bddOr(*dd, func, prod);
	    Cudd_Ref(tmp);
	    Cudd_IterDerefBdd(*dd, func);
	    func = tmp;
	}
	
	double num_minterm = Cudd_CountMinterm(*dd, func, numPI_ini);
	Cudd_IterDerefBdd(*dd, func);
	this_real_er = num_minterm/pow(2.0, numPI_ini);
	pattern_er.insert(make_pair(all_exdc_minterms, this_real_er));			
	ftime(&et);
    double rt_comp_er = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
//    cout << "@runtime for comp_er: " << rt_comp_er << endl;   

	return 0;
}





/*exdc_real_er_v2()*/
double exdc_real_er_v2(BnetNetwork *net, DdManager **dd, vector<char*> &cutnodes, string &sexdc, int num_digit)
{
	double this_real_er;
	BnetNode *nd, *auxnd;	
	DdNode *func;
	if(num_digit == 1)
	{
		char *cnode = cutnodes[0];
		st_lookup(net->hash, cnode, &nd);
		if(sexdc[0] == '1')
			this_real_er = nd->rp;
		else if(sexdc[0] == '0')
			this_real_er = 1 - nd->rp;
		return this_real_er;
	}
	else  
	{
		DdNode *prod = Cudd_ReadOne(*dd);
	    Cudd_Ref(prod);
//	    cout << "sexdc = " << sexdc << endl;
		for(int i = 0; i < sexdc.size(); i++)
		{
			if(sexdc[i] == '-')
				continue;
			char *cnode = cutnodes[i];
//			cout << "cnode: " << cnode << endl;
			if(!st_lookup((net)->hash, cnode, &auxnd))
			{
				cout << "current node doesn't exixt in st_table!" << endl;
				exit(1);
			}
			if(auxnd->dd == NULL)
			{
				cout << "auxnd->dd is NULL!" << endl;
				exit(1);
			}
			DdNode *var;
			if (sexdc[i] == '1')
				var = auxnd->dd;
		    else 
				var = Cudd_Not(auxnd->dd);
			if(var == NULL) 
			{
				cout << "var is NULL!" << endl;
				exit(1);
			}

			DdNode *tmp = Cudd_bddAnd((*dd), prod, var);
			if (tmp == NULL)
			{
				cout << "tmp is NULL!" << endl;
				exit(1);
			}			
			Cudd_Ref(tmp);
			Cudd_IterDerefBdd(*dd, prod);
//			Cudd_IterDerefBdd(*dd, var);
			prod = tmp;
		}

		double num_minterm = Cudd_CountMinterm(*dd, prod, numPI_ini);
//		cout << "num_minterm = " << num_minterm << endl;		
	    Cudd_IterDerefBdd(*dd, prod);
		double this_real_er = num_minterm/pow(2.0, numPI_ini);		
		return this_real_er;
	}
}



void build_tree_from_exp(string &str, btNode **root)
{
   *root = CreateInfixTree(str);
//   int numPT = 0;
//   cout << "compute numPT: " << endl;
//   compNumPT(root, numPT); 
//   cout << "inorder print: " << endl; 
   InOrderPrintTree(*root);
//   cout << endl; 
//   cout << "number of cubes: " << numPT << endl;
   
}


void add_space_star(string &str, string &new_str)
{
	string s;

	//add space around '(' and ')'
		string::size_type n = 0;
		string leftp = "( ";
		string rightp = " )";
		while((n = str.find('(', n)) != string::npos)
		{
			str.replace(n, 1, leftp);
			n += leftp.size();
		}
		n = 0;
		while((n = str.find(')', n)) != string::npos)
		{
			str.replace(n, 1, rightp);
			n += rightp.size();
		}
		//add '*' in the right place
//		cout << "current str: " << str << endl;
		string preStr("*");
		istringstream ss(str);
		string newStr;
		while(ss >> s)
		{
			if(s == "+" || s == ")")
			{
				newStr.append(s);
				newStr.append(" ");
			}
			else
			{
				if(preStr != "*" && preStr != "+" && preStr != "(")
					newStr.append("* ");
				newStr.append(s);
				newStr.append(" ");
			}
			preStr = s;
		}
		
	new_str = newStr;
}


void build_btree(vector<string> &lit_unit, vector<string> &sop_set, vector<string> &factor_set, map<string, btNode*> &factor_exp_trees)
{
	map<string, btNode*>::iterator itrm_sb1;
	
	string str, s;
	for(int i = 0; i < lit_unit.size(); i++)
	{
		str = lit_unit[i];
		string newStr;		
		add_space_star(str, newStr);
//		cout << "new str: " << newStr << endl;		
		if(newStr.find('(') != string::npos || newStr.find(')') != string::npos)
			factor_set.push_back(newStr);
		else
			sop_set.push_back(newStr);
	}
	
	//build a binary tree for each factor expression
//	cout << "factor_set: " << endl;
	for(int i = 0; i < factor_set.size(); i++)	
	{
//		cout << factor_set[i] << endl;
		btNode *root;			
		build_tree_from_exp(factor_set[i], &root);
		factor_exp_trees.insert(pair<string, btNode*>(factor_set[i], root));		
	}

}


int check_ignore_case(BnetNetwork *net, istringstream &ss1, string &cstr)
{
	string s;
	BnetNode *tmp;
	int num_node = 0;
  	while(ss1 >> s)
  	{
  		if(s == "*" || s == "(" || s == ")")
  			continue;
  		num_node++;
  	}
  	if(num_node == 1)  		
  	{
  		string ncstr;
  		string::size_type n = 0;
  		if((n = cstr.find('\'', n)) != string::npos)
  			ncstr = cstr.substr(0, n);
  		else
  			ncstr = cstr;
  		char *cname = new char[50];
  		strcpy(cname, ncstr.c_str());
  		st_lookup(net->hash, cname, &tmp);
  		if(tmp->type == BNET_INPUT_NODE)
  		{
  			cout << "the only node in this sop is an input node, ignore!" << endl;
  			delete []cname;
  			return 1;
  		}
  		else
  			delete []cname;
  	}
  	
  	return 0;

}

void get_true_name(string &s, string &ncstr, int &sign)
{
	//remove trailing spaces
  	int pos_space = -1;
  	for(int j = 0; j < s.size(); j++)
  		if(s[j] == ' ')
  		{
  			pos_space = j;
  			break;
  		}
  	if(pos_space != -1)
  		s = s.substr(0, pos_space);
  		
  	//remove last '\'' if exists
  	string::size_type n = 0;
  	sign = 1;  //positive literal
  	if((n = s.find('\'', n)) != string::npos)
  	{
  		ncstr = s.substr(0, n);
  		sign = 0;  //positive literal
  	}
  	else
  		ncstr = s;
}


void retrieve_cube(istringstream &ss, map<string, int> &name_pos, int num_input, string &ccube)
{
	string s;
	map<string, int>::iterator itrm_si;
	map<int, int> pos_sign;
  	map<int, int>::iterator itrmi;
  	
  	while(ss >> s)
  	{
  		if(s == "*" || s == "(" || s == ")")
  			continue;
  		string ncstr;
  		int sign;
  		get_true_name(s, ncstr, sign);
  		itrm_si = name_pos.find(ncstr);
  		pos_sign.insert(pair<int, int>(itrm_si->second, sign));
  	}
  	for(int i = 0; i < num_input; i++)
  	{
  		itrmi = pos_sign.find(i);
  		if(itrmi == pos_sign.end())
			ccube.append(1, '-');
		else
		{
			int sign = itrmi->second;
			if(sign == 1)
				ccube.append(1, '1');
			else
				ccube.append(1, '0');
		}
  	}
}


void update_index_inv_pla(btNode *root, map<int, btNode*> &leaf_set, map<string, int> &name_pos, int num_input, map<int, vector<string> > &index_inv_pla, multimap<string, int> &node_index, set<string> &current_inv_pla)
{
//	cout << "in update_index_inv_pla: "	<< endl;

	map<int, btNode*>::iterator itrm_ib;
	for(itrm_ib = leaf_set.begin(); itrm_ib != leaf_set.end(); itrm_ib++)
	{
		int index = itrm_ib->first;
		btNode *leaf = itrm_ib->second;
		node_index.insert(make_pair(leaf->data, index));
//		cout << "leaf node: " << leaf->data << endl;
		//obtain index_inv_pla
		vector<string> inv_cubes;
		get_involve_cubes(root, leaf, inv_cubes);
		vector<string> inv_pla;
//		cout << "inv_cubes: " << endl;
		for(int j = 0; j < inv_cubes.size(); j++)
		{
//			cout << inv_cubes[j] << endl;
			string this_inv_pla;
			get_involve_pla(leaf_set, name_pos, num_input, inv_cubes[j], this_inv_pla);
			current_inv_pla.insert(this_inv_pla);
			inv_pla.push_back(this_inv_pla);
		//	cout << this_inv_pla << endl;
		}
		index_inv_pla.insert(make_pair(index, inv_pla));
	}
}



//exdc_factor
void exdc_factor_new_v2(BnetNetwork *net, DdManager **dd, char *cnode, vector<char*> &unsort_cutnodes, map<string, struct score_pla> &sim_record, vector<string> &org_pla, vector<string> &dont_care, vector<string> &local_dc_set, set<char*> &po_set, map<string, vector<string> > &po_dont_care_map, map<string, int> &internal_index, vector<struct index_flag> &input_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int num_output, map<string, map<string, double> > &node_pattern_rate, struct score_pla &max_sp, double threshold_er, double threshold_em, double real_em, int &min_modified_po, int &max_modified_po, int iIndex)
{
    //iterators 
    struct timeb st, et, st1, et1, st2, et2;   
    multimap<double, string>::iterator itrmm_ds; 
    map<string, double>::iterator itrm_sd, itrm_sd1;  
    map<string, struct score_pla>::iterator itrm_ss; 
    multimap<double, vector<string> >::iterator itrm_dv;
    set<int>::iterator itrs0, itrs1, itrs2;
    set<string>::iterator itrss, itrss0, itrss1;
    set<char*>::iterator itrs_char;
    map<string, set<char*> >::iterator itrm_sc;
    map<string, set<string> >::iterator itrm_sss;
    map<string, vector<string> >::iterator itrm_sv;
    map<string, struct wi_pair>::iterator itrm_sw;
    map<int, int>::iterator itrmi;
    
    	//variables
	char com[100];
	string str, s;
	BnetNode *nd, *tmp, *auxnd;	
	max_sp.score = -(1e+9);
	max_sp.lit_save = 0;
	max_sp.real_er = 1;
	max_sp.ave_em = 1e+9;
	max_sp.status = 0;
	double max_score = -(1e+9);
    	int lit_save, num_lines_org_pla = 0, flag_within = 0, cur_em_ratio;
  	double score, this_real_er, max_error_mag = 1;
  	vector<string> final_pla;
  	map<set<string>, double> pattern_er;

	//print out current node's info
	cout << endl << "$current node: " << cnode << endl;
	int num_input = unsort_cutnodes.size();
  	st_lookup(net->hash, cnode, &nd);  	
  	cout << "org_pla: " << endl;
	for(int i = 0; i < org_pla.size(); i++)
	{
	    cout << org_pla[i] << endl;
	    num_lines_org_pla++;
	}
	cout << endl;
	
	//get the factored form of cnode
	cout << endl << "$factor form of bignode: " << endl;		
	write_bignode_pla(net, cnode);   		                        
	sprintf(com, "sis -t none -f ./script/print_factor_org.rug > factor.txt");
	system(com);
			
	//get name_pos
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
	
	//if current node is PO itself and mode = 2 (maximum error magnitude), then check whether its weight is larger than weight_limit
	int weight_limit = log(threshold_em)/log(2);
	
	//if current node has only one affected po, then first get its weight and compute its maximum error magnitude
	int flag_po = 0;
	int min_weight = 1000, max_weight = -1;
	int only_weight = -1;
	if (po_set.size() == 1) 
	{
		cout << "Current node has one affected po!" << endl;
		itrs_char = po_set.begin();
		string po(*itrs_char);
		if(iIndex > 0) po.append("sim");
		itrm_sw = sim_output_wi.find(po);
		int weight = itrm_sw->second.weight;
		min_weight = max_weight = only_weight = weight;
		max_error_mag = pow(2.0, weight);
		
		if (mode == 2)
		{
			st_lookup(net->hash, cnode, &nd);
			if (nd->nfo == 1)
			{
				char *fanout = nd->fanouts[0];
				st_lookup(net->hash, fanout, &tmp);
				if (tmp->type == BNET_OUTPUT_NODE)
				{
					cout << "weight = " << weight << ", weight_limit = " << weight_limit << endl;
					if (weight > weight_limit)
					{
						cout << "cnode is PO itself and its weight is larger than weight_limit!" << endl;
						return;
					}
					else
				}
			}		
		}		
	}
	else
	{
		cout << "Current node has multiple affected pos!" << endl;
		for(itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
		{
			string po(*itrs_char);
			if(iIndex > 0) po.append("sim");
			itrm_sw = sim_output_wi.find(po);
			int weight = itrm_sw->second.weight;
			if (weight < min_weight) min_weight = weight;
			if (weight > max_weight) max_weight = weight;
		}
		max_error_mag = pow(2.0, min_weight);
	}
	
//	cout << "max_error_mag = " << max_error_mag << endl;
	int pre_min_modified_po = min_modified_po;
	int pre_max_modified_po = max_modified_po;
	cout << "min_modified_po = " << min_modified_po << ", max_modified_po = " << max_modified_po << endl;
	cout << "min_weight = " << min_weight << ", max_weight = " << max_weight << endl;
	int this_min_modified_po = (min_weight < min_modified_po)? min_weight: min_modified_po;
	int this_max_modified_po = (max_weight > max_modified_po)? max_weight: max_modified_po;
	min_modified_po = this_min_modified_po;
	max_modified_po = this_max_modified_po;
	cout << "this_min_modified_po = " << this_min_modified_po << ", this_max_modified_po = " << this_max_modified_po << endl;
	if (this_min_modified_po > weight_limit) return;
	
	//read sub_abs.blif and its POs
	ifstream fin;
 	fin.open("./files/sub_abs/sub_abs_mtp.blif", ios::in);
 	string newline;
 	int start = 0;
 	vector<string> sub_abs_ckt;
 	while(getline(fin, str))
 	{
 		istringstream ss(str);
 		ss >> s;
 		if(s == ".end")
 			break;
 		if(s == ".names" && start == 0) start = 1;
		newline.clear();
 		if(start)
 		{
 			if (s == ".names")
	 		{
	 			newline.append(s);
	 			newline.append(" ");
	 			while (ss >> s)
	 			{
	 				s.append("_sub ");
	 				newline.append(s);
	 			}
	 		}
	 		else
	 			newline = str;
 		}	
 		sub_abs_ckt.push_back(newline);
 	}
 	vector<string> sub_abs_pi, sub_abs_po;
 	FILE *fp = fopen("./files/sub_abs/sub_abs_mtp.blif", "r");
    BnetNetwork *net_sub = Bnet_ReadNetwork(fp);
    for (int i = 0; i < net_sub->npis; i++)
    {
    	char *pi = net_sub->inputs[i];
    	string pi_str(pi);
    	pi_str.append("_sub");
    	sub_abs_pi.push_back(pi_str);
    }
    for (int i = 0; i < net_sub->npos; i++)
    {
    	char *po = net_sub->outputs[i];
    	string po_str(po);
    	po_str.append("_sub");
    	sub_abs_po.push_back(po_str);
    }
    fclose(fp);


#ifdef use_record	
	if(!sim_record.empty())
	{
		string cnode_str(cnode);
		itrm_ss = sim_record.find(cnode_str);
		if(itrm_ss != sim_record.end())
		{
			struct score_pla this_score_pla = itrm_ss->second;		  	
			cout << "found match in sim_record: " << endl;
			set<string> old_dc = this_score_pla.dc;
			set<string> new_dc;
			for(int i = 0; i < dont_care.size(); i++)
				new_dc.insert(dont_care[i]);
			for(int i = 0; i < local_dc_set.size(); i++)
				new_dc.insert(local_dc_set[i]);
			set<string> old_insig = this_score_pla.insig;
			final_pla = this_score_pla.pla;	
			string cdc(num_input, '-');
			if(old_dc != new_dc || old_insig != insig)
				cout << "don't-cares are changed or input nodes are changed!" << endl;
			else if(final_pla.empty() || final_pla[0] == cdc)
				cout << "this node is simplified to be 0 or 1!" << endl;
			else
			{
				cout << "don't-cares and input nodes are unchanged!" << endl;
				cout << "this_final_pla: " << endl;
				for(int i = 0; i < final_pla.size(); i++)
			  		cout << final_pla[i] << endl;
				max_sp = this_score_pla;
				return;
			}
		}
	}
#endif

	//obtain pattern_rate for cnode
  	map<string, map<string, double> >::iterator itrm_ssd;
  	string sname(cnode);
	itrm_ssd = node_pattern_rate.find(sname);
	map<string, double> pattern_rate = itrm_ssd->second;
//	cout << "patterns with error rate <= threshold_er: " << endl;
	//If all input patterns have probability larger than threshold_er, than directly return.
	int num_useful = 0;
	if (dont_care.empty())
	{
		for(itrm_sd = pattern_rate.begin(); itrm_sd != pattern_rate.end(); itrm_sd++)
		{
			double er = itrm_sd->second;
			if((er-threshold_er)/threshold_er <= inacc_ratio)
			{
				num_useful++;
				break;
			}
		}
		if(!num_useful)
		{
			max_sp.score = -1;
			max_sp.real_er = 1;
			cout << "no useful patterns!" << endl;
			return; 
		}
	}
  			
	//read the factored form of cnode from file factor.txt
	string ffe;
	read_factor_v2(ffe);
	
	//build a binary tree for the factored-form-expression
	btNode *ini_root;
	string ffe_space;
	add_space_star(ffe, ffe_space);
	build_tree_from_exp(ffe_space, &ini_root);
	cout << "ffe_space = " << ffe_space << endl;
  	
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
	 
	//start check every possible simplification	
	set<int> over_er_index;	
	set<string> exdc_cubes, all_exdc_minterms;
	vector<string> sim_org_pla, sim_org_pla_tmp, final_org_pla;
	double er_devi;
	int num_lit_org = num_lit;
	if(num_lit > 5) num_lit = 5;           //set the limit for literal save
	//ASPDAC modification: ignore the case where lit_save = 1
	for(int i = 2; i <= num_lit; i++)
	{
		if(i == num_lit)
		{
			cout << endl << "##lit_save = " << i << endl;
			cout << "case 1. removing this whole factor: " << endl;
	//		cout << "ini_current_inv_pla: " << endl;
			exdc_cubes.clear();
			for(itrss = ini_current_inv_pla.begin(); itrss != ini_current_inv_pla.end(); itrss++)
			{
	//			cout << *itrs << endl;
				exdc_cubes.insert(*itrss);
			}
			all_exdc_minterms.clear();
			for(itrss = exdc_cubes.begin(); itrss != exdc_cubes.end(); itrss++)
			{
				string cube = *itrss;
				exp_cube_set(cube, all_exdc_minterms);
			}

			double this_simu_er = 0;	
			//status = 1 means that there are non-dc minterms that will cause errors
			int status = simu_real_er(all_exdc_minterms, dont_care, local_dc_set, pattern_rate, num_deviation, this_simu_er);
			this_real_er = this_simu_er;
			cout << "this_simu_er = " << this_real_er << ", threshold_er = " << threshold_er << endl; 
			er_devi = (this_simu_er - threshold_er)/threshold_er;
			if(er_devi > inacc_ratio)
				cout << "!error rate is beyond threshold_er" << endl;	
			else
			{
				cout << "!error rate is within threshold_er" << endl;	
				//compute the averge error magnitude
				if (ignore_flag && this_real_er * max_error_mag > em_ratio_single * threshold_em)
					cout << "out of threshold_em!" << endl;	
			    else
			    {
			    	sim_org_pla.clear();
					for(int k = 0; k < org_pla.size(); k++)
					{
						itrss = ini_current_inv_pla.find(org_pla[k]);
						if(itrss == ini_current_inv_pla.end())
							sim_org_pla.push_back(org_pla[k]);
					}	
					final_org_pla.clear();			
					lit_save = get_save_new(net, cnode, 0, sim_org_pla, final_org_pla, iIndex);
					if (lit_save == 1)
					{
						cout << "lit_save = 1, ignore this benefit." << endl;
						continue;
					}						
					double cur_error_mag;
					if(this_simu_er == 0) 
					{
						cur_error_mag = real_em;
						cout << "po_set: " << endl;
						for(itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
						{
							cout << *itrs_char << " ";
							if(status)
							{
								string po(*itrs_char);
								if(iIndex > 0) po.append("sim");
								itrm_sw = sim_output_wi.find(po);
								int weight = itrm_sw->second.weight;
								cur_error_mag += pow(2.0, weight) * exp_ratio;
							}
						}
						cout << endl;
					}
					else
						cur_error_mag = estimate_error_mag_bdd(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, this_min_modified_po, this_max_modified_po, final_org_pla, threshold_em, weight_limit, mode);
					cout << "cur_error_mag = " << cur_error_mag << endl;
					cout << "threshold_em = " << threshold_em << endl;
					//if this node itself is a PO, then em_ratio = 1
					cur_em_ratio = em_ratio;
					if (flag_po == 1) cur_em_ratio = 1;
					if(cur_error_mag <= cur_em_ratio * threshold_em)
					{
						cout << "within threshold_em!" << endl;
						if (lit_save > 1)
						{
							if (flag_within == 0) flag_within = 1;
							if (score_mode == 1)
							{
								if(this_real_er == 0 && cur_error_mag == 0)
									score = lit_save * 100;
								else					
									score = lit_weight*lit_save/num_lit_org + er_weight*(threshold_er - this_real_er)/threshold_er + em_weight*(threshold_em - cur_error_mag)/threshold_em;
							}
							else if (score_mode == 2)
							{
								if (this_real_er == 0)
								{
									if (cur_error_mag == 0) score = lit_save * (1e+10);
									else score = lit_save / (exp_ratio * cur_error_mag);
								}
								else score = lit_save / (this_real_er*cur_error_mag);
							}

							cout << "lit_save = " << lit_save << endl; 										
							cout << "score = " << score << endl; 
							if(score > max_score)
							{
								max_sp.score = score;
								max_sp.lit_save = lit_save;
								max_sp.real_er = this_real_er;
								max_sp.ave_em = cur_error_mag;
								max_sp.pla = final_org_pla;
								if (status && this_real_er == 0)
									max_sp.status = 1;
								if (!status && this_real_er == 0)
								{
									if (pre_min_modified_po == 1000 && pre_max_modified_po == -1)
									{
										min_modified_po = 1000;
										max_modified_po	= -1;
									}
								}
								vector<string> all_dc = dont_care;
								all_dc.insert(all_dc.end(), local_dc_set.begin(), local_dc_set.end());
								set<string> dc_set;
								for(int k = 0; k < all_dc.size(); k++)
									dc_set.insert(all_dc[k]);
								max_sp.dc = dc_set;
								max_sp.insig = insig;
								max_score = score;
							}
						}
					}
				}
			}
			
			cout << "case 2. make this whole factor become 1: " << endl;
			string cdc(num_input, '-');
			vector<string> org_cubes;
			org_cubes.push_back(cdc);
			all_exdc_minterms.clear();
			minus_cubes(org_cubes, org_pla, all_exdc_minterms);

			this_simu_er = 0;	
			status = simu_real_er(all_exdc_minterms, dont_care, local_dc_set, pattern_rate, num_deviation, this_simu_er);
			this_real_er = this_simu_er;
			cout << "this_simu_er = " << this_real_er << ", threshold_er = " << threshold_er << endl; 
			er_devi = (this_simu_er - threshold_er)/threshold_er;
			if(er_devi > inacc_ratio)
				cout << "!error rate is beyond threshold_er" << endl;	
			else
			{
				cout << "!error rate is within threshold_er" << endl;	
				//compute the averge error magnitude
				if (ignore_flag && this_real_er * max_error_mag > em_ratio_single * threshold_em)
					cout << "out of threshold_em!" << endl;	
			    else
			    {
			    	sim_org_pla.clear();
					final_org_pla.clear();
					lit_save = get_save_new(net, cnode, 1, sim_org_pla, final_org_pla, iIndex);
					if (lit_save == 1)
					{
						cout << "lit_save = 1, ignore this benefit." << endl;
						continue;
					}		
					double cur_error_mag;
					if(this_simu_er == 0) 
					{
						cur_error_mag = real_em;
						cout << "po_set: " << endl;
						for(itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
						{
							cout << *itrs_char << " ";
							if(status)
							{
								string po(*itrs_char);
								if(iIndex > 0) po.append("sim");
								itrm_sw = sim_output_wi.find(po);
								int weight = itrm_sw->second.weight;
								cur_error_mag += pow(2.0, weight) * exp_ratio;
							}
						}
						cout << endl;
					}
					else
						cur_error_mag = estimate_error_mag_bdd(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, this_min_modified_po, this_max_modified_po, final_org_pla, threshold_em, weight_limit, mode);
					cout << "cur_error_mag = " << cur_error_mag << endl;
					cout << "threshold_em = " << threshold_em << endl;
					//if this node itself is a PO, then em_ratio = 1
					cur_em_ratio = em_ratio;
					if (flag_po == 1) cur_em_ratio = 1;
					if(cur_error_mag <= cur_em_ratio * threshold_em)
					{
						cout << "within threshold_em!" << endl;
						if (lit_save > 1)
						{
							if (flag_within == 0) flag_within = 1;	
							if (score_mode == 1)
							{
								if(this_real_er == 0 && cur_error_mag == 0)
									score = lit_save * 100;
								else					
									score = lit_weight*lit_save/num_lit_org + er_weight*(threshold_er - this_real_er)/threshold_er + em_weight*(threshold_em - cur_error_mag)/threshold_em;
							}
							else if (score_mode == 2)
							{
								if (this_real_er == 0)
								{
									if (cur_error_mag == 0) score = lit_save * (1e+10);
									else score = lit_save/(exp_ratio*cur_error_mag);
								}
								else score = lit_save/(this_real_er*cur_error_mag);
							}
							
							cout << "lit_save = " << lit_save << endl; 										
							cout << "score = " << score << endl; 
							if(score > max_score)
							{
								max_sp.score = score;
								max_sp.lit_save = lit_save;
								max_sp.real_er = this_real_er;
								max_sp.ave_em = cur_error_mag;
								max_sp.pla = final_org_pla;
								cout << "status = " << status << ", this_real_er = " << this_real_er << endl;
								if (status && this_real_er == 0)
									max_sp.status = 1;
								if (!status && this_real_er == 0)
								{
									if (pre_min_modified_po == 1000 && pre_max_modified_po == -1)
									{
										min_modified_po = 1000;
										max_modified_po	= -1;
									}
								}
								vector<string> all_dc = dont_care;
								all_dc.insert(all_dc.end(), local_dc_set.begin(), local_dc_set.end());
								set<string> dc_set;
								for(int k = 0; k < all_dc.size(); k++)
									dc_set.insert(all_dc[k]);
								max_sp.dc = dc_set;
								max_sp.insig = insig;
								max_score = score;
							}
						}
					}
				}
			}			
			break;
		}//if(i == num_lit)
		
		int lit_save = i;
		cout << endl << "##lit_save = " << i << endl;
		cout << "num_lit_org = " << num_lit_org << ", i = " << i << endl;
		double num_choice = (factorial(num_lit_org) / factorial(i) ) / factorial(num_lit_org - i);
		cout << "num_choice = " << num_choice << endl;
		int index_vec[num_lit];
		int index_limit[num_lit];
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
		/*	cout << endl << "index_limit: ";
			for (int k = 0; k < i; k++)
				cout << index_limit[k];
			cout << endl;
		*/
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
    			//	cout << "@runtime for update_tree: " << rt_tree << endl;    
    				
				index_inv_pla.clear();
				node_index.clear();
				update_index_inv_pla(root, new_leaf_set, name_pos, num_input, index_inv_pla, node_index, current_inv_pla); 
				
			}//for(int k = 0; k < i; k++)
			
			freeTree(root);
			ftime(&et1);
			double rt_step2 = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
    	//	cout << "@runtime for check_this_sim: " << rt_step2 << endl;      	
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
    //		cout << "@runtime for obtain_sim_pla: " << rt_step3 << endl;    	
				
			//step6. compute this_real_er & obtain lit_save, score
			if(exdc_cubes.empty())
				this_real_er = 0;
			else
			{
				all_exdc_minterms.clear();
				for(itrss = exdc_cubes.begin(); itrss != exdc_cubes.end(); itrss++)
				{
					string cube = *itrss;
					exp_cube_set(cube, all_exdc_minterms);
				}
			}	

			double this_simu_er = 0;	
			int status = simu_real_er(all_exdc_minterms, dont_care, local_dc_set, pattern_rate, num_deviation, this_simu_er);
			this_real_er = this_simu_er;
			cout << "this_simu_er = " << this_real_er << ", threshold_er = " << threshold_er << endl; 
			er_devi = (this_simu_er - threshold_er)/threshold_er;
			if(er_devi > inacc_ratio)
			{
				cout << "!error rate is beyond threshold_er" << endl;
				if(lit_save == 1) over_er_index.insert(index_vec[0]);    		
			}
			else
			{
				cout << "!error rate is within threshold_er" << endl;					
				if(sim_org_pla == org_pla)
				{
					cout << "sim_org_pla == org_pla" << endl;
					score = -1;
				}
				else
				{
					//compute the averge error magnitude
					if (ignore_flag && this_real_er * max_error_mag > em_ratio_single * threshold_em)
						cout << "out of threshold_em!" << endl;	
			    	else
			    	{
			    		final_org_pla.clear();
						lit_save = get_save_new(net, cnode, this_real_er, sim_org_pla, final_org_pla, iIndex);
						double cur_error_mag;
						if(this_simu_er == 0) 
						{
							cur_error_mag = real_em;
							cout << "po_set: " << endl;
							for(itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
							{
								cout << *itrs_char << " ";
								if(status)
								{
									string po(*itrs_char);
									if(iIndex > 0) po.append("sim");
									itrm_sw = sim_output_wi.find(po);
									int weight = itrm_sw->second.weight;
									cur_error_mag += pow(2.0, weight) * exp_ratio;
								}
							}
							cout << endl;
						}
						else
							cur_error_mag = estimate_error_mag_bdd(net, cnode, sub_abs_ckt, sub_abs_pi, sub_abs_po, this_min_modified_po, this_max_modified_po, final_org_pla, threshold_em, weight_limit, mode);
						cout << "cur_error_mag = " << cur_error_mag << endl;
						cout << "threshold_em = " << threshold_em << endl;
						//if this node itself is a PO, then em_ratio = 1
						cur_em_ratio = em_ratio;
						if (flag_po == 1) cur_em_ratio = 1;
						if(cur_error_mag <= cur_em_ratio * threshold_em)
						{		
							cout << "within threshold_em!" << endl;	
							if (lit_save > 1)
							{
								if (flag_within == 0) flag_within = 1;		
								if (score_mode == 1)
								{
									if(this_real_er == 0 && cur_error_mag == 0)
										score = lit_save * 100;
									else					
										score = lit_weight*lit_save/num_lit + er_weight*(threshold_er - this_real_er)/threshold_er + em_weight*(threshold_em - cur_error_mag)/threshold_em;
								}
								else if (score_mode == 2)
								{
									if (this_real_er == 0)
									{
										if (cur_error_mag == 0) score = lit_save * (1e+10);
										else score = lit_save/(exp_ratio*cur_error_mag);
									}
									else score = lit_save/(this_real_er*cur_error_mag);
								}	
								    
								cout << "lit_save = " << lit_save << endl; 									
								cout << "score = " << score << endl; 
								if(score > max_score)
								{
									max_sp.score = score;
									max_sp.lit_save = lit_save;
									max_sp.real_er = this_real_er;
									max_sp.ave_em = cur_error_mag;
									max_sp.pla = final_org_pla;
									if (status && this_real_er == 0)
										max_sp.status = 1;
									if (!status && this_real_er == 0)
									{
										if (pre_min_modified_po == 1000 && pre_max_modified_po == -1)
										{
											min_modified_po = 1000;
											max_modified_po	= -1;
										}
									}
									vector<string> all_dc = dont_care;
									all_dc.insert(all_dc.end(), local_dc_set.begin(), local_dc_set.end());
									set<string> dc_set;
									for(int k = 0; k < all_dc.size(); k++)
										dc_set.insert(all_dc[k]);
									max_sp.dc = dc_set;
									max_sp.insig = insig;
									max_score = score;
								}
							}
						}
					}
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
		}//while(num_total < num_choice)		    
	}//for(int i = 1; i <= num_lit; i++)	    
	
	//run simulation again to get more accurate average error magnitude
	cout << endl << "max_sp: " << cnode << endl;
	cout << "score : " << max_sp.score << ", real_er: " << max_sp.real_er << ", ave_em = " << max_sp.ave_em << ", status = " << max_sp.status << endl;
	
}

