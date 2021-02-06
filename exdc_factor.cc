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
#include "head/sim_new.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/cudd_build.h"
#include "cudd/cudd_comp.h"
#include "cudd/cudd_dst.h"

using namespace std;

extern int numPI_ini;

/*
functions in this file:

*/

//Global variables and external variables 


/*dc_real_er_v2()*/
void dc_real_er_v2(BnetNetwork *net, DdManager **dd, vector<char*> &cutnodes, map<string, double> &dc_include, double &total_er)
{
	set<string>::iterator itrs, itrs1;
	map<string, double>::iterator itrm_sd;
	
    /*Build bdds for don't cares in dc_include*/    	
    double total_dc_er = 0;
	char com[100];
	BnetNode *auxnd;
	DdNode *tmp;
    DdNode *func, *prod, *var, *var_tmp;
	for(itrm_sd = dc_include.begin(); itrm_sd != dc_include.end(); itrm_sd++)
	{
		prod = Cudd_ReadOne(*dd);
	    Cudd_Ref(prod);
		string dc = itrm_sd->first;	
//		cout << endl << "current dc = " << dc << endl;
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
//			cout << "auxnd: " << auxnd->name << ", type = " << auxnd->type << endl;

			if (dc[i] == '1')
				var = auxnd->dd;
		    else 
				var = Cudd_Not(auxnd->dd);

			if(var == NULL) 
			{
				cout << "var is NULL!" << endl;
				exit(1);
			}
			
			tmp = Cudd_bddAnd(*dd, prod, var);
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
//		func = prod;
		double num_minterm = Cudd_CountMinterm(*dd, prod, numPI_ini);
//		cout << "num_minterm = " << num_minterm << endl;
//		Cudd_IterDerefBdd(*dd, func);
	    Cudd_IterDerefBdd(*dd, prod);
		double this_real_er = num_minterm/pow(2.0, numPI_ini);
		itrm_sd->second = this_real_er;
		total_er += this_real_er;
	}
}


int exdc_cubes_real_er(BnetNetwork *net, DdManager **dd, vector<string> &dont_care, vector<string> &local_dc_set, vector<char*> &cutnodes, vector<string> &exdc_cubes, map<vector<string>, double> &pattern_er, double &this_real_er, double threshold)
{
	struct timeb st, et;
	set<string>::iterator itrs, itrs1;
	map<string, double>::iterator itrm_sd;
	map<vector<string>, double>::iterator itrm_vd;
	BnetNode *auxnd;
	
	itrm_vd = pattern_er.find(exdc_cubes);
	if(itrm_vd != pattern_er.end())
	{
		this_real_er = itrm_vd->second;
		cout << "found match in pattern_er!" << endl;
		return 0;
	}
	

	ftime(&st);
	ofstream fout;
	fout.open("./pla_files/exdc_cubes.pla", ios::out);
	string str = exdc_cubes[0];
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
	for(int i = 0; i < exdc_cubes.size(); i++)
	{
		if(isIncludeVec(dont_care, exdc_cubes[i]) || isIncludeVec(local_dc_set, exdc_cubes[i]))
			continue;
		fout << exdc_cubes[i] << " 1" << endl;
//		cout << exdc_cubes[i] << endl;
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

	//estimate the lower bound of the error rate
/*	for(int i = 0; i < exdc_cubes_sim.size(); i++)
	{	
		string exdc = exdc_cubes_sim[i];
		double lb = 1;
		for(int j = 0; j < exdc.size(); j++)
		{
			if(exdc[j] == '-')
				continue;
			char *cnode = cutnodes[j];
			st_lookup((net)->hash, cnode, &auxnd);
			double this_sp;
			if(exdc[j] == '1')
				this_sp = auxnd->rp;
			else
				this_sp = 1 - auxnd->rp;
			lb = lb * this_sp;
		}
//		cout << "estimated lb = " << lb << endl;
		if(lb > threshold)
		{
			this_real_er = lb;
			cout << "final, lb = " << lb << endl;
			return 1;
		}
	
	}	
*/	
	
    /*Build bdds for don't cares in dc_include*/ 
    int num_exdc_cubes = exdc_cubes_sim.size(); 	
	ftime(&st);
	DdNode *tmp;
    DdNode *func, *prod, *var, *var_tmp;
    func = Cudd_ReadLogicZero(*dd);
	Cudd_Ref(func);
//	cout << "exdc_cubes_sim: " << endl;
	for(int i = 0; i < exdc_cubes_sim.size(); i++)
	{		
		prod = Cudd_ReadOne(*dd);
	    Cudd_Ref(prod);
		string dc = exdc_cubes_sim[i];
//		cout << dc << endl;	
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
	pattern_er.insert(make_pair(exdc_cubes, this_real_er));			
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



void build_btree(vector<string> &lit_unit, vector<string> &sop_set, vector<string> &factor_set, map<string, btNode*> &factor_exp_trees)
{
	map<string, btNode*>::iterator itrm_sb1;
	
	string str, s;	
	for(int i = 0; i < lit_unit.size(); i++)
	{
		str = lit_unit[i];
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




//exdc_factor
void exdc_factor(BnetNetwork *net, DdManager **dd, char *cnode, vector<char*> &unsort_cutnodes, map<string, struct score_pla> &sim_record, vector<string> &org_pla, vector<string> &dont_care, vector<string> &local_dc_set, struct score_pla &max_sp, double threshold, int iIndex)
{
//    cout << endl << "Start exdc_factor: " << endl;
    //iterators and variables
    struct timeb st, et, st1, et1, st2, et2;   
    multimap<double, string>::iterator itrmm_ds; 
    map<string, double>::iterator itrm_sd, itrm_sd1;  
    map<string, struct score_pla>::iterator itrm_ss; 
    multimap<double, vector<string> >::iterator itrm_dv;
    
    char com[100];
    string str, s;
    BnetNode *nd, *tmp, *auxnd;	
	max_sp.score = 0;
	max_sp.real_er = 1;
	double max_score;
    int lit_save;
  	double score, this_real_er;
  	vector<string> final_pla;

	//print out current node's info
	cout << endl << "$current node: " << cnode << endl;
	int num_input = unsort_cutnodes.size();
  	st_lookup(net->hash, cnode, &nd);  	
  	cout << "org_pla: " << endl;
	for(int i = 0; i < org_pla.size(); i++)
	    cout << org_pla[i] << endl;
	cout << endl;
	
/*	if(!sim_record.empty())
	{
		string cnode_str(cnode);
		itrm_ss = sim_record.find(cnode_str);
		if(itrm_ss != sim_record.end())
		{
			struct score_pla this_score_pla = itrm_ss->second;		  	
			cout << "searching sim_record: " << endl;
			final_pla = this_score_pla.pla;			  	
			cout << "this_final_pla: " << endl;
			for(int i = 0; i < final_pla.size(); i++)
			  cout << final_pla[i] << endl;
			max_sp = this_score_pla;
			return;
		}
	}
*/
    cout << "unsort_cutnodes: " << endl;
    map<string, int> name_pos;
    map<string, int>::iterator itrm_si;
	for(int i = 0; i < unsort_cutnodes.size(); i++)
	{
	    cout << unsort_cutnodes[i] << " ";
	    st_lookup(net->hash, unsort_cutnodes[i], &tmp);
	    cout << tmp->rp << endl;
	    string name(unsort_cutnodes[i]);
	    name_pos.insert(pair<string, int>(name, i));
	}
	cout << endl;

	//get the factored form of cnode
	cout << endl << "$factor form of bignode: " << endl;		
	write_bignode_pla(net, cnode);   			                        
	sprintf(com, "sis -t none -f ./script/print_factor_org.rug > factor.txt");
	system(com);
	
	
	//read the factored form of cnode from file factor.txt
	vector<string> lit_unit;
	read_factor(lit_unit);
	
	//build a binary tree for each of the element in the factored form
//	cout << endl << "$build binary tree for elements of the factored form: " << endl;
	vector<string> sop_set;
	vector<string> factor_set;
	map<string, btNode *> factor_exp_trees; 
	map<string, btNode *>::iterator itrm_sb, itrm_sb1;
  	build_btree(lit_unit, sop_set, factor_set, factor_exp_trees);
  	
  	//start checking all possible exdcs
  	cout << endl << "**************************************" << endl;
  	cout << "$pick exdcs for sop_set: " << endl;  	
  	map<vector<string>, double > pattern_er;
  	ftime(&st);
  	for(int i = 0; i < sop_set.size(); i++)
  	{  		
  		string cstr = sop_set[i];
  		cout << endl << "current sop: " << cstr << endl;
  		//check cases that can be ignored  		
  		istringstream ss(cstr);
  		istringstream ss1(cstr);
  		int ig = check_ignore_case(net, ss1, cstr);
  		if(ig) continue;
  			 		
  		//retrieve the cube from the expression:  ccube
  		string ccube;
  		retrieve_cube(ss, name_pos, num_input, ccube);
  		cout << "ccube = " << ccube << endl;
  		
		//find unique minterm for current cube and check if its error rate is within threshold
		ftime(&st1);
		find_cand_reduce(net, dd, cnode, org_pla, unsort_cutnodes, dont_care, local_dc_set, max_sp, ccube, pattern_er, threshold, iIndex);
		ftime(&et1);
		double rt_reduce = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
    	cout << "@runtime for reduce: " << rt_reduce << endl;    

		//find MSEOPs of current cube and make qualified ones as candidates
		ftime(&st2);
    	find_cand_exp(net, dd, cnode, org_pla, unsort_cutnodes, dont_care, local_dc_set, max_sp, ccube, pattern_er, threshold, iIndex);    
    	ftime(&et2);
		double rt_expand = ((et2.time - st2.time)*1000 + (et2.millitm - st2.millitm))/1000.0;
    	cout << "@runtime for expand: " << rt_expand << endl;    
    	        	
  	}
  	max_score = max_sp.score;
  	ftime(&et);
    double rt_sop_pick = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    cout << "@runtime for sop pick: " << rt_sop_pick << endl;    

	cout << endl << "**************************************" << endl;
	cout << "$pick exdcs for factor_set: " << endl;  	
  	for(itrm_sb = factor_exp_trees.begin(); itrm_sb != factor_exp_trees.end(); itrm_sb++)
  	{
  		ftime(&st);
  		string cstr = itrm_sb->first;
  		cout << endl << "current factor: " << cstr << endl; 
  		btNode *root = itrm_sb->second;
  		//find all leaf nodes and set their indexes: leaf_set
  		map<int, btNode*> leaf_set;
  		map<int, btNode*>::iterator itrm_ib;
  		int start_ind = 0;
	    visitleaf(root, leaf_set, start_ind);	
	    cout << "leaf_set: " << endl;
	    for(itrm_ib = leaf_set.begin(); itrm_ib != leaf_set.end(); itrm_ib++)
	    	cout << itrm_ib->first << ". " << itrm_ib->second->data << endl;		
  		//compute expression at each node  		
		comp_exp(&root);
		
	    for(itrm_ib = leaf_set.begin(); itrm_ib != leaf_set.end(); itrm_ib++)
	    {
	    	cout << endl << "leaf node: " << itrm_ib->second->data << endl;
	    	btNode *leaf = itrm_ib->second;
	    	vector<string> inv_cubes;
	    	get_involve_cubes(root, leaf, inv_cubes);
	    	vector<string> inv_pla;
	//    	cout << "inv_pla: " << endl;
	    	int num_lit_inv = 0;
	    	for(int j = 0; j < inv_cubes.size(); j++)
	    	{
	    		string this_inv_pla;
	    		get_involve_pla(leaf_set, name_pos, num_input, inv_cubes[j], this_inv_pla);
	    		inv_pla.push_back(this_inv_pla);
	    		for(int k = 0; k < this_inv_pla.size(); k++)
	    			if(this_inv_pla[k] == '0' || this_inv_pla[k] == '1')
	    				num_lit_inv++;
	//    		cout << this_inv_pla << endl;
	    	}
	    	vector<string> exp_cubes;
	    	vector<string> exdc_cubes;
	    	//case 1 : make this variable change to 0
	    	cout << endl << "case 1: make this variable change to 0" << endl;
	    	exp_cubes = inv_pla;
	    	cout << "0. exp_cubes: " << exp_cubes.size() << endl;
			for(int q = 0; q < exp_cubes.size(); q++)
				cout << exp_cubes[q] << endl;
	    	int flag_type = 0;
	    	find_unique_minterm_exdc(org_pla, dont_care, local_dc_set, exp_cubes, exdc_cubes, flag_type);
	    	cout << "0. exdc_cubes: " << exdc_cubes.size() << endl;
			for(int q = 0; q < exdc_cubes.size(); q++)
				cout << exdc_cubes[q] << endl;
			if(exdc_cubes.empty())
				this_real_er = 0;
	    	else
			{
				int res = exdc_cubes_real_er(net, dd, dont_care, local_dc_set, unsort_cutnodes, exdc_cubes, pattern_er, this_real_er, threshold);
				if(res)
					cout << "!lower bound is beyond threshold" << endl;
			}
	     	cout << "this_real_er = " << this_real_er << endl; 			
			if(this_real_er <= threshold)
			{	
				cout << "!error rate is within threshold" << endl;
				vector<string> sim_org_pla, final_org_pla;
				for(int k = 0; k < org_pla.size(); k++)
					if(!search_in_vec(inv_pla, org_pla[k])) 
						sim_org_pla.push_back(org_pla[k]);
				lit_save = get_save_new(net, cnode, sim_org_pla, final_org_pla, iIndex);
				if(this_real_er == 0)
					score = lit_save * 10000000;
				else
					score = lit_save / this_real_er;
				cout << "lit_save = " << lit_save << endl;   					
				cout << "score = " << score << endl; 
				if(score > max_score)
				{
					max_sp.score = score;
					max_sp.real_er = this_real_er;
					max_sp.pla = final_org_pla;
					max_score = score;
				}
			}
			else
				cout << "!error rate is beyond threshold" << endl;
				
				
	    	//case 2: make this variable change to 1
	    	cout << endl << "case 2: make this variable change to 1" << endl;
	    	string name = itrm_ib->second->data;
	    	string true_name;
	    	int sign;
	    	get_true_name(name, true_name, sign);
	    	itrm_si = name_pos.find(true_name);
	    	int pos = itrm_si->second;
	    	vector<string> exp_inv_pla;
	    	for(int j = 0; j < inv_pla.size(); j++)
	    	{
	    		string str(inv_pla[j]);
	    		str[pos] = '-';
	    		exp_inv_pla.push_back(str);
	    	}
	    	exp_cubes = exp_inv_pla;
	    	flag_type = 1;
	    	exdc_cubes.clear();
	    	find_unique_minterm_exdc(org_pla, dont_care, local_dc_set, exp_cubes, exdc_cubes, flag_type);
	    	cout << "1. exdc_cubes: " << exdc_cubes.size() << endl;
			for(int q = 0; q < exdc_cubes.size(); q++)
				cout << exdc_cubes[q] << endl;
			if(exdc_cubes.empty())
				this_real_er = 0;
	    	else
			{
				int res = exdc_cubes_real_er(net, dd, dont_care, local_dc_set, unsort_cutnodes, exdc_cubes, pattern_er, this_real_er, threshold);
				if(res)
					cout << "!lower bound is beyond threshold" << endl;
			}
	     	cout << "this_real_er = " << this_real_er << endl; 
			if(this_real_er <= threshold)
			{	
				cout << "!error rate is within threshold" << endl;
				vector<string> sim_org_pla, final_org_pla;
				for(int k = 0; k < org_pla.size(); k++)
					if(!search_in_vec(inv_pla, org_pla[k])) 
						sim_org_pla.push_back(org_pla[k]);
				for(int k = 0; k < exp_inv_pla.size(); k++)
					sim_org_pla.push_back(exp_inv_pla[k]);

				lit_save = get_save_new(net, cnode, sim_org_pla, final_org_pla, iIndex);
				if(this_real_er == 0)
					score = lit_save * 10000000;
				else
					score = lit_save / this_real_er;
				cout << "lit_save = " << lit_save << endl;   
				cout << "this_real_er = " << this_real_er << endl; 					
				cout << "score = " << score << endl; 
				if(score > max_score)
				{
					max_sp.score = score;
					max_sp.real_er = this_real_er;
					max_sp.pla = final_org_pla;
					max_score = score;
				}
			}
			else
				cout << "!error rate is beyond threshold" << endl;	
	    	
	    }		    	
  	}
  	

/*  if(new_score_map.empty())
  		return -1;
  	itrm_dv = new_score_map.end();
  	itrm_dv--;
  	score = itrm_dv->first;
  	final_pla = itrm_dv->second;
*/
	final_pla = max_sp.pla;
  	cout << "this_final_pla: " << endl;
  	for(int i = 0; i < final_pla.size(); i++)
  		cout << final_pla[i] << endl;
  		
  	
  		
  		
//	return score;

}

