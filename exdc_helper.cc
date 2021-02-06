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
#include "head/helper.h"
#include "head/read_file.h"
#include "head/write_func.h"
#include "head/helper.h"
#include "head/exdc_helper.h"
#include "head/loc_sim_main_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"

using namespace std;

extern int numPI_ini;
extern double total_deviation;
extern int num_deviation;

/*
functions in this file:

*/

//Global variables and external variables 


int measure_dist(string cube, string sexdc)
{
	int dist = 0;
	for(int i = 0; i < cube.size(); i++)
	{
		char c1 = cube[i];
		char c2 = sexdc[i];
		if(c1 == '1' && c2 == '0')
			dist += 1;
		else if(c1 == '0' && c2 == '1') 
			dist += 1;
	}
	return dist;
}


int get_save_new(BnetNetwork *net, char *cnode, double sp, vector<string> &sim_org_pla, vector<string> &final_org_pla, int iIndex)
{
	BnetNode *nd;
	multimap<double, string>::iterator itrm_ds;
	struct timeb st, et;
	ftime(&st);
	
/*	cout << "sim_org_pla: " << endl;
	for(int i = 0; i < sim_org_pla.size(); i++)
		cout << sim_org_pla[i] << endl;
*/

	char com[100];   
	//write bigNode.pla & bigNode_sim.pla
	write_bignode_pla_sim(net, cnode, sp, sim_org_pla);
	
	ifstream fin;
    string str, s;
    string filename = "./pla_files/bigNode_sim.pla";
  

	//simplify bigNode_sim.pla			                        
	sprintf(com, "sis -t none -f ./script/sim_pla_bignode.rug > ./pla_files/pla.txt");
	system(com);
	
	//read the simplified pla	
    filename = "./pla_files/bigNode_sim.pla";
    fin.open(filename.c_str(), iostream::in);
//    cout << "bigNode_sim.pla: " << endl;
    while(getline(fin, str))
    {
    //	cout << str << endl;
    	istringstream ss(str);
    	ss >> s;
    	if(s[0] != '.')
    		final_org_pla.push_back(s);
    }
    fin.close();

	//print literal count of bigNode.pla & bigNode_sim.pla
	sprintf(com, "sis -t none -f ./script/print_fac.rug > sis.txt");
	system(com);
	sprintf(com, "sis -t none -f ./script/print_factor.rug");
	system(com);
	int this_area_save;	
	int res = read_sis_result(this_area_save);
	if(res)
	{
		int pol = -1;
	//	st_lookup(net->hash, cnode, &nd);
		if(sp < 0.5)
			pol = 0;
		else
			pol = 1;
		write_ckt_sim_const(net, cnode, pol);		
		if(iIndex == 0)
		{
			//sweep ckt_org.blif & ckt_sim_tmp.blif
			sprintf(com, "sis -t none -f ./script/sim_ckt_sweep_v1.rug > sis_sim.txt");
			system(com);
			sprintf(com, "sis -t none -f ./script/print_fac_whole_v1.rug > sis.txt");
		}
		else
		{
			//sweep ckt_sim.blif & ckt_sim_tmp.blif
			sprintf(com, "sis -t none -f ./script/sim_ckt_sweep_v2.rug > sis_sim.txt");
			system(com);
			sprintf(com, "sis -t none -f ./script/print_fac_whole_v2.rug > sis.txt");
		}
		system(com);
		int res1 = read_sis_result(this_area_save);
	}  
//	cout << "end of get_save_new!" << endl;
	
	ftime(&et);
	double rt = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	cout << "runtime for get_save_new: " << rt << endl;
	return this_area_save;              
}




void insert_lit(string &str, set<int> &pos_x, string &new_str)
{
	set<int>::iterator itrs;
	
	new_str = str;
	for(itrs = pos_x. begin(); itrs != pos_x.end(); itrs++)
	{
		int loc = *itrs;
		new_str.insert(loc, 1, '-');
	}
}


//find all the cubes that cover the given minterm
void find_cover_cubes(string &minterm, int num_input, vector<string> &MSEOP)
{
//	cout << "find_cover_cubes: " << endl;
	set<int> pos_x;
	string new_minterm, old_minterm;
	old_minterm = minterm;
	for(int i = 0; i < minterm.size(); i++)
	{
		if(minterm[i] == '-')
		{
			pos_x.insert(i);
			num_input--;
		}
		else
			new_minterm.append(1, minterm[i]);
	}
	minterm = new_minterm;
	
	int i = 1;	
	while(i <= num_input)
	{
		//check for MSEOPs that save i literals
//		cout << endl << "i = " << i << endl;
		int quotient = factorial(num_input - i);
	    int num_MSEOP_i = (factorial(num_input) / factorial(i) ) / quotient;
		int index_vec[num_input];
		int index_limit[num_input];
		int q = 0, p = 0;
	    for(int j = 0; j < i; j++)
	    {	
	    	index_vec[q++] = j;
	    	index_limit[p++] = num_input - i + j;
	    }

	    int num_total = 0;
	    while(num_total < num_MSEOP_i)
	    {	    		    	
/*	    	cout << "index_vec: " << endl;
		    for(int k = 0; k < i; k++)
		    	cout << index_vec[k];
		    cout << endl;
*/
	    	string nstr;
	    	nstr.assign(minterm);
	    	for(int k = 0; k < i; k++)
	    		nstr[index_vec[k]] = '-';
//	    	cout << "nstr = " << nstr << endl; 
			string new_nstr;
			insert_lit(nstr, pos_x, new_nstr);	  
			MSEOP.push_back(new_nstr);
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
	    }
	    i = i + 1;
	}//end of for loop
	
	minterm = old_minterm;
	return;
}

void exp_cube(string &org_cube, vector<string> &exp_cubes)
{
	vector<int> index_X;
	int numX = 0;
	for(int i = 0; i < org_cube.size(); i++)
		if(org_cube[i] == '-')
		{
			index_X.push_back(i);
			numX++;
		}
	int num_total = (int)pow(2.0, numX);
	for(int i = 0; i < num_total; i++)
	{
		vector<int> bin;
		int2bin(i, bin, numX);
		string str(org_cube);
		for(int j = 0; j < numX; j++)
			str[index_X[j]] = bin[j] + 48;

		exp_cubes.push_back(str);
	}
}

void exp_cube_set(string &org_cube, set<string> &exp_cubes)
{
	vector<int> index_X;
	int numX = 0;
	for(int i = 0; i < org_cube.size(); i++)
		if(org_cube[i] == '-')
		{
			index_X.push_back(i);
			numX++;
		}
	int num_total = (int)pow(2.0, numX);
	for(int i = 0; i < num_total; i++)
	{
		vector<int> bin;
		int2bin(i, bin, numX);
		string str(org_cube);
		for(int j = 0; j < numX; j++)
			str[index_X[j]] = bin[j] + 48;

		exp_cubes.insert(str);
	}
}

//minus the given minterm from the cube org_cube and store the left minterms in new_minterms 
void minus_minterm(string &org_cube, string &minterm, vector<string> &new_minterms)
{
//	cout << "In minus_cube: " << endl;
	//expand the org_cube and remove the minus_cube
	vector<int> index_X;
	int numX = 0;
//	cout << "index_X:" << endl;
	for(int i = 0; i < org_cube.size(); i++)
		if(org_cube[i] == '-')
		{
			index_X.push_back(i);
//			cout << i;
			numX++;
		}
	int num_total = (int)pow(2.0, numX);
	for(int i = 0; i < num_total; i++)
	{
		vector<int> bin;
		int2bin(i, bin, numX);
		string str(org_cube);
		for(int j = 0; j < numX; j++)
			str[index_X[j]] = bin[j] + 48;
			
		if(str != minterm)
			new_minterms.push_back(str);
	}
}


//minus the given cubes from the original org_cubes and store the left minterms in new_minterms 
void minus_cubes(vector<string> &org_cubes, vector<string> &mcubes, set<string> &new_minterms)
{
	struct timeb st, et;
	ftime(&st);
	//expand the org_cube and remove the minterms included in mcubes
	for(int j = 0; j < org_cubes.size(); j++)
	{
		string org_cube = org_cubes[j];
		vector<int> index_X;
		int numX = 0;
		for(int i = 0; i < org_cube.size(); i++)
			if(org_cube[i] == '-')
			{
				index_X.push_back(i);
				numX++;
			}
		int num_total = (int)pow(2.0, numX);
		for(int i = 0; i < num_total; i++)
		{
			vector<int> bin;
			int2bin(i, bin, numX);
			string str(org_cube);
			for(int j = 0; j < numX; j++)
				str[index_X[j]] = bin[j] + 48;
				
			if(!isIncludeVec(mcubes, str))
			{
//				cout << "not included in intersect cubes!" << endl;
				new_minterms.insert(str);
			}
		}
	}
	
	ftime(&et);
    double rt_minus_cubes = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
//	cout << "@runtime for minus_cubes: " << rt_minus_cubes << endl;   
}


int isConflict(string &cube1, string &cube2)
{
	for(int i = 0; i < cube1.size(); i++)
	{
		if((cube1[i] == '0' && cube2[i] == '1') || (cube1[i] == '1' && cube2[i] == '0'))
			return 1;
	}
	return 0;
}


void find_unique_minterm_one(multimap<int, string> &sorted_pla, string &cube, string &unique_minterm)
{
	multimap<int, string>::iterator itrm_is, itrm_is1;
	
	string cstr = cube;
//	cout << endl << "cstr: " << cstr << endl;
	//Find the cubes that are not conflicted with current cube
	vector<string> non_conf_pla;
//	cout << "non_conf_pla: " << endl;
	for(itrm_is1 = sorted_pla.begin(); itrm_is1 != sorted_pla.end(); itrm_is1++)
	{
	    if(itrm_is1->second == cstr)
	    	continue;
	    if(isConflict(cstr, itrm_is1->second))
	    	continue;
	    non_conf_pla.push_back(itrm_is1->second);
//	    cout << itrm_is1->second << endl;
	}
	if(non_conf_pla.size() == 0)
	{
	    int numX = 0;
	    for(int j = 0; j < cstr.size(); j++)
	    	if(cstr[j] == '-')
	    		numX++;
	    if(numX == 0)
	    {
	    	unique_minterm = cstr;
	    	cout << "found um: " << unique_minterm << endl;
	    }	
	}
	//if the number of non-conflict cubes is 1
	else if(non_conf_pla.size() == 1)
	{
	    string ncstr = non_conf_pla[0];
	    int flag_more = 0;
	    int numX = 0;
	    for(int j = 0; j < cstr.size(); j++)
	    {
	    	if(cstr[j] == '-')
	    	{
	    		numX++;
	    		if(ncstr[j] == '1')
	    			cstr[j] = '0';
	    		else if(ncstr[j] == '0')
	    			cstr[j] = '1';
	    	}
	    	if(numX > 1)
	    	{
		    	flag_more = 1;   // there are more than 1 unique minterms for this pla
		    	break;
		    }
	    }
	    if(!flag_more)
	    {	    		
	    	unique_minterm = cstr;
	    	cout << "found um: " << unique_minterm << endl;
	    }
	}
	//if the number of non-conflict cubes > 1
	else
	{
		int flag_more = 0;
	    int numX = 0;
	    for(int j = 0; j < cstr.size(); j++)
	    {
	    	if(cstr[j] == '-')
	    	{
	    		numX++;
	    		for(int i = 0; i < non_conf_pla.size(); i++)
	    		{
	    			string ncstr = non_conf_pla[i];
	    			if(ncstr[j] == '0')
	    			{
	    				cstr[j] = '1';
	    				break;
	    			}
	    			else if(ncstr[j] == '1')
	    			{
	    				cstr[j] = '0';
	    				break;
	    			}
	    		}
	    	}
	    	if(numX > 1)
	    	{
		    	flag_more = 1;   // there are more than 1 unique minterms for this pla
		    	break;
		    }	    			
	   	}
	    if(!flag_more)
	    {
	    	unique_minterm = cstr;
	    	cout << "found um: " << unique_minterm << endl;
	    }
    }
}

void find_unique_minterm_helper(string &cube, string &ncstr, set<string> &unique_minterms)
{
	int flag_more = 0;
	int numX = 0;
	vector<int> pos_x;
	string diff, unique_minterm;
	
	if(isInclude(ncstr, cube))
		return;
	    
	string cstr = cube;
	for(int j = 0; j < cstr.size(); j++)
	{
	    	if(cstr[j] == '-')
	    	{	   
	    		if(ncstr[j] == '0' || ncstr[j] == '1') 		
		    	{
		    		diff.append(1, ncstr[j]);
		    		pos_x.push_back(j);
		    		numX++;
		    		if(ncstr[j] == '1')
		    			cstr[j] = '0';
		    		else if(ncstr[j] == '0')
		    			cstr[j] = '1';
		    	}
	    	}
	    	if(numX > 1)
	    	{
		    	flag_more = 1;   // there are more than 1 unique minterms for this pla
		    	break;
		    }
	}
	    if(!flag_more)
	    {	    		
	    	unique_minterm = cstr;
	    	unique_minterms.insert(unique_minterm);
	    }
	    else
	    {
	    	int num_total = (int)pow(2.0, numX);
	    	cout << "diff = " << diff << endl;
			for(int i = 0; i < num_total; i++)
			{
				vector<int> bin;
				int2bin(i, bin, numX);
				string str;
				for(int j = 0; j < numX; j++)
				{
					char c = bin[j] + 48;
					str.append(1, c);
				}
				cout << "str = " << str << endl;
				if(diff != str)
				{
					string um(cube);
					for(int j = 0; j < diff.size(); j++)
						um[pos_x[j]] = str[j];
					cout << "um = " << um << endl;
					unique_minterms.insert(um);
				}
			}
	    }
}


void find_unique_minterm_multi(multimap<int, string> &sorted_pla, string &cube, set<string> &unique_minterms)
{
	multimap<int, string>::iterator itrm_is, itrm_is1;
	string unique_minterm;
	string cstr = cube;

	//Find the cubes that are not conflicted with current cube
	vector<string> non_conf_pla;
	for(itrm_is1 = sorted_pla.begin(); itrm_is1 != sorted_pla.end(); itrm_is1++)
	{
	    if(itrm_is1->second == cstr)
	    	continue;
	    if(isConflict(cstr, itrm_is1->second))
	    	continue;
	    non_conf_pla.push_back(itrm_is1->second);
	}
	
	if(non_conf_pla.size() == 0)
	    unique_minterms.insert(cube);
	else if(non_conf_pla.size() == 1)
	{
		string ncstr = non_conf_pla[0];
		find_unique_minterm_helper(cube, ncstr, unique_minterms);	    
	}	
	else //if the number of non-conflict cubes > 1
	{
		vector<string> inter_cubes;
		for(int i = 0; i < non_conf_pla.size(); i++)
		{
			string intersect_cube;
			intersect(non_conf_pla[i], cube, intersect_cube);
			inter_cubes.push_back(intersect_cube);
		}
		vector<string> org_cubes;
		org_cubes.push_back(cube);
		minus_cubes(org_cubes, inter_cubes, unique_minterms);	
    }
}


void find_unique_minterms(multimap<int, string> &sorted_pla, vector<string> &unique_minterms)
{
	
//	cout << endl <<  "Coming into find_unique_minterms! " << endl;
	multimap<int, string>::iterator itrm_is, itrm_is1;
	vector<string> non_conf_pla;
	
	string unique_minterm;
	for(itrm_is = sorted_pla.begin(); itrm_is != sorted_pla.end(); itrm_is++)
	{
	    string cstr = itrm_is->second;
	    find_unique_minterm_one(sorted_pla, cstr, unique_minterm);
	    unique_minterms.push_back(unique_minterm);
	}
/*	cout << "unique_minterms: " << endl;
	for(int j = 0; j < unique_minterms.size(); j++)
		cout << unique_minterms[j] << endl;
*/
}




void find_adj_minterms(string &cube, vector<string> &adj_set)
{
	vector<string> tmp_adj_set;
	for(int i = 0; i < cube.size(); i++)
	{
		if(cube[i] != '-')
		{
			string adj_min(cube);
			if(cube[i] == '1')
				adj_min[i] = '0';
			else
				adj_min[i] = '1';
			tmp_adj_set.push_back(adj_min);
		}
	}	
	for(int i = 0; i < tmp_adj_set.size(); i++)
		exp_cube(tmp_adj_set[i], adj_set);
}


void get_exp_cube_helper(string &cpi, string &ccube, vector<string> &exp_cubes)
{
//	cout << "in get_exp_cube_helper, cpi = " << cpi << ", ccube = " << ccube << endl;
	string exp_cube;
	for(int i = 0; i < cpi.size(); i++)
	{
		if(cpi[i] == ccube[i])
			exp_cube.append(1, cpi[i]);
		else
		{
			if(ccube[i] == '0')
				exp_cube.append(1, '1');
			else if(ccube[i] == '1')
				exp_cube.append(1, '0');			
		}
	}
	exp_cubes.push_back(exp_cube);
}

void get_exp_cube(string &cpi, string &ccube, vector<string> &exp_cubes)
{
//	cout << "in get_exp_cube, cpi = " << cpi << ", ccube = " << ccube << endl;
	int num_exp_lit = 0;
	vector<int> loc_x;
	for(int i = 0; i < cpi.size(); i++)
	{
		if(cpi[i] == '-' && ccube[i] != '-')
		{
			num_exp_lit++;	
			loc_x.push_back(i);	
		}
	}
	string last_str = ccube;
	for(int i = 0; i < num_exp_lit; i++)
	{
		int pos = loc_x[i];
		string new_str(last_str);
		new_str[pos] = '-';
		get_exp_cube_helper(new_str, last_str, exp_cubes);
		last_str = new_str;
	}	
}


void intersect(string &pi1, string &pi2, string &inter_cube)
{
	for(int i = 0; i < pi1.size(); i++)
	{
		if(pi1[i] == pi2[i])
			inter_cube.append(1, pi1[i]);
		else 
		{
			if(pi1[i] == '-' && pi2[i] != '-')
				inter_cube.append(1, pi2[i]);
			else if(pi1[i] != '-' && pi2[i] == '-')
				inter_cube.append(1, pi1[i]);		
			else if(pi1[i] != '-' && pi2[i] != '-')	
			{
				inter_cube.clear();
				return;
			}
		}
	}
}


void find_diff_cubes(vector<string> &ini_pla, vector<string> &sim_pla, set<string> &exdc_cubes)
{
	string inter_cube;
	vector<string> intersect_cubes;
	int flag_not_include = 0;
	for(int i = 0; i < ini_pla.size(); i++)
		for(int j = 0; j < sim_pla.size(); j++)
		{
			string pla = sim_pla[j];
			if(!isIncludeVec(ini_pla, pla))
				flag_not_include = 1;
			inter_cube.clear();
			intersect(pla, ini_pla[i], inter_cube);
			if(!inter_cube.empty())
				intersect_cubes.push_back(inter_cube);
		}
	//minus from exp_cube the intersected cubes	
	if(!flag_not_include)
		minus_cubes(ini_pla, sim_pla, exdc_cubes);
	else	
	{
		minus_cubes(ini_pla, intersect_cubes, exdc_cubes);
		minus_cubes(sim_pla, intersect_cubes, exdc_cubes);
	}
}


void find_unique_minterm_exdc(vector<string> &org_pla, vector<string> &exp_cubes, set<string> &exdc_cubes, int flag_type)
{
	//find intersected cubes from exp_cube
	multimap<int, string>::iterator itrm_is;
	vector<string> intersect_cubes;
	string inter_cube;
	vector<string> other_org_pla;
	
	if(flag_type == 0)   
		for(int i = 0; i < org_pla.size(); i++)
		{
			if(search_in_vec(exp_cubes, org_pla[i]))
				continue;
			other_org_pla.push_back(org_pla[i]);
		}
	else if(flag_type == 1)
		other_org_pla = org_pla;
	
//	cout << "inter_cube: " << endl;
	for(int i = 0; i < exp_cubes.size(); i++)
		for(int j = 0; j < other_org_pla.size(); j++)
		{
			string pla = other_org_pla[j];
			inter_cube.clear();
			intersect(pla, exp_cubes[i], inter_cube);
			if(!inter_cube.empty())
			{
				intersect_cubes.push_back(inter_cube);
//				cout << inter_cube << endl;
			}
		}
	//minus from exp_cube the intersected cubes
	minus_cubes(exp_cubes, intersect_cubes, exdc_cubes);
}


/*
void find_cand_reduce(BnetNetwork *net, DdManager **dd, char *cnode, vector<string> &org_pla, vector<char*> &unsort_cutnodes, vector<string> &dont_care,  struct score_pla &this_sp, string &ccube, map<set<string>, double> &pattern_er, map<string, double> &pattern_rate, double threshold, int iIndex)
{
	double this_real_er, score;
	int lit_save;
	BnetNode *nd;
	int num_input = unsort_cutnodes.size();
	
	//sort org_pla
	multimap<int, string> sorted_pla;	
	for(int i = 0; i < org_pla.size(); i++)
	{
	    string str = org_pla[i];
	    int num_digit = 0;
	    for(int j = 0; j < str.size(); j++)
	    	if(str[j] != '-')
	    		num_digit++;
	    sorted_pla.insert(pair<int, string>(num_digit, str));
	}       
	
	st_lookup(net->hash, cnode, &nd);
	set<string> unique_minterms;
//	find_unique_minterm_one(sorted_pla, ccube, unique_minterm);
	find_unique_minterm_multi(sorted_pla, ccube, unique_minterms);

	set<string> all_exdc_minterms;
	set<string>::iterator itrs; 
	for(itrs = unique_minterms.begin(); itrs != unique_minterms.end(); itrs++)
	{
		string cube = *itrs;
		exp_cube_set(cube, all_exdc_minterms);	
	}
	cout << "all_exdc_minterms: " << endl;
	for(itrs = all_exdc_minterms.begin(); itrs != all_exdc_minterms.end(); itrs++)
		cout << *itrs << endl;
	cout << endl;
	
#ifdef use_bdd
	int res = exdc_cubes_real_er(net, dd, dont_care, unsort_cutnodes, all_exdc_minterms, pattern_er, this_real_er, threshold);
	cout << "this_real_er = " << this_real_er << endl;
#endif 

#ifdef use_simu	
	double this_simu_er = 0;	
	simu_real_er(all_exdc_minterms, dont_care, pattern_rate, num_deviation, this_simu_er);
	this_real_er = this_simu_er;
	cout << "this_simu_er = " << this_simu_er << endl; 	
#endif
	
//	double diff = abs(this_simu_er - this_real_er);				
//	total_deviation += diff;
//	cout << "deviation = " << diff << endl;
//	cout << "total_devi = " << total_deviation << endl;

	pattern_er.insert(make_pair(unique_minterms, this_real_er));
	if(this_real_er <= threshold)
	{
		vector<string> sim_org_pla, final_org_pla;
		for(int k = 0; k < org_pla.size(); k++)
		{
			if(org_pla[k] != ccube)
				sim_org_pla.push_back(org_pla[k]);
		}
		cout << "sim_org_pla: " << endl;
		for(int k = 0; k < sim_org_pla.size(); k++)
			cout << sim_org_pla[k] << endl;
		lit_save = get_save_new(net, cnode, 0, sim_org_pla, final_org_pla, iIndex);
		if(this_real_er == 0)
			score = lit_save * 10000000;
		else
			score = lit_save / (this_real_er*(1-this_real_er));
		
		cout << "lit_save = " << lit_save << endl;  
		cout << "this_real_er = " << this_real_er << endl;  					
		cout << "score = " << score << endl;
		if(score > this_sp.score) 
		{
			this_sp.score = score;
			this_sp.real_er = this_real_er;
			this_sp.pla = final_org_pla;
		}
	}
}
*/

/*
void find_cand_exp(BnetNetwork *net, DdManager **dd, char *cnode, vector<string> &org_pla, vector<char*> &unsort_cutnodes, vector<string> &dont_care,  struct score_pla &this_sp, string &ccube, map<set<string>, double> &pattern_er, map<string, double> &pattern_rate, double threshold, int iIndex)
{
	struct timeb st, et, st1, et1, st2, et2;   
	 map<string, double>::iterator itrm_sd;
	 set<string>::iterator itrs;
	 int lit_save;
	 double score, this_real_er;
	 BnetNode *nd;
	 double max_score = this_sp.score;	
	st_lookup(net->hash, cnode, &nd);
	int num_input = ccube.size();
	 	 
	 //check all the adjacent minterms of only minterm
//	 cout << "$check all adjacent minterms" << endl;
	 set<string> over_er_minterm;
	 vector<string> adj_set;	 
//	 find_adj_minterms(ccube, adj_set);
//	 map<string, double> adj_er_set;
//	 for(int i = 0; i < adj_set.size(); i++)
//	 {
//	 	string amin = adj_set[i];    
//	 	if(isIncludeVec(dont_care, amin))
//	 		this_real_er = 0;
//	    else
//	    	this_real_er = exdc_real_er_v2(net, dd, unsort_cutnodes, amin, num_input);
//	    cout << "adjacent minterm: " << amin << ", sp = " << this_real_er << endl;
//	    adj_er_set.insert(make_pair(amin, this_real_er));
//	    if(this_real_er > threshold)
//	    	over_er_minterm.insert(amin);				    			
//	 }

	 //check all MSEOPs of this only minterm: find_cover_cubes
	 cout << "$check all MSEOPs" << endl;
	 vector<string> MSEOP;
	 ftime(&st1);
	 find_cover_cubes(ccube, num_input, MSEOP);
	 ftime(&et1);
	 double rt_cover = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
	 set<string> over_er_cubes;
//     cout << "@runtime for find_cover_cubes: " << rt_cover << endl;

	 set<string> exdc_cubes;    
	 for(int i = 0; i < MSEOP.size(); i++)
	 {
	 	string cpi = MSEOP[i];
	 	cout << endl << "CPI: " << cpi << endl;
	 	//get the cube of the expanded portion: exp_cube
	 	ftime(&st1);
		vector<string> exp_cubes;
		get_exp_cube(cpi, ccube, exp_cubes);
		ftime(&et1);
	 	double rt_get_exp_cube = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
		//cout << "@runtime for get_exp_cube: " << rt_get_exp_cube << endl;    

		//get the cubes that are to be specified as exdcs included in exp_cubes: exdc_cubes
		ftime(&st1);		
		int flag_type = 1;
		exdc_cubes.clear();
		find_unique_minterm_exdc(org_pla, exp_cubes, exdc_cubes, flag_type);
		ftime(&et1);
	 	double rt_find_unique_minterm_exdc = ((et1.time - st1.time)*1000 + (et1.millitm - st1.millitm))/1000.0;
		//cout << "@runtime for find_unique_minterm_exdc: " << rt_find_unique_minterm_exdc << endl;    

		int flag_include = 0;	
		int flag_continue = 0;	    
		if(!flag_include)
		{
			for(itrs = over_er_cubes.begin(); itrs != over_er_cubes.end(); itrs++)
			{
				string cube = *itrs;
				if(isInclude(cpi, cube))
				{
					cout << "this CPI covers some cube in over_er_cubes" << endl;
					flag_continue = 1;
					break;
				}
			}
			if(flag_continue)
				continue;
				
			set<string> all_exdc_minterms;
			if(exdc_cubes.empty())
				this_real_er = 0;
	    	else
	    	{	    		
				for(itrs = exdc_cubes.begin(); itrs != exdc_cubes.end(); itrs++)
				{
					string cube = *itrs;
					exp_cube_set(cube, all_exdc_minterms);
				}
				cout << "all_exdc_minterms: " << endl;
				for(itrs = all_exdc_minterms.begin(); itrs != all_exdc_minterms.end(); itrs++)
					cout << *itrs << endl;
				cout << endl;
#ifdef use_bdd	
				int res	= exdc_cubes_real_er(net, dd, dont_care, unsort_cutnodes, all_exdc_minterms, pattern_er, this_real_er, threshold);
				cout << "this_real_er = " << this_real_er << endl; 
#endif
			}
 
#ifdef use_simu    		
			double this_simu_er = 0;	
			simu_real_er(all_exdc_minterms, dont_care, pattern_rate, num_deviation, this_simu_er);
			this_real_er = this_simu_er;	
			cout << "this_simu_er = " << this_simu_er << endl; 					
#endif	
					
	//		cout << "this_real_er = " << this_real_er << endl; 
			if(this_real_er <= threshold)
			{	
				cout << "!error rate is within threshold" << endl;
				vector<string> sim_org_pla, final_org_pla;
				for(int k = 0; k < org_pla.size(); k++)
				{
					if(org_pla[k] != ccube) 
						sim_org_pla.push_back(org_pla[k]);
				}
				sim_org_pla.push_back(cpi);
				lit_save = get_save_new(net, cnode, 1, sim_org_pla, final_org_pla, iIndex);
				if(this_real_er == 0)
					score = lit_save * 10000000;
				else
					score = lit_save / (this_real_er*(1-this_real_er));
		
				cout << "lit_save = " << lit_save << endl;   
				cout << "this_real_er = " << this_real_er << endl; 					
				cout << "score = " << score << endl; 
				if(score > max_score)
				{
					this_sp.score = score;
					this_sp.real_er = this_real_er;
					this_sp.pla = final_org_pla;
					max_score = score;
				}
			}
			else
			{
				cout << "!error rate is beyond threshold" << endl;
				over_er_cubes.insert(cpi);
			}
		}
	}
}
*/




//function: exdc_cubes_real_er

int exdc_cubes_real_er(BnetNetwork *net, DdManager **dd, vector<string> &dont_care, vector<char*> &cutnodes, set<string> &all_exdc_minterms, map<set<string>, double> &pattern_er, double &this_real_er, double threshold)
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
		if(isIncludeVec(dont_care, minterm))
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





//function: /*exdc_real_er_v2()*/
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
