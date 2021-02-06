#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <cmath>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include "head/graph.h"
#include "head/node.h"
#include "head/edge.h"
#include "head/HashTable.h"
#include "head/call_abc.h"
#include "head/basics.h"
#include "head/helper.h"
#include "head/write_func.h"
#include "head/read_file.h"
#include "head/sim.h"
#include "head/exdc_helper.h"
#include "head/loc_sim_main_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/cudd_build_v2.h"
#include "cudd/cudd_comp.h"
#include "cudd/cudd_dst.h"
#include "cudd/bnet.h"
#include "cudd/ntr.h"

using namespace std;

//Global variables and external variables
extern int numPI_ini, numPO_ini;
extern int sample_num;
extern int sample_size;
extern vector<string> ckt_org_po;
extern map<string, set<char*> > po_tfi;
extern map<string, set<string> > po_inputs;
extern map<string, vector<string> > po_cone_string;

/*
functions in this file
1. int extNum(string &s, int flag, int &flag_hyphen)
2. bool remove_vector(vector<Edge*> &NumVec, vector<int> & FlagVec)
3. void find_tranfanin_hash(map<int, Node*> &nodes, HashTable &tfi, int node)
4. void find_tranfanout_hash(map<int, Node*> &nodes, HashTable &tranFanOut, int lineNumber)
*/

//round to the nearest hundredth" or "round to 2 decimal places"
void round_error_rate(double &er, int precision)
{
	if (er < 1e-8) return;
	double er_round;
	double er_tmp = er;
	int num = 0;
	while (er_tmp < pow(10.0, precision-1))
	{
		er_tmp *= 10;
		num++;
	}	
	int er_tmp_int = er_tmp;
	if (er_tmp - er_tmp_int < 0.5)
		er_round = er_tmp_int/pow(10.0, num);
	else
		er_round = (er_tmp_int+1)/pow(10.0, num);	

	er = er_round;
}


void pick_max_sp(multimap<double, struct score_pla> &sim_record_top, double T_em, struct score_pla &max_sp, double &max_score_second)
{
	multimap<double, struct score_pla>::iterator itrmm_ds;

	//step1 update sim_record_top
	if (sim_record_top.empty())
	{
		cout << "sim_record_top is empty!" << endl;
		sim_record_top.insert(make_pair(max_sp.score, max_sp));
		max_score_second = max_sp.score;
		return;
	}
	if (sim_record_top.size() == 1)
	{
		cout << "only one candidate in sim_record_top!" << endl;
		sim_record_top.insert(make_pair(max_sp.score, max_sp));
	}
	else if (sim_record_top.size() == 2)
	{
		cout << "two candidates in sim_record_top!" << endl;
		itrmm_ds = sim_record_top.begin();
		cout << "node1: " << itrmm_ds->second.node << ", score = " << itrmm_ds->first << endl; 
		itrmm_ds = sim_record_top.end();
		itrmm_ds--;
		cout << "node2: " << itrmm_ds->second.node << ", score = " << itrmm_ds->first << endl; 
		itrmm_ds = sim_record_top.begin();
		if (max_sp.score > itrmm_ds->second.score)
		{
			sim_record_top.erase(itrmm_ds);
			sim_record_top.insert(make_pair(max_sp.score, max_sp));
		}	
	}

	//step2 pick the better one from the two candidates
	itrmm_ds = sim_record_top.begin();
	struct score_pla sp1 = itrmm_ds->second;
	itrmm_ds = sim_record_top.end();
	itrmm_ds--;
	struct score_pla sp2 = itrmm_ds->second;
	cout << "sp1: " << sp1.node << endl;
	cout << "score = " << sp1.score << ", lit_save = " << sp1.lit_save << ", real_er = " << sp1.real_er << ", max_weight = " << sp1.max_weight << endl;
	cout << "sp2: " << sp2.node << endl;
	cout << "score = " << sp2.score << ", lit_save = " << sp2.lit_save << ", real_er = " << sp2.real_er << ", max_weight = " << sp2.max_weight << endl;
	if (sp2.score > 1.5 * sp1.score)
	{
		cout << "case1. sp2.score > 1.5 * sp1.score" << endl;
		max_sp = sp2;
	}
	else
	{	
	/*	double score_1_vp = sp1.lit_save/(sp1.real_er+0.003) * pow(2.718, -sp1.max_weight/T_em);	
		double score_2_vp = sp2.lit_save/(sp2.real_er+0.003) * pow(2.718, -sp2.max_weight/T_em);	
		if (score_1_vp <= score_2_vp)
			max_sp = sp2;
		else
		{
			if (sp1.real_er <= 0.003 && sp2.real_er <= 0.003)
			{
				if (sp2.lit_save >= sp1.lit_save) max_sp = sp2;
				else max_sp = sp1;
			}
			else if (sp1.real_er > 0.003 && sp2.real_er < 0.003)
				max_sp = sp2;
			else if (sp1.real_er < 0.003 && sp2.real_er > 0.003)
				max_sp = sp1;
			else
			{
				double score_1_vm = sp1.lit_save/(sp1.real_er-0.003) * pow(2.718, -sp1.max_weight/T_em);
				double score_2_vm = sp2.lit_save/(sp2.real_er-0.003) * pow(2.718, -sp2.max_weight/T_em);	
				if (score_1_vm > score_2_vm)
					max_sp = sp1;
				else 
					max_sp = sp2;
			}
		}
	*/
	
		if (sp2.real_er > 1.5 * sp1.real_er) max_sp = sp1;
		else max_sp = sp2;
		if (max_sp.score == sp1.score) cout << "flipped!" << endl;
	}

	if (max_sp.score == sp1.score) max_score_second = sp2.score;
	else max_score_second = sp1.score;
}




/*find_tranfanout*/
void find_tranfanout(BnetNetwork *net, char *cnode, set<char*> &tfo, set<char*> &po_set, int &affected_po_size)
{
//	cout << "find tfo for node " << cnode << endl;

	BnetNode *nd, *tmp;
	st_lookup(net->hash, cnode, &nd);
	tfo.insert(nd->name);
	
	//If corrent line is a PO, then return;
	if (nd->type == BNET_OUTPUT_NODE) 
        {
		po_set.insert(nd->name);
		affected_po_size++;
		return;
	}
	
//	cout << "nfo = " << nd->nfo << endl;
	for(int i = 0; i < nd->nfo; i++)
	{	
		char *outnode = nd->fanouts[i];
//		cout << "outnode: " << outnode << endl;
		st_lookup(net->hash, outnode, &tmp);
		if(!tmp->visited)
		{
			tfo.insert(tmp->name);
			tmp->visited = 1;
			find_tranfanout(net, outnode, tfo, po_set, affected_po_size);
		}
	}
}


/*find_tranfanout_level*/
void find_tranfanout_level(BnetNetwork *net, char *cnode, int cur_level, set<char*> &tfo, map<char*, int> &po_level, int &affected_po_size)
{
//	cout << "find tfo for node " << cnode << ", cur_level = " << cur_level << endl;

	BnetNode *nd, *tmp;
	st_lookup(net->hash, cnode, &nd);
	tfo.insert(nd->name);
	
	//If corrent line is a PO, then return;
	if (nd->type == BNET_OUTPUT_NODE) 
        {
	//	po_set.insert(nd->name);
		po_level.insert(make_pair(nd->name, cur_level));
		affected_po_size++;
		return;
	}
	
//	cout << "nfo = " << nd->nfo << endl;
	for(int i = 0; i < nd->nfo; i++)
	{	
		char *outnode = nd->fanouts[i];
//		cout << "outnode: " << outnode << endl;
		st_lookup(net->hash, outnode, &tmp);
		if(!tmp->visited)
		{
			tfo.insert(tmp->name);
			tmp->visited = 1;
			find_tranfanout_level(net, outnode, cur_level+1, tfo, po_level, affected_po_size);
		}
	}
}


/*find_tranfanin*/
void find_tranfanin_supp(BnetNetwork *net, char *cnode, set<char*> &tfi, set<char*> &pi_supp)
{
//	cout << "find tfo for node " << cnode << endl;

	BnetNode *nd, *tmp;
	st_lookup(net->hash, cnode, &nd);
	tfi.insert(nd->name);
	
	//If corrent line is a PI, then return;
	if (nd->type == BNET_INPUT_NODE) 
	{
		pi_supp.insert(nd->name);
		return;
	}
	
	for(int i = 0; i < nd->ninp; i++)
	{	
		char *innode = nd->inputs[i];
		st_lookup(net->hash, innode, &tmp);
		if(!tmp->visited)
		{
			tfi.insert(tmp->name);
			tmp->visited = 1;
			find_tranfanin_supp(net, innode, tfi, pi_supp);
		}
	}
}


void find_tranfanin(BnetNetwork *net, char *cnode, set<char*> &tfi)
{
//	cout << "find tfo for node " << cnode << endl;

	BnetNode *nd, *tmp;
	st_lookup(net->hash, cnode, &nd);
	tfi.insert(nd->name);
	
	//If corrent line is a PI, then return;
	if (nd->type == BNET_INPUT_NODE) 
	{
		tfi.insert(nd->name);
		return;
	}
	
	for(int i = 0; i < nd->ninp; i++)
	{	
		char *innode = nd->inputs[i];
		st_lookup(net->hash, innode, &tmp);
		if(!tmp->visited)
		{
			tfi.insert(tmp->name);
			tmp->visited = 1;
			find_tranfanin(net, innode, tfi);
		}
	}
}


void find_tranfanout_po(BnetNetwork *net, char *cnode, set<char*> &po_set)
{
	BnetNode *nd, *tmp;
	st_lookup(net->hash, cnode, &nd);
	
	if (nd->type == BNET_OUTPUT_NODE) 
	{
		po_set.insert(nd->name);
		return;
	}
	
	for(int i = 0; i < nd->nfo; i++)
	{	
		char *outnode = nd->fanouts[i];
		st_lookup(net->hash, outnode, &tmp);
		if(!tmp->visited)
		{
			tmp->visited = 1;
			find_tranfanout_po(net, outnode, po_set);
		}
	}
}



void int2bin(int d, vector<int> &bin, int numBit)
{
	int mod = 0;
	vector<int> tmpstr;
	
	if(d == 0)
	{
		for(int i = 0; i < numBit; i++)
			bin.push_back(0);
		return;	
	}

	while(d > 0)
	{
		mod = d%2;
		d /= 2;
		tmpstr.push_back(mod);
	}
	unsigned int len = tmpstr.size();
	int minus = numBit - len;
	if(minus > 0)
		for(int i = 0; i < minus; i++)
			bin.push_back(0);	
	for(int i = len - 1; i >= 0; i--)
		bin.push_back(tmpstr[i]);
}


void int2bin_reverse(int d, vector<int> &bin, int numBit)
{
	int mod = 0;
	
	if(d == 0)
	{
		for(int i = 0; i < numBit; i++)
			bin.push_back(0);
		return;	
	}

	while(d > 0)
	{
		mod = d % 2;
		d /= 2;
		bin.push_back(mod);
	}
	
	unsigned int len = bin.size();
	int minus = numBit - len;
	if(minus > 0)
		for(int i = 0; i < minus; i++)
			bin.push_back(0);	
}




string int2grey(int indexHor, int nVarsHor)
{
//	cout << "indexHor = " << indexHor << ", nVarsHor = " << nVarsHor << endl;
	vector<int> bin;
	int2bin(indexHor, bin, nVarsHor);
	vector<int> grey;
	for(int i = 0; i < bin.size(); i++)
	{
		int bit;
		if(i > 0) 
			bit = bin[i] ^ bin[i-1];
		else
			bit = bin[i];
		grey.push_back(bit);
	}

	string grey_str;
	char s[1];
	for(int i = 0; i < grey.size(); i++)
	{
		sprintf(s, "%d", grey[i]);
		string str(s);
		grey_str.append(str);
	}
//	cout << "grey_str = " << grey_str << endl;
	return grey_str;
}



int call_vcs()
{
    //delete the old comp.txt first
    char com[100];
    sprintf(com, "rm -rf comp.txt");
    system(com);

    char command[100];
    sprintf(command, "vcs +v2k  ./verilog_files/ckt_org.v ./verilog_files/ckt_sim.v ./verilog_files/ckt_tb.v");
    system(command);
    
    sprintf(command, "./simv +ntb_random_seed=1");
    system(command);

    ifstream fin;
    fin.open("comp.txt", ifstream::in);
    string str;
    int num_diff = 0;
    cout << "comp.txt: " << endl;
    while(getline(fin, str))
    {
		cout << str << endl;
		num_diff++;
    }
    fin.close();
    return num_diff;
}


void call_vcs_org(vector<string> &rand, vector<string> &simu_res)
{
    char command[100];
    sprintf(command, "vcs +v2k  ./verilog_files/ckt_org_simu.v ./verilog_files/ckt_tb_org.v");
//    sprintf(command, "vcs +v2k  ./verilog_files/ckt_org_simu.v ./verilog_files/ckt_tb_org.v +ntb_random_seed=6");
    system(command);
    sprintf(command, "./simv");
//    sprintf(command, "./simv +ntb_random_seed=10");
    system(command);

    ifstream fin;
    fin.open("simu_res.txt", ifstream::in);
    string str;
    int i = 0;
    while(getline(fin, str))
    {
	//	simu_res.insert(pair<int, string>(i, str));		
	//	i++;
		simu_res.push_back(str);
    }
    fin.close();
    
    fin.open("rand.txt", ifstream::in);
    i = 0;
    while(getline(fin, str))
    {
	//	rand.insert(pair<int, string>(i, str));
	//	i++;
		rand.push_back(str);
    }
    fin.close();

    return;
}


void call_vcs_org_after_modify(vector<string> &simu_res)
{
    char command[100];
    sprintf(command, "vcs +v2k  ./verilog_files/ckt_org_simu_copy.v ./verilog_files/ckt_tb_org_after_modify.v");
    system(command);
    sprintf(command, "./simv");
    system(command);

    ifstream fin;
    fin.open("simu_res.txt", ifstream::in);
    string str;
    int i = 0;
    while(getline(fin, str))
    {
		simu_res.push_back(str);
    }
    fin.close();
    
    cout << "end of call_vcs_org_after_modify!" << endl;
    
    return;
}



void call_vcs_sim(vector<string> &rand, vector<string> &simu_res, string &po_sig, string &dut)
{
    char command[200];
//    sprintf(command, "vcs +v2k  ./verilog_files/ckt_org_sim_simu.v ./verilog_files/ckt_tb_sim.v");

	char *cdut = new char[dut.size()+1];
	for(int i = 0; i < dut.size(); i++)
		cdut[i] = dut[i];
	cdut[dut.size()] = '\0';
    sprintf(command, "vcs +v2k  ./verilog_files/%s.v ./verilog_files/ckt_tb_sim.v", cdut);
    delete []cdut;
    
//    sprintf(command, "vcs +v2k  ./verilog_files/ckt_org_sim_simu.v ./verilog_files/ckt_tb_sim.v +ntb_random_seed=8");
    system(command);
    sprintf(command, "./simv");
//    sprintf(command, "./simv +ntb_random_seed=10");
    system(command);

    ifstream fin;
    fin.open("simu_res.txt", ifstream::in);
    string str;
    int i = 0;
    while(getline(fin, str))
    {
	//	simu_res.insert(pair<int, string>(i, str));		
	//	i++;
		simu_res.push_back(str);
    }
    fin.close();
    
    fin.open("rand.txt", ifstream::in);
    i = 0;
    while(getline(fin, str))
    {
	//	rand.insert(pair<int, string>(i, str));
	//	i++;
		rand.push_back(str);
    }
    fin.close();
    
    fin.open("simu_po.txt", ifstream::in);
    while(getline(fin, str))
    {
		po_sig = str;
    }
    fin.close();

    return;
}


void call_vcs_sim_after_modify(vector<string> &simu_res, string &dut)
{
    char command[200];
	char *cdut = new char[dut.size()+1];
	for(int i = 0; i < dut.size(); i++)
		cdut[i] = dut[i];
	cdut[dut.size()] = '\0';
    sprintf(command, "vcs +v2k  ./verilog_files/%s.v ./verilog_files/ckt_tb_sim_after_modify.v", cdut);
    delete []cdut;    
    system(command);
    
    sprintf(command, "./simv");
    system(command);

    ifstream fin;
    fin.open("simu_res.txt", ifstream::in);
    string str;
    int i = 0;
    while(getline(fin, str))
		simu_res.push_back(str);
    fin.close();

    return;
}




void exp_cube(vector<string> &new_pla, string &s)
{
	set<string>::iterator itrs;

	cout << "in exp_cube: s = " << s << endl;
	int num_dc_bit = 0;
    for(int i = 0; i < s.size(); i++)
    {
    	if(s[i] == '-')
        	num_dc_bit++;
    }
    int num_dc = pow(2.0, num_dc_bit);    
	for(int i = 0; i < num_dc; i++)
    {
    	int index_dc_bit = 0;
    	string exp_this_cube;
        for(int j = 0; j < s.size(); j++)
        {
        	if(s[j] != '-')
            	exp_this_cube.append(1, s[j]);
            else
            {
            	vector<int> bin;
            	int2bin(i, bin, num_dc_bit);
            	char cc = bin[index_dc_bit++] + 48;
            	exp_this_cube.append(1, cc);
            }
        }
        new_pla.push_back(exp_this_cube);
    }  
    
}


void build_ddnode_cube(BnetNetwork *net, DdManager **dd, DdNode *prod, vector<char*> &cutnodes, string &sexdc)
{
	BnetNode *auxnd;
	
		for(int i = 0; i < sexdc.size(); i++)
		{
			if(sexdc[i] == '-')
				continue;
			char *cnode = cutnodes[i];
			if(!st_lookup((net)->hash, cnode, &auxnd))
			{
				cout << "current node doesn't exixt in st_table!" << endl;
				exit(1);
			}
			DdNode *var;
			if (sexdc[i] == '1')
				var = auxnd->dd;
		    else 
				var = Cudd_Not(auxnd->dd);
			DdNode *tmp = Cudd_bddAnd((*dd), prod, var);	
			Cudd_Ref(tmp);
			Cudd_IterDerefBdd(*dd, prod);
			prod = tmp;
		}
}


//exdc_minterm
void build_ddnode_sop(BnetNetwork *net, DdManager **dd, DdNode **func, vector<char*> &cutnodes, vector<string> &exdc_set)
{
	BnetNode *auxnd;	
	DdNode *prod, *tmp;
	set<string>::iterator itrs;

	*func = Cudd_ReadLogicZero(*dd);
	Cudd_Ref(*func);
	
	for(int i = 0; i < exdc_set.size(); i++)
	{
//		cout << "build ddnode for cube " << exdc_set[i] << endl;
		prod = Cudd_ReadOne(*dd);
		Cudd_Ref(prod);
		string sexdc = exdc_set[i];
		for(int j = 0; j < sexdc.size(); j++)
		{
			if(sexdc[j] == '-')
				continue;
			char *cnode = cutnodes[j];
//			cout << "current node: " << cnode << endl;
			if(!st_lookup((net)->hash, cnode, &auxnd))
			{
				cout << "current node doesn't exixt in st_table!" << endl;
				exit(1);
			}
			DdNode *var;
			if (sexdc[i] == '1')
				var = auxnd->dd;
		    else 
				var = Cudd_Not(auxnd->dd);
			if(var == NULL)
			{
				cout << "var = NULL!" << endl;
				exit(1);
			}
			DdNode *tmp = Cudd_bddAnd((*dd), prod, var);	
			Cudd_Ref(tmp);
			Cudd_IterDerefBdd(*dd, prod);
			prod = tmp;
		}		
		tmp = Cudd_bddOr(*dd, *func, prod);
	    if (tmp == NULL)
	    {
			cout << "0. tmp is NULL!" << endl;
			exit(1);
		}	
	    Cudd_Ref(tmp);
	    Cudd_IterDerefBdd(*dd, *func);
	    Cudd_IterDerefBdd(*dd, prod);
	    *func = tmp;
	}
	return;
}



//exdc_minterm
double exdc_minterm(BnetNetwork *net, DdManager **dd, vector<char*> &cutnodes, set<string> &diff_pla)
{
	double this_real_er;
	BnetNode *auxnd;	
	DdNode *prod, *func, *tmp;
	set<string>::iterator itrs;

	func = Cudd_ReadLogicZero(*dd);
	Cudd_Ref(func);
	
	for(itrs = diff_pla.begin(); itrs != diff_pla.end(); itrs++)
	{
		string sexdc = *itrs;
		build_ddnode_cube(net, dd, prod, cutnodes, sexdc);		
		tmp = Cudd_bddOr(*dd, func, prod);
	    if (tmp == NULL)
	    {
			cout << "1. tmp is NULL!" << endl;
			exit(1);
		}	
	    Cudd_Ref(tmp);
	    Cudd_IterDerefBdd(*dd,func);
	    Cudd_IterDerefBdd(*dd,prod);
	    func = tmp;
	}
	
	double num_minterm = Cudd_CountMinterm(*dd, func, numPI_ini);
	this_real_er = num_minterm/pow(2.0, numPI_ini);
	return this_real_er;
}




double comp_each_real_er(BnetNetwork *net, DdManager **dd, vector<char*> &cutnodes, vector<string> &input_nodes, vector<string> &ori_pla_v0)
{
	vector<string> new_pla_v0;
	set<string> diff_pla;
	set<string>::iterator itrs;
	map<string, int> new_pla, ori_pla;
	map<string, int>::iterator itrm_si, itrm_si1;  

	ifstream fin;
	fin.open("./pla_files/bigNode_sim.pla", ifstream::in);
	string str;
	int flag_start = 0;
	cout << "read new_pla: " << endl;
	while(getline(fin, str))
	{
		cout << str << endl;
		istringstream ss(str);
	    string s, s_first;
	    ss >> s;
	    s_first = s;
	    if(s == ".p")
	    {
	    	flag_start = 1;
	    	continue;
	    }
	    if(flag_start == 1)
	    {
//	    	cout << "0. s = " << s << endl;
	    	while(ss >> s);
	    	if(s == "2")
	    		break;
//	    	cout << "1. s = " << s << endl;
//	    	cout << "s_first: " << s_first << endl;
	    	exp_cube(new_pla_v0, s_first);    
	    }
    }   	    		
       
    //sort new_pla_v0
    permute(new_pla_v0, input_nodes, cutnodes);
       
    //Find out diff_pla         
    for(int i = 0; i < new_pla_v0.size(); i++)
    	new_pla.insert(pair<string, int>(new_pla_v0[i], -1));
    for(int i = 0; i < ori_pla_v0.size(); i++)
    	ori_pla.insert(pair<string, int>(ori_pla_v0[i], -1));
        
	cout << "new_pla: " << endl;
    for(itrm_si = new_pla.begin(); itrm_si != new_pla.end(); itrm_si++)
    	cout << itrm_si->first << endl;    
    
    cout << endl << "ori_pla: " << endl;

    for(itrm_si = ori_pla.begin(); itrm_si != ori_pla.end(); itrm_si++)
    {
    	cout << itrm_si->first << endl;
        itrm_si1 = new_pla.find(itrm_si->first);
        if(itrm_si1 != new_pla.end())
        {
        	itrm_si->second = 1;
        	itrm_si1->second = 1;	
        }
    }
       
    //Obtain diff_pla 
	for(itrm_si = ori_pla.begin(); itrm_si != ori_pla.end(); itrm_si++)
    {
        if(itrm_si->second == -1)
        	diff_pla.insert(itrm_si->first);
    }
    for(itrm_si = new_pla.begin(); itrm_si != new_pla.end(); itrm_si++)
    {
    	if(itrm_si->second == -1)
        	diff_pla.insert(itrm_si->first);
    }
                	
	//Build bdd to find the minterms corresponding to cubes in diff_pla     	
    cout << "diff_pla: " << endl;
    for(itrs = diff_pla.begin(); itrs != diff_pla.end(); itrs++)
     	cout << *itrs << endl;
     		
    double each_real_er;
    if(diff_pla.empty())
    	each_real_er = 0;
    else
		each_real_er = exdc_minterm(net, dd, cutnodes, diff_pla);

    return each_real_er;
}


void image_computation(BnetNetwork *net, DdManager **dd, DdNode *global_care, char *cnode, vector<string> &local_dc_set)
{
	map<int, int>::iterator itrmi;
	set<string> local_care_set;
	vector<int> index_set;	
	BnetNode *nd, *auxnd;
	st_lookup(net->hash, cnode, &nd);
	
//	cout << "In image_computation: " << endl;
//	cout << "current node: " << nd->name << endl;
/*	cout << "global_care: " << endl;
	Cudd_PrintDebug(*dd, global_care, numPI_ini, 2);
*/	
	int num_var = nd->ninp;
	DdNode **array_funcs = ALLOC(DdNode*, num_var);
		
	for(int i = 0; i < num_var; i++)
	{
		index_set.push_back(i);
		char *innode = nd->inputs[i];
//		cout << "input node: " << innode << endl;
		if(!st_lookup(net->hash, innode, &auxnd))
		{
			cout << "current input is not in hash!" << endl;
			exit(1);
		}		
		if(auxnd->dd == NULL)
		{
			cout << "current input's dd is NULL!" << endl;
			exit(2);
		}			
/*		cout << "input node's dd: " << endl;
		Cudd_PrintDebug(*dd,auxnd->dd, numPI_ini, 2);
		cout << "global_care: " << endl;
		Cudd_PrintDebug(*dd, global_care, numPI_ini, 2);
*/
		DdNode *local_care = Cudd_bddConstrain(*dd, auxnd->dd, global_care);	
		Cudd_Ref(local_care);
//		cout << "constrained node " << i << endl;
//		Cudd_PrintDebug(*dd, local_care, numPI_ini, 2);
		array_funcs[i] = local_care;				
	}
	vector<map<int, int> > lit_cubes;
	image_computation_recur(dd, array_funcs, num_var, index_set, lit_cubes);
//	cout << endl << "finally, lit_cubes: " << endl;
	for(int i = 0; i < lit_cubes.size(); i++)
	{
		string care_str;
		map<int, int> lit_phase = lit_cubes[i];
/*		for(itrmi = lit_phase.begin(); itrmi != lit_phase.end(); itrmi++)
			cout << "(" << itrmi->first << ", " << itrmi->second << ") ";
		cout << endl;
*/
		for(int j = 0; j < nd->ninp; j++)
		{
			itrmi = lit_phase.find(j);
			if(itrmi == lit_phase.end())
				care_str.append(1, '-');
			else if(itrmi->second == 1)
				care_str.append(1, '1');
			else if(itrmi->second == 0)
				care_str.append(1, '0');
		}
//		cout << "care_str: " << care_str << endl;
		local_care_set.insert(care_str);
	}
	
//	cout << "complement_cover: " << endl;
	complement_cover(local_care_set, local_dc_set, nd->ninp);	
	
//	FREE(array_funcs);
}

void image_computation_recur(DdManager **dd, DdNode **array_funcs, int &num_var, vector<int> &index_set, vector<map<int, int> >&lit_cubes)
{

//	cout << endl << "Come into image_computation_recur!" << endl;
	
/*	cout << "index_set: " << endl;
	for(int i = 0; i < index_set.size(); i++)
		cout << index_set[i] << " ";
	cout << endl;
*/
	set<int> const_index;
	set<int>::iterator itrs;
	vector<int> new_index_set;
	map<int, int> const_lit_phase;
	int thisone = -1;
	int count = 0;
		
	//Check constant bdd nodes	
	for(int i = 0; i < num_var; i++)
	{
		DdNode *this_func = array_funcs[i];
		if(this_func == Cudd_ReadLogicZero(*dd) || this_func == Cudd_ReadOne(*dd))
		{
			if(this_func == Cudd_ReadLogicZero(*dd))
			{
				const_lit_phase.insert(pair<int, int>(index_set[i], 0));
//				cout << "node " << index_set[i] << "'s bdd node is constant 0!" << endl;	
			}
			else if(this_func == Cudd_ReadOne(*dd))
			{
				const_lit_phase.insert(pair<int, int>(index_set[i], 1));
//				cout << "node " << index_set[i]<< "'s bdd node is constant 1!" << endl;
			}
			this_func = NULL;
			const_index.insert(index_set[i]);
			continue;
		}
		count++;
		if (thisone == -1 ) 
            thisone = i;         
	}
	lit_cubes.push_back(const_lit_phase);
	
	for(int i = 0; i < index_set.size(); i++)
	{
		itrs = const_index.find(index_set[i]);
		if(itrs == const_index.end())
			new_index_set.push_back(index_set[i]);
	}
	index_set = new_index_set;
	
/*	cout << "new_index_set: " << endl;
	for(int i = 0; i < index_set.size(); i++)
		cout << index_set[i] << " ";
	cout << endl;
*/
	
	//Update array_funcs to new_array_funcs
	int new_num_var = count;
	DdNode **pos_array_funcs = ALLOC(DdNode*, new_num_var);
	DdNode **neg_array_funcs = ALLOC(DdNode*, new_num_var);
//	DdNode **new_array_funcs = ALLOC(DdNode*, count);
	int j = 0;
	for(int i = 0; i < num_var; i++)
	{
/*		itrs = const_index.find(i);
		if(itrs == const_index.end())
		{
			new_array_funcs[j] = array_funcs[i];
			j++;
		}
*/
		if((array_funcs[i] != Cudd_ReadLogicZero(*dd)) && (array_funcs[i] != Cudd_ReadOne(*dd)))
		{
		//	new_array_funcs[j] = array_funcs[i];
			pos_array_funcs[j] = array_funcs[i];
			neg_array_funcs[j] = array_funcs[i];
			j++;
		}
	}
	num_var = new_num_var;	

	//If only one function is non-null, return;
	if (count <= 1) 
	{
		FREE(pos_array_funcs);
		FREE(neg_array_funcs);
		return;
	}
		

	//Output cofactoring: first positive cofactoring then negative cofactoring	
/*	for(int i = 0; i < num_var; i++)
	{
		pos_array_funcs[i] = new_array_funcs[i];
		neg_array_funcs[i] = new_array_funcs[i];
	}
*/
	int pos_num_var = num_var;
	int neg_num_var = num_var;
	vector<int> pos_index_set = index_set;
	vector<int> neg_index_set = index_set;
	vector<map<int, int> > pos_lit_cubes, neg_lit_cubes;
	map<int, int>::iterator itrmi;	
	DdNode *this_func, *new_func;
	//1. positive cofactor
//	cout << "-------------------------------" << endl;
//	cout << "Positive cofactoring over node " << pos_index_set[0] << ": " << endl;
	DdNode *co_func = pos_array_funcs[0];
	Cudd_Ref(co_func);	
	for(int i = 0; i < pos_num_var; i++)
	{
		this_func = pos_array_funcs[i];
//		cout << "node " << pos_index_set[i] << ":" << endl;
//		Cudd_PrintDebug(*dd, this_func, numPI_ini, 3);
		new_func = Cudd_bddConstrain(*dd, this_func, co_func);
		Cudd_Ref(new_func);
//		Cudd_IterDerefBdd(*dd, pos_array_funcs[i]);
		pos_array_funcs[i] = new_func;	
	}
//	pos_array_funcs[0] = Cudd_ReadOne(*dd);
	image_computation_recur(dd, pos_array_funcs, pos_num_var, pos_index_set, pos_lit_cubes);

	//2. negtative cofactor
//	cout << "-------------------------------" << endl;
//	cout << "Negative cofactoring over node " << neg_index_set[0] << ": " << endl;
	DdNode *co_func_neg = Cudd_Not(neg_array_funcs[0]);
	Cudd_Ref(co_func_neg);	
	for(int i = 0; i < neg_num_var; i++)
	{
		this_func = neg_array_funcs[i];
//		cout << "node " << neg_index_set[i] << ":" << endl;
//		Cudd_PrintDebug(*dd, this_func, numPI_ini, 3);
		new_func = Cudd_bddConstrain(*dd, this_func, co_func_neg);
		Cudd_Ref(new_func);
//		Cudd_IterDerefBdd(*dd, neg_array_funcs[i]);
		neg_array_funcs[i] = new_func;	
	}
//	neg_array_funcs[0] = Cudd_ReadLogicZero(*dd);
	image_computation_recur(dd, neg_array_funcs, neg_num_var, neg_index_set, neg_lit_cubes);
	
	FREE(pos_array_funcs);
	FREE(neg_array_funcs);
	

/*	cout << "after this recursion, const_lit_phase: " << endl;
	for(itrmi = const_lit_phase.begin(); itrmi != const_lit_phase.end(); itrmi++)
		cout << "(" << itrmi->first << ", " << itrmi->second << ")" << endl;
	cout << endl;
*/
	lit_cubes.clear();
//	cout << "after this recursion, pos_lit_cubes: " << endl;
	for(int i = 0; i < pos_lit_cubes.size(); i++)
	{
		map<int, int> lit_phase = pos_lit_cubes[i];
/*		for(itrmi = lit_phase.begin(); itrmi != lit_phase.end(); itrmi++)
			cout << "(" << itrmi->first << ", " << itrmi->second << ") ";
		cout << endl;
*/
		map<int, int> new_phase = lit_phase;
		for(itrmi = const_lit_phase.begin(); itrmi != const_lit_phase.end(); itrmi++)
			new_phase.insert(pair<int, int>(itrmi->first, itrmi->second));
		lit_cubes.push_back(new_phase);
	}
//	cout << "after this recursion, neg_lit_cubes: " << endl;
	for(int i = 0; i < neg_lit_cubes.size(); i++)
	{
		map<int, int> lit_phase = neg_lit_cubes[i];
/*		for(itrmi = lit_phase.begin(); itrmi != lit_phase.end(); itrmi++)
			cout << "(" << itrmi->first << ", " << itrmi->second << ") ";
		cout << endl;
*/
		map<int, int> new_phase = lit_phase;
		for(itrmi = const_lit_phase.begin(); itrmi != const_lit_phase.end(); itrmi++)
			new_phase.insert(pair<int, int>(itrmi->first, itrmi->second));
		lit_cubes.push_back(new_phase);
	}
//	cout << "after this recursion, lit_cubes: " << endl;
	for(int i = 0; i < lit_cubes.size(); i++)
	{
		map<int, int> lit_phase = lit_cubes[i];
/*		for(itrmi = lit_phase.begin(); itrmi != lit_phase.end(); itrmi++)
			cout << "(" << itrmi->first << ", " << itrmi->second << ") ";
		cout << endl;
*/
	}
	
	

	
}	


void complement_cover(set<string> &local_care_set, vector<string> &local_dc_set, int num_input)
{
	set<string>::iterator itrs;

	ofstream fout;
	fout.open("./blif_files/care_set.blif", ios::out);
	fout << ".model care_set" << endl;
	fout << ".inputs ";
	for(int i = 0; i < num_input; i++)
		fout << "i" << i << " ";
	fout << endl;
	fout << ".outputs nn" << endl;
	fout << ".names ";
	for(int i = 0; i < num_input; i++)
		fout << "i" << i << " ";
	fout << " nn" << endl;
	for(itrs = local_care_set.begin(); itrs != local_care_set.end(); itrs++)
		fout << *itrs << " 1" << endl;
	fout << ".end" << endl;
	fout.close();
	
	FILE *fp;
	fp = fopen("./blif_files/care_set.blif", "r");
	BnetNetwork *net = Bnet_ReadNetwork(fp);
	fclose(fp);
//	Bnet_PrintNetwork(net);
	
	DdManager *dd = NULL;
    dd = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    if (dd == NULL) exit(1);   
    cudd_build_v2(net, &dd, "./blif_files/care_set.blif", BNET_GLOBAL_DD);
    
    BnetNode *nd;
    char *output = net->outputs[0];
//    cout << "output node: " << output << endl;
    st_lookup(net->hash, output, &nd);
//    Cudd_PrintDebug(dd, nd->dd, num_input, 2);

    fp = fopen("dc_set.log", "w");
    dd->out = fp;
    DdNode *dc = Cudd_Not(nd->dd);
    Cudd_Ref(dc);
    Cudd_PrintMinterm(dd, dc);
    fclose(fp);
    
    Bnet_FreeNetwork(net);
    Cudd_Quit(dd);
    
    ifstream fin;
    fin.open("dc_set.log", ios::in);
    string str, s;
//    cout << "dc_set: " << endl;
    while(getline(fin, str))
    {
//    	cout << str << endl;
    	istringstream ss(str);
    	ss >> s;
    	local_dc_set.push_back(s);
    }

}




/*double sim_pick(BnetNetwork *net, char *cnode, vector<char*> &cutnodes, multimap<double, string> &exdc_set, multimap<double, string> &exdc_record, int iIndex)
{
	multimap<double, string>::iterator n_iter;
	BnetNode *nd;
    
    double this_score = 0, max_score = 0, include_er = 0, max_include_er = 0;
    int this_area_save;
    multimap<double, string> tmp_exdc_set;
    multimap<double, string> max_exdc_set;
    ifstream fin;
    
    st_lookup(net->hash, cnode, &nd);
    int nfo = nd->nfo;
    
    int count = 0;
    for(n_iter = exdc_set.begin(); n_iter != exdc_set.end(); n_iter++)
    {
    	max_include_er += n_iter->first;
    	if(n_iter->first > 0)
    		count++;
    }
//    cout << "max_include_er = " << max_include_er << endl;
//    cout << "count = " << count << endl;
    int j = 0, flag_index = exdc_set.size()-1;
    int flag_index_limit = exdc_set.size()- count;
//    cout << "limit for flag_index = " << flag_index_limit << endl;
    while(1)
    {
    	cout << endl << "j = " << j << endl;
    	include_er = 0;
    	for(n_iter = exdc_set.begin(); n_iter != exdc_set.end(); n_iter++)
    	{
    		include_er += n_iter->first;
    		if(n_iter->first == 0)
    			continue;
    		cout << n_iter->second << ": " << n_iter->first << endl;
    	}    	
    	cout << "include_er = " << include_er << endl;

		//write bigNode.pla
		write_bignode_pla(net, cnode, exdc_set);	
		char com[100];   	                        
		sprintf(com, "espresso ./pla_files/bigNode.pla > ./pla_files/bigNode_sim.pla");
		system(com); 

		sprintf(com, "sis -t none -f ./script/print_fac.rug > sis.txt");
	    system(com);
	//	this_area_save = read_sis_result();   
		read_sis_result(this_area_save);                 
	    cout << "this_area_save = " << this_area_save << endl;
	    
	    if(this_area_save <= 0 && include_er == max_include_er)
	    {
			if(!exdc_record.empty())
	    	{
	    		cout << "try exdc_record: " << endl;
	    		exdc_set = exdc_record;
	    		flag_index = exdc_set.size()-1;
	    		write_bignode_pla(net, cnode, exdc_set);	  	                        
				sprintf(com, "espresso ./pla_files/bigNode.pla > ./pla_files/bigNode_sim.pla");
				system(com); 
				sprintf(com, "sis -t none -f ./script/print_fac.rug > sis.txt");
			    system(com);
		//		this_area_save = read_sis_result(); 
				read_sis_result(this_area_save);               
			    cout << "this_area_save = " << this_area_save << endl;
			    this_score = this_area_save/include_er;
	    	}
	    	else
	    		break;
	    }	    	
	    else if(include_er > 0)
	    	this_score = this_area_save/include_er;
	    else if(include_er == 0 && this_area_save > 0)
	    	this_score = 10000000;
	    else
	    	this_score = 0;

	    	
	    cout << "this_score = " << this_score << endl;
//	    cout << "max_score = " << max_score << endl;
	    
	    if(this_score < max_score)
	    {
//	    	cout << " < max_score !" << endl;
	    	if(include_er == 0)
	    		break;
	    	n_iter = tmp_exdc_set.begin();
//		    cout << "retrieve exdc: " << n_iter->second << endl;
		    exdc_set.insert(pair<double, string>(n_iter->first, n_iter->second));
		    tmp_exdc_set.erase(n_iter);    	
	    		
	    	if(flag_index < flag_index_limit)
	    		break;	    		
		    int k = 0;
	    	for(n_iter = exdc_set.begin(); n_iter != exdc_set.end(); n_iter++, k++)
	    	{
	    		if(k == flag_index)
	    		{
	    			exdc_set.erase(n_iter);
	    			flag_index -= 1;
	    			break;
	    		}
	    	}	    
//	    	cout << "remove exdc: " << n_iter->second << endl;
	    	tmp_exdc_set.insert(pair<double, string>(n_iter->first, n_iter->second));    		    
	    }
	    else if(this_score >= max_score)
	    {
//	    	cout << " >= max_score !" << endl;
	    	max_score = this_score;
	    	max_exdc_set = exdc_set;
	    	if(flag_index < flag_index_limit)
	    		break;
	    	int k = 0;
	    	for(n_iter = exdc_set.begin(); n_iter != exdc_set.end(); n_iter++, k++)
	    	{
	    		if(k == flag_index)
	    		{
	    			exdc_set.erase(n_iter);
	    			flag_index -= 1;
	    			break;
	    		}
	    	} 	    
//	    	cout << "remove exdc: " << n_iter->second << endl;
	    	tmp_exdc_set.insert(pair<double, string>(n_iter->first, n_iter->second));
	    }
	    
	    j++;
	    if(j == pow(2.0, count))
	    	break;
    }
    exdc_set = max_exdc_set;
    cout << "after sim_pick, exdc_set: " << endl;
    for(n_iter = exdc_set.begin(); n_iter != exdc_set.end(); n_iter++)
    {
    	if(n_iter->first == 0)
    		continue;
    	cout << n_iter->second << ": " << n_iter->first << endl;
    }  
    cout << endl;  	
    return max_score;
		        
}
*/



bool isInclude(string &dc_sim, string &dc_org)
{
	for(int i = 0; i < dc_sim.size(); i++)
	{
		if(dc_sim[i] != dc_org[i] && dc_sim[i] != '-')
			return 0;
	}
	return 1;
}


bool isIncludeVec(vector<string> &dc_set, string &dc)
{
	for(int i = 0; i < dc_set.size(); i++)
	{
		if(isInclude(dc_set[i], dc))
			return 1;
	}
	return 0;
}

bool isIncludeSet(set<string> &dc_set, string &dc)
{
	set<string>::iterator itrs;
	for(itrs = dc_set.begin(); itrs != dc_set.end(); itrs++)
	{
		string str = *itrs;
		if(isInclude(str, dc))
			return 1;
	}
	return 0;
}


double factorial(int x)
{
	if(x == 0)
		return 1;
	else
		return (x == 1? x : x * factorial(x-1));
}


int search_in_vec(vector<string> &org_pla, string str)
{
	for(int i = 0; i < org_pla.size(); i++)
		if(str == org_pla[i])
			return 1;
			
	return 0;
}


int search_array_int(int *vec, int num, int d)
{
	for(int i = 0; i < num; i++)
		if(vec[i] == d)
			return 1;
			
	return 0;
}



int simu_real_er(set<string> &all_exdc_minterms, set<string> &dont_care_set, map<string, double> &pattern_rate, double &this_simu_er)
{
	map<string, double>::iterator itrm_sd;
	set<string>::iterator itrs;
	int num_non_dc = 0;
	this_simu_er = 0;
		
//	cout << "all_exdc_minterms: " << all_exdc_minterms.size() << endl;
	for(itrs = all_exdc_minterms.begin(); itrs != all_exdc_minterms.end(); itrs++)
	{
		string pattern = *itrs;
       	//	cout << pattern << endl;
	//	if(isIncludeVec(dont_care, pattern))
	        if (dont_care_set.find(pattern) != dont_care_set.end())
		{
	//		cout << "minterm " << pattern << " is dc." << endl;
			continue;
		}

		num_non_dc++;
		itrm_sd = pattern_rate.find(pattern);
		double rate = 0;
		if(itrm_sd == pattern_rate.end())
			rate = 0;
		else
			rate = itrm_sd->second;
		this_simu_er += rate;
	//	cout << "rate = " << rate << endl;
	//	cout << "pattern = " << pattern << ", rate = " << rate << endl;
	}
	if(num_non_dc > 0)	return 1;
	else return 0;
}


int check_status(set<string> &all_exdc_minterms, set<string> &dont_care_set, map<string, double> &pattern_rate)
{
	map<string, double>::iterator itrm_sd;
	set<string>::iterator itrs;
	int num_non_dc = 0;
		
//	cout << "all_exdc_minterms: " << all_exdc_minterms.size() << endl;
	for(itrs = all_exdc_minterms.begin(); itrs != all_exdc_minterms.end(); itrs++)
	{
		string pattern = *itrs;
	    if (dont_care_set.find(pattern) != dont_care_set.end())
		{
	//		cout << "minterm " << pattern << " is dc." << endl;
			continue;
		}
		num_non_dc++;
		break;
	}
	if(num_non_dc > 0)	return 1;
	else return 0;
}




int read_file(string &fileName)
{
	cout << "read file " << fileName << endl;
	ifstream fin;
	fin.open(fileName.c_str(), ios::in);
	if(!fin)
	{
		cout << "couldn't open the file " << fileName << "!" << endl;
		return 0;
	}
    cout << fileName << ": " << endl;
    string str;
	while(getline(fin, str))
		cout << str << endl;
	cout << endl;
	fin.close();
	
	return 1;
}


double estimate_ave_error_mag(set<char*> &po_set, map<string, vector<string> > &po_dont_care_map, map<string, struct wi_pair> &sim_output_wi, vector<struct index_flag> &input_index, vector<string> &rand, vector<string> &simu_res, vector<string> &dont_care, vector<string> &local_dc_set, set<string> &all_exdc_minterms, int iIndex)
{
	//iterators
	set<string>::iterator itrss;
	set<char*>::iterator itrs_char;
	map<string, vector<string> >::iterator itrm_sv;
	map<string, struct wi_pair>::iterator itrm_sw;
    map<int, int>::iterator itrmi;
    map<string, int>::iterator itrm_si;
    double cur_ave_error_mag = 0;

	cout << endl << "enter into estimate_ave_error_mag!" << endl;
/*	cout << "po_dont_care_map: " << endl;
	for (itrm_sv = po_dont_care_map.begin(); itrm_sv != po_dont_care_map.end(); itrm_sv++)
	{
		string po = itrm_sv->first;
		vector<string> dc = itrm_sv->second;
		cout << "po = " << po << endl;
		for(int i = 0; i < dc.size(); i++)
			cout << dc[i] << endl;
	}
*/
	for(itrss = all_exdc_minterms.begin(); itrss != all_exdc_minterms.end(); itrss++)
	{
		string vec = *itrss;
		cout << "$check vector " << vec << endl;
		vector<string> aff_po_vec; //store the POs that are affected by this vector
		if(isIncludeVec(dont_care, vec) || isIncludeVec(local_dc_set, vec))
			continue;
		for(itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
		{
			char *po = *itrs_char;
			cout << "$affected po: " << po << endl;
			string str_po(po);
			itrm_sv = po_dont_care_map.find(str_po);
			vector<string> dont_care_po = itrm_sv->second;
			if(dont_care_po.size() > 0) cout << "non empty!" << endl;
			for(int i = 0; i < dont_care_po.size(); i++)
				cout << dont_care_po[i] << endl;
			if(!isIncludeVec(dont_care_po, vec))
				aff_po_vec.push_back(str_po);
		}
		//compute the probabilities of the combination of the affected POs when current vector happens (2016.8.17)
		//sim_output_wi, simu_res
		vector<int> aff_po_index, aff_po_weight;
		map<string, int> aff_po_pat_map;
		for(int i = 0; i < aff_po_vec.size(); i++)
		{
			string po = aff_po_vec[i];
			if(iIndex > 0) po.append("sim");
			itrm_sw = sim_output_wi.find(po);
			if(itrm_sw == sim_output_wi.end())
			{
				cout << "no such po " << po << " in sim_output_si!" << endl;
				exit(3);
			}
			int index = itrm_sw->second.index;
			int weight = itrm_sw->second.weight;
			aff_po_index.push_back(index);
			aff_po_weight.push_back(weight);
		}
		cout << "aff_po_weight: " << endl;
		for(int i = 0; i < aff_po_weight.size(); i++)
			cout << aff_po_weight[i] << " ";
		cout << endl;
	/*	cout << "aff_po_index: " << endl;
		for(int i = 0; i < aff_po_index.size(); i++)
			cout << aff_po_index[i] << " ";
		cout << endl;
		cout << "input_index: " << endl;
		for(int i = 0; i < input_index.size(); i++)
		{
			int index = input_index[i].index;
			int flag = input_index[i].flag;
			cout << "(" << index << ", " << flag << ")" << endl;
		}
	*/
		for(int i = 0; i < simu_res.size(); i++)
		{
			string output_pat = simu_res[i];
			string input_pat = rand[i];
			string local_input_pat;
			char c;
			for(int i = 0; i < input_index.size(); i++)
			{
				int index = input_index[i].index;
				int flag = input_index[i].flag;
				if(flag == 1) c = input_pat[index];
				else c = output_pat[index];
				local_input_pat.append(1, c);
			}
			if(local_input_pat == vec)
			{
				string aff_po_pat;
				for(int j = 0; j < aff_po_index.size(); j++)
				{
					char c = output_pat[aff_po_index[j]];
					aff_po_pat.append(1, c);
				}
			//	cout << "aff_po_pat = " << aff_po_pat << endl;
				itrm_si = aff_po_pat_map.find(aff_po_pat);
				if(itrm_si == aff_po_pat_map.end()) 
					aff_po_pat_map.insert(make_pair(aff_po_pat, 1));
				else
				{
					int &this_num = itrm_si->second;
					this_num++;
				}
			}
		}
		double max_value = 0;
		for(int i = 0; i < aff_po_weight.size(); i++)
			max_value += pow(2.0, aff_po_weight[i]);
		cout << "max_value = " << max_value << endl;
		for(itrm_si = aff_po_pat_map.begin(); itrm_si != aff_po_pat_map.end(); itrm_si++)
		{
			int num = itrm_si->second;
			double prob = (double)num/sample_num;
			string pat = itrm_si->first;
			double value = 0;
			for(int i = 0; i < pat.size(); i++)
				if(pat[i] == '1')
					value += pow(2.0, aff_po_weight[i]);
		//	cout << "pat = " << pat << ", value = " << value << endl;
			double mag;
			if(value > max_value/2)	mag = value;
			else mag = max_value - value;
		//	cout << "prob = " << prob << ", mag = " << mag << endl;
			cur_ave_error_mag += prob * mag;
		}
	}
	
	return cur_ave_error_mag;
}



double get_ave_error_mag_after_modify(vector<string> &simu_res, vector<string> &simu_res_after_modify, struct po_index_set &pis, int iIndex)
{
	double real_ave_error_mag = 0;
	
	if (iIndex == 0)
	{
		for(int i = 0; i < simu_res.size(); i++)
		{
			string org_vec = simu_res[i];
			string sim_vec = simu_res_after_modify[i];
			double sim_value = 0, org_value = 0;
			string sim_sig, org_sig;
			int k = 0;
			for(int j = pis.org_index_start; j <= pis.org_index_end; j++, k++)
			{
				char c_org = org_vec[j];
				char c_sim = sim_vec[j];
				if (c_org == '1')
					org_value += pow(2.0, k);
				if (c_sim == '1')
					sim_value += pow(2.0, k);
			/*	if (c_org != c_sim) 
					real_ave_error_mag += pow(2.0, k);	
			*/
			}
			real_ave_error_mag += abs(org_value - sim_value);
		}
	}
	else
	{
		for(int i = 0; i < simu_res_after_modify.size(); i++)
		{
			string vec = simu_res_after_modify[i];
			double sim_value = 0, org_value = 0;
			string sim_sig, org_sig;
			int k = 0;
			for(int j = pis.sim_index_start; j <= pis.sim_index_end; j++, k++)
			{
				char c = vec[j];
				sim_sig.append(1, c);
				if(c == '1') sim_value += pow(2.0, k);
			}
			k = 0;
			for(int j = pis.org_index_start; j <= pis.org_index_end; j++, k++)
			{
				char c = vec[j];
				org_sig.append(1, c);
				if(c == '1') org_value += pow(2.0, k);
			}
			real_ave_error_mag += abs(sim_value - org_value);
		}		
	}

	real_ave_error_mag /= sample_num;	
	return real_ave_error_mag;
}



/*void get_conf_interval(vector<string> &simu_res, vector<string> &simu_res_after_modify, struct po_index_set &pis, double conf_prob, double &lb_ave_error_mag, double &ub_ave_error_mag, int iIndex)
{
//	int sample_size = 500;
	int sample_times = sample_num/sample_size;
	double *sampling = new double[sample_times];
	cout << "sample_size = " << sample_size << ", sample_times = " << sample_times << endl;
	
	if (iIndex == 0)
	{		
		for(int i = 0; i < sample_times; i++)
		{
			sampling[i] = 0;
			int start = i * sample_size;
			for(int j = start; j < start + sample_size; j++)
			{
				string org_vec = simu_res[j];
				string sim_vec = simu_res_after_modify[j];
				double sim_value = 0, org_value = 0;
				string sim_sig, org_sig;
				int k = 0;
				for(int j = pis.org_index_start; j <= pis.org_index_end; j++, k++)
				{
					char c_org = org_vec[j];
					char c_sim = sim_vec[j];
					if (c_org == '1')
						org_value += pow(2.0, k);
					if (c_sim == '1')
						sim_value += pow(2.0, k);
				}
				sampling[i] += abs(org_value - sim_value);
			}
			sampling[i] /= sample_size;
		}
	}
	else
	{
		for(int i = 0; i < sample_times; i++)
		{
			sampling[i] = 0;
			int start = i * sample_size;
			for(int j = start; j < start + sample_size; j++)
			{
				string vec = simu_res_after_modify[j];
				double sim_value = 0, org_value = 0;
				string sim_sig, org_sig;
				int q = 0;
				for(int k = pis.sim_index_start; k <= pis.sim_index_end; k++, q++)
				{
					char c = vec[k];
					sim_sig.append(1, c);
					if(c == '1') sim_value += pow(2.0, q);
				}
				q = 0;
				for(int k = pis.org_index_start; k <= pis.org_index_end; k++, q++)
				{
					char c = vec[k];
					org_sig.append(1, c);
					if(c == '1') org_value += pow(2.0, q);
				}
				sampling[i] += abs(sim_value - org_value);
			}
			sampling[i] /= sample_size;
		}		
	}

	double sampling_mean = 0;
	cout << "sampling: " << endl;
	for(int i = 0; i < sample_times; i++)
	{
		cout << sampling[i] << " ";
		sampling_mean += sampling[i];
	}
	cout << endl;
	sampling_mean /= sample_times;
	
	double S_variance = 0, S_std;
	for(int i = 0; i < sample_times; i++)
		S_variance += pow((sampling[i] - sampling_mean), 2);
	S_variance /= (sample_times - 1);
	S_std = sqrt(S_variance);
	
	double sampling_std = S_std/sqrt(sample_times);
//	double sampling_std = S_std;
	cout << "sampling_mean = " << sampling_mean << ", sampling_std = " << sampling_std << endl;
	
	double diff_times_sampling_std;
	if (conf_prob - 0.95 < 1e-8)
	{
		if (sample_times == 10)
			diff_times_sampling_std = 2.262;
		else if (sample_times == 20)
			diff_times_sampling_std = 2.093;
		else if (sample_times == 5)
			diff_times_sampling_std = 2.776;
		else
			diff_times_sampling_std = 1.96;
	}
	else if (conf_prob - 0.99 < 1e-8)
	{
		if (sample_times == 10)
			diff_times_sampling_std = 3.250;
		else if(sample_times == 20)
			diff_times_sampling_std = 2.861;
		else if (sample_times == 5)
			diff_times_sampling_std = 4.604;
		else
			diff_times_sampling_std = 2.576;
	}
	
	lb_ave_error_mag = sampling_mean - diff_times_sampling_std * sampling_std;
	ub_ave_error_mag = sampling_mean + diff_times_sampling_std * sampling_std;
	
	cout << "lb_ave_error_mag = " << lb_ave_error_mag << ", ub_ave_error_mag = " << ub_ave_error_mag << endl;
	
	delete []sampling;
}
*/


void get_conf_interval(vector<string> &simu_res, vector<string> &simu_res_after_modify, struct po_index_set &pis, double conf_prob, double & sampling_mean, double &lb_ave_error_mag, double &ub_ave_error_mag, int iIndex)
{
	double *sampling = new double[sample_num];
	double *org_sampling = new double[sample_num];
	double *sim_sampling = new double[sample_num];
	double sampling_sum = 0;
	double org_sampling_sum = 0, sim_sampling_sum = 0;
	
	if (iIndex == 0)
	{		
		for(int i = 0; i < simu_res.size(); i++)
		{
			string org_vec = simu_res[i];
			string sim_vec = simu_res_after_modify[i];
			double sim_value = 0, org_value = 0;
			string sim_sig, org_sig;
			int k = 0;
			for(int j = pis.org_index_start; j <= pis.org_index_end; j++, k++)
			{
				char c_org = org_vec[j];
				char c_sim = sim_vec[j];
				org_sig.append(1, c_org);
				sim_sig.append(1, c_sim);
				if (c_org == '1')
					org_value += pow(2.0, k);
				if (c_sim == '1')
					sim_value += pow(2.0, k);
			}
			sampling[i] = abs(org_value - sim_value);
			sampling_sum += sampling[i];
			org_sampling[i] = org_value;
			org_sampling_sum += org_sampling[i];
			sim_sampling[i] = sim_value;
			sim_sampling_sum += sim_sampling[i];
			
		/*	if (sim_value != org_value)
			{
				cout << endl << "ii = " << i << endl;
				cout << "sim_sig = " << sim_sig << endl;
				cout << "org_sig = " << org_sig << endl;
				cout << "diff_index: ";
				for (int k = 0; k < sim_sig.size(); k++)
					if (sim_sig[k] != org_sig[k])
						cout << k << " ";
				cout << endl;
			}
		*/
		}
	}
	else
	{
		for(int i = 0; i < simu_res.size(); i++)
		{
			sampling[i] = 0;
			string vec = simu_res_after_modify[i];
			double sim_value = 0, org_value = 0;
			string sim_sig, org_sig;
			int q = 0;
			for(int k = pis.sim_index_start; k <= pis.sim_index_end; k++, q++)
			{
				char c = vec[k];
				sim_sig.append(1, c);
				if(c == '1') sim_value += pow(2.0, q);
			}
			q = 0;
			for(int k = pis.org_index_start; k <= pis.org_index_end; k++, q++)
			{
				char c = vec[k];
				org_sig.append(1, c);
				if(c == '1') org_value += pow(2.0, q);
			}
			sampling[i] = abs(sim_value - org_value);
			sampling_sum += sampling[i];
			org_sampling[i] = org_value;
			org_sampling_sum += org_sampling[i];
			sim_sampling[i] = sim_value;
			sim_sampling_sum += sim_sampling[i];
			
		/*	if (sim_value != org_value)
			{
				cout << endl << "ii = " << i << endl;
				cout << "sim_sig = " << sim_sig << endl;
				cout << "org_sig = " << org_sig << endl;
				cout << "diff_index: ";
				for (int k = 0; k < sim_sig.size(); k++)
					if (sim_sig[k] != org_sig[k])
						cout << k << " ";
				cout << endl;
			}
		*/
		}		
	}

	double org_sampling_mean = org_sampling_sum / sample_num;
	cout << "org_sampling_mean = " << org_sampling_mean << endl;
	double sim_sampling_mean = sim_sampling_sum / sample_num;
	cout << "sim_sampling_mean = " << sim_sampling_mean << endl;

//	double sampling_mean = sampling_sum / sample_num;
	sampling_mean = sampling_sum / sample_num;
	double S_variance = 0, S_std;
	cout << "sampling: " << endl;
	for(int i = 0; i < sample_num; i++)
	{
		if (sampling[i] > 0)
			cout << sampling[i] << " ";
		S_variance += pow((sampling[i]-sampling_mean), 2);
	}
	cout << endl;
	S_variance /= (sample_num-1);
	S_std = sqrt(S_variance);
	S_std = S_std/sqrt(sample_num);

	cout << "sampling_mean = " << sampling_mean << ", sampling_std = " << S_std << endl;
	
	double diff_times_sampling_std;
	if ((conf_prob - 0.9) < 1e-8)
		diff_times_sampling_std = 1.28;
	else if ((conf_prob - 0.95) < 1e-8)
		diff_times_sampling_std = 1.96;
	else if ((conf_prob - 0.99) < 1e-8)
		diff_times_sampling_std = 2.58;
	
	lb_ave_error_mag = sampling_mean - diff_times_sampling_std * S_std;
	ub_ave_error_mag = sampling_mean + diff_times_sampling_std * S_std;
	
	cout << "lb_ave_error_mag = " << lb_ave_error_mag << ", ub_ave_error_mag = " << ub_ave_error_mag << endl;
	
	delete []sampling;
	delete []org_sampling;
	delete []sim_sampling;
}



double estimate_ave_error_mag_bdd(BnetNetwork *net, char *cnode, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, int this_min_modified_po, int this_max_modified_po, vector<string> &final_pla, double threshold_em, int &real_aff_po_min, int &real_aff_po_max, double &real_er_whole)
{
	//iterators
	set<string>::iterator itrss;
	set<char*>::iterator itrs_char;
    
    	//variables
    	struct timeb st, et, st1, et1, st2, et2; 
    	BnetNode *nd; 
    	double cur_ave_error_mag = 0;
    	int weight;
    
    	cout << endl << "enter into estimate_ave_error_mag_bdd!" << endl;
    	cout << "cnode = " << cnode << endl;
    	//Write the comparing circuit
    	string filename = "./blif_files/compare_ckt_ave.blif";
	write_compare_circuit_ave(net, cnode, ckt_org_po, sub_abs_ckt, sub_abs_pi, sub_abs_po, final_pla, filename);
	//Read the comparing circuit into a BnetNetwork
	FILE *fp = fopen(filename.c_str(), "r");
	BnetNetwork *net_comp = Bnet_ReadNetwork(fp);
	fclose(fp);
	const char *com = "sis -t none -f ./script/sim_compare_ave.rug > sis.out";
	system(com);
	//Build bdd to compute the real average error	
	ftime(&st);         
    	DdManager *dd_comp = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    	int numPI = net_comp->npis;	
    	cur_ave_error_mag = cudd_build_v4(net_comp, &dd_comp, filename.c_str(), BNET_GLOBAL_DD, this_min_modified_po, this_max_modified_po, threshold_em, numPI, real_aff_po_min, real_aff_po_max, real_er_whole);
    	ftime(&et);
    	double rt_build_bdd = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    	cout << "rt_build_bdd = " << rt_build_bdd << endl;
    	Bnet_FreeNetwork_Bdd(net_comp, dd_comp);
    	Cudd_Quit(dd_comp);

    	return cur_ave_error_mag;
}


void get_real_er_acc(BnetNetwork *net, char *cnode, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, vector<string> &comparator_ckt, vector<string> &comparator_pi, vector<string> &comparator_po, vector<int> &comp_number, vector<string> &final_pla, double &real_er_whole)
{
	//iterators
	set<string>::iterator itrss;
	set<char*>::iterator itrs_char;
    
    	//variables
    	struct timeb st, et, st1, et1, st2, et2; 
    	BnetNode *nd; 
    	double cur_ave_error_mag = 0;
    	int weight;
    
    	cout << endl << "enter into get_real_er_acc!" << endl;
    	cout << "cnode = " << cnode << endl;
    	//Write the comparing circuit
    	string filename = "./blif_files/compare_ckt_max_bdd.blif";
    	write_compare_circuit_max(net, cnode, ckt_org_po, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number, final_pla, filename, 1);
	//Read the comparing circuit into a BnetNetwork
	FILE *fp = fopen(filename.c_str(), "r");
	BnetNetwork *net_comp = Bnet_ReadNetwork(fp);
	fclose(fp);
//	const char *com = "sis -t none -f ./script/sim_compare_max_bdd.rug > sis.out";
//	system(com);
	//Build bdd to compute the real average error	
	ftime(&st);         
    	DdManager *dd_comp = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    	int numPI = net_comp->npis;	
    	cudd_build_v5(net_comp, &dd_comp, filename.c_str(), BNET_GLOBAL_DD, numPI, real_er_whole);
    	ftime(&et);
    	double rt_build_bdd = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    	cout << "rt_build_bdd = " << rt_build_bdd << endl;
    	Bnet_FreeNetwork_Bdd(net_comp, dd_comp);
    	Cudd_Quit(dd_comp);

}

/*
double estimate_worst_error_bdd(BnetNetwork *net, char *cnode, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, int this_min_modified_po, int this_max_modified_po, vector<string> &final_pla, double threshold_em, int weight_limit, int mode)
{
    //iterators
    set<string>::iterator itrss;
    set<char*>::iterator itrs_char;
    
    //variables
    struct timeb st, et; 
    BnetNode *nd; 
    
    cout << endl << "enter into estimate_worst_error_bdd!" << endl;
    cout << "cnode = " << cnode << endl;
    //Write the comparing circuit
    string filename = "./blif_files/compare_ckt.blif";
    write_compare_circuit(net, cnode, ckt_org_po, sub_abs_ckt, sub_abs_pi, sub_abs_po, final_pla, filename);
    //Read the comparing circuit into a BnetNetwork
    FILE *fp = fopen(filename.c_str(), "r");
    BnetNetwork *net_comp = Bnet_ReadNetwork(fp);
    fclose(fp);
    const char *com = "sis -t none -f ./script/sim_compare.rug > sis.out";
    system(com);
    //Build bdd to compute the real error rate
    ftime(&st);	        
    DdManager *dd_comp = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    int numPI = net->npis; 
	int weight_limit_plus;
	int flag_power;
	if (abs(pow(2.0, weight_limit) - threshold_em) < 1e-9)
	{
		flag_power = 1;
		weight_limit_plus = weight_limit - 1;
	}
	else
	{
		flag_power = 0;
		weight_limit_plus = weight_limit;
	}
    double cur_error_mag = cudd_build_v3(net_comp, &dd_comp, filename.c_str(), BNET_GLOBAL_DD, this_min_modified_po, this_max_modified_po, weight_limit_plus, flag_power, mode);
    ftime(&et);
    double rt_build_bdd = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    cout << "rt_build_bdd = " << rt_build_bdd << endl;
    Bnet_FreeNetwork_Bdd(net_comp, dd_comp);
    Cudd_Quit(dd_comp);

    return cur_error_mag;
}
*/


int check_valid_error_mag(BnetNetwork *net, char *cnode, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, vector<string> &comparator_ckt,  vector<string> &comparator_pi, vector<string> &comparator_po, vector<int> &comp_number, vector<string> &final_pla, int BDD_SAT_mode)
{
    //iterators
    set<string>::iterator itrss;
    set<char*>::iterator itrs_char;
    
    //variables
    struct timeb st, et; 
    BnetNode *nd; 
    int result;
    
    cout << endl << "enter into check_valid_error_mag!" << endl;
    cout << "cnode = " << cnode << endl;
	
    ftime(&st);
    //Write the comparing circuit
    string filename = "./blif_files/compare_ckt_max_sat.blif";
    write_compare_circuit_max(net, cnode, ckt_org_po, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number, final_pla, filename, 0);
    ftime(&et);    
    double rt_write = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    cout << "runtime for write: " << rt_write << endl;
    
    //Read the comparing circuit into a BnetNetwork
    ftime(&st);
    const char *com = "sis -t none -f ./script/sim_compare_max_sat.rug > sis.out";
    system(com);
    ftime(&et);    
    double rt_sim = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    cout << "runtime for sim: " << rt_sim << endl;

/*    FILE *fp = fopen(filename.c_str(), "r");
    BnetNetwork *net_comp = Bnet_ReadNetwork(fp);
    fclose(fp);
    if (BDD_SAT_mode == 1) //use BDD
    {
	    //Build bdd to check
	    ftime(&st);	        
	    DdManager *dd_comp = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
	    cudd_build_v2(net_comp, &dd_comp, filename.c_str(), BNET_GLOBAL_DD);
	    ftime(&et);
	    double rt_build_bdd = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	    cout << "rt_build_bdd = " << rt_build_bdd << endl;
	    char *out_node = net_comp->outputs[0];
		st_lookup(net_comp->hash, out_node, &nd);
		int numPI = net_comp->npis;	
		double num_diff = Cudd_CountMinterm(dd_comp, nd->dd, numPI);
		Bnet_FreeNetwork_Bdd(net_comp, dd_comp);
	    	Cudd_Quit(dd_comp);
		if (num_diff > 0) return false;
		else return true;
	    
	}
	else
*/
	{
		//Use SAT to check
	    ftime(&st);	        
	    call_abc(6); //call abc to convert blif to cnf file
	    char com[100];
	    sprintf(com, "minisat_static -no-luby -rinc=1.5 -phase-saving=0 -rnd-freq=0.02 compare_ckt_max_sat.cnf sat.out");
	    system(com);
	    ftime(&et);
	    double rt_sat = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	    cout << "rt_sat = " << rt_sat << endl;

            ifstream fin;
            fin.open("sat.out");
            string str;
            getline(fin, str);
            if (str.size() == 3) return false;
            else if (str.size() == 5) return true; 
            else  
            {   
                 cout << "error in sat.out!" << endl;
                 exit(-1);
            } 
            fin.close();
	}

    return result;
}




/*
double estimate_po_sp_IS(BnetNetwork *net, char *cnode, set<char*> &po_set, map<string, struct wi_pair> &sim_output_wi, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &final_pla, int weight_bound, int iIndex)
{
	//iterators
	set<string>::iterator itrss;
	set<char*>::iterator itrs_char;
	map<string, set<char*> >::iterator itrm_sc;
	map<string, set<string> >::iterator itrm_sss;
	map<string, vector<string> >::iterator itrm_sv;
	map<string, struct wi_pair>::iterator itrm_sw;
    
    //variables
    struct timeb st, et, st1, et1, st2, et2; 
    BnetNode *nd; 
    
    for (itrs_char = po_set.begin(); itrs_char != po_set.end(); itrs_char++)
    {
    	char *cur_po = *itrs_char;
    	string cur_po_str(cur_po);
    	if(iIndex > 0) cur_po_str.append("sim");
		itrm_sw = sim_output_wi.find(cur_po_str);
		int weight = itrm_sw->second.weight;
		if (weight < weight_bound) continue;
		
		//write the circuit 
		string filename = "./blif_files/cone_po_sub.blif";
		itrm_sc = po_tfi.find(cur_po_str);
		set<char*> &tfi = itrm_sc->second;
		itrm_sss = po_inputs.find(cur_po_str);
		set<string> &inputs_set = itrm_sss->second;
		itrm_sv = po_cone_string.find(cur_po_str);
		vector<string> cone_string = itrm_sv->second;
		write_cone_po_sub(net, tfi, inputs_set, cone_string, cnode, cur_po_str, final_pla, filename);
		
		//estimate the signal probability using IS (importance sampling)
		double cur_sp = estimate_po_sp_each(net, cnode, sub_abs_ckt, sub_abs_pi, final_pla, filename, iIndex);
		cout << "cur_sp = " << cur_sp << endl;	
    }
}



double estimate_po_sp_each(BnetNetwork *net, char *cnode, string filename)
{
	//iterators
	set<string>::iterator itrss;
	set<char*>::iterator itrs_char;
    
    //variables
    struct timeb st, et, st1, et1, st2, et2; 
    BnetNode *nd; 
    double cur_ave_error_mag = 0;
    int weight;
    
    cout << "cnode = " << cnode << endl;
   
   	ftime(&st); 
   	double sp = is_simu(net, cnode, filename)
    ftime(&et);
    
	double rt_is_each = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	cout << "rt_is_each = " << rt_is_each << endl;
	
	return sp;
}

*/



//find_ave_sp()
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
