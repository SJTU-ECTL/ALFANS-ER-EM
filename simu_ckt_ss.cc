#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
//#include <unordered_map>
#include <cmath>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include "head/queue.h"
#include "head/call_abc.h"
#include "head/write_func.h"
#include "head/helper.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/bnet.h"

using namespace std;
extern int numPI_ini, numPO_ini;
extern int sample_num;


//Global variables and external variables

void simu_ckt_ss(BnetNetwork *net, map<string, string> &node_signatures, map<string, double> &node_sp, double &real_er, int iIndex, int flag)
{
	struct timeb startTime, endTime;
	BnetNode *nd, *tmp;
	BnetTabline *f;

	if(iIndex == 0 && flag == 0)
	{	
		//get node_set: all nodes excluding PI nodes
		int num_input = net->npis;
		int num_output = 0;
		vector<char*> node_set;
		nd = net->nodes;
		cout << "node_set: " << endl;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			node_set.push_back(nd->name);
			cout << nd->name << " ";
			num_output++;
			nd = nd->next;
		}
		cout << endl;

		//write ckt_org_simu.blif
		string filename = "./blif_files/ckt_org_simu.blif";
	    ofstream fout;
	    fout.open(filename.c_str(), iostream::out);
	    fout << ".model ckt_org_simu" << endl;
	    fout << ".inputs ";
	    for(int i = 0; i < net->npis; i++)
	        fout << "n" << net->inputs[i] << " ";
	    fout  << endl;    
	    fout << ".outputs ";
	    for(int i = 0; i < node_set.size(); i++)
	    	fout << "n" << node_set[i] << " ";
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
	            fout << "n" << nd->inputs[i] << " ";
	        fout << "n" << nd->name << endl;
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
	
		//convert ckt_org_simu.blif to ckt_org_simu.v
		int res = call_abc(0); //ckt_org_simu.blif -> ckt_org_simu.v		
		//write testbench file
    	gen_testbench_org(num_input, num_output, sample_num);    	
    	//call vcs to run the simulation
    	vector<string> rand, simu_res;
    	rand.reserve(sample_num);
    	simu_res.reserve(sample_num);
		call_vcs_org(rand, simu_res);
		
		//obtain input_signatures
		cout << "obtain node_signatures of PI nodes: " << endl;
		for(int i = 0; i < net->npis; i++)
		{
			char *input = net->inputs[i];
			string sig;
			double num_one = 0;
			for(int j = 0; j < rand.size(); j++)
			{
				string vec = rand[j];
				char c = vec[i];
				sig.append(1, c);
				if(c == '1')
					num_one++;
			}
			double sp = num_one/sample_num;
			cout << "input node: " << input << ", sp = " << sp << endl;
			string str(input);
			node_sp.insert(make_pair(str, sp));
			node_signatures.insert(make_pair(str, sig));
		}
				
		//obtain node_signatures
		cout << "obtain node_signatures of internal nodes: " << endl;			
		int index = 0;
		nd = net->nodes;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			char *cnode = nd->name;
			string str(cnode);
			string sig;
			double num_one = 0;
			for(int j = 0; j < simu_res.size(); j++)
			{
				string vec = simu_res[j];
				char c = vec[index];
				sig.append(1, c);
				if(c == '1')
					num_one++;
			}
			double sp = num_one/sample_num;			
			node_sp.insert(make_pair(str, sp));
			node_signatures.insert(pair<string, string>(str, sig));
			nd = nd->next;
			index++;
		}
	}//if(iIndex == 0 && flag == 0)
	
	else 
	{
		//step1. get node_set: all nodes excluding PI and nodes in ckt_org
		int num_input = net->npis;
		int num_output = 0;
		vector<char*> node_set;
		nd = net->nodes;
		cout << "node_set: " << endl;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			string sname(nd->name);
			if((sname.find("sim") == string::npos) || nd->type == BNET_OUTPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			node_set.push_back(nd->name);
			num_output++;
			nd = nd->next;
		}
		char *outnode = net->outputs[0];
		st_lookup(net->hash, outnode, &nd);
		node_set.push_back(nd->name);
		num_output++;
		
		//write ckt_org_sim_simu.blif
		string filename = "./blif_files/ckt_org_sim_simu.blif";
	    ofstream fout;
	    fout.open(filename.c_str(), iostream::out);
	    fout << ".model ckt_org_sim_simu" << endl;
	    fout << ".inputs ";
	    for(int i = 0; i < net->npis; i++)
	        fout << "n" << net->inputs[i] << " ";
	    fout  << endl;    
	    fout << ".outputs ";
	    for(int i = 0; i < node_set.size(); i++)
	    	fout << "n" << node_set[i] << " ";
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
	            fout << "n" << nd->inputs[i] << " ";
	        fout << "n" << nd->name << endl;
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
		
		
		//step2.1 convert ckt_org_sim_simu.blif to ckt_org_sim_simu.v
		int res = call_abc(1); 	
		//step2.2 write testbench file
		string dut = "ckt_org_sim_simu";
    	gen_testbench_sim(num_input, num_output, sample_num, dut);    	
    	//step2.3 call vcs to run the simulation
    	vector<string> rand, simu_res;
    	string po_sig;
		call_vcs_sim(rand, simu_res, po_sig, dut);
		int num_diff = 0;
		for(int j = 0; j < po_sig.size(); j++)
			if(po_sig[j] == '1')
				num_diff++;
		cout << "num_diff = " << num_diff << endl;
		real_er = num_diff/(double)sample_num;		
		
		//step3. obtain node_signatures
		cout << "obtain node_signatures of PI nodes: " << endl;
		for(int i = 0; i < net->npis; i++)
		{
			char *input = net->inputs[i];
			string snode(input);
			string sig;
			double num_one = 0;
			for(int j = 0; j < rand.size(); j++)
			{
				string vec = rand[j];
				char c = vec[i];
				sig.append(1, c);
				if(c == '1')
					num_one++;
			}
			double sp = num_one/sample_num;
	//		cout << "input node: " << input << ", sp = " << sp << endl;
			node_sp.insert(make_pair(snode, sp));
			node_signatures.insert(make_pair(snode, sig));
		}
				
		//step4. obtain node_signatures of internal nodes
		cout << "obtain node_signatures of internal nodes: " << endl;			
		int index = 0;
		nd = net->nodes;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			char *cnode = nd->name;
			string snode(cnode);
			if((snode.find("sim") == string::npos) || nd->type == BNET_OUTPUT_NODE)
			{
				nd = nd->next;
				continue;
			}			
			string sig;
			double num_one = 0;
			for(int j = 0; j < simu_res.size(); j++)
			{
				string vec = simu_res[j];
				char c = vec[index];
				sig.append(1, c);
				if(c == '1')
					num_one++;
			}
			double sp = num_one/sample_num;
	//		cout << "internal node: " << snode << ", sp = " << sp << endl;
			
			string snode_org;
			if(snode.find("sim") != string::npos)
				snode_org = snode.substr(0, snode.length()-3);
			node_sp.insert(make_pair(snode_org, sp));
			node_signatures.insert(make_pair(snode_org, sig));
			nd = nd->next;
			index++;
		}

	}
    
}

