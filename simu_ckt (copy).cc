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
#include "head/loc_sim_main.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/bnet.h"

using namespace std;
extern int numPI_ini, numPO_ini;
int sample_num = 10000;


//Global variables and external variables

void simu_ckt(BnetNetwork *net, map<string, map<string, double> > &node_pattern_rate, map<string, double> &node_sp, double &real_er, int iIndex, int flag)
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
			num_output++;
			nd = nd->next;
		}

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
		if (res != 0)
		{
			cout << "error in call_abc(0)!" << endl;
			exit(1);
		}		
		//write testbench file
    	gen_testbench_org(num_input, num_output, sample_num);    	
    	//call vcs to run the simulation
    	vector<string> rand, simu_res;
    	rand.reserve(sample_num);
    	simu_res.reserve(sample_num);
		call_vcs_org(rand, simu_res);
		
		//obtain input_signatures
		cout << "obtain input_signatures: " << endl;
		map<char*, string> input_signatures;
		map<char*, string>::iterator itrm_cs;
		for(int i = 0; i < net->npis; i++)
		{
			char *input = net->inputs[i];
			st_lookup(net->hash, input, &tmp);
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
			input_signatures.insert(make_pair(tmp->name, sig));
		}
				
		//obtain node_signatures
		cout << "obtain node_signatures: " << endl;
		map<char*, string> node_signatures;				
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
			string signature;
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
	//		cout << "internal node: " << cnode << ", sp = " << sp << endl;
			string str(cnode);
			node_sp.insert(make_pair(str, sp));
			node_signatures.insert(pair<char*, string>(cnode, sig));
			nd = nd->next;
			index++;
		}
		
		//obtain node_pattern_rate from simulation result
		struct timeb st, et;		
		cout << "obtain node_pattern_rate: " << endl;		
		index = 0;
		nd = net->nodes;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			char *cnode = nd->name;
	//		cout << "cnode: " << cnode << endl;
			
			ftime(&st);
			vector<string> in_signatures;
			in_signatures.reserve(20);
			for(int j = 0; j < nd->ninp; j++)
			{
				char *innode = nd->inputs[j];
				st_lookup(net->hash, innode, &tmp);
				if(tmp->type == BNET_INPUT_NODE)
					itrm_cs = input_signatures.find(tmp->name);
				else
					itrm_cs = node_signatures.find(tmp->name);
				string signature = itrm_cs->second;
				in_signatures.push_back(signature);
			}
			ftime(&et);
    		double rt_npe1 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	 //   	cout << "@runtime for npe1: " << rt_npe1 << endl;  
			
			ftime(&st);
			map<string, double> pattern_count;
			map<string, double>::iterator itrm_sd;
			for(int j = 0; j < sample_num; j++)
			{
				string pattern;
				for(int k = 0; k < in_signatures.size(); k++)
				{
					string sig = in_signatures[k];
					char c = sig[j];
					pattern.append(1, c);
				}
				itrm_sd = pattern_count.find(pattern);
	//			cout << "this_pattern = " << pattern << ", cnode = " << cnode << endl;
				if(itrm_sd == pattern_count.end())
					pattern_count.insert(make_pair(pattern, 1));
				else
					itrm_sd->second += 1;
			}
			ftime(&et);
    		double rt_npe2 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	    	cout << "@runtime for npe2: " << rt_npe2 << endl;  
	    	
	    	ftime(&st);
			map<string, double> pattern_rate;
			for(itrm_sd = pattern_count.begin(); itrm_sd != pattern_count.end(); itrm_sd++)
			{
				double count = itrm_sd->second;
				double rate = count/(double)sample_num;
				itrm_sd->second = rate;
		//		cout << "pattern: " << itrm_sd->first << ", rate: " << rate << endl;				
			}
			string node_str(cnode);
			node_pattern_rate.insert(make_pair(node_str, pattern_count));
			ftime(&et);
    		double rt_npe3 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	    	cout << "@runtime for npe3: " << rt_npe3 << endl;  
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
	    fout << ".end" << endl;
	    fout.close();
				
		//step2.1 convert ckt_org_sim_simu.blif to ckt_org_sim_simu.v
		int res = call_abc(1); 	
		if (res != 0)
		{
			cout << "error in call_abc(1)!" << endl;
			exit(1);
		}
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
		
		//step3. obtain input_signatures
		cout << "obtain input_signatures: " << endl;
		map<char*, string> input_signatures;
		map<char*, string>::iterator itrm_cs;
		for(int i = 0; i < net->npis; i++)
		{
			char *input = net->inputs[i];
			st_lookup(net->hash, input, &tmp);
			string sig;
			for(int j = 0; j < rand.size(); j++)
			{
				string vec = rand[j];
				char c = vec[i];
				sig.append(1, c);
			}
			input_signatures.insert(make_pair(tmp->name, sig));
		}
				
		//step4. obtain node_signatures
		cout << "obtain node_signatures: " << endl;
		map<char*, string> node_signatures;				
		int index = 0;
		nd = net->nodes;
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
			char *cnode = nd->name;
			string signature;
			string sig;
			for(int j = 0; j < simu_res.size(); j++)
			{
				string vec = simu_res[j];
				char c = vec[index];
				sig.append(1, c);
			}
			node_signatures.insert(pair<char*, string>(cnode, sig));
			nd = nd->next;
			index++;
		}

		//step5. obtain node_pattern_rate from simulation result
		struct timeb st, et;		
		cout << "obtain node_pattern_rate: " << endl;		
		index = 0;
		nd = net->nodes;
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
			char *cnode = nd->name;
		//	cout << "cnode: " << cnode << " ";
			vector<string> in_signatures;
			in_signatures.reserve(20);
			for(int j = 0; j < nd->ninp; j++)
			{
				char *innode = nd->inputs[j];
				st_lookup(net->hash, innode, &tmp);
				if(tmp->type == BNET_INPUT_NODE)
					itrm_cs = input_signatures.find(tmp->name);
				else
					itrm_cs = node_signatures.find(tmp->name);
				string signature = itrm_cs->second;
				in_signatures.push_back(signature);
			}
			
			ftime(&st);
			map<string, double> pattern_count;
			map<string, double>::iterator itrm_sd;
			for(int j = 0; j < sample_num; j++)
			{
				string pattern;
				for(int k = 0; k < in_signatures.size(); k++)
				{
					string sig = in_signatures[k];
					char c = sig[j];
					pattern.append(1, c);
				}
				itrm_sd = pattern_count.find(pattern);
				if(itrm_sd == pattern_count.end())
					pattern_count.insert(make_pair(pattern, 1));
				else
					itrm_sd->second += 1;
			}
			ftime(&et);
    		double rt_npe2 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	  //  	cout << "@runtime for npe2: " << rt_npe2 << endl;  
	    	
			map<string, double> pattern_rate;
			for(itrm_sd = pattern_count.begin(); itrm_sd != pattern_count.end(); itrm_sd++)
			{
				double count = itrm_sd->second;
			//	cout << "pattern: " << itrm_sd->first << ", count: " << count << endl;				
				double rate = count/(double)sample_num;
				itrm_sd->second = rate;
			}

			string node_str;
			if(sname.find("sim") != string::npos)
				node_str = sname.substr(0, sname.length()-3);
			node_pattern_rate.insert(make_pair(node_str, pattern_count));
			nd = nd->next;
			index++;
		}
	}
}



void simu_ckt_both(BnetNetwork *net, map<string, map<string, double> > &node_pattern_rate, map<string, double> &node_sp, map<string, int> &internal_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int &num_output, double &ave_error_mag, double &real_er, int iIndex, int flag)
{
	struct timeb startTime, endTime;
	BnetNode *nd, *tmp;
	BnetTabline *f;
	map<string, struct wi_pair>::iterator itrm_sw;

	if(iIndex == 0 && flag == 0)
	{	
		cout << "in simu_ckt_both, iIndex = 0: " << endl;
		//get node_set: all nodes excluding PI nodes
		int num_input = net->npis;
		num_output = 0;
		vector<char*> node_set;
		nd = net->nodes;
		cout << "node_set: " << endl;
		struct wi_pair wi;
		int input_index = 0, output_index = 0;
		while(nd != NULL)
		{
			string str(nd->name);
			if(nd->type == BNET_INPUT_NODE)
			{
				internal_index.insert(make_pair(str, input_index));
				input_index++;
				nd = nd->next;
				continue;
			}
			if(nd->type == BNET_OUTPUT_NODE)
			{
				wi.weight = output_index++;
				wi.index = num_output;
				sim_output_wi.insert(make_pair(str, wi));
			}
			else
				internal_index.insert(make_pair(str, num_output));
			node_set.push_back(nd->name);
			num_output++;
			nd = nd->next;
		}
		
		//set pis
		pis.org_index_end = num_output - 1;
		pis.org_index_start = num_output - net->npos;

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
	    cout << "0.outputs " << endl;
	    int num = 0;
	    for(int i = 0; i < node_set.size(); i++)
	    {
	    	string str(node_set[i]);
	    	if(sim_output_wi.find(str) != sim_output_wi.end())
	    		continue;
	    	cout << node_set[i] << " ";
	    	num++;
	    }
	    cout  << endl;
	    cout << "num = " << num << endl;	    
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
	    fout << ".end" << endl;	
	    fout.close();
	
		//convert ckt_org_simu.blif to ckt_org_simu.v
		int res = call_abc(0); //ckt_org_simu.blif -> ckt_org_simu.v	
		if (res != 0)
		{
			cout << "error in call_abc(0)!" << endl;
			exit(1);
		}	
		//write testbench file
    	gen_testbench_org(num_input, num_output, sample_num);    	
    	//call vcs to run the simulation
    //	simu_res.reserve(sample_num);
		call_vcs_org(rand, simu_res);
		
		//obtain input_signatures
		cout << "obtain input_signatures: " << endl;
		map<char*, string> input_signatures;
		map<char*, string>::iterator itrm_cs;
		for(int i = 0; i < net->npis; i++)
		{
			char *input = net->inputs[i];
			st_lookup(net->hash, input, &tmp);
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
			input_signatures.insert(make_pair(tmp->name, sig));
		}
				
		//obtain node_signatures
		cout << "obtain node_signatures: " << endl;
		map<char*, string> node_signatures;				
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
			string signature;
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
	//		cout << "internal node: " << cnode << ", sp = " << sp << endl;
			string str(cnode);
			node_sp.insert(make_pair(str, sp));
			node_signatures.insert(pair<char*, string>(cnode, sig));
			nd = nd->next;
			index++;
		}
		
		//obtain node_pattern_rate from simulation result
		struct timeb st, et;		
		cout << "obtain node_pattern_rate: " << endl;		
		index = 0;
		nd = net->nodes;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			char *cnode = nd->name;
	//		cout << "cnode: " << cnode << endl;
			
			ftime(&st);
			vector<string> in_signatures;
			in_signatures.reserve(20);
			for(int j = 0; j < nd->ninp; j++)
			{
				char *innode = nd->inputs[j];
	//			cout << "innode: " << innode << " ";
				st_lookup(net->hash, innode, &tmp);
				if(tmp->type == BNET_INPUT_NODE)
					itrm_cs = input_signatures.find(tmp->name);
				else
					itrm_cs = node_signatures.find(tmp->name);
				string signature = itrm_cs->second;
			//	cout << "insig: " << signature << endl;
				in_signatures.push_back(signature);
			}
			ftime(&et);
    		double rt_npe1 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	 //   	cout << "@runtime for npe1: " << rt_npe1 << endl;  
			
			ftime(&st);
			map<string, double> pattern_count;
			map<string, double>::iterator itrm_sd;
			for(int j = 0; j < sample_num; j++)
			{
				string pattern;
				for(int k = 0; k < in_signatures.size(); k++)
				{
					string sig = in_signatures[k];
					char c = sig[j];
					pattern.append(1, c);
				}
				itrm_sd = pattern_count.find(pattern);
	//			cout << "this_pattern = " << pattern << ", cnode = " << cnode << endl;
				if(itrm_sd == pattern_count.end())
					pattern_count.insert(make_pair(pattern, 1));
				else
					itrm_sd->second += 1;
			}
			ftime(&et);
    		double rt_npe2 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	    //	cout << "@runtime for npe2: " << rt_npe2 << endl;  
	    	
	    	ftime(&st);
			map<string, double> pattern_rate;
			for(itrm_sd = pattern_count.begin(); itrm_sd != pattern_count.end(); itrm_sd++)
			{
				double count = itrm_sd->second;
				double rate = count/(double)sample_num;
				itrm_sd->second = rate;
		//		cout << "pattern: " << itrm_sd->first << ", rate: " << rate << endl;				
			}
			string node_str(cnode);
			node_pattern_rate.insert(make_pair(node_str, pattern_count));
			ftime(&et);
    		double rt_npe3 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	    //	cout << "@runtime for npe3: " << rt_npe3 << endl;  
			nd = nd->next;
			index++;
		}		
	}//if(iIndex == 0 && flag == 0)
	
	else 
	{
		//step1. get node_set: all nodes excluding PI and nodes in ckt_org
		int num_input = net->npis;
		num_output = 0;
		vector<char*> node_set;
		vector<string> org_outputs, sim_outputs;
		set<string> org_outputs_set;
		pis.org_index_start = -1;
		pis.org_index_end = -1;
		pis.sim_index_start = -1;
		pis.sim_index_end = -1;
		nd = net->nodes;
		cout << "node_set: " << endl;
		int input_index = 0, output_index = 0;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE || nd->type == BNET_OUTPUT_NODE)
			{				
				if(nd->type == BNET_INPUT_NODE)
				{
					string str(nd->name);
					internal_index.insert(make_pair(str, input_index));
					input_index++;
				}
				nd = nd->next;
				continue;
			}
			string sname(nd->name);
			if((sname.find("sim") != string::npos) && (sname.find("_out") == string::npos))
			{
				node_set.push_back(nd->name);
				internal_index.insert(make_pair(sname, num_output));
				num_output++;
				nd = nd->next;
				continue;
			}
			//put into node_set the outputs of original network and of simplified network
			if(sname.find("inter_") != string::npos)
			{
				string sin0(nd->inputs[0]);
				org_outputs.push_back(sin0);
				if(pis.org_index_start == -1) pis.org_index_start = num_output;				
				string sin1(nd->inputs[1]);
				sim_outputs.push_back(sin1);
				struct wi_pair wi;
				wi.weight = output_index++;
				wi.index = -1;
				sim_output_wi.insert(make_pair(sin1, wi));
				
			}
			nd = nd->next;
		}		
		char *outnode = net->outputs[0];
		st_lookup(net->hash, outnode, &nd);
		node_set.push_back(nd->name);
		num_output++;
		cout << "num_output = " << num_output << endl;
		
		num_output += org_outputs.size() * 2;
		pis.org_index_end = pis.org_index_start + org_outputs.size() - 1;
		pis.sim_index_start = pis.org_index_end + 1;
		pis.sim_index_end = pis.sim_index_start + sim_outputs.size() - 1;
		cout << "org_index_start = " << pis.org_index_start << ", org_index_end = " << pis.org_index_end << endl;
		cout << "sim_index_start = " << pis.sim_index_start << ", sim_index_end = " << pis.sim_index_end << endl;
		
		cout << "org_outputs: " << endl;
		for(int i = 0; i < org_outputs.size(); i++)
			cout << org_outputs[i] << " ";
		cout << endl;
		cout << endl << "sim_outputs: " << endl;
		int start_index = pis.sim_index_start;
		for(int i = 0; i < sim_outputs.size(); i++)
		{
			cout << sim_outputs[i] << " ";
			itrm_sw = sim_output_wi.find(sim_outputs[i]);
			itrm_sw->second.index = start_index++;
		//	sim_output_index.insert(make_pair(sim_outputs[i], start_index++));	
		}
		cout << endl;
				
		//step2. write ckt_org_sim_simu.blif: add the output nodes of the original network as outputs (for error magnitude metric)
		string filename = "./blif_files/ckt_org_sim_simu.blif";
	    ofstream fout;
	    fout.open(filename.c_str(), iostream::out);
	    fout << ".model ckt_org_sim_simu" << endl;
	    fout << ".inputs ";
	    for(int i = 0; i < net->npis; i++)
	        fout << "n" << net->inputs[i] << " ";
	    fout  << endl;    
	    fout << ".outputs ";
	    for(int i = 0; i < node_set.size() - 1; i++)
	    	fout << "n" << node_set[i] << " ";
	    cout << endl << "1.outputs "<< node_set.size() << endl;
	    for(int i = 0; i < node_set.size() - 1; i++)
	    	cout << node_set[i] << " ";
	    cout << endl;
	    for(int i = 0; i < org_outputs.size(); i++)
	    	fout << "n" << org_outputs[i] << " ";
	    for(int i = 0; i < sim_outputs.size(); i++)
	    	fout << "n" << sim_outputs[i] << " ";
	    fout << "n" << node_set[node_set.size()-1] << endl;	    
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
	    fout << ".end" << endl;
	    fout.close();
				
		//step2.1 convert ckt_org_sim_simu.blif to ckt_org_sim_simu.v
		int res = call_abc(1); 	
		if (res != 0)
		{
			cout << "error in call_abc(1)!" << endl;
			exit(1);
		}
		//step2.2 write testbench file
		string dut = "ckt_org_sim_simu";
    	gen_testbench_sim(num_input, num_output, sample_num, dut);    	
    	//step2.3 call vcs to run the simulation
    //	vector<string> rand;
    	string po_sig;
		call_vcs_sim(rand, simu_res, po_sig, dut);
		int num_diff = 0;
		for(int j = 0; j < po_sig.size(); j++)
			if(po_sig[j] == '1')
				num_diff++;
		cout << "num_diff = " << num_diff << endl;
		real_er = num_diff/(double)sample_num;		
		
		//step3. obtain input_signatures
		cout << "obtain input_signatures: " << endl;
		map<char*, string> input_signatures;
		map<char*, string>::iterator itrm_cs;
		for(int i = 0; i < net->npis; i++)
		{
			char *input = net->inputs[i];
			st_lookup(net->hash, input, &tmp);
			string sig;
			for(int j = 0; j < rand.size(); j++)
			{
				string vec = rand[j];
				char c = vec[i];
				sig.append(1, c);
			}
			input_signatures.insert(make_pair(tmp->name, sig));
		}
				
		//step4. obtain node_signatures
		cout << "obtain node_signatures: " << endl;
		map<char*, string> node_signatures;				
		int index = 0;
		nd = net->nodes;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE || nd->type == BNET_OUTPUT_NODE)
			{
				nd = nd->next;
				continue;
			}								
			string sname(nd->name);	
			if((sname.find("sim") != string::npos) && (sname.find("_out") == string::npos))
			{
				char *cnode = nd->name;
				string signature;
				string sig;
				for(int j = 0; j < simu_res.size(); j++)
				{
					string vec = simu_res[j];
					char c = vec[index];
					sig.append(1, c);
				}
				node_signatures.insert(pair<char*, string>(cnode, sig));
				index++;
			}
			nd = nd->next;
		}
		
		//obtain ave_error_mag		
/*	    vector<int> org_each_bit_sum, sim_each_bit_sum;
		int k = 0;
		for(int i = sim_index_start; i <= sim_index_end; i++, k++)
		{
			string sig;
			int cur_sum = 0;
			for(int j = 0; j < simu_res.size(); j++)
			{
				string vec = simu_res[j];
				char c = vec[i];
				sig.append(1, c);
				if(c == '1') cur_sum += 1;
			}
			sim_each_bit_sum.push_back(cur_sum);
		}				
		k = 0;
		for(int i = org_index_start; i <= org_index_end; i++, k++)
		{
			string sig;
			int cur_sum = 0;
			for(int j = 0; j < simu_res.size(); j++)
			{
				string vec = simu_res[j];
				char c = vec[i];
				sig.append(1, c);
				if(c == '1') cur_sum += 1;
			}
			org_each_bit_sum.push_back(cur_sum);
		}
		
		for(int i = 0; i < org_outputs.size(); i++)
		{
			int each_bit_sum_sim = sim_each_bit_sum[i];
			int each_bit_sum_org = org_each_bit_sum[i];
			ave_error_mag += abs(each_bit_sum_sim - each_bit_sum_org) * pow(2.0,i);
		}
		ave_error_mag /= sample_num;
		cout << "org_each_bit_sum: ";
		for(int i = 0; i < org_each_bit_sum.size(); i++) cout << org_each_bit_sum[i] << " ";
		cout << endl;
		cout << "sim_each_bit_sum: ";
		for(int i = 0; i < sim_each_bit_sum.size(); i++) cout << sim_each_bit_sum[i] << " ";
		cout << endl;
*/
		cout << "in simu_ckt_both: " << endl;
		for(int i = 0; i < simu_res.size(); i++)
		{
			string vec = simu_res[i];
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
			ave_error_mag += abs(sim_value - org_value);
			if (abs(sim_value - org_value) > 0)
				cout << abs(sim_value - org_value) << " ";
		}
		cout << endl;
		ave_error_mag /= sample_num;
		cout << "after simulation: ave_error_mag = " << ave_error_mag << endl;

		//step5. obtain node_pattern_rate from simulation result
		struct timeb st, et;		
		cout << "obtain node_pattern_rate: " << endl;		
		nd = net->nodes;
		while(nd != NULL)
		{
			if(nd->type == BNET_INPUT_NODE  || nd->type == BNET_OUTPUT_NODE)
			{
				nd = nd->next;
				continue;
			}
			string sname(nd->name);
			if((sname.find("sim") != string::npos) && (sname.find("_out") == string::npos))
			{				
				char *cnode = nd->name;
				vector<string> in_signatures;
				in_signatures.reserve(20);
				for(int j = 0; j < nd->ninp; j++)
				{
					char *innode = nd->inputs[j];
					st_lookup(net->hash, innode, &tmp);
					if(tmp->type == BNET_INPUT_NODE)
						itrm_cs = input_signatures.find(tmp->name);
					else
						itrm_cs = node_signatures.find(tmp->name);
					string signature = itrm_cs->second;
					in_signatures.push_back(signature);
				}
				
				ftime(&st);
				map<string, double> pattern_count;
				map<string, double>::iterator itrm_sd;
				for(int j = 0; j < sample_num; j++)
				{
					string pattern;
					for(int k = 0; k < in_signatures.size(); k++)
					{
						string sig = in_signatures[k];
						char c = sig[j];
						pattern.append(1, c);
					}
					itrm_sd = pattern_count.find(pattern);
					if(itrm_sd == pattern_count.end())
						pattern_count.insert(make_pair(pattern, 1));
					else
						itrm_sd->second += 1;
				}
				ftime(&et);
	    		double rt_npe2 = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
		  //  	cout << "@runtime for npe2: " << rt_npe2 << endl;  
		    	
				map<string, double> pattern_rate;
				for(itrm_sd = pattern_count.begin(); itrm_sd != pattern_count.end(); itrm_sd++)
				{
					double count = itrm_sd->second;
				//	cout << "pattern: " << itrm_sd->first << ", count: " << count << endl;				
					double rate = count/(double)sample_num;
					itrm_sd->second = rate;
				}

				string node_str;
				node_str = sname.substr(0, sname.length()-3);
				node_pattern_rate.insert(make_pair(node_str, pattern_count));
			}
			nd = nd->next;
		}
	
	}
    
}



void simu_ckt_both_after_modify(BnetNetwork *net, int &num_output, string cnode, struct score_pla &max_sp, int num_lines_org_pla, vector<string> &simu_res_after_modify, int iIndex)
{
	struct timeb startTime, endTime;
	BnetNode *nd, *tmp;
	BnetTabline *f;
	char *char_node = new char[cnode.size()+1];
	for(int i = 0; i < cnode.size(); i++)
		char_node[i] = cnode[i];
	char_node[cnode.size()] = '\0';
	cout << "char_node = " << char_node << endl;
	
	if (iIndex == 0)
	{
		/* step1. modify ckt_org_simu.blif */
		string snode = "n";
		snode.append(cnode);
		cout << "snode = " << snode << endl;
				
		//copy "./blif_files/ckt_org_simu.blif" to "./blif_files/ckt_org_simu_copy.blif"
		char com0[100];
		sprintf(com0, "rm -rf ./blif_files/ckt_org_simu_copy.blif");
		system(com0);
		sprintf(com0, "cp ./blif_files/ckt_org_simu.blif ./blif_files/ckt_org_simu_copy.blif");
		system(com0);
		
		//write script using sed command to delete lines for current node and the last line
		ofstream fout;
		fout.open("./script/del_1.rug", ios::out); 
		fout << "sed -i 's/ckt_org_simu/ckt_org_simu_copy/' ./blif_files/ckt_org_simu_copy.blif" << endl;
		fout << "sed -i '/";
		for (int i = 0; i < snode.size(); i++)
		{
			if (snode[i] == '[' || snode[i] == ']')
				fout << "\\";
			fout << snode[i];
		}
		fout << "$/,+" << num_lines_org_pla << "d' ./blif_files/ckt_org_simu_copy.blif" << endl;
		fout.close();
		char com[100];
		sprintf(com, "source ./script/del_1.rug");
		system(com);
		fout.open("./script/del_2.rug", ios::out);
		fout << "sed -i '$d' ./blif_files/ckt_org_simu_copy.blif" << endl;
		fout.close();
		sprintf(com, "source ./script/del_2.rug");
		system(com);
		
		//write new pla into ckt_org_simu_copy.blif
		st_lookup(net->hash, char_node, &nd);
		fout.open("./blif_files/ckt_org_simu_copy.blif", ios::app);
		vector<string> final_pla = max_sp.pla;
		if (final_pla.empty())
			fout << endl << ".names n" << cnode << endl;
		else
		{
			fout << ".names ";
			for(int i = 0; i < nd->ninp; i++)
	        	fout << "n" << nd->inputs[i] << " ";
			fout << "n" << cnode << endl;
			for(int i = 0; i < final_pla.size(); i++)
				fout << final_pla[i] << " 1" << endl;
	/*		if (final_pla.empty())
			{
				for (int i = 0; i < nd->ninp; i++) fout << "-";
				fout << " 0" << endl;
			}
	*/
		}
		fout << ".end" << endl;
		
		sprintf(com, "sis -t none -f ./script/sim_ckt_org_copy.rug > sis_ckt.txt");
    	system(com); 
		
		/* step2. convert ckt_org_simu_copy.blif to ckt_org_simu_copy.v */
		int res = call_abc(3);
		if (res != 0)
		{
			cout << "error in call_abc(3)!" << endl;
			exit(1);
		}
		
		/* step3. generate testbench file: reading rand.txt to get original input sets */
		int num_input = net->npis;
		gen_testbench_org_after_modify(num_input, num_output, sample_num);
		
		/* step4. run logic simulation on ckt_org_simu.v */
		call_vcs_org_after_modify(simu_res_after_modify);
		delete []char_node;
		return;
	}
	else
	{
		/* step1. modify ckt_org_sim_simu.blif */
		string snode = "n";
		snode.append(cnode);
		snode.append("sim");
		cout << "snode = " << snode << endl;
		
		//copy "./blif_files/ckt_org_sim_simu.blif" to "./blif_files/ckt_org_sim_simu_copy.blif"
		char com0[100];
		sprintf(com0, "rm -rf ./blif_files/ckt_org_sim_simu_copy.blif");
		system(com0);
		sprintf(com0, "cp ./blif_files/ckt_org_sim_simu.blif ./blif_files/ckt_org_sim_simu_copy.blif");
		system(com0);
		
		//write script using sed command to delete lines for current node and the last line
		ofstream fout;
		fout.open("./script/del_1.rug", ios::out);
		fout << "sed -i 's/ckt_org_sim_simu/ckt_org_sim_simu_copy/' ./blif_files/ckt_org_sim_simu_copy.blif" << endl;
		fout << "sed -i '/";
		for (int i = 0; i < snode.size(); i++)
		{
			if (snode[i] == '[' || snode[i] == ']')
				fout << "\\";
			fout << snode[i];
		}
		fout << "$/,+" << num_lines_org_pla << "d' ./blif_files/ckt_org_sim_simu_copy.blif" << endl;
		fout.close();
		char com[100];
		sprintf(com, "source ./script/del_1.rug");
		system(com); 
		fout.open("./script/del_2.rug", ios::out);
		fout << "sed -i '$d' ./blif_files/ckt_org_sim_simu_copy.blif" << endl;
		fout.close();
		sprintf(com, "source ./script/del_2.rug");
		system(com);
		
		//write new pla into ckt_org_sim_simu_copy.blif
		st_lookup(net->hash, char_node, &nd);
		fout.open("./blif_files/ckt_org_sim_simu_copy.blif", ios::app);
		vector<string> final_pla = max_sp.pla;
		if (final_pla.empty())
			fout << endl << ".names n" << cnode << "sim" << endl;
		else
		{
			fout << ".names ";
			for(int i = 0; i < nd->ninp; i++)
	        {
	        	st_lookup(net->hash, nd->inputs[i], &tmp);
	        	if(tmp->type == BNET_INPUT_NODE)
	        		fout << "n" << nd->inputs[i] << " ";
	        	else
	        		fout << "n" << nd->inputs[i] << "sim ";
	        }
			fout << "n" << cnode << "sim" << endl;			
			for(int i = 0; i < final_pla.size(); i++)
				fout << final_pla[i] << " 1" << endl;
		/*	if (final_pla.empty())
			{
				for (int i = 0; i < nd->ninp; i++) fout << "-";
				fout << " 0" << endl;
			}
		*/
		}
		fout << ".end" << endl;
		
		sprintf(com, "sis -t none -f ./script/sim_ckt_sim_copy.rug > sis_ckt.txt");
    	system(com); 
		
		/* step2. convert ckt_org_sim_simu_copy.blif to ckt_org_sim_simu_copy.v */
		int res = call_abc(4);
		if (res != 0)
		{
			cout << "error in call_abc(4)!" << endl;
			exit(1);
		}
		
		/* step3. generate testbench file: reading rand.txt to get original input sets */
		int num_input = net->npis;
		string dut = "ckt_org_sim_simu_copy";
		gen_testbench_sim_after_modify(num_input, num_output, sample_num, dut);
		
		/* step4. run logic simulation on ckt_org_sim_simu_copy.v */
		call_vcs_sim_after_modify(simu_res_after_modify, dut);
		delete []char_node;
		return;
	
	}
}


void simu_assure(BnetNetwork *net_org, BnetNetwork *net, string &snode, vector<string> &sim_pla, double &real_er)
{
	BnetNode *nd, *tmp;
	BnetTabline *tl, *f;
	FILE *fp;

    string filename = "./blif_files/ckt_assure.blif";
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_assure" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs outnode" << endl;
 
 	//write ckt_org.blif 	
 	ifstream fin;
 	fin.open("./blif_files/ckt_org.blif", ios::in);
 	string str, s;
 	int flag_start = 0;
 	while(getline(fin, str))
 	{
 		istringstream ss(str);
 		ss >> s;
 		if(s == ".end")
 			break;
 		if(s == ".names")
 			flag_start = 1;
 		if(flag_start)
 			fout << str << endl; 
 	}
 	fin.close(); 	
 	fout << endl << "# end of ckt_org.blif" << endl;
 	
 	//write ckt_sim.blif
	nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE) 
		{
			string cnode(nd->name);
			if(cnode == snode)
			{
				nd = nd->next;
				continue;
			}
			fout << ".names ";
			for(int i = 0; i < nd->ninp; i++)
			{
				st_lookup(net->hash, nd->inputs[i], &tmp);
				if(tmp->type == BNET_INPUT_NODE)
					fout << nd->inputs[i] << " ";
				else
					fout << nd->inputs[i] << "sim ";
			}
			fout << nd->name << "sim" << endl;
		    tl = nd->f;
		    while (tl != NULL) 
		    {
				if (tl->values != NULL) 
					fout << tl->values << " " << 1 - nd->polarity << endl;
				else 
					fout << 1 - nd->polarity << endl;
				tl = tl->next;
		    }
		}
		nd = nd->next;
    }
    fout << endl;
        
    if(sim_pla.empty())
	{
    	fout << ".names " << snode << "sim" << endl;
    	fout << "0" << endl;
    }
    else
    {
    	string first = sim_pla[0];
    	int nchar = first.size();
    	string dc(nchar, '-');
    	if(sim_pla.size() == 1 && first == dc)
    	{
    		fout << ".names " << snode << "sim" << endl;
    		fout << "1" << endl;
    	}
    	else
    	{
    		char *mn = new char[snode.size()+1];
    		strcpy(mn, snode.c_str());    
    		st_lookup(net->hash, mn, &nd);
    		fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
			{
				st_lookup(net->hash, nd->inputs[i], &tmp);
				if(tmp->type == BNET_INPUT_NODE)
					fout << nd->inputs[i] << " ";
				else
					fout << nd->inputs[i] << "sim ";
			}
			fout << nd->name << "sim" << endl;
			for(int i = 0; i < sim_pla.size(); i++)
				fout << sim_pla[i] << " 1" << endl;
			delete []mn;
		}
	}
 	 	
 	int i = 0;
 	fout << endl;
 	for(i = 0; i < net->npos; i++)
 	{
 		fout << ".names ";
 		char *outnode = net->outputs[i];
 		string str(outnode);
 		string str_sim(str);
 		str_sim.append("sim");
 		string str_newout("inter_n");
 		char tmp[100];
 		sprintf(tmp, "%d", i);
 		string index(tmp);
 		str_newout.append(index);
 		fout << str << " " << str_sim << " " << str_newout << endl;
 		fout << "10 1" << endl;
 		fout << "01 1" << endl;
 	}
 	fout << endl << ".names ";
 	for(i = 0; i < net->npos; i++)
 	{ 		
 		string str_newout("inter_n"); 		
 		char tmp[100];
 		sprintf(tmp, "%d", i);
 		string index(tmp);
 		str_newout.append(index);
 		fout << str_newout << " ";
 	}
 	fout << "outnode" << endl;
 	for(i = 0; i < net->npos; i++)
 		fout << "0";
 	fout << " 0" << endl;
 	fout << ".end" << endl;
 	fout.close();	
 	
 	//sweep ckt_assure.blif
/* 	char com[200];
	sprintf(com, "sis -t none -f ./script/sim_ckt_sweep_assure.rug > ./script/sis.txt"); 
	system(com);
*/ 
	//write ckt_assure_simu.blif
    fp = fopen("./blif_files/ckt_assure.blif", "r");
    BnetNetwork *net_assure = Bnet_ReadNetwork(fp);
    fclose(fp);
//    cout << "ckt_assure.blif: " << endl;
//    Bnet_PrintNetwork(net_assure);
    filename = "./blif_files/ckt_assure_simu.blif";
	fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_assure_simu" << endl;
	fout << ".inputs ";
	for(int i = 0; i < net_assure->npis; i++)
		fout << "n" << net_assure->inputs[i] << " ";
	fout  << endl;    
	fout << ".outputs ";
	for(int i = 0; i < net_assure->npos; i++)
		fout << "n" << net_assure->outputs[i] << " ";
	fout  << endl;	    
	nd = net_assure->nodes;
	while(nd != NULL)
	{
		if(nd->type != BNET_INPUT_NODE)
	    {
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
	    }	
	    nd = nd->next;
	}
	fout << ".end" << endl;
	fout.close();
 			
	//step2.1 convert ckt_assure_simu.blif to ckt_assure_simu.
	int res = call_abc(2); 
	if(res)
	{
		cout << "call_abc(2) fails!" << endl;
		exit(1);
	}	
	//step2.2 write testbench file
	int num_input = net->npis;
	int num_output = 1;
	string dut = "ckt_assure_simu";
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
		
}
