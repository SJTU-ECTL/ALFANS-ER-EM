#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <deque>
#include <list>
#include <cmath>
#include <cassert>
#include "head/helper.h"
#include "cudd/bnet.h"

using namespace std;

//external variables
extern int sample_num, numPI_ini, numPO_ini;


/*
functions in this files:
1. void Graph::gen_ckt_tb()
2. void Graph::gen_ckt_tb_bn(int input_size, int output_size)
3. void Graph::gen_ckt_tb_saif(int iIndex)
4. void Graph::gen_ckt_tb_simu(vector<string> &error_invec, vector<string> &all_nodes)
5. void Graph::write_ckt_bignode(Node *cnode, map<int, int> &in_sig, vector<string> &exdc_set)  with exdcs
6. void Graph::write_ckt_bignode_mv(Node *cnode, map<int, int> &in_sig, vector<string> &exdc_set)
7. void Graph::write_ckt_org()
8. void Graph::write_ckt_sim(Node *node, vector<string> &sim_node_nodes, map<int, int> &const_PO)
9. void Graph::write_ckt_bignode(Node *cnode, map<int, int> &in_sig, vector<int> &sort_insig)
10. void Graph::add_exdc_pla(vector<string> &exdc_set)
11. void Graph::cluster_ckt(int max_save_node, map<int, int> &in_sig, int iIndex)
*/


/*7. write_ckt_org()*/
void write_ckt_org(BnetNetwork *net)
{
	BnetNode *nd;
	BnetTabline *f;
	set<char*>::iterator itrs;
	
//	Bnet_PrintNetwork(net);

    /*write circuit for the original circuit*/
    string filename = "./blif_files/ckt_org.blif";
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_org" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    for(int i = 0; i < net->npos; i++)
        fout << net->outputs[i] <<  "_out ";
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
            fout << nd->inputs[i] << " ";
        fout << nd->name << endl;
        f = nd->f;
        while(f != NULL)
        {
        	if (f->values != NULL) 
        	{
		    	fout << f->values << " " << 1 - nd->polarity << endl;
			} 
			else 
			{
			    fout <<  1 - nd->polarity << endl;
			}
        	f = f->next;
        }
        nd = nd->next;
    }
    
    fout  << endl;
    for(int i = 0; i < net->npos; i++)
    {
    	fout << ".names " << net->outputs[i] <<  " " << net->outputs[i] << "_out" << endl;
    	fout << "1 1" << endl;
    }
    
    
    fout << ".end" << endl;
    fout.close();
   
}


void write_mvsis_rug_po(char *po, char *node)
{  
    ofstream fout;
    char filename[50] = "./script/mvsis_dc_po.rug";
    fout.open(filename, ios::out);
    
    string fn = "./blif_files/cone_po_";
    string po_str(po);
    fn.append(po_str);
    fn.append(".blif");
    fout << "read_blif " << fn << endl;
    fout << "mfs -s -w 33 -k -N " << node << endl;
//    cout << "read_blif " << fn << endl;
//    cout << "mfs -s -w 22 -k -N " << node << endl;
    fout.close();
}

void write_mvsis_rug(int iIndex, char *node)
{
	ofstream fout;
	char filename[100] = "./script/mvsis_dc.rug";
	fout.open(filename, ios::out);
	if(iIndex == 0)
		fout << "read_blif ./blif_files/ckt_org.blif" << endl;
	else
		fout << "read_blif ./blif_files/ckt_sim.blif" << endl;
	fout << "mfs -s -w 33 -k -N " << node << endl;
//	fout << "mfs -s -w 55 -k -N " << node << endl;
	fout.close();
}

void write_mvsis_rug_MFFC(char *node)
{
	ofstream fout;
	char filename[100] = "./script/mvsis_dc.rug";
	fout.open(filename, ios::out);
	fout << "read_blif ./blif_files/ckt_cluster.blif" << endl;
	fout << "mfs -s -w 22 -k -N " << node << endl;
//	fout << "mfs -s -w 55 -k -N " << node << endl;
	fout.close();
}


void write_abc_rug(int iIndex, char *node)
{
	ofstream fout;
	char filename[100] = "./script/abc_dc.rug";
	fout.open(filename, ios::out);
	if(iIndex == 0)
		fout << "read_blif ./blif_files/ckt_org.blif" << endl;
	else
		fout << "read_blif ./blif_files/ckt_sim.blif" << endl;
	fout << "mfs -N " << node << " -w -g" << endl;
	fout.close();
}


void write_abc_rug_po(char *po, char *node)
{  
    ofstream fout;
    char filename[50] = "./script/abc_dc_po.rug";
    fout.open(filename, ios::out);
    
    string fn = "./blif_files/cone_po_";
    string po_str(po);
    fn.append(po_str);
    fn.append(".blif");
    fout << "read_blif " << fn << endl;
    fout << "mfs -N " << node << " -w" << endl;
    fout.close();
}



/*5. write_ckt_bignode() without exdc*/
void write_ckt_bignode(BnetNetwork *net, char *cnode)
{
//    cout << "Coming into write_ckt_bignode!" << endl;
	multimap<double, string>::iterator itrmm_ds;


	BnetNode *nd;
	st_lookup(net->hash, cnode, &nd);

    string dir = "./blif_files/";
    string filename = "bigNode.blif";
    dir.append(filename);
    filename = dir;
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model bigNode" << endl;
    fout << ".inputs ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
	fout << endl;
    fout << ".outputs " << cnode << endl;
    
    fout << ".names ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
    fout << cnode << endl;

    BnetTabline *t = nd->f;
	while(t != NULL) 
	{
		if(t->values != NULL) 
			fout << t->values << " " << 1 - nd->polarity << endl;
		else
			 fout << 1 - nd->polarity << endl;
	
		t = t->next;
	}
   
    
    fout << ".end" << endl;
    fout.close();
}




/*5. write_ckt_bignode() with exdc*/
void write_ckt_bignode(BnetNetwork *net, char *cnode, vector<char*> &cutnodes, multimap<double, string> &exdc_set)
{
//    cout << "Coming into write_ckt_bignode!" << endl;
	multimap<double, string>::iterator itrmm_ds;


	BnetNode *nd;
	st_lookup(net->hash, cnode, &nd);

    string dir = "./blif_files/";
    string filename = "bigNode.blif";
    dir.append(filename);
    filename = dir;
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model bigNode" << endl;
    fout << ".inputs ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
	fout << endl;
    fout << ".outputs " << cnode << endl;
    
    fout << ".names ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
    fout << cnode << endl;

    BnetTabline *t = nd->f;
	while(t != NULL) 
	{
		if(t->values != NULL) 
			fout << t->values << " " << 1 - nd->polarity << endl;
		else
			 fout << 1 - nd->polarity << endl;
	
		t = t->next;
	}
    
    fout << ".exdc " << endl;
    fout << ".inputs ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
	fout << endl;
    fout << ".outputs " << cnode << endl;
    fout << ".names ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
    fout << cnode << endl;

    for(itrmm_ds = exdc_set.begin(); itrmm_ds != exdc_set.end(); itrmm_ds++)
        fout << itrmm_ds->second << " 1" << endl;
    
    fout << ".end" << endl;
    fout.close();
}


/*5. write_ckt_bignode() without exdc*/
void write_ckt_bignode_MFFC(BnetNetwork *net, char *cnode, set<char*> &cMFFC, map<char*, char*> &in_sig)
{
	set<char*>::iterator itrs;
	map<char*, char*>::iterator itrm_cc;

	BnetNode *nd, *tmp;
	st_lookup(net->hash, cnode, &nd);

    string filename = "./blif_files/bigNode.blif";
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model bigNode" << endl;
    fout << ".inputs ";
    for(itrm_cc = in_sig.begin(); itrm_cc != in_sig.end(); itrm_cc++)
        fout << itrm_cc->first << " ";
	fout << endl;
    fout << ".outputs " << cnode << endl;
    
    for(itrs = cMFFC.begin(); itrs != cMFFC.end(); itrs++)
    {
    	char *node = *itrs;
        st_lookup(net->hash, node, &tmp);
        if(tmp->type == BNET_INPUT_NODE)
            continue;
        fout << ".names "; 
        for(int i = 0 ; i < tmp->ninp; i++)
        	fout << tmp->inputs[i] << " ";
        fout << node << endl;       
 		BnetTabline *t = tmp->f;
		while(t != NULL) 
		{
			if(t->values != NULL) 
				fout << t->values << " " << 1 - nd->polarity << endl;
			else
				 fout << 1 - nd->polarity << endl;
		
			t = t->next;
		}               
    }    
    
    fout << ".end" << endl;
    fout.close();
}



void write_bignode_pla(BnetNetwork *net, char *cnode)
{
	multimap<double, string>::iterator itrmm_ds;

	BnetNode *nd;
	st_lookup(net->hash, cnode, &nd);

    string filename = "./pla_files/bigNode.pla";
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".i " << nd->ninp << endl;;
    fout << ".o 1" << endl;
    fout << ".ilb ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
	fout << endl;
	fout << ".ob " << cnode << endl;
    
    BnetTabline *t = nd->f;
	while(t != NULL) 
	{
		if(t->values != NULL) 
			fout << t->values << " " << 1 - nd->polarity << endl;
		else
			 fout << 1 - nd->polarity << endl;
	
		t = t->next;
	}
        
    fout << ".end" << endl;
    fout.close();
}



void write_bignode_pla_sim(BnetNetwork *net, char *cnode, double sp, vector<string> &sim_org_pla)
{
	multimap<double, string>::iterator itrmm_ds;

	BnetNode *nd;
	st_lookup(net->hash, cnode, &nd);

    //write ./pla_files/bigNode_sim.pla
    string dir = "./pla_files/";
    string filename = "bigNode_sim.pla";
    dir.append(filename);
    filename = dir;
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".i " << nd->ninp << endl;;
    fout << ".o 1" << endl;
    fout << ".ilb ";
    for(int i = 0; i < nd->ninp; i++)
        fout << nd->inputs[i] << " ";
	fout << endl;
	fout << ".ob " << cnode << endl;
	
	if(!sim_org_pla.empty())
    	for(int i = 0; i < sim_org_pla.size(); i++)
    		fout << sim_org_pla[i] << " 1" << endl;
    else
    {
    	for(int i = 0; i < nd->ninp; i++)
    			fout << "-";
    /*	if(nd->rp < 0.5)	
    		fout << " 0" << endl;
    	else
    		fout << " 1" << endl;
    */
    	if(sp < 0.5)
   			fout << " 0" << endl;
    	else
    		fout << " 1" << endl; 	
    		
    }
    
    fout << ".end" << endl;
    fout.close();
    
}





/*8. write_ckt_sim()*/
void write_ckt_sim(BnetNetwork *net, char *cnode, vector<string> &final_pla)
{
	set<char*>::iterator itrs;
	BnetNode *nd;
	BnetTabline *tl;
	

    /*write circuit for the simplified circuit*/
    string path = "./blif_files/";
    string filename = "ckt_sim.blif";
    path.append(filename);
    filename = path;
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_sim" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    for(int i = 0; i < net->npos; i++)
        fout << net->outputs[i] <<  " ";
    fout  << endl;
    
    
    nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, cnode)) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
				fout << nd->inputs[i] << " ";
			fout << nd->name << endl;
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

	cout << "simplified bignode: " << endl;
	st_lookup(net->hash, cnode, &nd);	
	if(final_pla.empty())
	{
		fout << endl << ".names " << cnode << endl;
		cout << endl << ".names " << cnode << endl;
	}
	else
	{
		fout << endl << ".names ";
		for(int i = 0; i < nd->ninp; i++)
			fout << nd->inputs[i] << " ";
		fout << nd->name << endl;
		for(int i = 0; i < final_pla.size(); i++)
		    fout << final_pla[i] << " 1" << endl;
		    
		for(int i = 0; i < nd->ninp; i++)
			cout << nd->inputs[i] << " ";
		cout << nd->name << endl;
		for(int i = 0; i < final_pla.size(); i++)
		    cout << final_pla[i] << " 1" << endl;
	}
	cout << endl;
    fout << ".end" << endl;
    fout.close();
    
}


/*8. write_ckt_sim_const()*/
void write_ckt_sim_const(BnetNetwork *net, char *cnode, int pol)
{
	set<char*>::iterator itrs;
	BnetNode *nd;
	BnetTabline *tl;
	

    /*write circuit for the simplified circuit*/
    string path = "./blif_files/";
    string filename = "ckt_sim_tmp.blif";
    path.append(filename);
    filename = path;
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_sim_tmp" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    for(int i = 0; i < net->npos; i++)
        fout << net->outputs[i] <<  " ";
    fout  << endl;
    
    
    nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, cnode)) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
				fout << nd->inputs[i] << " ";
			fout << nd->name << endl;
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

	fout << ".names " << cnode << endl;
	fout << pol << endl;
	
    fout << ".end" << endl;
    fout.close();
    
}




/*8. write_ckt_sim()*/
void write_ckt_sim_file(BnetNetwork *net, char *cnode, vector<string> &sim_node_nodes, string &filename)
{
	set<char*>::iterator itrs;
	BnetNode *nd;
	BnetTabline *tl;
	

    /*write circuit for the simplified circuit*/
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_sim_after" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    for(int i = 0; i < net->npos; i++)
        fout << net->outputs[i] <<  " ";
    fout  << endl;
    
    
    nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, cnode)) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
				fout << nd->inputs[i] << " ";
			fout << nd->name << endl;
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

	for(int i = 0; i < sim_node_nodes.size(); i++)
	    fout << sim_node_nodes[i] << endl;
    fout << ".end" << endl;
    fout.close();
    
}






void write_ckt_comb(BnetNetwork *net)
{	
	set<char*>::iterator itrs;
	BnetNode *nd, *tmp;
	BnetTabline *tl;

	/*write circuit for the simplified circuit*/
    string filename = "./blif_files/ckt_org_sim.blif";
    ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model ckt_org_sim" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs outnode" << endl;
 
 	//write ckt_org.blif
/* 	string filename_org = "./blif_files/ckt_org.blif";
    FILE *fp = fopen(filename_org.c_str(), "r");
    BnetNetwork *net_org = Bnet_ReadNetwork(fp); 
    fclose(fp);
    nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
				fout << "n" << nd->inputs[i] << " ";
			fout << "n" << nd->name << endl;
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
*/
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
 	
 	//write ckt_sim.blif
/* 	nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
			{
				st_lookup(net->hash, nd->inputs[i], &tmp);
				if(tmp->type == BNET_INPUT_NODE)
					fout << "n" << nd->inputs[i] << " ";
				else
					fout << "n" << nd->inputs[i] << "sim ";
			}
			fout << "n" << nd->name << "sim" << endl;
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
*/ 	
 	
 	fin.open("./blif_files/ckt_sim.blif", ios::in);
 	flag_start = 0;
 	int flag_continue = 0;
 	while(getline(fin, str))
 	{
 		istringstream ss(str);
 		ss >> s;
 		if(s == ".end")
 			break;
 		if(s == ".names")
 			flag_start = 1;
 		if(flag_start)
 		{
 			if(s == ".names" || flag_continue == 1)
 			{
 				string str_new;
 				str_new.append(s);
 				str_new.append(" ");
 				string s_sim;
 				while(ss >> s)
 				{
 					char *name_ch = const_cast<char*>(s.c_str());
 					st_lookup(net->hash, name_ch, &nd);
 					if(!st_lookup(net->hash, name_ch, &nd))
					{
 						s_sim = s;
 						s_sim.append("sim");
 						str_new.append(s_sim);
 					} 	
 					else if(nd->type == BNET_INPUT_NODE) 						
 						str_new.append(s);
 					else if(s == "\\")
 					{ 						
 						str_new.append(s);
 						flag_continue = 1;
 					}
 					else
 					{
 						s_sim = s;
 						s_sim.append("sim");
 						str_new.append(s_sim);
 					} 
 					str_new.append(" ");						
 				}
 				if(s != "\\")
 					flag_continue = 0;
 				str = str_new;
 			}
 			fout << str << endl; 
 		}
 	}
 	fin.close();
 	
 	int i = 0;
 	fout << endl;
 	for(i = 0; i < net->npos; i++)
 	{
 		fout << ".names ";
 		char *outnode = net->outputs[i];
/* 		string str = "n";
 		string ou(outnode);
 		str.append(ou);
*/
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
}



void gen_testbench_org(int num_input, int num_output, int sample_num)
{
	srand((unsigned)time(NULL));
	unsigned int rand_number1 = rand();
	cout << "rand_number1 = " << rand_number1 << endl;

    ofstream fout;
    string filename = "./verilog_files/ckt_tb_org.v";
    fout.open(filename.c_str(), ofstream::out);
//    fout << "`timescale 1ps/1ps" << endl;	
    fout << "module ckt_tb_org();" << endl;
    fout << "  " << "parameter M = " << num_input << ", N = " << num_output << ";" << endl;
    fout << "  " << "parameter sample_num = " << sample_num << ";" << endl;
    fout << "  " << "parameter H1 = M/2 , H2 = M - H1;" << endl;
    fout << "  " << "reg [M-1:0] rand;" << endl;
    fout << "  " << "reg bit;" << endl;
    fout << "  " << "reg [H1-1:0] rand1;" << endl;
    fout << "  " << "reg [H2-1:0] rand2;" << endl;
    fout << "  " << "wire [N-1:0] y;" << endl;
    fout << "  " << "integer i, j, seed, fp1, fp2, fp3;" << endl;
    fout << "  " << "integer seed1 = " << rand_number1 << ";" << endl;
    fout << "  " << "ckt_org_simu dut(";
    for(int i = num_input - 1; i >= 0; i--)
        fout << "rand[" << i << "], ";
    for(int j = num_output - 1; j >= 0; j--)
    {
    	if(j != 0)
        	fout << "y[" << j << "], ";	
        else
        	fout << "y[" << j << "] ";	
    }
    fout << ");" << endl;
    
    fout << "  " << "initial begin " << endl;
    fout << "  	" << "fp1 = $fopen(\"rand.txt\", \"w\");" << endl;
    fout << "  	" << "fp2 = $fopen(\"simu_res.txt\", \"w\");" << endl;

    fout << "  	" << "for(i=0; i < sample_num; i=i+1) begin" << endl;    
	fout << "  		" << "for(j=0; j < M; j=j+1) begin" << endl;
	fout << "     		" << "bit = $random(seed1)%2;" << endl;
	fout << "     		" << "rand = {{rand}, bit};" << endl;
	fout << "   	" << "end" << endl;
    fout << "       " << "#1 $fwrite(fp1, \"%b\\n\", rand);" << endl;
    fout << "       " << "$fwrite(fp2, \"%b\\n\", y);" << endl;
    fout << "   " << "end" << endl;
    fout << "  " << "$fclose(fp1);" <<endl;
    fout << "  " << "$fclose(fp2);" <<endl;
    fout << "  " << "end" << endl;
    fout << "endmodule" << endl;
    fout.close();
}


void gen_testbench_org_after_modify(int num_input, int num_output, int sample_num, int status)
{
	srand((unsigned)time(NULL));
	unsigned int rand_number1 = rand();
	cout << "rand_number1 = " << rand_number1 << endl;

    ofstream fout;
    string filename = "./verilog_files/ckt_tb_org_after_modify.v";
    fout.open(filename.c_str(), ofstream::out);
    fout << "module ckt_tb_org_after_modify();" << endl;
    fout << "  " << "parameter M = " << num_input << ", N = " << num_output << ";" << endl;
    fout << "  " << "parameter sample_num = " << sample_num << ";" << endl;
    fout << "  " << "parameter H1 = M/2 , H2 = M - H1;" << endl;
    fout << "  " << "reg [M-1:0] rand;" << endl;
    fout << "  " << "reg bit;" << endl;
    fout << "  " << "wire [N-1:0] y;" << endl;
    fout << "  " << "integer i, j, seed, cnt, fp_r, fp_w, fp1, fp2;" << endl;
    fout << "  " << "integer seed1 = " << rand_number1 << ";" << endl;
    fout << "  " << "ckt_org_simu_copy dut(";
    for(int i = num_input - 1; i >= 0; i--)
        fout << "rand[" << i << "], ";
    for(int j = num_output - 1; j >= 0; j--)
    {
    	if(j != 0)
        	fout << "y[" << j << "], ";	
        else
        	fout << "y[" << j << "] ";	
    }
    fout << ");" << endl;

	if (!status)    
	{
	    fout << "  " << "initial begin " << endl;
	    fout << "  	" << "fp_w = $fopen(\"simu_res.txt\", \"w\");" << endl;
	    fout << "  	" << "fp_r = $fopen(\"rand.txt\", \"r\");" << endl;
	    
	//    fout << "   " << "while (!$feof(fp_r)) begin" << endl;
	    fout << "  	" << "for(i=0; i < sample_num; i=i+1) begin" << endl;   
	    fout << "     " << "cnt = $fscanf(fp_r, \"%b\", rand);" << endl;
	    fout << "     " << "#1 $fwrite(fp_w, \"%b\\n\", y);" << endl;
	    fout << "   " << "end" << endl;
	    fout << "  " << "$fclose(fp_r);" <<endl;
	    fout << "  " << "$fclose(fp_w);" <<endl;
	    fout << "  " << "end" << endl;
	    fout << "endmodule" << endl;
	}
	else
	{
		fout << "  " << "initial begin " << endl;
	    fout << "  	" << "fp1 = $fopen(\"rand.txt\", \"w\");" << endl;
	    fout << "  	" << "fp2 = $fopen(\"simu_res.txt\", \"w\");" << endl;

	    fout << "  	" << "for(i=0; i < sample_num; i=i+1) begin" << endl;    
		fout << "  		" << "for(j=0; j < M; j=j+1) begin" << endl;
		fout << "     		" << "bit = $random(seed1)%2;" << endl;
		fout << "     		" << "rand = {{rand}, bit};" << endl;
		fout << "   	" << "end" << endl;
	    fout << "       " << "#1 $fwrite(fp1, \"%b\\n\", rand);" << endl;
	    fout << "       " << "$fwrite(fp2, \"%b\\n\", y);" << endl;
	    fout << "   " << "end" << endl;
	    fout << "  " << "$fclose(fp1);" <<endl;
	    fout << "  " << "$fclose(fp2);" <<endl;
	    fout << "  " << "end" << endl;
	    fout << "endmodule" << endl;
	}
    fout.close();
}


void gen_testbench_sim(int num_input, int num_output, int sample_num, string &dut)
{
    ofstream fout;
    string filename = "./verilog_files/ckt_tb_sim.v";
    fout.open(filename.c_str(), ofstream::out);
    
    srand((unsigned)time(0));
	unsigned int rand_number1 = rand();
	cout << "rand_number1 = " << rand_number1 << endl;

    fout << "module ckt_tb_sim();" << endl;
    fout << "  " << "parameter M = " << num_input << ", N = " << num_output << ";" << endl;
    fout << "  " << "parameter sample_num = " << sample_num << ";" << endl;
    fout << "  " << "reg [M-1:0] rand;" << endl;
    fout << "  " << "reg bit;" << endl;
    fout << "  " << "wire [N-1:0] y;" << endl;
    fout << "  " << "integer i, j, fp1, fp2, fp3;" << endl;
    fout << "  " << "integer seed1 = " << rand_number1 << ";" << endl;
//    fout << "  " << "ckt_org_sim_simu dut(";
	fout << "  " << dut << " dut(";
    for(int i = num_input - 1; i >= 0; i--)
        fout << "rand[" << i << "], ";
    for(int j = num_output - 1; j >= 0; j--)
    {
    	if(j != 0)
        	fout << "y[" << j << "], ";	
        else
        	fout << "y[" << j << "] ";	
    }
    fout << ");" << endl;
    
    fout << "  " << "initial begin " << endl;
    fout << "  " << "fp1 = $fopen(\"rand.txt\", \"w\");" << endl;
    fout << "  " << "fp2 = $fopen(\"simu_res.txt\", \"w\");" << endl;
    fout << "  " << "fp3 = $fopen(\"simu_po.txt\", \"w\");" << endl;
    
    fout << "  " << "for(i=0; i < sample_num; i=i+1) begin" << endl;    
//	fout << "     " << "rand1 = {$random(seed1)};" << endl;
//    fout << "     " << "rand2 = {$random(seed2)};" << endl;
//    fout << "     " << "rand = {{rand1}, {rand2}};" << endl;
	fout << "  	" << "for(j=0; j < M; j=j+1) begin" << endl;
	fout << "     " << "bit = $random(seed1)%2;" << endl;
	fout << "     " << "rand = {{rand}, bit};" << endl;
	fout << "   " << "end" << endl;
    fout << "     " << "#1 $fwrite(fp1, \"%b\\n\", rand);" << endl;
    fout << "     " << "$fwrite(fp2, \"%b\\n\", y);" << endl;
    fout << "     " << "$fwrite(fp3, \"%d\", y[0]);" << endl;
//    fout << "     " << "$display(\"%d, %d\\n\", rand1, rand2);" << endl;
    
    fout << "  " << "end" << endl;
    fout << "  " << "$fclose(fp1);" <<endl;
    fout << "  " << "$fclose(fp2);" <<endl;
    fout << "  " << "$fclose(fp3);" <<endl;
    fout << "  " << "end" << endl;
    fout << "endmodule" << endl;
    fout.close();
}


void gen_testbench_sim_after_modify(int num_input, int num_output, int sample_num, int status, string &dut)
{
    ofstream fout;
    string filename = "./verilog_files/ckt_tb_sim_after_modify.v";
    fout.open(filename.c_str(), ofstream::out);
    
    srand((unsigned)time(0));
	unsigned int rand_number1 = rand();
	cout << "rand_number1 = " << rand_number1 << endl;

    fout << "module ckt_tb_sim_after_modify();" << endl;
    fout << "  " << "parameter M = " << num_input << ", N = " << num_output << ";" << endl;
    fout << "  " << "parameter sample_num = " << sample_num << ";" << endl;
    fout << "  " << "reg [M-1:0] rand;" << endl;
    fout << "  " << "reg bit;" << endl;
    fout << "  " << "wire [N-1:0] y;" << endl;
    fout << "  " << "integer i, j, cnt, fp1, fp2, fp3;" << endl;
    fout << "  " << "integer seed1 = " << rand_number1 << ";" << endl;
	fout << "  " << dut << " dut(";
    for(int i = num_input - 1; i >= 0; i--)
        fout << "rand[" << i << "], ";
    for(int j = num_output - 1; j >= 0; j--)
    {
    	if(j != 0)
        	fout << "y[" << j << "], ";	
        else
        	fout << "y[" << j << "] ";	
    }
    fout << ");" << endl;
    
    if (!status)
    {
	    fout << "  " << "initial begin " << endl;
	    fout << "   " << "fp1 = $fopen(\"rand.txt\", \"r\");" << endl;
	    fout << "   " << "fp2 = $fopen(\"simu_res.txt\", \"w\");" << endl;
	    fout << "   " << "fp3 = $fopen(\"simu_po.txt\", \"w\");" << endl;
	       
		fout << "  	" << "for(i=0; i < sample_num; i=i+1) begin" << endl;   
	    fout << "    " << "cnt = $fscanf(fp1, \"%b\", rand);" << endl;
	    fout << "    " << "#1 $fwrite(fp2, \"%b\\n\", y);" << endl;
	    fout << "    " << "$fwrite(fp3, \"%d\", y[0]);" << endl;
	    fout << "   " << "end" << endl;
	    
	    fout << "   " << "$fclose(fp1);" <<endl;
	    fout << "   " << "$fclose(fp2);" <<endl;
	    fout << "   " << "$fclose(fp3);" <<endl;
	    fout << "  " << "end" << endl;
	    fout << "endmodule" << endl;
	}
	else
	{
		fout << "  " << "initial begin " << endl;
	    fout << "  " << "fp1 = $fopen(\"rand.txt\", \"w\");" << endl;
	    fout << "  " << "fp2 = $fopen(\"simu_res.txt\", \"w\");" << endl;
	    fout << "  " << "fp3 = $fopen(\"simu_po.txt\", \"w\");" << endl;
	    
	    fout << "  " << "for(i=0; i < sample_num; i=i+1) begin" << endl;    
		fout << "  	" << "for(j=0; j < M; j=j+1) begin" << endl;
		fout << "     " << "bit = $random(seed1)%2;" << endl;
		fout << "     " << "rand = {{rand}, bit};" << endl;
		fout << "   " << "end" << endl;
	    fout << "     " << "#1 $fwrite(fp1, \"%b\\n\", rand);" << endl;
	    fout << "     " << "$fwrite(fp2, \"%b\\n\", y);" << endl;
	    fout << "     " << "$fwrite(fp3, \"%d\", y[0]);" << endl;
	    
	    fout << "  " << "end" << endl;
	    fout << "  " << "$fclose(fp1);" <<endl;
	    fout << "  " << "$fclose(fp2);" <<endl;
	    fout << "  " << "$fclose(fp3);" <<endl;
	    fout << "  " << "end" << endl;
	    fout << "endmodule" << endl;
	}
    fout.close();
}






/*1. gen_ckt_tb()*/
void gen_ckt_tb(int num_input, int num_output)
{
    ofstream fout;
    string path = "./verilog_files/";
    string filename = "ckt_tb.v";
    path.append(filename);
    filename = path;
    fout.open(filename.c_str(), ofstream::out);
//    fout << "`timescale 1ps/1ps" << endl;
    fout << "module ckt_tb();" << endl;
    fout << "  " << "parameter M = " << num_input << ", N = " << num_output << ";" << endl;
    fout << "  " << "parameter snum = " << sample_num << ";" << endl;
    fout << "  " << "parameter H1 = M/2 , H2 = M - H1;" << endl;
    fout << "  " << "reg [M-1:0] rand;" << endl;
    fout << "  " << "reg [H1-1:0] rand1;" << endl;
    fout << "  " << "reg [H2-2:0] rand2;" << endl;
    fout << "  " << "wire [N-1:0] y_ini, y_sim;" << endl;
    fout << "  " << "integer i, fp;" << endl;
    fout << "  " << "ckt_org dut1(";
    for(int i = numPI_ini - 1; i >= 0; i--)
        fout << "rand[" << i << "], ";
    for(int j = numPO_ini - 1; j >= 0; j--)
    {
        if(j == 0)
            fout << "y_ini[" << j << "] ";
        else
            fout << "y_ini[" << j << "], ";
    }
    fout << ");" << endl;
    fout << "  " << "ckt_sim dut2(";
    for(int i = numPI_ini - 1; i >= 0; i--)
        fout << "rand[" << i << "], ";
    for(int j = numPO_ini - 1; j >= 0; j--)
    {
        if(j == 0)
            fout << "y_sim[" << j << "] ";
        else
            fout << "y_sim[" << j << "], ";
    }
    fout << ");" << endl;
    fout << "  " << "initial begin " << endl;
    fout << "  " << "fp = $fopen(\"comp.txt\", \"w\");" << endl;
    if(num_input <= 20)
    {
    	fout << "  " << "for(i=0; i < 2**(M); i=i+1) begin" << endl;
	fout << "    " << "rand = i;" << endl;
    }
    else
    {
        fout << "  " << "for(i=0; i < 2**snum; i=i+1) begin" << endl;
		fout << "     " << "rand1 = {$random};" << endl;
        fout << "     " << "rand2 = {$random};" << endl;
        fout << "     " << "rand = {rand1, rand2};" << endl;
     }   
//    fout << "     " << "#1 $fwrite(fpp, \"rand = %b, y_ini = %b, y_sim = %b\\n\", rand, y_ini, y_sim);" << endl;
//    fout << "     " << "#1 if(y_ini != y_sim)" << endl;
//    fout << "      " << "$display(\"%b, y_ini = %b, y_sim = %b\\n\", rand, y_ini, y_sim);" << endl;
    fout << "     " << "#1 if(y_ini != y_sim)" << endl;
    fout << "      " << "$fwrite(fp, \"%b, y_ini = %b, y_sim = %b\\n\", rand, y_ini, y_sim);" << endl;
    fout << "  " << "end" << endl;
    fout << "  " << "$fclose(fp);" <<endl;
//    fout << "  " << "$fclose(fpp);" <<endl;
    fout << "  " << "end" << endl;
    fout << "endmodule" << endl;
    fout.close();
}



void rewrite_bigNode_sim(char *cnode)
{
	ifstream fin;
	fin.open("./blif_files/bigNode_sim.blif", ios::in);
	ofstream fout;
	fout.open("./blif_files/bigNode_sim_rew.blif", ios::out);
	string str, s;
//	cout << "bigNode_sim.blif: " << endl;
	while(getline(fin, str))
	{
//		cout << str << endl;
		istringstream ss(str);
		ss >> s;
		if(s == ".names" || s == ".exdc")
			break;
		fout << str << endl;
	}
	fin.close();

	char com[100];
	sprintf(com, "sis -t none -f ./script/blif2pla.rug > sim_bignode.txt");
    system(com);  
    fin.open("./pla_files/bigNode_sim.pla", ios::in);
//    cout << "bigNode_sim.pla: " << endl;
    fout << ".names ";
    while(getline(fin, str))
	{
//		cout << str << endl;
		istringstream ss(str);
		ss >> s;
		if(s == ".ilb")
		{
			while(ss >> s)
				fout << s << " ";
			fout << cnode << endl;
		}
		else if(s[0] == '.')
			continue;
		else
		{
			while(ss >> s);
			if(s == "2")
				break;
			fout << str << endl;
		}
	}
	fout << ".end" << endl;
	fout.close();
	fin.close();

/*	char com1[100];	
	sprintf(com1, "rm -rf ./blif_files/bigNode_sim.blif");
    system(com1);  
    char com2[200];	
    sprintf(com2, "cp ./blif_files/bigNode_sim_rew.blif ./blif_files/bigNode_sim.blif");
    system(com2); 
*/    
/*    fin.open("./blif_files/bigNode_sim_rew.blif", ios::in);
    cout << "after rewrite, bigNode_sim_rew.blif: " << endl;
	while(getline(fin, str))
	{
		cout << str << endl;
	}
	fin.close();
*/    
}


void write_ckt_whole(BnetNetwork *net)
{
	char com[100];
	sprintf(com, "cp ./blif_files/ckt_sim.blif ./blif_files/ckt_sim_whole.blif");
	system(com);
	sprintf(com, "sed -i '/.end/d' ./blif_files/ckt_sim_whole.blif");
	system(com);
	ifstream fin;
	fin.open("minterm.log", ios::in);
	string str;
	vector<string> minterm;
	while(getline(fin, str))
		minterm.push_back(str);
	fin.close();
	ofstream fout;
	fout.open("./blif_files/ckt_sim_whole.blif", ios::app);
	fout << endl << ".exdc" << endl;
	fout << ".inputs ";
	for(int i = 0; i < net->npis; i++)
		fout << net->inputs[i] << " ";
	fout << endl;
	fout << ".outputs ";
	for(int i = 0; i < net->npos; i++)
		fout << net->outputs[i] << " ";
	fout << endl;
	for(int i = 0; i < net->npos; i++)
	{
		fout << ".names ";
		for(int j = 0; j < net->npis; j++)
			fout << net->inputs[j] << " ";
		fout << net->outputs[i] << endl;
		for(int j = 0; j < minterm.size(); j++)
			fout << minterm[j] << endl;
	}
		
	fout << ".end" << endl;
	fout.close();
}



void write_cone_circuit(BnetNetwork *net, int iIndex, char *po, set<char*> &tfi, string &fn, set<string> &inputs_set, vector<string> &cone_string)
{
	BnetNode *nd, *tmp;
	BnetTabline *tl;
	set<char*>::iterator itrs;

	ofstream fout;
	fout.open(fn.c_str(), ios::out);
	fout << ".model cone_" << po << endl;
	fout << ".inputs ";
	for(int i = 0; i < net->npis; i++)
	{
		itrs = tfi.find(net->inputs[i]);
		if(itrs == tfi.end()) continue;
		fout << net->inputs[i] << " ";
		
		if (iIndex == 0)
		{	
			string cin(net->inputs[i]);
			inputs_set.insert(cin);
		}
	}
	fout << endl;
	fout << ".outputs " << po << endl;
	
	nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type == BNET_INPUT_NODE) 
		{
			nd = nd->next;
			continue;
		} 
		itrs = tfi.find(nd->name);
		if(itrs == tfi.end())
		{
			nd = nd->next;
			continue;
		} 
		
		fout << ".names ";
		for (int i = 0; i < nd->ninp; i++)
			fout << nd->inputs[i] << " ";
		fout << nd->name << endl;
	    tl = nd->f;
	    while (tl != NULL) 
	    {
			if (tl->values != NULL) 
				fout << tl->values << " " << 1 - nd->polarity << endl;
			else 
				fout << 1 - nd->polarity << endl;
			tl = tl->next;
	    }
	    
		if (iIndex == 0)
		{
		    //fill cone_string
		    string str(".names");
		    str.append(" ");
		    for (int i = 0; i < nd->ninp; i++)
			{
				string fin(nd->inputs[i]);
				str.append(fin);
				str.append(" ");
			}
			string snode(nd->name);
			str.append(snode);
			cone_string.push_back(str);
			tl = nd->f;
		    while (tl != NULL) 
		    {
		    	str.clear();
				if (tl->values != NULL) 
				{
					string value(tl->values);
					value.append(" ");
					if (1-nd->polarity) value.append("1");
					else 			    value.append("0");
					str = value;
				}
				else 
				{
					string value;
					if (1-nd->polarity) value.append("1");
					else 			    value.append("0");
					str = value;
				}
				cone_string.push_back(str);
				
				tl = tl->next;
		    }
		}
	    
		nd = nd->next;
    }
    fout << ".end" << endl;
    
    fout.close();
}


void write_cone_po_xor(BnetNetwork *net, set<char*> &tfi, set<string> &inputs_set, vector<string> &cone_string, char *cnode, string &po, vector<string> &new_pla, string &filename)
{
	BnetNode *nd, *tmp;
	BnetTabline *tl;
	set<char*>::iterator itrs;
	set<string>::iterator itrss;

	ofstream fout;
	fout.open(filename.c_str(), ios::out);
	fout << ".model cone_xor_" << po << endl;
	fout << ".inputs ";
	for(itrss = inputs_set.begin(); itrss != inputs_set.end(); itrss++)
		fout << *itrss << " ";
	fout << endl;
	fout << ".outputs final_po" << endl;
	
	//print the lines for the original cone
	for(int i = 0; i < cone_string.size(); i++)
		fout << cone_string[i] << endl;
	fout << endl;	
	
	//print the lines for the current cone
	nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type == BNET_INPUT_NODE || !strcmp(nd->name, cnode)) 
		{
			nd = nd->next;
			continue;
		} 
		itrs = tfi.find(nd->name);
		if(itrs == tfi.end())
		{
			nd = nd->next;
			continue;
		} 
		
		fout << ".names ";
		for (int i = 0; i < nd->ninp; i++)
		{
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if (tmp->type == BNET_INPUT_NODE)
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
	    
	    nd = nd->next;
	}
	    
    //write the lines for cnode
    st_lookup(net->hash, cnode, &nd);	
	if(new_pla.empty())
		fout << endl << ".names " << cnode << "sim" << endl;
	else
	{
		fout << endl << ".names ";
		for(int i = 0; i < nd->ninp; i++)
		{
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if (tmp->type == BNET_INPUT_NODE)
				fout << nd->inputs[i] << " ";
			else
				fout << nd->inputs[i] << "sim ";
		}
		fout << nd->name << "sim" << endl;
		for(int i = 0; i < new_pla.size(); i++)
		    fout << new_pla[i] << " 1" << endl;
	}
	
	//write the XOR gate
	fout << endl << ".names " << po << " " << po << "sim final_po" << endl;
	fout << "10 1" << endl;
	fout << "01 1" << endl; 
	
    fout << ".end" << endl;
    fout.close();
}



void write_cone_po_sub(BnetNetwork *net, set<char*> &tfi, set<string> &inputs_set, vector<string> &cone_string, char *cnode, string &po, vector<string> &new_pla, string &filename)
{
	BnetNode *nd, *tmp;
	BnetTabline *tl;
	set<char*>::iterator itrs;
	set<string>::iterator itrss;

	ofstream fout;
	fout.open(filename.c_str(), ios::out);
	fout << ".model cone_sub_" << po << endl;
	fout << ".inputs ";
	for(itrss = inputs_set.begin(); itrss != inputs_set.end(); itrss++)
		fout << *itrss << " ";
	fout << endl;
	fout << ".outputs final_po" << endl;
	
	//print the lines for the original cone
	for(int i = 0; i < cone_string.size(); i++)
		fout << cone_string[i] << endl;
	fout << endl;	
	
	//print the lines for the current cone
	nd = net->nodes;    
    while (nd != NULL) 
    {
		if(nd->type == BNET_INPUT_NODE || !strcmp(nd->name, cnode)) 
		{
			nd = nd->next;
			continue;
		} 
		itrs = tfi.find(nd->name);
		if(itrs == tfi.end())
		{
			nd = nd->next;
			continue;
		} 
		
		fout << ".names ";
		for (int i = 0; i < nd->ninp; i++)
		{
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if (tmp->type == BNET_INPUT_NODE)
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
	    
	    nd = nd->next;
	}
	    
    //write the lines for cnode
    st_lookup(net->hash, cnode, &nd);	
	if(new_pla.empty())
		fout << endl << ".names " << cnode << "sim" << endl;
	else
	{
		fout << endl << ".names ";
		for(int i = 0; i < nd->ninp; i++)
		{
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if (tmp->type == BNET_INPUT_NODE)
				fout << nd->inputs[i] << " ";
			else
				fout << nd->inputs[i] << "sim ";
		}
		fout << nd->name << "sim" << endl;
		for(int i = 0; i < new_pla.size(); i++)
		    fout << new_pla[i] << " 1" << endl;
	}
	
	//write the XOR gate
	fout << endl << ".names " << po << " " << po << "sim final_po" << endl;
	fout << "10 1" << endl;
	fout << "01 1" << endl; 
	
    fout << ".end" << endl;
    fout.close();
}




void write_compare_circuit_ave(BnetNetwork *net, char *cnode, vector<string> &ckt_org_po, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, vector<string> &final_pla, string &filename)
{
	ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model compare_circuit" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    for (int i = 0; i < sub_abs_po.size(); i++)
    	fout << sub_abs_po[i] << " ";
    fout << endl;

	//step1. write the original correct circuit
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
 	
 	//step2. write the current circuit with cnode changed
 	fout << "#simplified circuit: " << endl;
 	BnetTabline *tl;
 	BnetNode *nd = net->nodes; 
 	BnetNode *tmp;   
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, cnode)) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
			{
				st_lookup(net->hash, nd->inputs[i], &tmp);
				if (tmp->type == BNET_INPUT_NODE)
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

	cout << "simplified bignode: " << endl;
	st_lookup(net->hash, cnode, &nd);	
	if(final_pla.empty())
		fout << endl << ".names " << cnode << "sim" << endl;
	else
	{
		fout << endl << ".names ";
		for(int i = 0; i < nd->ninp; i++)
		{
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if (tmp->type == BNET_INPUT_NODE)
				fout << nd->inputs[i] << " ";
			else
				fout << nd->inputs[i] << "sim ";
		}
		fout << nd->name << "sim" << endl;
		for(int i = 0; i < final_pla.size(); i++)
		    fout << final_pla[i] << " 1" << endl;
	}
 	
 	
 	//step3. add the buffers that connect the two circuits with the sub_abs circuit
// 	fout << endl << "# buffers_1: " << endl;
 	int index = 0;
 	for (int i = 0; i < ckt_org_po.size(); i++)
 	{
 		fout << ".names " << ckt_org_po[i] << " " << sub_abs_pi[i] << endl;
 		fout << "1 1" << endl;
 		index++; 		
 	}
 	fout << endl << "# buffers_2: " << endl;
 	for (int j = 0; j < net->npos; j++)
 	{
 		fout << ".names " << net->outputs[j] << "sim " << sub_abs_pi[index++] << endl;
 		fout << "1 1" << endl;
 	}
 	
 	//step4. write the sub_abs circuit
 	fout << endl << "# sub_abs circuit" << endl;
 	for (int i = 0; i < sub_abs_ckt.size(); i++)
 		fout << sub_abs_ckt[i] << endl;

 		
 	fout << ".end" << endl;
 	fout.close();
}


void write_compare_circuit_max(BnetNetwork *net, char *cnode, vector<string> &ckt_org_po, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, vector<string> &comparator_ckt,  vector<string> &comparator_pi, vector<string> &comparator_po, vector<int> &comp_number, vector<string> &final_pla, string &filename, int flag_or)
{
	ofstream fout;
    fout.open(filename.c_str(), iostream::out);
    fout << ".model compare_circuit_new" << endl;
    fout << ".inputs ";
    for(int i = 0; i < net->npis; i++)
        fout << net->inputs[i] <<  " ";
    fout  << endl;    
    fout << ".outputs ";
    if (flag_or == 0)
    {
	for (int i = 0; i < comparator_po.size(); i++)
    		fout << comparator_po[i] << " ";
    }
    else if (flag_or == 1)
    	fout << "or_output" << endl;
   fout << endl;
	
	//step1. write the original correct circuit
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
 	
 	//step2. write the current circuit with cnode changed
 	fout << "#simplified circuit: " << endl;
 	BnetTabline *tl;
 	BnetNode *nd = net->nodes; 
 	BnetNode *tmp;   
    while (nd != NULL) 
    {
		if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, cnode)) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
			{
				st_lookup(net->hash, nd->inputs[i], &tmp);
				if (tmp->type == BNET_INPUT_NODE)
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

	cout << "simplified bignode: " << endl;
	st_lookup(net->hash, cnode, &nd);	
	if(final_pla.empty())
		fout << endl << ".names " << cnode << "sim" << endl;
	else
	{
		fout << endl << ".names ";
		for(int i = 0; i < nd->ninp; i++)
		{
			st_lookup(net->hash, nd->inputs[i], &tmp);
			if (tmp->type == BNET_INPUT_NODE)
				fout << nd->inputs[i] << " ";
			else
				fout << nd->inputs[i] << "sim ";
		}
		fout << nd->name << "sim" << endl;
		for(int i = 0; i < final_pla.size(); i++)
		    fout << final_pla[i] << " 1" << endl;
	}
 	
 	
 	//step3. add the buffers that connect the two circuits with the sub_abs circuit
// 	fout << endl << "# buffers_1: " << endl;
 	int index = 0;
 	for (int i = 0; i < ckt_org_po.size(); i++)
 	{
 		fout << ".names " << ckt_org_po[i] << " " << sub_abs_pi[i] << endl;
 		fout << "1 1" << endl;
 		index++; 		
 	}
// 	fout << endl << "# buffers_2: " << endl;
 	for (int j = 0; j < net->npos; j++)
 	{
 		fout << ".names " << net->outputs[j] << "sim " << sub_abs_pi[index++] << endl;
 		fout << "1 1" << endl;
 	}
 	
 	//step4. write the sub_abs circuit
 	fout << endl << "# sub_abs circuit" << endl;
 	for (int i = 0; i < sub_abs_ckt.size(); i++)
 		fout << sub_abs_ckt[i] << endl;
 		
 	//step5. add the buffers that connect the sub_abs circuit with the comparator circuit
 	//and the comparator circuit
 	index = 0;
 	if (flag_or == 0)
	{
		for (int i = 0; i < sub_abs_po.size(); i++)
		{
			fout << ".names " << sub_abs_po[i] << " " << comparator_pi[i] << endl;
			fout << "1 1" << endl;
			index++; 		
		}
		for (int i = 0; i < comparator_pi.size()/2; i++)
		{
			fout << ".names " << comparator_pi[index++] << endl;
			fout << comp_number[i] << endl;
		}
 		
 		fout << endl << "# comparator circuit" << endl;
 		for (int i = 0; i < comparator_ckt.size(); i++)
 			fout << comparator_ckt[i] << endl;
	}	
	else if (flag_or == 1)
	{
 		fout << endl << "# the big OR gate" << endl;
		fout << ".names ";
		for (int i = 0; i < sub_abs_po.size(); i++)
			fout << sub_abs_po[i] << " ";
		fout << "or_output" << endl;
		for (int i = 0; i < sub_abs_po.size(); i++)
			fout << "0";
		fout << " 0" << endl;
	}
 		
 	fout << ".end" << endl;
 	fout.close();
}
