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
#include "head/basics.h"
#include "cudd/bnet.h"
#include "cudd/cudd_build_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"


using namespace std;

extern int numPI_ini;

/*
functions in this file:

*/

//Global variables and external variables


void sub_sim_ckt(BnetNetwork *net, set<char*> &ts_MFFC, set<char*> &ss_TFI, char *ts_best, char *ss_best, int &pol_best, int iIndex)
{
//	cout << "in sub_sim_ckt, ts_best = " << ts_best << ", ss_best = " << ss_best << endl;
	
	set<char*>::iterator itrs, itrs1;
	BnetNode *nd, *tmp;
	BnetTabline *tl;
	
/*	cout << "ts_MFFC:" << endl;
	for(itrs = ts_MFFC.begin(); itrs != ts_MFFC.end(); itrs++)
		cout << *itrs << " ";
	cout << endl;
	
	cout << "ss_TFI:" << endl;
	for(itrs = ss_TFI.begin(); itrs != ss_TFI.end(); itrs++)
		cout << *itrs << " ";
	cout << endl;
*/	
    string filename = "./blif_files/ckt_sim.blif";
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
    
    if(!strcmp(ss_best, "c0") || !strcmp(ss_best, "c1"))
    {
    	nd = net->nodes;    
	    while (nd != NULL) 
	    {
	    	itrs = ts_MFFC.find(nd->name);
		//	if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, ts_best)  && (itrs == ts_MFFC.end())) //ignore PI nodes and ts_best node
			if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, ts_best)) //ignore PI nodes and ts_best node
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
	    st_lookup(net->hash, ts_best, &nd);
	    fout << ".names ";
		for (int i = 0; i < nd->ninp; i++)
			fout << nd->inputs[i] << " ";
		fout << nd->name << endl;
		for (int i = 0; i < nd->ninp; i++)
			fout << "-";
	    if(!strcmp(ss_best, "c0"))
			fout << " 0" << endl;
	    else if(!strcmp(ss_best, "c1"))
	    	fout << " 1" << endl;
	    fout << ".end" << endl;
    	fout.close();    	
    	return;
    }
    
    char *new_ts_name = new char[30];
    strcpy(new_ts_name, ss_best);
    strcat(new_ts_name, "_ss");
    char *index = new char[5]; 
    sprintf(index, "%d", iIndex);
    strcat(new_ts_name, index);   
 /*   cout << "ts_MFFC: " << endl;
    for(itrs = ts_MFFC.begin(); itrs != ts_MFFC.end(); itrs++)
    {
    	char *node = *itrs;
    	cout << node << " ";
    }
 */
    cout << endl;
    nd = net->nodes;    
    while (nd != NULL) 
    {
    	itrs = ts_MFFC.find(nd->name);
		if(nd->type != BNET_INPUT_NODE && strcmp(nd->name, ts_best) && (itrs == ts_MFFC.end())) 
		{
			fout << ".names ";
			for (int i = 0; i < nd->ninp; i++)
			{				
				st_lookup(net->hash, nd->inputs[i], &tmp);
				if(!strcmp(ts_best, tmp->name))
					fout << new_ts_name << " ";
				else
					fout << nd->inputs[i] << " ";
			}
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
    
    st_lookup(net->hash, ss_best, &nd);
    itrs = ts_MFFC.find(nd->name);
    int flag_ignore;
    if(itrs != ts_MFFC.end())
    {
    	cout << "ss is in tfi of ts!" << endl;
    	flag_ignore = 0;
    }
    else 
    	flag_ignore = 1;
    
    if(!flag_ignore)
    {
    	for(itrs = ss_TFI.begin(); itrs != ss_TFI.end(); itrs++)
	    {
	    	char *cnode = *itrs;
	    	st_lookup(net->hash, cnode, &nd);
			itrs1 = ts_MFFC.find(nd->name);
			if(itrs1 == ts_MFFC.end())
				continue;
	    	if(nd->type == BNET_INPUT_NODE)
	    		continue;
	    	fout << endl << ".names ";
	    	for(int i = 0; i < nd->ninp; i++)			
				fout << nd->inputs[i] << " ";
			fout << cnode << endl;
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
    }	
    
    if(pol_best == 1)
	{
		fout << endl << ".names " << ss_best << " " << new_ts_name << endl;
		fout << "1 1" << endl;
	}
	else if(pol_best == 0)
	{
		fout << endl << ".names " << ss_best << " " << new_ts_name << endl;
		fout << "0 1" << endl;
	}
	
	st_lookup(net->hash, ts_best, &nd);
	if(nd->type == BNET_OUTPUT_NODE)
	{
		fout << ".names " << new_ts_name << " " << ts_best << endl;
		fout << "1 1" << endl;
	}
	
	fout << ".end" << endl;
	fout.close();
	
	delete []new_ts_name;
	delete []index;
		
	cout << "after substitution, ckt_sim.blif: " << endl;
	ifstream fin;
	fin.open("./blif_files/ckt_sim.blif", ios::in);
	string str;
	while(getline(fin, str))
		cout << str << endl;
	cout << endl;
	fin.close();

    
}
