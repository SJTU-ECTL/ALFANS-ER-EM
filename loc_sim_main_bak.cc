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
#include "head/write_func.h"
#include "head/comp_real_er.h"
#include "head/call_abc.h"
#include "head/helper.h"
#include "head/sim.h"
#include "cudd/bnet.h"
#include "cudd/cudd_build_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"


using namespace std;

extern int numPI_ini;
int sample_num;

/*
functions in this file:

*/

//Global variables and external variables



/*1. loc_sim_ite*/
int loc_sim_ite(BnetNetwork *net, BnetNetwork *&net_comb, DdManager **dd_comb, char *last_node, vector<char*> &last_inputs, double threshold, double &real_er, int &iIndex)
{
    //iterators and variables
    struct timeb startTime, endTime; 
    BnetNode *nd;
    multimap<double, char*> ave_sp;
    multimap<double, char*>::iterator itrm_dc;
    FILE *fp;
    
    cout <<"step0. net && net_comb: " << endl;
/*    Bnet_PrintNetwork(net);
    if(iIndex > 0)
    	Bnet_PrintNetwork(net_comb);
*/
    
    /*Build BDD for the entire circuit: ckt_org.blif or ckt_sim.blif and compute signal probabilities*/
    cout << "****************************" << endl;
    cout << "step1. Build BDD for current network and compute signal probabilities: " << endl;
    ftime(&startTime);
    string path_bdd = "./blif_files/";
    string filename_bdd;
    if(iIndex == 0)
		filename_bdd = "ckt_org.blif";
    else
		filename_bdd = "ckt_sim.blif";
    path_bdd.append(filename_bdd);
    filename_bdd = path_bdd;
    DdManager *dd = NULL;
    dd = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    if (dd == NULL) return -1;      
    cudd_build_v2(net, &dd, filename_bdd.c_str(), BNET_GLOBAL_DD);
    nd = net->nodes;
    while(nd != NULL)
    {
//    	cout << "Current node: " << nd->name << endl;
    	double num_minterm = Cudd_CountMinterm(dd, nd->dd, numPI_ini);
    	double sp = num_minterm/pow(2.0, numPI_ini);
    	nd->rp = sp;
    	nd->p = findmin(sp);
//    	cout << "node " << nd->name << ", sp: " << sp << endl;
    	nd = nd->next;
    }
    ftime(&endTime);
    double rt_step1 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for compute signal probabilities: " << rt_step1 << endl;
    
    //Image computation: utilize dd_comb
/*    if(iIndex > 0) 
	{
	    cout << "****************************" << endl;
	    cout << "step2. Read ckt_sim_last.blif and build BDD for ckt_sim_last.blif: " << endl;
	    cout << "last_node: " << last_node << endl;
	    FILE *fp;
	    fp = fopen("./blif_files/ckt_sim_last.blif", "r");
	    net_last = Bnet_ReadNetwork(fp);
	    Bnet_PrintNetwork(net_last);
	    fclose(fp);
	    dd_last = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
	    if (dd_last == NULL) return -1;      
	    cudd_build_v2(net_last, &dd_last, "./blif_files/ckt_sim_last.blif", BNET_GLOBAL_DD);
    }
*/    
       
    //Obtain the average value of input node signal probabilities for each two-level nodes
    cout << "****************************" << endl;
    cout << "step3. find_ave_sp: " << endl;
    ftime(&startTime);
    find_ave_sp(net, ave_sp);
    ftime(&endTime);
    double rt_step3 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for step3: " << rt_step3 << endl;
    
        
    //Find exdcs for each big node and simplify them using these exdcs. Obtain 
    //the one with maximum node save
    cout << "****************************" << endl;
    cout << "step4. find_exdc_sim: " << endl;	
    ftime(&startTime);
    //store the files for the simplifed big node
    map<char*, vector<string> > node_sim_files; 
    map<char*, vector<string> >::iterator itrm_cv;
    double max_score = 0;
    char *max_save_node = find_exdc_sim(net, &dd, net_comb, dd_comb, ave_sp, node_sim_files, last_node, last_inputs, max_score, threshold, iIndex);     
    ftime(&endTime);
    double rt_step4 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for step4: " << rt_step4 << endl;     
    Cudd_Quit(dd);
   	
    if(iIndex > 0)  
    {
    	Bnet_FreeNetwork(net_comb);
    	Cudd_Quit(*dd_comb);
    } 
     
    if(max_score <= 0)
    {
/*    	if(iIndex > 0)
	    {
	    	string path = "./blif_files/";
		    string filename = "ckt_org_sim.blif";
		    path.append(filename);
		    filename = path;
		    write_ckt_comb(net);
		    fp = fopen(filename.c_str(), "r");
		    cout << "read combined network: " << endl;
		    BnetNetwork *net_comb = Bnet_ReadNetwork(fp); 
		    fclose(fp);    
			DdManager *dd_comb = NULL;
			dd_comb = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
		    if (dd_comb == NULL) return -1;      
			cudd_build_v2(net_comb, &dd_comb, filename.c_str());
			char *out_node = net_comb->outputs[0];
			st_lookup(net_comb->hash, out_node, &nd);
			FILE *fp_out;
			fp_out = fopen("minterm.log", "w");
			dd_comb->out = fp_out;
			Cudd_PrintMinterm(dd_comb, nd->dd);
			fclose(fp_out);
			write_ckt_whole(net);
			char com[100];
			sprintf(com, "sis -t none -f ./script/sim_ckt_whole.rug > sis_whole.txt");
			system(com);
		}
*/
	    return -1;
    }
    

    //Extract the structure for the simplified max_node
    cout << "****************************" << endl;
    cout << "step5. extract_sim_nodes: " << endl;	
    ftime(&startTime);    
    itrm_cv = node_sim_files.find(max_save_node);
    vector<string> sim_node_lines = itrm_cv->second;
    vector<string> sim_node_nodes;
    extract_sim_nodes(sim_node_lines, sim_node_nodes);  
    cout << "sim_node_nodes:" << endl;
    for(int i = 0; i < sim_node_nodes.size(); i++)
    	cout << sim_node_nodes[i] << endl;
    cout << endl;
    double rt_step5 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for step5: " << rt_step5 << endl;   

    //Write the whole simplified circuit into ckt_sim.blif
    cout << "****************************" << endl;
    cout << "step6. write_ckt_sim && sweep && update it to add buffer: " << endl;
    ftime(&startTime);
    write_ckt_sim(net, max_save_node, sim_node_nodes);
    cout << "call SIS to sweep ckt_sim further" << endl;  
	char com[100];      
    sprintf(com, "sis -t none -f ./script/sim_ckt_v2.rug > sis_ckt.txt");
    system(com); 
    Bnet_FreeNetwork(net);
    update_ckt(net);
    ftime(&endTime);
    double rt_step6 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for step6: " << rt_step6 << endl;  

    //Build bdd to compute the real error rate
	cout << "****************************" << endl;
    cout << "step7. Write ckt_org_blif and build bdd to compute real error rate: " << endl;
    ftime(&startTime);
    string path = "./blif_files/";
    string filename = "ckt_org_sim.blif";
    path.append(filename);
    filename = path;
    write_ckt_comb(net);
    fp = fopen(filename.c_str(), "r");
    cout << "read combined network: " << endl;
    net_comb = Bnet_ReadNetwork(fp); 
//    cout << "print combined network: " << endl;
//    Bnet_PrintNetwork(net_comb); 
    fclose(fp);  
    *dd_comb = NULL;
    *dd_comb = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);    
	cudd_build_v2(net_comb, dd_comb, filename.c_str(), BNET_GLOBAL_DD);
	char *out_node = net_comb->outputs[0];
	st_lookup(net_comb->hash, out_node, &nd);
	double num_diff = Cudd_CountMinterm(*dd_comb, nd->dd, numPI_ini);
	real_er = num_diff/pow(2.0, numPI_ini);
//	cout << "diff_pla: " << endl;
//	Cudd_PrintDebug(*dd_comb, nd->dd, numPI_ini, 2);
	ftime(&endTime);
    double rt_step7 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for step7: " << rt_step7 << endl;  
	
    return 0;
}





/*2. loc_sim_main()*/
void loc_sim_main(BnetNetwork *net, string &filename, double threshold)
{
    struct timeb st, et; 
	double threshold_org = threshold;
    int iIndex = 0;
    char com[100];
    string str;
    
    write_ckt_org(net);   //add buffers to output nodes
    ifstream fin;
    fin.open("./blif_files/ckt_org.blif", ios::in);
    cout << "modified ckt_org: " << endl;
    while(getline(fin, str))
    	cout << str << endl;
    cout << endl;    
    fin.close();
    
    Bnet_FreeNetwork(net);
    FILE *fp;
    fp = fopen("./blif_files/ckt_org.blif", "r");
    net = Bnet_ReadNetwork(fp);
    fclose(fp);
    
    //If this is the first iteration, write the inital circuit into ckt_org.blif
/*    char new_filename[100] = "./blif_files/ckt_org.blif";    
    sprintf(com, "cp %s %s", filename.c_str(), new_filename);
    system(com);
    sprintf(com, "sed -i '/^.model/d' ./blif_files/ckt_org.blif");
    system(com);
    sprintf(com, "sed -i '1i .model ckt_org' ./blif_files/ckt_org.blif");
    system(com);
*/
    
    
    //start the iteration of local simplification
    double real_er = 0;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "i" << iIndex << ". iterations start: " << endl;
    cout << "threshold = " << threshold << endl;
	char *last_node = new char[50];
	vector<char*> last_inputs;
	BnetNetwork *net_comb;
	DdManager *dd_comb;
    ftime(&st);  
    int res = loc_sim_ite(net, net_comb, &dd_comb, last_node, last_inputs, threshold, real_er, iIndex);
    ftime(&et);
    double rt_loc_sim_ite = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    cout << "i" << iIndex << ". runtime for loc_sim_ite: " << rt_loc_sim_ite << endl;

    if(res == -1)
    {
    	cout << "max_save = 0. No more optimization!" << endl;
    	delete []last_node;
    	for(int i = 0; i < last_inputs.size(); i++)
    		delete []last_inputs[i];
    	return;
    }
    cout << endl << "##############################################" << endl;
    cout << "report: " << endl;
    cout << "1. real_er = " << real_er << endl;    
	cout <<  "2. current stats by SIS: " << endl;        
    sprintf(com, "sis -t none -f ./script/map.rug");
    system(com);
    
    iIndex++;
    threshold =  threshold_org -  real_er; 
	cout << "updated threshold = " << threshold << endl;
    
    		
    while(real_er < threshold_org)
    {
        cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    	cout << "i" << iIndex << ". iterations start: " << endl;
        ftime(&st);  
        res = loc_sim_ite(net, net_comb, &dd_comb, last_node, last_inputs, threshold, real_er, iIndex);
	    ftime(&et);
    	rt_loc_sim_ite = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    	cout << "i" << iIndex << ". runtime for loc_sim_ite: " << rt_loc_sim_ite << endl;
    	    
    	if(res == -1)
	    {
			cout << "max_save = 0. No more optimization!" << endl;
			delete []last_node;
			for(int i = 0; i < last_inputs.size(); i++)
    			delete []last_inputs[i];
			return;
	    }	
	       
        cout << endl << "##############################################" << endl;
		cout << "report: " << endl;
		cout << "1. real_er = " << real_er << endl;    
		cout <<  "2. current stats by SIS: " << endl;        
		sprintf(com, "sis -t none -f ./script/map.rug");
		system(com);
	                        	    
	    iIndex++;
	    threshold = threshold_org -  real_er; 
	    cout << "updated threshold = " << threshold << endl;
	    	    	    		    	
    }
        
    delete []last_node;
    return;

}
