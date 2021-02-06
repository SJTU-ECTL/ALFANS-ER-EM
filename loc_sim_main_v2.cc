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
#include "head/loc_sim_main_v2.h"
#include "head/sim_new_abc_max.h"
#include "head/sim_new_abc_ave.h"
#include "head/simu_ckt.h"
#include "head/read_file.h"
#include "cudd/bnet.h"
#include "cudd/cudd_build_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"


using namespace std;

//Global variables
vector<string> ckt_org_po;
map<string, set<char*> > po_tfi;
map<string, set<string> > po_inputs;
map<string, vector<string> > po_cone_string;
vector<string> sub_abs_ckt, sub_abs_pi, sub_abs_po;
vector<string> comparator_ckt, comparator_pi, comparator_po;
vector<int> comp_number;

int flag_dyn = 0;
double er_weight = 1, em_weight = 1;
double ini_ratio = er_weight/em_weight;
double ave_error_max, ave_error_max_cur, max_weight_max, max_weight_max_cur, real_er_max, real_er_max_cur;
int lit_save_pre = 0;

//external variables
extern int em_class; //1: max, 2: ave
extern int numPI_ini, numPO_ini;
extern double ini_threshold_er, ini_threshold_em;
extern string diff_file, comparator_file;


/*
functions in this file:
int loc_sim_ite();

int loc_sim_main();
*/




/*1. loc_sim_ite*/
int loc_sim_ite(BnetNetwork *net, BnetNetwork *&net_comb, DdManager **dd_comb, map<string, struct score_pla> &sim_record, multimap<double, struct score_pla> &sim_record_top, map<string, map<string, double> > &node_pattern_rate, map<string, int> &internal_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int &num_output, double threshold_er, double threshold_em, double &real_er, double &real_em, int &min_modified_po, int &max_modified_po, double &T_er, double T_er1, double &T_em, double T_em1, int &iIndex)
{
    //variables
    struct timeb startTime, endTime; 
    BnetNode *nd, *nd1;
    multimap<double, char*> ave_sp;
    FILE *fp;
    DdManager *dd = NULL; 
    cout << "at the start of loc_sim_ite, real_er = " << real_er << endl;    
    cout << "min_modified_po = " << min_modified_po << ", max_modified_po = " << max_modified_po << endl;

    //iterators
    multimap<double, char*>::iterator itrm_dc;    
    set<char*>::iterator itrs;
    map<string, struct wi_pair>::iterator itrm_sw;
    
	//print out the affected primary outputs of each node
    	nd = net->nodes;
	while(nd != NULL)
	{
		nd1 = net->nodes;
		while(nd1 != NULL)
		{
			nd1->visited = 0;
			nd1 = nd1->next;
		}
		set<char*> po_set;
		find_tranfanout_po(net, nd->name, po_set);
		cout << endl << "current node " << nd->name << "'s affected pos: " << po_set.size() << endl;
		for(itrs = po_set.begin(); itrs != po_set.end(); itrs++)
			cout << *itrs << " ";
		cout << endl;
		nd = nd->next;
	}
	
	map<string, double> node_sp_simu;
	map<string, double>::iterator itrm_sd;
	double ave_error_mag = 0;
	
	if(iIndex == 0)
	{
		cout << "running simulation for the first iteration: " << endl;
		simu_ckt_both(net, node_pattern_rate, node_sp_simu, internal_index, sim_output_wi, pis, rand, simu_res, num_output, ave_error_mag, real_er, iIndex, 0);
	}
    
    	//Find exdcs for each big node and simplify them using these exdcs. Obtain 
    	//the one with maximum node save
    	cout << "****************************" << endl;
    	cout << "step2. find_exdc_sim: " << endl;	
    	ftime(&startTime);    
    	vector<string> final_pla;  //store the files for the simplifed big node
    	double max_score = -1;
	char *max_save_node;
    	cout << "iteration = " << iIndex << ", T_em = " << T_em << ", T_er = " << T_er << endl;
	if (em_class == 1)
	{
    		max_save_node = find_exdc_sim_max(net, &dd, net_comb, dd_comb, final_pla, sim_record, sim_record_top, internal_index, sim_output_wi, pis, rand, simu_res, num_output, node_pattern_rate, max_score, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_em, T_em1, iIndex);     
	}
	else
	{
		if (flag_dyn && iIndex > 0)
		{
			double real_er_tmp;
			if (real_er == 0) real_er_tmp = 0.005;
			else real_er_tmp = real_er;
			double er_margin = ini_threshold_er/real_er_tmp;
			double em_margin = ini_threshold_em/real_em;
			double ratio = ini_ratio * (em_margin/er_margin); 
			em_weight = 1/(1+ratio);
			er_weight = 1 - em_weight;
		}
    		max_save_node = find_exdc_sim_ave(net, &dd, net_comb, dd_comb, final_pla, sim_record, sim_record_top, internal_index, sim_output_wi, pis, rand, simu_res, num_output, node_pattern_rate, max_score, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_em, T_em1, T_er, T_er1, iIndex);     
	}
		
    	ftime(&endTime);
    	double rt_step2 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    	cout << "runtime for step2 find_exdc_sim: " << rt_step2 << endl;     
	
    	if(iIndex > 0)  
    	{
    		Bnet_FreeNetwork_Bdd(net_comb, *dd_comb);
		#ifdef use_bdd
    		Cudd_Quit(*dd_comb);
		#endif
    	} 
    	if(max_score <= 0)
    	{
    		Bnet_FreeNetwork_Bdd(net, dd);
		#ifdef use_bdd
    		Cudd_Quit(dd);
		#endif
		return -1;
    	}    
    
    	//Write the whole simplified circuit into ckt_sim.blif
    	cout << "****************************" << endl;
    	cout << "step3. write_ckt_sim && sweep && map: " << endl;
    	ftime(&startTime);
    	write_ckt_sim(net, max_save_node, final_pla);
    	char com[100];
    	cout << "call SIS to sweep ckt_sim further" << endl;  	      
    	sprintf(com, "sis -t none -f ./script/sim_ckt_v2.rug > sis_ckt.txt");
    	system(com); 
    	cout <<  "Current stats by SIS: " << endl;        
    	sprintf(com, "sis -t none -f ./script/map0.rug");
    	system(com);    
    	Bnet_FreeNetwork_Bdd(net, dd);
    	update_ckt(net, dd);      
    	ftime(&endTime);
    	double rt_step3 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    	cout << "runtime for step3: " << rt_step3 << endl;  
    
    	cout << "****************************" << endl;
    	cout << "step4. Write ckt_org_sim.blif and read the combined network: " << endl;    
    	write_ckt_comb(net);
    	string filename = "./blif_files/ckt_org_sim.blif";
    	fp = fopen(filename.c_str(), "r");
    	net_comb = Bnet_ReadNetwork(fp); 
    	fclose(fp); 
    
    	cout << "****************************" << endl;
    	cout << "step5. run simulation: " << endl;    
    	write_ckt_comb(net);
    	filename = "./blif_files/ckt_org_sim.blif";
    	fp = fopen(filename.c_str(), "r");
    	ftime(&startTime); 
	cout << "running simulation for checking real error rate (em_class = 1) and obtain the error rates of local input patterns: " << endl;
	node_pattern_rate.clear();
	node_sp_simu.clear();
	internal_index.clear();
	sim_output_wi.clear();
	rand.clear();
	simu_res.clear();
	double real_er_simu = 0;
	double real_er_prev = real_er;
//	simu_ckt(net_comb, node_pattern_rate, node_sp_simu, real_er_simu, iIndex, 1);
	simu_ckt_both(net_comb, node_pattern_rate, node_sp_simu, internal_index, sim_output_wi, pis, rand, simu_res, num_output, ave_error_mag, real_er_simu, iIndex, 1);	
	cout << "real_er_simu = " << real_er_simu << endl;
	if (em_class == 1) real_er = real_er_simu;	
 	ftime(&endTime);
    	double rt_step5 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    	cout << "runtime for step5: " << rt_step5 << endl;  
	
        return 0;
}



/*2. loc_sim_main()*/
void loc_sim_main_v2(BnetNetwork *net, string &filename, double threshold_er, double threshold_em, double T_em0, double T_em1, double T_er0, double T_er1)
{
	//variables
    struct timeb st, et; 
    ini_threshold_er = threshold_er; 
    ini_threshold_em = threshold_em;
    int iIndex = 0;
    char com[100];
    string str;
    BnetNode *nd;
    
    //iterators
    set<char*>::iterator itrs;
        
    write_ckt_org(net);   //add buffers to output nodes
    sprintf(com, "cp ./blif_files/ckt_org.blif ./blif_files/ckt_org0.blif");
    system(com);
    
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
    //get ckt_org_po
    for (int i = 0; i < net->npos; i++)
    {
    	char *po = net->outputs[i];
    	string po_str(po);
    	ckt_org_po.push_back(po_str);
    }
    
    //Read sub_abs file and comparator file
    read_diff_comparator(diff_file, comparator_file, threshold_em, sub_abs_ckt, sub_abs_pi, sub_abs_po, comparator_ckt, comparator_pi, comparator_po, comp_number);         

        //start the iteration of local simplification
    	double real_er = 0, real_em = 0;
    	int min_modified_po = 1000, max_modified_po = -1;
 	ave_error_max = threshold_em;
	max_weight_max = net->npos - 1;
    	real_er_max = 0.005;
    	ave_error_max_cur = 0;
    	real_er_max_cur = 0;
	max_weight_max_cur = 0;
    	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    	cout << "i" << iIndex << ". iterations start: " << endl;
    	cout << "threshold_er = " << threshold_er << endl;
    	cout << "threshold_em = " << threshold_em << endl;
	BnetNetwork *net_comb;
	DdManager *dd_comb;
	map<string, struct score_pla> sim_record;
	map<string, struct score_pla>::iterator itrm_ss;
	multimap<double, struct score_pla> sim_record_top;
	map<string, map<string, double> > node_pattern_rate;
	map<string, int> internal_index;
	map<string, struct wi_pair> sim_output_wi;
	struct po_index_set pis = {-1, -1, -1, -1};
	vector<string> rand, simu_res;
	int num_output;
    	ftime(&st);  
    	sim_record_top.clear();
	double T_em, T_er;
	if (iIndex == 0)
	{
		T_er = T_er0;
		T_em = T_em0;
	}
    	int res = loc_sim_ite(net, net_comb, &dd_comb, sim_record, sim_record_top, node_pattern_rate, internal_index, sim_output_wi, pis, rand, simu_res, num_output, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_er, T_er1, T_em, T_em1, iIndex);
    	ftime(&et);
    	double rt_loc_sim_ite = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    	cout << "i" << iIndex << ". runtime for loc_sim_ite: " << rt_loc_sim_ite << endl;

	if(res == -1)
	{
		cout << "max_save = 0. No more optimization!" << endl;
		return;
        }
    	cout << endl << "##############################################" << endl;
    	cout << "report: " << endl;
    	cout << "real_er = " << real_er << endl;  
    	if (em_class == 2)
    		cout << "real_em = " << real_em << endl;    
	
    	iIndex++;
    	threshold_er =  ini_threshold_er -  real_er;
    		
    while(real_er <= ini_threshold_er && real_em <= ini_threshold_em)
    {
        cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    	cout << "i" << iIndex << ". iterations start: " << endl;
    	cout << "threshold_er = " << threshold_er << endl;
    	cout << "threshold_em = " << threshold_em << endl;
        ftime(&st);  
	cout << "T_em = " << T_em << endl;
    	sim_record_top.clear();
        res = loc_sim_ite(net, net_comb, &dd_comb, sim_record, sim_record_top, node_pattern_rate, internal_index, sim_output_wi, pis, rand, simu_res, num_output, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_er, T_er1, T_em, T_em1, iIndex);
	ftime(&et);
    	rt_loc_sim_ite = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
    	cout << "i" << iIndex << ". runtime for loc_sim_ite: " << rt_loc_sim_ite << endl;
    	    
    	if(res == -1)
	{
		cout << "max_save = 0. No more optimization!" << endl;
		cout << "final real_er = " << real_er << endl;    
		if (em_class == 2)
			cout << "final real_em = " << real_em << endl;    
		return;
	}	
	       
        cout << endl << "##############################################" << endl;
	cout << "report: " << endl;
	cout << "real_er = " << real_er << endl;    
	if (em_class == 2)
		cout << "real_em = " << real_em << endl;    
	                        	    
	iIndex++;
	threshold_er = ini_threshold_er -  real_er;
    }
        
    return;

}
