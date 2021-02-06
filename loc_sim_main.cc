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
#include "head/sim_new_abc.h"
#include "head/simu_ckt.h"
#include "head/loc_sim_main.h"
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

int flag_round = 0;
double threshold_org_er;
double lit_weight = 0.1, er_weight = 0.5, em_weight = 0.5;
double ini_ratio = em_weight/er_weight;
double ave_error_max, ave_error_max_cur, max_weight_max, max_weight_max_cur, real_er_max, real_er_max_cur;
int lit_save_max, lit_save_max_cur;

//external variables
extern int em_class; //1: max, 2: ave
extern int numPI_ini, numPO_ini;
extern double ini_threshold_er, ini_threshold_em;
extern string diff_file, comparator_file;


/*
functions in this file:

*/




/*1. loc_sim_ite*/
int loc_sim_ite(BnetNetwork *net, BnetNetwork *&net_comb, DdManager **dd_comb, map<string, struct score_pla> &sim_record, multimap<double, struct score_pla> &sim_record_top, map<string, map<string, double> > &node_pattern_rate, map<string, int> &internal_index, map<string, struct wi_pair> &sim_output_wi, struct po_index_set &pis, vector<string> &rand, vector<string> &simu_res, int &num_output, double threshold_er, double threshold_em, double &real_er, double &real_em, int &min_modified_po, int &max_modified_po, double &T_er, double T_em0, double T_em1, int &iIndex)
{
    //variables
    struct timeb startTime, endTime; 
    BnetNode *nd, *nd1;
    multimap<double, char*> ave_sp;
    FILE *fp;
    DdManager *dd = NULL; 
    cout << "at the start of loc_sim_ite, real_er = " << real_er << endl;    

    //iterators
    multimap<double, char*>::iterator itrm_dc;    
    set<char*>::iterator itrs;
    map<string, struct wi_pair>::iterator itrm_sw;
    
    //write circuits for transitive fanin cone of each primary output
/*    cout <<"step0. write circuits for transitive fanin cone: " << endl;
    po_tfi.clear();
    for(int i = 0; i < net->npos; i++)
    {
    	nd = net->nodes;
    	while(nd != NULL)
    	{
    		nd->visited = 0;
    		nd = nd->next;
    	}
    	char *po = net->outputs[i];
    	st_lookup(net->hash, po, &nd);
    	set<char*> tfi;
    	find_tranfanin(net, po, tfi);
    	string po_str(po);
    	string fn = "./blif_files/cone_po_";    	
    	fn.append(po_str);
    	fn.append(".blif");
    	char com[100];
    	sprintf(com, "rm -rf %s", fn.c_str());
    	system(com);
    	//update po_tfi
	po_tfi.insert(make_pair(po_str, tfi));
    	//write the cone circuit for current po
    	set<string> inputs_set;
    	vector<string> cone_string;
    	write_cone_circuit(net, iIndex, po, tfi, fn, inputs_set, cone_string);
		//if iIndex == 0, initialize po_inputs and po_cone_string
    	if (iIndex == 0)
    	{
    		po_inputs.insert(make_pair(po_str, inputs_set));
    		po_cone_string.insert(make_pair(po_str, cone_string));
    	}
    }
*/
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
	
#ifdef use_simu
	if(iIndex == 0)
	{
		cout << "running simulation for the first iteration: " << endl;
//		simu_ckt(net, node_pattern_rate, node_sp_simu, real_er, iIndex, 0);
		simu_ckt_both(net, node_pattern_rate, node_sp_simu, internal_index, sim_output_wi, pis, rand, simu_res, num_output, ave_error_mag, real_er, iIndex, 0);
	/*	cout << "0. sim_output_wi: " << endl;
		for(itrm_sw = sim_output_wi.begin(); itrm_sw != sim_output_wi.end(); itrm_sw++)
			cout << itrm_sw->first << ": " << itrm_sw->second.weight << ", " << itrm_sw->second.index << endl;	
	*/
	}
#endif
    
    	//Find exdcs for each big node and simplify them using these exdcs. Obtain 
    	//the one with maximum node save
    	cout << "****************************" << endl;
    	cout << "step2. find_exdc_sim: " << endl;	
    	ftime(&startTime);    
    	vector<string> final_pla;  //store the files for the simplifed big node
    	double max_score = -1;
    	cout << "threshold_em = " << threshold_em << endl;
//      if (T > 7 && iIndex > 0 && iIndex % 3 == 0)
//		T = T/4;
//    	if (iIndex > 0 && T > 67)
//    		T = T * 0.9;
	double T_em;
	if (iIndex == 0) T_em = T_em0;
/*	else
	{
		cout << "real_er = " << real_er << ", threshold_er = " << threshold_er << ", T_em0 = " << T_em0 << endl;
	//	T_em = T_em0 * pow((1 - real_er/ini_threshold_er), pow(iIndex, 0.5));
		if (real_er <= 0.5 * ini_threshold_er) T_em = T_em * pow(0.8, pow(iIndex, 0.5));
		else  T_em = T_em * pow(0.4, pow(iIndex, 0.5));
		if (T_em < T_em1) T_em = T_em1;
	}
*/

    	cout << "iteration = " << iIndex << ", T_er = " << T_er << ", T_em = " << T_em << endl;
	if (em_class == 2)
    	{
		if (iIndex > 0)
		{
		/*
			double margin_em = 1-real_em/ini_threshold_em;
			double margin_er = 1-real_er/ini_threshold_er;
			double margin_ratio = margin_er/margin_em;      	
			cout << "margin_em = " << margin_em << ", margin_er = " << margin_er << ", margin_ratio = " << margin_ratio << endl;
			double new_ratio = ini_ratio * margin_ratio;
			cout << "new_ratio = " << new_ratio << endl;
			er_weight = 1/(1 + new_ratio);
			em_weight = 1 - er_weight;
		*/
		}
	/*	er_weight = 0.8 * pow(0.95, sqrt(iIndex));
		if (er_weight < 0.7) er_weight = 0.7;
		em_weight = 1 - er_weight;
	*/
    		cout << "iteration = " << iIndex << ", er_weight = " << er_weight << ", em_weight = " << em_weight << endl;
    	}
    	char *max_save_node = find_exdc_sim(net, &dd, net_comb, dd_comb, final_pla, sim_record, sim_record_top, internal_index, sim_output_wi, pis, rand, simu_res, num_output, node_pattern_rate, max_score, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_er, T_em, T_em1, iIndex);     
	cout << "after find_exdc_sim, real_er = " << real_er << endl;
    	ftime(&endTime);
    	double rt_step2 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    	cout << "runtime for find_exdc_sim: " << rt_step2 << endl;     
	
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
//    	sprintf(com, "sis -t none -f ./script/sim_ckt_v2.rug > sis_ckt.txt");  //sweep
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
    
    	ftime(&startTime); 
    
    	//Build bdd to compute the real error rate and the real average error magnitude	        
/*  	DdManager *dd_final = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);    
	cudd_build_v2(net_comb, dd_final, filename.c_str(), BNET_GLOBAL_DD);
	char *out_node = net_comb->outputs[0];
	st_lookup(net_comb->hash, out_node, &nd);
	double num_diff = Cudd_CountMinterm(*dd_comb, nd->dd, numPI_ini);
	double real_er_bdd = num_diff/pow(2.0, numPI_ini);
	cout << "real_er_bdd = " << real_er_bdd << endl;
	real_er = real_er_bdd;
	//compute the ave_error_mag from the BDD	
	real_em = ave_error_mag;
*/	
	
#ifdef use_simu	
	cout << "running simulation for checking real error rate: " << endl;
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
 #endif   

 
 	ftime(&endTime);
    	double rt_step5 = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    	cout << "runtime for step5: " << rt_step5 << endl;  
	
        return 0;
}



/*2. loc_sim_main()*/
void loc_sim_main(BnetNetwork *net, string &filename, double threshold_er, double threshold_em, double T_em0, double T_em1)
{
	//variables
    struct timeb st, et; 
    threshold_org_er = threshold_er;
    double threshold_org_em = threshold_em;
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
    double T_er = 0.2;
    double T_em = T_em0;
 	ave_error_max = threshold_em;
	max_weight_max = net->npos - 1;
    	real_er_max = 0.005;
        lit_save_max = 1;
    	ave_error_max_cur = 0;
    	real_er_max_cur = 0;
	max_weight_max_cur = 0;
	lit_save_max_cur = 0;
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
    int res = loc_sim_ite(net, net_comb, &dd_comb, sim_record, sim_record_top, node_pattern_rate, internal_index, sim_output_wi, pis, rand, simu_res, num_output, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_er, T_em0, T_em1, iIndex);
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
    threshold_er =  threshold_org_er -  real_er;
//    threshold_em =  threshold_org_em -  real_em;
    		
    while(real_er <= threshold_org_er && real_em <= threshold_org_em)
    {
        cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    	cout << "i" << iIndex << ". iterations start: " << endl;
    	cout << "threshold_er = " << threshold_er << endl;
    	cout << "threshold_em = " << threshold_em << endl;
        ftime(&st);  
	cout << "T_em = " << T_em << endl;
    	sim_record_top.clear();
        res = loc_sim_ite(net, net_comb, &dd_comb, sim_record, sim_record_top, node_pattern_rate, internal_index, sim_output_wi, pis, rand, simu_res, num_output, threshold_er, threshold_em, real_er, real_em, min_modified_po, max_modified_po, T_er, T_em0, T_em1, iIndex);
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
	    threshold_er = threshold_org_er -  real_er;
   // 	    threshold_em =  threshold_org_em - real_em;
    }
        
    return;

}
