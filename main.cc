#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include "cudd/bnet.h"
#include "head/red_same.h"
#include "head/loc_sim_main_v2.h"
#include "head/basics.h"
#include "cudd/cudd_build_v2.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"

using namespace std;

//global variables
int numPI_ini, numPO_ini;
double ini_threshold_er, ini_threshold_em;
int num_one_po_equal = 0, num_one_po_unequal = 0;
string diff_file, comparator_file;
int em_class;
double T_em0;



int main(int argc, char **argv)
{
    struct timeb startTime, endTime;                         //Record the computing time.
    
    if(argc < 7)
    {
        cout << "Correct usage: ./main input_file em_class threshold_er threshold_em T_em0 T_em1 T_er T_er1 > output_file " << endl;
        exit(1);
    }

	ftime(&startTime);

    //*************************************//
    //*****Initialize the circuit**********//
    //*************************************//
    cout << endl << "**************************************************************" << endl;
    cout << "Read the Boolean network: " << endl;
    char com[100];
    sprintf(com, "rm -rf ./blif_files/ckt_sim.blif");
    system(com);
    ftime(&startTime);    
    string filename = argv[1]; 
    FILE *fp;
    fp = fopen(filename.c_str(), "r");
    BnetNetwork *net = Bnet_ReadNetwork(fp);
//    Bnet_PrintNetwork(net);
    numPI_ini = net->npis;
    numPO_ini = net->npos;
    fclose(fp);
    ftime(&endTime);
    double runtime_init = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "runtime for initialize: " << runtime_init << endl;

//	red_same(net);

//    return 0;

	std::size_t pos = filename.find("cla");
    if (pos != std::string::npos) 
    {
		diff_file = "./files/sub_abs/sub_abs_17.blif";
		comparator_file = "./files/comparator/comparator_17.blif";
    }
	else 
	{
		pos = filename.find("ksa");	
		if (pos != std::string::npos) 
		{
			diff_file = "./files/sub_abs/sub_abs_17.blif";
			comparator_file = "./files/comparator/comparator_17.blif";
		}
		else 
		{
			pos = filename.find("rca");	
			if (pos != std::string::npos) 
			{
				diff_file = "./files/sub_abs/sub_abs_33.blif";
				comparator_file = "./files/comparator/comparator_33.blif";
			}
			else
			{
				pos = filename.find("wal");
				if (pos != std::string::npos)
				{
					diff_file = "./files/sub_abs/sub_abs_mtp.blif";
					comparator_file = "./files/comparator/comparator_mtp.blif";
				}
				else
				{
					pos = filename.find("sin");	
					if (pos != std::string::npos)
					{
						diff_file = "./files/sub_abs/sub_abs_sin.blif";
						comparator_file = "./files/comparator/comparator_sin.blif";
					}
				}
			}
		}
	}
	


    
    if (net == NULL) {
        cout << "Syntax error in " << filename << endl;
		exit(2);
    }
    
    em_class = atoi(argv[2]);
    double threshold_er = atof(argv[3]);
    ini_threshold_er = threshold_er;
    double threshold_em = atof(argv[4]);
    ini_threshold_em = threshold_em;
    double cur_max_em = 0;
//    double T_em0 = atof(argv[5]);
    T_em0 = atof(argv[5]);
    double T_em1 = atof(argv[6]);
//    double T_er0 = atof(argv[7]);
//    double T_er1 = atof(argv[8]);
    double T_er0 = 1, T_er1 = 1;
    loc_sim_main_v2(net, filename, threshold_er, threshold_em, T_em0, T_em1, T_er0, T_er1);
    
    cout << "call SIS to sweep ckt_sim further" << endl;  	      
    sprintf(com, "sis -t none -f ./script/sim_ckt_v2.rug > sis_ckt.txt");  //sweep
    system(com); 
    
    cout << endl << "##############################################" << endl;
	cout <<  "Final result: " << endl;
	sprintf(com, "sis -t none -f ./script/map.rug");
	system(com);
	
	ftime(&endTime);
	
	double rt_locsim = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
	cout << "total runtime for locsim: " << rt_locsim << endl;
	
	cout << "num_one_po_equal = " << num_one_po_equal << endl;
	cout << "num_one_po_unequal = " << num_one_po_unequal << endl;
    
    return 0;
}
