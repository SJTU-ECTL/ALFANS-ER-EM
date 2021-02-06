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
#include "head/graph.h"
#include "head/edge.h"
#include "head/node.h"
#include "head/queue.h"
#include "head/HashTable.h"
#include "head/call_abc.h"
#include "head/helper.h"
#include "head/read_file.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/cudd_build.h"
#include "cudd/cudd_comp.h"
#include "cudd/cudd_dst.h"

using namespace std;

extern int numPI_ini;

/*
functions in this file:
*/

//Global variables and external variables 


double comp_real_er(BnetNetwork *net, BnetNetwork *net_sim, DdManager **dd)
{
	BnetNode *node1, *node2;
    char *oname;
    DdNode *odd1, *odd2, *tmp, *func, *x, *y;
    int i;
    int pr;
    
	func = Cudd_ReadLogicZero(*dd);
	Cudd_Ref(func);
	
	for(int i = 0; i < net->noutputs; i++) 
	{
		oname = net->outputs[i];
		printf("output: %s\n", oname);
		if (!st_lookup(net->hash, oname, &node1)) 
		{
			cout << "error: output node " << oname << " is not in hash table!" << endl;
		    exit(1);
		}
		odd1 = node1->dd;
		if(odd1 == NULL)
		{
			cout << "odd1 = NULL!" << endl;
			exit(1);
		}
		if (!st_lookup(net_sim->hash, oname, &node2)) 
		{
	    	cout << "error: output node " << oname << " is not in net_sim!" << endl;
		    exit(1);
		} 
	    odd2 = node2->dd;
	    if(odd2 == NULL)
	    {
			cout << "odd2 = NULL!" << endl;
			exit(1);
		}
		
		
		double num_odd1 = Cudd_CountMinterm(*dd, odd1, numPI_ini);
    	cout << "num_odd1 = " << num_odd1 << endl;
    	Cudd_PrintDebug(*dd, odd1, numPI_ini, 2);
    	
    	double num_odd2 = Cudd_CountMinterm(*dd, odd2, numPI_ini);
    	cout << "num_odd2 = " << num_odd2 << endl;
    	Cudd_PrintDebug(*dd, odd2, numPI_ini, 2);
		
		DdNode *inter = Cudd_bddXor(*dd, odd1, odd2);
//		Cudd_PrintDebug(*dd, inter, numPI_ini, 2);
		Cudd_Ref(inter);		
    	
		tmp = Cudd_bddOr(*dd, func, inter);
//		Cudd_PrintDebug(*dd, tmp, numPI_ini, 2);
		Cudd_Ref(tmp);
		
		Cudd_IterDerefBdd(*dd, func);
		Cudd_IterDerefBdd(*dd, inter);
		func = tmp;
    }
    
    double num_diff = Cudd_CountMinterm(*dd, func, numPI_ini);
    cout << "In comp_real_er, num_diff = " << num_diff << endl;
    double real_er = num_diff/pow(2.0, numPI_ini);
    
    Cudd_IterDerefBdd(*dd, func);

    return real_er;
}

