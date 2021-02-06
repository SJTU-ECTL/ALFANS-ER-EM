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
#include "head/exdc_factor_new_v2.h"
#include "head/helper.h"
#include "head/read_file.h"
#include "head/write_func.h"
#include "head/loc_sim_main.h"
#include "head/simu_ckt.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"
#include "cudd/cudd_build.h"
#include "cudd/cudd_comp.h"
#include "cudd/cudd_dst.h"

using namespace std;

#define M 100
#define N 10000

extern int numPI_ini;
extern double ini_threshold_er, ini_threshold_em;

/*
functions in this file:

*/

//Global variables and external variables 


//


//starting_sel: starting points selection
void starting_sel(BnetNetwork *net, char *cnode, vector<int> &start_points)
{
	
}


//algorithm 2: one-dim sampling
void create_new_sampling_point()
{

}


//algorithm 1: Gibbs sampling
void gibbs_sampling(int dim)
{
	//step 1. select an initial starting point
	vector<int> start_points;
	starting_sel(start_points);
	
/*	for (int t = 0; t < N; t++)
	{
		for (int i = 0; i < dim; i++)
		{
			//draw from the conditional pdf to create a new sampling point
			create_new_sampling_point();
		}
	}
*/

}


//importance sampling
void is_simu(BnetNetwork *net, char *cnode, string filename)
{
	gibbs_sampling(dim);
	
	//obtain_q_opt();
	
	//simu_from_q_opt();
}




