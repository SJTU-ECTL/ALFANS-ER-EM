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
#include "head/HashTable.h"
#include "head/helper.h"
#include "head/stack.h"

using namespace std;



/*int read_sis_result( )
{
//    cout << "In read_sis_result: " << endl;
    ifstream fin;
    string str;
    string filename = "sis.txt";
    fin.open(filename.c_str(), iostream::in);
    int num[2];
    int i = 0;
    cout << "area result: " << endl;
    while(getline(fin, str))
    {
    	cout << str << endl;
        if(i >= 2)
            break;
        istringstream ss(str);
        string s;
        ss >> s;
//		if(s == "Total")
//		if(s == "lits(sop)=")
		if(s == "lits(fac)=")
	    {
			while(ss >> s);
//			cout << "area: " << s << endl;
			num[i] = atoi(s.c_str());
        	i++;
		}
    }
    if(i == 1)
        num[1] = 0;
    int area_save = num[0] - num[1];
    cout << "num[0] = " << num[0] << ", num[1] = " << num[1] << endl;
    fin.close();
    return area_save; 
}   
*/


int read_espresso_result(int num_lit_org)
{
    ifstream fin;
    string str;
    string filename = "./pla_files/bigNode_sim.pla";
    fin.open(filename.c_str(), iostream::in);
    int num[2];
    int i = 0;
//    cout << "area result: " << endl;
    int num_lit_sim = 0;
    int flag_start = 0;
    while(getline(fin, str))
    {
//    	cout << str << endl;
        istringstream ss(str);
        string s, s1;
        ss >> s;
        if(s[0] == '.')
        	continue;
	    ss >> s1;
	    if(s1 == "2")
	    	break;
		for(int i = 0; i < s.size(); i++)
		{
			if(s[i] != '-')
				num_lit_sim++;
		}
    }

    int area_save = num_lit_org - num_lit_sim;
    cout << "num_lit_org = " << num_lit_org << ", num_lit_sim = " << num_lit_sim << endl;
    fin.close();
    return area_save; 
}   


//read lits(fac)
int read_sis_result(int &area_save)
{
    cout << "In read_sis_result: " << endl;
    ifstream fin;
    string str;
    string filename = "sis.txt";
    fin.open(filename.c_str(), iostream::in);
    int num[2];
    int i = 0;
//    cout << "area result: " << endl;
    while(getline(fin, str))
    {
//    	cout << str << endl;
        if(i >= 2)
            break;
 /*       int flag = 0;
        istringstream ss(str);
        istringstream ss1(str);
        string s;
        int num_s = 0;
 */
        string::size_type n = 0;
		n = str.find("lits(fac)=");
		if(n != string::npos)
		{
			string rev_num_str;
			int len = str.size();
			for(int j = len-1; j >= 0; j--)
			{
				if(isdigit(str[j]))
					rev_num_str.append(1, str[j]);
				if(str[j] == '=')
					break;
			}
	//		cout << "rev_num_str = " << rev_num_str << endl;
			string num_str;
			for(int j = rev_num_str.size()-1; j >= 0; j--)
				num_str.append(1, rev_num_str[j]);
	//		cout << "num_str = " << num_str << endl;
			num[i] = atoi(num_str.c_str());
			i++;
		}
/*        while(ss >> s)
        	num_s++;
        if(num_s == 4)
	    {
	        while(ss1 >> s)
				if(s == "lits(fac)=")
			    {
			    	flag = 1;
					break;				
				}
			if(flag == 1)
			{
				ss1 >> s;
				num[i] = atoi(s.c_str());
			    i++;
			}
		}
		else if(num_s == 2)
		{
			while(ss1 >> s);
			string::size_type n = 0;
			n = s.find("lits(fac)=");
			int size_num = s.size() - 10; 
			string num_sim_lit = s.substr(n+10, size_num);
			num[i] = atoi(num_sim_lit.c_str());
			i++;
		}
*/
    }
    if(i == 1)
        num[1] = 0;
    area_save = num[0] - num[1];
    cout << "num[0] = " << num[0] << ", num[1] = " << num[1] << endl;
    fin.close();
    
    if(num[1] == 0)
    	return 1;
    else
    	return 0; 
}    
   
   
//read lits(fac)
int read_num_lit()
{
    ifstream fin;
    string str, s;
    string filename = "./output/lit.txt";
    fin.open(filename.c_str(), iostream::in);
    int num;
//    cout << "lit result: " << endl;
    while(getline(fin, str))
    {
//    	cout << str << endl;
    	if(str.find("lits(fac)=") == string::npos)
    		continue;
        int flag = 0;
        istringstream ss(str);
        istringstream ss1(str);
        int num_s = 0;
        while(ss >> s)
        	num_s++;
        if(num_s == 4)
	    {
	        while(ss1 >> s)
				if(s == "lits(fac)=")
			    {
			    	flag = 1;
					break;				
				}
			if(flag == 1)
			{
				ss1 >> s;
				num = atoi(s.c_str());
			}
		}
		else if(num_s == 2)
		{
			while(ss1 >> s);
			string::size_type n = 0;
			n = s.find("lits(fac)=");
			int size_num = s.size() - 10; 
			string num_sim_lit = s.substr(n+10, size_num);
			num = atoi(num_sim_lit.c_str());
		}
    }
	
	return num;

}
   
   
   

/*read_mvsis_result*/
void read_mvsis_result(vector<string> &dont_care,  vector<string> &insig_string, string &filename)
{
//	cout << endl << "Coming into read_mvsis_result!" << endl;
	ifstream fin;
    string str;
//    string filename = "mvsis.txt";
    fin.open(filename.c_str(), iostream::in);
	vector<string> varVerStr, varHorStr;
	int nVarsVer = 0, nVarsHor = 0;
	int flag_start = 0, first_doll = 0;
	while(getline(fin, str))
	{
//		cout << str << endl;
		if(str.empty())
			continue;
//		cout << "flag_start = " << flag_start << endl;
//		cout << "str = " << str << endl;
		istringstream ss(str);
		string s;
		ss >> s;
		if(s == "$$")
		{		
			first_doll++;
			if(first_doll == 1)
			{				
				int flag = 0;
				while(ss >> s)
				{
					if(flag == 0)
					{	
						if(s != "\\")
							varVerStr.push_back(s);
						else
							flag = 1;
					}
					else
						varHorStr.push_back(s);
				}
				nVarsVer = varVerStr.size();
				nVarsHor = varHorStr.size();
			}		
		}
		
		if(first_doll == 1)
		{
			string::size_type position = str.find('|');
			if(position == string::npos)
				continue;
			else
			{
				//cout << "| is here!" << endl;
				string ver = s;
				vector<int> indexHor;
				int ind = 0;
				while(ss >> s)
				{
					if(s == "1")
						indexHor.push_back(ind);
					if(s == "0" || s == "1" || s == "-")
						ind++;
				}
				for(int i = 0; i < indexHor.size(); i++)
				{
					string hor = int2grey(indexHor[i], nVarsHor);
					string mint = ver;
					mint.append(hor);
				}
			}
		}				
		else if(first_doll == 2)
		{
			string::size_type position = str.find('|');
			//cout << "position = " << position << endl;
			if(position == string::npos)
				continue;
			else
			{
				//cout << "| is here!" << endl;
				string ver = s;
				//cout << "ver = " << ver << endl;
				vector<int> indexHor;
				int ind = 0;
				while(ss >> s)
				{
					if(s == "-")
						indexHor.push_back(ind);
					if(s== "0" || s == "1" || s == "-")
						ind++;
				}
				for(int i = 0; i < indexHor.size(); i++)
				{
					string hor = int2grey(indexHor[i], nVarsHor);
					string exdc = ver;
					exdc.append(hor);
					dont_care.push_back(exdc);
//					cout << "exdc = " << exdc << endl;
				}
			}
		}
	}
	fin.close();

	insig_string = varVerStr;
	vector<string>::iterator it = insig_string.end();
	insig_string.insert(it, varHorStr.begin(), varHorStr.end());
/*	cout << "insig_string: " << endl;
	for(int i = 0; i < insig_string.size(); i++)
		cout << insig_string[i] << " ";
	cout << endl;	
*/
}


void read_abc_result(vector<string> &dont_care, string &filename)
{
	ifstream fin;
    string str, s;
    fin.open(filename.c_str(), iostream::in);
	int flag_start = 0;
	while(getline(fin, str))
	{
		istringstream ss(str);
		ss >> s;
		if (s == ".names")
		{
			flag_start = 1;
			continue;
		}
		if (flag_start)
		{
			if (str[0] != '0' && str[0] != '1') break;
			dont_care.push_back(str);
		}
	}
	fin.close();
}


int read_abc_level(string &filename)
{
	int level;
	ifstream fin;
	fin.open(filename.c_str(), ios::in);
	string str, s;
	while (getline(fin, str))
	{
		istringstream ss(str);
		ss >> s;
		if (s == "Level")
		{
			ss >> s;
			ss >> s;
			stringstream tmp(s);
			tmp >> level;
			break;
		}
	}	
	fin.close();
	return level;
}



void read_factor(vector<string> &lit_unit)
{
	ifstream fin;
	fin.open("factor.txt", ios::in);
	string str, s;
	
	string lit_set;
	int start_store = 0;
	int flag_start = 0;
	int line = 0;
	int flag_quit = 0;
	stack factor_stack(200);
	cout << "factor.txt: " << endl;
	while(getline(fin, str))
	{
		istringstream ss(str);
		cout << "str = " << str << endl;
		if(factor_stack.empty())
			lit_set.clear();
		if(line == 0)
			flag_start = 0;
		else
			flag_start = 1;
		while(ss >> s)
		{		
			if(s == ".model")
			{
				flag_quit = 1;
				break;
			}	
			if(s == " ")
				continue;
			if(line == 0 && s == "=")
				flag_start = 1;			
			if(flag_start)
			{
//				cout << "s: " << s << endl;
				if(s == "+" && factor_stack.empty())
				{					
//					cout << "lit_set: " << lit_set << endl;	
					if(!lit_set.empty())				
						lit_unit.push_back(lit_set);
					lit_set.clear();
				}
				else
				{
					string lit = s;
					if(lit[0] == '(')
						factor_stack.push(1);
					else if(lit[lit.size()-1] == ')')
						int right = factor_stack.pop();
					if(s == "=")
						continue;
					lit_set.append(s);
					lit_set.append(" ");
//					cout << "lit_set: " << lit_set << endl;	
				}
			}
		}
		if(!lit_set.empty() && factor_stack.empty())	
			lit_unit.push_back(lit_set);
		line++;
		if(flag_quit)
			break;
	}
	
	cout << "lit_unit: " << endl;
	for(int i = 0; i < lit_unit.size(); i++)
		cout << lit_unit[i] << endl;
	cout << endl;
}	



void read_factor_v2(string &factor)
{
	ifstream fin;
	fin.open("factor.txt", ios::in);
	string str, s;
	
	int flag_start = 0;
	int line = 0;
	int flag_quit = 0;
	while(getline(fin, str))
	{
		istringstream ss(str);
		if(line == 0)
			flag_start = 0;
		else
			flag_start = 1;
		while(ss >> s)
		{		
			if(s == ".model")
			{
				flag_quit = 1;
				break;
			}	
			if(s == " ")
				continue;
			if(line == 0 && s == "=")
				flag_start = 1;			
			if(flag_start)
			{
				if(s == "=")
					continue;
				factor.append(s);
				factor.append(" ");
			}
		}
		line++;
		if(flag_quit)
			break;
	}
	
	cout << "factor = " << factor << endl;
}


void read_diff_comparator(string &diff_file, string &comparator_file, int threshold_em, vector<string> &sub_abs_ckt, vector<string> &sub_abs_pi, vector<string> &sub_abs_po, vector<string> &comparator_ckt, vector<string> &comparator_pi, vector<string> &comparator_po, vector<int> &comp_number) 
{
	/* step1. read from the sub_abs file to get sub_abs_ckt and sub_abs_pi, sub_abs_po */
	//read diff file and its POs
	ifstream fin;
 	fin.open(diff_file.c_str(), ios::in);
 	string newline;
 	int start = 0;
        string str, s;
 	while(getline(fin, str))
 	{
 		istringstream ss(str);
 		ss >> s;
 		if(s == ".end")
 			break;
 		if(s == ".names" && start == 0) start = 1;
		newline.clear();
 		if(start)
 		{
 			if (s == ".names")
	 		{
	 			newline.append(s);
	 			newline.append(" ");
	 			while (ss >> s)
	 			{
	 				s.append("_sub ");
	 				newline.append(s);
	 			}
	 		}
	 		else
	 			newline = str;
 		}	
 		sub_abs_ckt.push_back(newline);
 	}
 	fin.close();

        cout << "diff_file: " << diff_file << endl;
 	FILE *fp = fopen(diff_file.c_str(), "r");
	if (fp == NULL)
		cout << "opening diff_file fails!" << endl;
    BnetNetwork *net_sub = Bnet_ReadNetwork(fp);
    for (int i = 0; i < net_sub->npis; i++)
    {
    	char *pi = net_sub->inputs[i];
    	string pi_str(pi);
    	pi_str.append("_sub");
    	sub_abs_pi.push_back(pi_str);
    }
    for (int i = 0; i < net_sub->npos; i++)
    {
    	char *po = net_sub->outputs[i];
    	string po_str(po);
    	po_str.append("_sub");
    	sub_abs_po.push_back(po_str);
    }
    fclose(fp);
    
    
	/* step2. read from the comparator file to get comparator_ckt and comparator_pi, comparator_po */
	//read diff file and its POs
 	fin.open(comparator_file.c_str(), ios::in);
 	start = 0;
 	while(getline(fin, str))
 	{
 		istringstream ss(str);
 		ss >> s;
 		if(s == ".end")
 			break;
 		if(s == ".names" && start == 0) start = 1;
		newline.clear();
 		if(start)
 		{
 			if (s == ".names")
	 		{
	 			newline.append(s);
	 			newline.append(" ");
	 			while (ss >> s)
	 			{
	 				s.append("_comp ");
	 				newline.append(s);
	 			}
	 		}
	 		else
	 			newline = str;
 		}	
 		comparator_ckt.push_back(newline);
 	}
 	fin.close();
 	fp = fopen(comparator_file.c_str(), "r");
    BnetNetwork *net_comparator = Bnet_ReadNetwork(fp);
    for (int i = 0; i < net_comparator->npis; i++)
    {
    	char *pi = net_comparator->inputs[i];
    	string pi_str(pi);
    	pi_str.append("_comp");
    	comparator_pi.push_back(pi_str);
    }
    for (int i = 0; i < net_comparator->npos; i++)
    {
    	char *po = net_comparator->outputs[i];
    	string po_str(po);
    	po_str.append("_comp");
    	comparator_po.push_back(po_str);
    }
    fclose(fp);
    
    
	/* step3. get comp_number from threshold_em */
	int numBit = comparator_pi.size()/2;
	int2bin_reverse(threshold_em, comp_number, numBit);
}
