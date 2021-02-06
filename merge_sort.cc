#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <cassert>
#include "head/helper.h"
#include "head/graph.h"
#include "head/edge.h"
#include "head/node.h"
#include "head/queue.h"


void Graph::merge(vector<Node*> &ini_list, int left, int mid, int right)
{
 //   cout << "Coming into merge!" << endl;
 //   cout << "left = " << left << ", mid = " << mid << ", right = " << right << endl;
    vector<Node*>::iterator itrv1, itrv2;
    
    vector<Node*> left_list, right_list;
 //   cout << "left_list: " << endl;
    for(int i = left; i <= mid; i++)
    {
        left_list.push_back(ini_list[i]);
  //      cout << ini_list[i]->index << ", " << ini_list[i]->gate.p << endl;;
    }
 //   cout << "right_list: " << endl;
    for(int i = mid+1; i <= right; i++)
    {
        right_list.push_back(ini_list[i]);
 //       cout << ini_list[i]->index << ", " << ini_list[i]->gate.p << endl;
    }

    vector<Node*> sort_list;
    itrv1 = left_list.begin();
    itrv2 = right_list.begin();
    while(itrv1 != left_list.end() && itrv2 != right_list.end())     
    {
        if((*itrv1)->gate.p <= (*itrv2)->gate.p)
        {
            sort_list.push_back(*itrv1);
            itrv1++;
        }
        else
        {
            sort_list.push_back(*itrv2);
            itrv2++;
        }        
    }
    if(itrv1 == left_list.end())
       sort_list.insert(sort_list.end(), itrv2, right_list.end());
    else if(itrv2 == right_list.end())
       sort_list.insert(sort_list.end(), itrv1, left_list.end());
    
    for(int i = left; i <= right; i++)
        ini_list[i] = sort_list[i-left];
       
}

void Graph::mergeSort_helper(vector<Node*> &ini_list, int left, int right)
{
 //   cout << "left = " << left << ", right = " << right << endl;
    if(left >= right)
        return;
    int mid = (left + right)/2;
    mergeSort_helper(ini_list, left, mid);
    mergeSort_helper(ini_list, mid+1, right);
    merge(ini_list, left, mid, right);
}


double Graph::mergeSort(vector<int> &merge_sort_list)
{
    map<int, Node*>::iterator itrm_in;
    vector<Node*> ini_list;
        
    for(itrm_in = nodes.begin(); itrm_in != nodes.end(); itrm_in++)
    {        
    //    cout << "current node: " << itrm_in->second->index << ", rp = " << itrm_in->second->gate.rp;
        double mp = findmin(itrm_in->second->gate.rp);
        itrm_in->second->gate.p = mp;
   //     cout << ", p = " << itrm_in->second->gate.p << endl;
        ini_list.push_back(itrm_in->second);
    }
    
    int left = 0, right = ini_list.size()-1;
    mergeSort_helper(ini_list, left, right);
    
//    cout << "merge_sort: " << endl;
    double minp;
    for(int i = 0; i < ini_list.size(); i++)
    {
        if (i == 0)
            minp = ini_list[i]->gate.p;
        merge_sort_list.push_back(ini_list[i]->index);
//      cout << i << ", node = " << ini_list[i]->index << ", p = " << ini_list[i]->gate.p << endl;
    } 
    return minp;
}


