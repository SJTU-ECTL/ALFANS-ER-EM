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
#include "head/stack.h"
#include "head/basics.h"
#include "head/helper.h"
#include "head/read_file.h"
#include "head/write_func.h"
#include "head/exdc_helper.h"
#include "head/exdc_factor.h"
#include "head/btree.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cudd.h"
#include "/home/wuyi/usr/CUDD/cudd-2.5.0/cudd/cuddInt.h"

using namespace std;

extern int numPI_ini;

/*
functions in this file:

*/

//Global variables and external variables 

btNode *CreateNode(const string& x)
{
   btNode * p = new btNode;
   p->parent = p->left = p->right = NULL;
   p->data = x;
   p->numPT = 0;
   return p;
}

bool IsOperator(const string& x)
{
   // Since the only impact of parentheses () is on precedence, 
   // they are not considered as operators here
   return ((x.length() == 1) && (x[0] == '*' || x[0] == '+'));

}

bool IsLeftParenthesis(const string& x)
{
   return x == "(";
}

bool IsRightParenthesis(const string& x)
{
   return x == ")";
}

bool IsOperand(const string& x)
{
/*   int y;
   stringstream ss(x);
   if (ss >> y) return true;
   else return false;
*/
	if(x != "*" && x != "+" && x != "(" && x != ")")
		return true;
	return false;
}

int GetPrecedence(const string& x)
{
   assert(IsOperator(x));
   if (x[0] == '*') return 2;
   else return 1;
}

btNode * CreateInfixTree(const string& exp)
{
   // create a dummy root with minimal precedence
   // its content is trivial
   btNode * root = CreateNode("0");
   root->precedence = INT_MIN;

   // the previous operand of current operator
   btNode * preOperand = NULL;
   // the previous operator of current operator
   btNode * preOperator = root;
   // the impact of preceding parenthesis, if any
   int correction = 0;

   string token;
   stringstream ss(exp);

   while (ss >> token)
   {
      if (IsOperand(token))
      {
         preOperand = CreateNode(token);
      }
      else if (IsOperator(token))
      {
         btNode * p = CreateNode(token);
         p->precedence = GetPrecedence(token) + correction;
         if (p->precedence > preOperator->precedence)
         {
            p->left = preOperand;
            preOperator->right = p;
            p->parent = preOperator;
         }
         else
         {
            preOperator->right = preOperand;
            btNode * q = preOperator->parent;
            while (p->precedence <= q->precedence) q = q->parent;

            p->left = q->right;
            q->right = p;
            p->parent = q;
         }
         preOperand = NULL;
         preOperator = p;        
      }//else if (IsOperator(token)
      else if (IsLeftParenthesis(token))
      {
         correction += 2;
      }
      else if (IsRightParenthesis(token))
      {
         correction -= 2;
      }
      else
      {
         cout << "illegal token found: " << token << endl;
         break;
      }
   }//while

   if (preOperand == NULL)
       cout << "illegal expression: cannot end with operator: "
            << preOperator->data << endl;
   else preOperator->right = preOperand;

   // delete dummy root
   btNode * realRoot = root->right;
   delete root;
   if (realRoot) realRoot->parent = NULL;
   

   return realRoot;
   
   
}

void PostOrderPrintTree(btNode * node)
{
   if (node)
   {
      PostOrderPrintTree(node->left);
      PostOrderPrintTree(node->right);
      cout << node->data << " ";
   }
}


void InOrderPrintTree(btNode * node)
{
   if (node)
   {
   	 
   	  btNode *left_child = node->left;
   	  if(left_child != NULL)
   	  	left_child->parent = node;
      InOrderPrintTree(node->left);
      
//      cout << node->data << " ";
      
      btNode *right_child = node->right;
   	  if(right_child != NULL)
   	  	right_child->parent = node;
      InOrderPrintTree(node->right);
      
   }
}

void copyTree(btNode *root, btNode *new_root)
{
	if(root)
	{
/*		cout << endl << "node: " << root->data << endl;
		if(root->left != NULL)
			cout << "left child: " << root->left->data << endl;
		if(root->right != NULL)
			cout << "right child: " << root->right->data << endl;
*/			
		new_root->data = root->data;
		new_root->precedence = root->precedence;
		new_root->numPT = root->numPT;
		new_root->ind = root->ind;
		new_root->exp = root->exp;
	
		if(root->left == NULL)
			new_root->left = NULL;
		else
		{
			btNode *left = new btNode;
			copyTree(root->left, left);
			new_root->left = left;
			left->parent = new_root;
		}
		
		if(root->right == NULL)
			new_root->right = NULL;
		else
		{
			btNode *right = new btNode;
			copyTree(root->right, right);
			new_root->right = right;
			right->parent = new_root;
		}
/*		cout << "current node: " << new_root->data << endl;
		if(new_root->left != NULL)
			cout << "left child: " << new_root->left->data << endl;
		if(new_root->right != NULL)
			cout << "right child: " << new_root->right->data << endl;
*/
	}
}


void freeTree(btNode *root)
{
	if(root)
	{		
		if(root->left == NULL && root->right == NULL)
		{
//			cout << "free node " << root->data << endl;
			btNode *parent = root->parent;
			if(parent == NULL)
			{
				delete root;
				return;
			}
			else
			{
				if(parent->left == root)
					parent->left = NULL;
				else if(parent->right == root)
					parent->right = NULL;
				delete root;
			}
		}
		else 
		{
			if(root->left != NULL)
				freeTree(root->left);	
			if(root->right != NULL)	
				freeTree(root->right);			
			freeTree(root);
		}
	}
}


void removeLeaf(btNode *root, btNode *leaf)
{
	if(leaf == root)
	{
		root = NULL;
		return;
	}
	else
	{
		btNode *other_child, *ppnode;
		btNode *pnode = leaf->parent;
		if(pnode == root)
		{
			if(pnode->left == leaf)
				other_child = pnode->right;
			else if(pnode->right == leaf)
				other_child = pnode->left;
			delete leaf;	
			root = other_child;
			root->parent = NULL;
			delete pnode;
		}
		else
		{
			if(pnode->left == leaf)
				other_child = pnode->right;
			else if(pnode->right == leaf)
				other_child = pnode->left;
			delete leaf;	
			
			ppnode = pnode->parent;
			other_child->parent = ppnode;
			if(ppnode->left == pnode)
				ppnode->left = other_child;
			else if(ppnode->right = pnode)
				ppnode->right = other_child;
			delete pnode;
		}
		return;
	}
}


int removeLeafNode(btNode **root, string &leaf_node)
{
	multimap<string, btNode*> node_leaf_set;
	multimap<string, btNode*>::iterator itrm_sb;
	
	visitleafnode(*root, node_leaf_set);
	int num = node_leaf_set.count(leaf_node);
	if(num > 1)
	{
	 	cout << "same node appears more than once!" << endl;
	 	return 1;
	}
	else
	{
		itrm_sb = node_leaf_set.find(leaf_node);
		btNode *leaf = itrm_sb->second;
	//	cout << "in removeleafnode, leaf = " << leaf->data << endl;

		if(leaf == *root)
		{
			*root = NULL;
			return 0;
		}
		else
		{
			btNode *other_child, *ppnode;
			btNode *pnode = leaf->parent;
			if(pnode == *root)
			{
				if(pnode->left == leaf)
					other_child = pnode->right;
				else if(pnode->right == leaf)
					other_child = pnode->left;
	//			cout << "other_child = " << other_child->data << endl;
				delete leaf;	
				other_child->parent = NULL;
				*root = other_child;
	//			cout << "root node: " << (*root)->data << endl;
			//	delete pnode;
			}
			else
			{
				if(pnode->left == leaf)
					other_child = pnode->right;
				else if(pnode->right == leaf)
					other_child = pnode->left;
				delete leaf;	
				
				ppnode = pnode->parent;
				other_child->parent = ppnode;
				if(ppnode->left == pnode)
					ppnode->left = other_child;
				else if(ppnode->right = pnode)
					ppnode->right = other_child;
				delete pnode;
			}
			return 0;
		}
		
	}
}


void InOrderPrintExp0(btNode * node)
{
   if (node)
   {
      InOrderPrintExp0(node->left);
      
      cout << node->data << " ";
      
      InOrderPrintExp0(node->right);
      
   }
}



void InOrderPrintExp(btNode * node)
{
   if (node)
   {
      InOrderPrintExp(node->left);
      
      cout << "node " << node->data << ": " << endl;
      
      vector<string> exp = node->exp;
      for(int i = 0; i < exp.size(); i++)
      	cout << exp[i] << endl;
      
      InOrderPrintExp(node->right);
      
   }
}


void compNumPT(btNode **node, int &num)
{
   if (*node)
   {	
   	  if((*node)->left == NULL && (*node)->left == NULL)
   	  {
   	  	num = 1;
   	  	(*node)->numPT = num;
   	  	return;
   	  }
   	  int num_left = 0;
      compNumPT(&(*node)->left, num_left);
      int num_right = 0;
      compNumPT(&(*node)->right, num_right);
      if((*node)->data == "+")  //+
      	num = num_left + num_right;
      else if((*node)->data == "*")  //*
      	num = num_left * num_right;
      (*node)->numPT = num;
   }
}


void visitleaf(btNode *root, map<int, btNode*> &leaf_set, int &ind)
{
	if(root != NULL)
		root->exp.clear();
	if(root->left != NULL)
		visitleaf(root->left, leaf_set, ind);
	if(root->right != NULL)
		visitleaf(root->right, leaf_set, ind);
	if(root->left == NULL && root->right == NULL)
	{
		leaf_set.insert(pair<int, btNode*>(ind, root));
		root->ind = ind;
		ind++;
	//	node_leaf_set.insert(pair<string, btNode*>(root->data, root));
	}
}

void visitleafnode(btNode *root, multimap<string, btNode*> &node_leaf_set)
{
	if(root->left != NULL)
		visitleafnode(root->left, node_leaf_set);
	if(root->right != NULL)
		visitleafnode(root->right, node_leaf_set);
	if(root->left == NULL && root->right == NULL)
	{
		node_leaf_set.insert(pair<string, btNode*>(root->data, root));
	}
}



void compInvCubeNum(btNode *root, btNode *leaf, int &num_cube)
{
	int num_sibling;
	btNode *pnode = leaf;
	btNode *cnode;
	
	while(1)
	{		
		cnode = pnode;
		pnode = pnode->parent;
		
		if(pnode->data == "*")
		{
			if(pnode->left == cnode)
				num_sibling = pnode->right->numPT;
			else if(pnode->right == cnode)
				num_sibling = pnode->left->numPT;		
			num_cube = num_cube * num_sibling;
		}
		
		if(pnode == root)
			break;
	}
}


void comp_exp(btNode **root)
{
	if (*root)
    {
    	if((*root)->left == NULL && (*root)->right == NULL)
    	{
    		ostringstream ss;
    		ss << (*root)->ind;	
    		string str = ss.str();
    		(*root)->exp.push_back(str);
    		return;
    	}
    	
        comp_exp(&(*root)->left);
        comp_exp(&(*root)->right);
        
        vector<string> exp_left = (*root)->left->exp;
        vector<string> exp_right = (*root)->right->exp;
        vector<string> new_exp;
        
        if((*root)->data == "*")  
        {
        	for(int i = 0; i < exp_left.size(); i++)
        		for(int j = 0; j < exp_right.size(); j++)
        		{
        			string str;
        			str.append(exp_left[i]);
        			str.append(exp_right[j]);
        			new_exp.push_back(str);
        		}
        }
        else if((*root)->data == "+")
        {
        	new_exp.insert(new_exp.end(), exp_left.begin(), exp_left.end());
        	new_exp.insert(new_exp.end(), exp_right.begin(), exp_right.end());
        }
        
        (*root)->exp = new_exp;
    }
}

void printVec(vector<string> &cubes)
{
//	cout << "size of this vector: " << cubes.size() << endl;
	for(int i = 0; i < cubes.size(); i++)
		cout << cubes[i] << endl;
}


void get_involve_cubes(btNode *root, btNode *leaf, vector<string> &inv_cubes)
{
	inv_cubes = leaf->exp;
	btNode *cnode = leaf;
	btNode *pnode;
	vector<string> exp_sibling;
	
	if(leaf == root)
	{
		inv_cubes = leaf->exp;
		return;
	}
	
	while(1)
	{
		pnode = cnode->parent;
		if(pnode->data == "*")
		{
			if(pnode->left == cnode)
				exp_sibling = pnode->right->exp;
			else if(pnode->right == cnode)
				exp_sibling = pnode->left->exp;
			
			vector<string> new_exp;
			for(int i = 0; i < exp_sibling.size(); i++)
        		for(int j = 0; j < inv_cubes.size(); j++)
        		{
        			string str;
        			str.append(exp_sibling[i]);
        			str.append(inv_cubes[j]);
        			new_exp.push_back(str);
        		}
        	inv_cubes = new_exp;
		}
		
		cnode = pnode;
		
		if(pnode == root)
			break;		
					
	}
}


void get_involve_pla(map<int, btNode*> &leaf_set, map<string, int> &name_pos, int len, string &this_inv_cube, string &this_inv_pla)
{
	map<int, btNode*>::iterator itrm_ib;
	map<string, int>::iterator itrm_si;
	map<int, int> pos_sign;
	map<int, int>::iterator itrmi;
	string name, true_name;
	int sign;
	
	for(int i = 0; i < this_inv_cube.size(); i++)
	{
		char c = this_inv_cube[i];
		int p = c - 48;
		itrm_ib = leaf_set.find(p);
		name = itrm_ib->second->data;		
		get_true_name(name, true_name, sign);
		itrm_si = name_pos.find(true_name);
		pos_sign.insert(pair<int, int>(itrm_si->second, sign));
	}
	
	string str;
	str.append(len, '-');
	for(itrmi = pos_sign.begin(); itrmi != pos_sign.end(); itrmi++)
	{
		int pos = itrmi->first;
		int sign = itrmi->second;
		if(sign)
			str[pos] = '1';
		else if(!sign)
			str[pos] = '0';
	}
	this_inv_pla = str;
}


void build_tree_from_exp(string &str, btNode **root)
{
   *root = CreateInfixTree(str);
//   int numPT = 0;
//   cout << "compute numPT: " << endl;
//   compNumPT(root, numPT); 
//   cout << "inorder print: " << endl; 
   InOrderPrintTree(*root);
//   cout << endl; 
//   cout << "number of cubes: " << numPT << endl;
   
}


/*
void build_btree(vector<string> &lit_unit, vector<string> &sop_set, vector<string> &factor_set, map<string, btNode*> &factor_exp_trees)
{
	map<string, btNode*>::iterator itrm_sb1;
	
	string str, s;
	for(int i = 0; i < lit_unit.size(); i++)
	{
		str = lit_unit[i];
		string newStr;		
		add_space_star(str, newStr);
//		cout << "new str: " << newStr << endl;		
		if(newStr.find('(') != string::npos || newStr.find(')') != string::npos)
			factor_set.push_back(newStr);
		else
			sop_set.push_back(newStr);
	}
	
	//build a binary tree for each factor expression
//	cout << "factor_set: " << endl;
	for(int i = 0; i < factor_set.size(); i++)	
	{
//		cout << factor_set[i] << endl;
		btNode *root;			
		build_tree_from_exp(factor_set[i], &root);
		factor_exp_trees.insert(pair<string, btNode*>(factor_set[i], root));		
	}
}
*/



