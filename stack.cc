#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include "./head/stack.h"
#include "cudd/bnet.h"

using namespace std; 
	
//Constructor
stack::stack(int length)                          
{  
   data = new int[length];
   size = length;
   len = 0;
}  

stack::~stack()                          
{  
	delete []data;
}  
	  
//Push into the stack
void stack::push(int d)  
{  
	if(len < size)
		data[len] = d;
	else
	{
		size = size * 2;
		int *new_data = new int[size];
		int i = 0;
		while(data[i] != '\0')
		{
			new_data[i] = data[i];
			i++;
		}
		delete []data;
		data = new_data;
		new_data = NULL;		
	}
	len = len + 1;
	
}  
	  
//Pop from the stack 
int stack::pop()  
{  
	if(len <= 0)
		return -1;
	int pop =  data[len-1];
	data[len-1] = -1;
	len = len - 1;
	return pop;
}  
	  

//Get the size of the stack  
int stack::get_num()  
{  
	return len;
}  
	      
//Judge if the stack is empty  
bool stack::empty()  
{  
	if(len == 0)
		return 1;
	return 0;
}  

void stack::traverse()
{
	int i = 0;
	while(data[i] != '\0')
		cout << data[i++];
}
	
	
