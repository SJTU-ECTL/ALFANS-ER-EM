#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <utility>

#include "head/HashTable.h"


using namespace std;


int hashsize[]={11, 29, 59, 113, 223, 489, 997, 1499, 2999, 4999, 9973, 19997, 40009, 104393, 113381, 120011, 130003, 140009, 150001, 160001, 170003, 180001, 190027, 200003, 210011, 220009, 230003, 240007, 250007, 260003, 270001, 280001, 290011, 300007, 717419, 1253621, 1930801, 3125051, 4351433, 6866927, 8150399, 9458311}; 

HashTable::HashTable()
{}

int HashTable::Init_HashTable()
{  
    int i;  
    count = 0;  
    size = 0;
    c = 0;
    Hash_length = hashsize[0];  
    elem = new Elemtype[Hash_length];  
    if(!elem)  
    {  
        cout << "malloc fails when initialize hash table!"<<endl;
        exit(1);  
    }  
    for(i = 0; i < Hash_length; i++)  
        elem[i] = NULLKEY;  
    return 1;  
}  
  
void HashTable::Destroy_HashTable()  
{  
    delete []elem;  
    elem = NULL;  
    count = 0;  
    size = 0;  
}  
  
unsigned HashTable::Hash(KeyType k) 
{  
    return k % Hash_length;  
}  
  
void HashTable::Collision(int &p,int d) 
{  
	//	cout<<"In collision, p = "<<p<<", d = "<<d<<", Hash_length = "<<Hash_length<<endl;
        p = (int)((p + pow(d,2))) % Hash_length;
       // p = (int)((p + d)) % Hash_length;
    //    cout<<"In collision, p = "<<p<<endl;
}  
  
bool HashTable::Search_Hash(KeyType k,int &p)
{  
    //cout << "coming into search hash !"<<endl;
        
        c=0;  
        p=Hash(k); 
    //    cout <<"hash length = "<<Hash_length<<", In search, p = "<<p;
    //    cout << ", elem["<<p<<"] = "<<elem[p]<<endl;
        while(elem[p]!=NULLKEY && elem[p]!=k)    //elem[p] is occupied by other element.
        {  
            c++;  
            p = Hash(k);
            if(c < Hash_length)  
            {           	
                Collision(p,c);  
                if( p < 0 )
                    return 0;
            }
            else  
                return 0; 
        }  
        if(elem[p]==k)  
            return 1;  
        else  
            return 0;  
}  

int HashTable::Insert_Hash(Elemtype e) 
{  
  
        int p;  
        if(Search_Hash(e,p) == 1)    //The element has been in the table.
            return -1;

        else if(c<hashsize[size]/2) 
        {  
        	int i = 0; 
		  	while(i < Hash_length)
		  	{
		  		p = Hash(e);
		  		Collision(p,i);  
		  		if(elem[p] == NULLKEY || elem[p] == DELKEY)
		  		{
		  			elem[p] = e;
		  			count++;  
		  			return 1;
		  		}
		  		else i++;
			}
        }  
        else  
        {
            ReCreate_HashTable(); 
            return 0;
        }
        
}  

int HashTable::Get_Data(int p)  
{  
    return elem[p];
        //cout<<"number："<<elem[p]<<endl;  
}

void HashTable::Del_Data(int p)
{
	if(elem[p] != NULLKEY)
		elem[p] = DELKEY;
}  
  
void HashTable::ReCreate_HashTable() 
{  

        int i, count2 = count;  
        Elemtype *p, *elem2 = new Elemtype[count];  
        p = elem2;  
 
        for(i=0;i<Hash_length;i++) 
            if(elem[i]!=NULLKEY)  
                *p++=*(elem+i);  
        count=0;
        
        size++; 
        Hash_length=hashsize[size];  
        p=new Elemtype[Hash_length];  
        if(!p)  
        {  
            cout<<""<<endl;  
            exit(0);  
        }  
        delete []elem;
        elem=p;  
        for(i=0;i<Hash_length;i++)  
            elem[i]=NULLKEY;  
        for(p=elem2;p<elem2+count2;p++) 
            int result = Insert_Hash(*p); 
        delete []elem2;
    }  
  
void HashTable::Traverse_HashTable()  
{  
    //	cout << "Traverse the hashtable: "<<Hash_length<<endl; 
        for(int i=0;i<Hash_length;i++)  
        {
        	if(elem[i] == NULLKEY || elem[i] == DELKEY)  
        		continue;
                //cout<<i<<","<<NULLKEY<<endl;  
            else if(elem[i]!=NULLKEY)  
                cout <<elem[i]<<" ";  
        }
		cout << endl;  
}  
  
    
int HashTable::Get_Size()  
{  
    	//int size = 0;
    	//for(int i = 0; i < Hash_length; i++)
    	//	if(elem[i] != NULLKEY)
    	//		size++;
    	return Hash_length;
}  

int HashTable::Data_Size()  
{  
    int data_size = 0;
    for(int i = 0; i < Hash_length; i++)
    	if(elem[i] != NULLKEY && elem[i] != DELKEY)
    		data_size++;
    return data_size;
}  
    
bool HashTable::empty()  
{  
    int size = 0;
    for(int i = 0; i < Hash_length; i++)
        if(elem[i] != NULLKEY &&  elem[i] != DELKEY)
    	    size++;
    	if(size == 0)
    		return 1;
    	else
    		return 0;
}  
      
