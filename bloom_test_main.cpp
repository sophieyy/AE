#include <stdio.h>
#include <stdlib.h>
#include "std_bf.h"
#include <fstream>
#include <windows.h>
#include <math.h>
#include <iostream>
using namespace std;


#define PREFIX_NUM	20
#define PREFIX_LEN 4
/*
void test_bf_query(StandardBF *sBF1)
{

	ifstream fin;
	fin.open("E:\\Dev-Cpp\\router\\test_ipprefix.4bytes");
	int j = 0;
	unsigned int x;
	unsigned char prefix[PREFIX_NUM][PREFIX_LEN];
	if (!fin.is_open())
	{
		cout << "fail" << endl;
		return;
	}

	while (!fin.eof() && j<PREFIX_NUM)
	{

		for (int i = 0; i<PREFIX_LEN; i++)
		{
			fin >> x;

			//	cout<<x<<" ";
			prefix[j][i] = (char)x;
		}
		j++;
		//cout<<endl;
	}
	fin.close();

	for (int i = 0; i<PREFIX_NUM; i++)
	{
		for (j = 0; j<PREFIX_LEN; j++)
		{
			cout << (int)prefix[i][j] << ends;
		}
		if (sBF1->query(prefix[i], PREFIX_LEN))
			cout << "the prefix is found.\n";
		else
			cout << "the prefix is not found.\n";
		cout << endl;
	}


}*/

void _main(int argc,char* argv[])
{
	cout << _pgmptr << endl;
	int x=0,j=0;
	unsigned char prefix[PREFIX_NUM][PREFIX_LEN];

	/*char * fileName;
	if (argc>1)
		fileName = argv[1];
	else
	{
		string fileN = "ipprefix.4bytes";
		fileName = (char *)fileN.c_str();
	}

	ifstream fin(fileName);*/
	ifstream fin("ipprefix.4bytes");

	while (!fin.eof() && j<PREFIX_NUM)
	{
		for (int i=0;i<PREFIX_LEN;i++)
		{
			fin>>x;
			prefix[j][i]=(char)x;
		}
		j++;
		//cout<<endl;
	}
	fin.close();

	for(int i=0;i<PREFIX_NUM;i++)
	{
		for(j=0;j<PREFIX_LEN;j++)
			cout<<(int)prefix[i][j]<<ends;
		cout<<endl;
	}

	uint k = 10;
	uint m = (int)(k * PREFIX_NUM / log(2)); 

	StandardBF *sBF1=new StandardBF(m,k);

	
	for (uint i=0;i<PREFIX_NUM;i++)
		sBF1->insert(prefix[i],PREFIX_LEN);

	sBF1->outputOHABF("sbf1.bf");
	cout << "k = " << sBF1->Get_bf_k() <<endl;
	cout << "m = " << sBF1->Get_bf_m() <<endl;
	cout << "n = " << sBF1->Get_bf_n() <<endl;

	uchar p[PREFIX_LEN];

	p[0]=10;
	p[1]=1;
	p[2]=0;
	p[3]=2;
	cout<<(int)p[0]<<(int)p[1]<<(int)p[2]<<(int)p[3]<<endl;
	if(sBF1->query(p,PREFIX_LEN))
		cout << "the prefix is found.\n";
	else
		cout << "the prefix is not found.\n";

	for(j=0;j<PREFIX_NUM;j++)
	{
		if(sBF1->query(prefix[j],PREFIX_LEN))
			cout<<(int)prefix[j][0]<<ends<<(int)prefix[j][1]<<ends<<(int)prefix[j][2]<<ends<<(int)prefix[j][3]<<ends<<"the prefix is found.\n";
		else
			cout<<(int)prefix[j][0]<<ends<<(int)prefix[j][1]<<ends<<(int)prefix[j][2]<<ends<<(int)prefix[j][3]<<ends<<"the prefix is not found.\n";
	}
	//test_bf_query(sBF1);
	system("pause");
}



