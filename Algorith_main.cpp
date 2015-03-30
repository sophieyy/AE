#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <windows.h>
#include <math.h>

#include <iostream>
#include <string>
#include "bloom_filter.h"

//#define PREFIX_NUM	20
//#define PREFIX_LEN 4
using namespace std;
int _main()
{ 
	/*
	rrctest10000k.txt
rrctest100000k.txt
rrctest1000000k.txt
	*/
	//char *filename = "rand_32_2000K_15738.tr";
	//char *filename = "rand_32_10000K_5526.tr";
	//char *filename = "rand_32_1000K_10093.tr";
	//	char *filename = "rrctest.txt"; //1000k
//	char *filename = "rrctest10000k.txt";
//	char *filename = "rrctest100000k.txt";
	 //char *filename = "rand_32_100000K_13713.tr";
	
	char filename[50];
	char filelist[] = "filelist.txt";
	ifstream fin(filelist);
	while (!fin.eof())
	{
		fin >> filename;
		cout << filename << endl;
		false_positive_std_bf(filename);
		//q_test_false_positive_qwstdbf(filename);
		
	}
	fin.close();
	
	
	//bloom_filter(filename);
	//bloom_filter_false_positive_ratio();

	//false_positive_std_qw_bf(filename); //change sBF

	//q_test_false_positive_qwstdbf(filename);

	//test_false_positive_qwstdbf(filename);
	//test_false_positive_stdbf(filename);
	
	//cout << "hello world" << endl;
	system("pause");
	return 0;
}