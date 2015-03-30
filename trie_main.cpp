#include "Fib.h"
#include <iostream>
#include <windows.h>	
#include <time.h>
#include "bloom_filter.h"
#include <set>
//#include <string>
#include <iomanip>
using namespace std;

set<string> routerSet;
bitset<65536> bit_array(0);

bool	bIpFormat=true;						
#define LINESIZE	66
#define BF_NUM 16
string DupString(char cElement,int iLen,bool bNewLine){
	string strRet;
	for (int i=1;i<=iLen;i++)
	{
		if(i==iLen && bNewLine==true)
			strRet += "\n";
		else
			strRet += cElement;
	}
	return strRet;
}

bool Exist(char * FileName)
{
	FILE*   fp;
	if((fp=fopen(FileName, "r "))==NULL)
	{
		return false;
	}
	else
	{
		fclose(fp);
		return true;
	}
}

void ShowTipInfo()
{

	printf("Trie Merging Algorithm\n");
	printf("%s",DupString('=',LINESIZE,true).c_str());
	printf("RIB Format  Binary(B) or IP(I):");

	char cResponse=getchar();
	//printf("%c",cResponse);
	while(cResponse!='B'||cResponse!='b'||cResponse!='I'||cResponse!='i')
	{
		if(cResponse=='B'||cResponse=='b')
		{
			bIpFormat=false;
			break;
		}
		if(cResponse=='I'||cResponse=='i')
		{
			break;
		}

		cResponse=getchar();
		//printf("%c",cResponse);
	}
}


int count_ip_cpr(char *ribFile2)
{
	static int fib_number = 1;
	CFib tFib2 = CFib();
	ofstream fout("result725",ios::app);
	tFib2.BuildFibFromFile(ribFile2);
	tFib2.sub_CountOutputMergedTrie(tFib2.m_pTrie,0,0);
	int count = tFib2.cpr_count;
	tFib2.GetNodeCounts();
	cout << "Fib " << fib_number << ":\t" << count << " \t pre_fix : \t"<<tFib2.solidNodeCount<<endl;
	fout << "Fib " << fib_number << ":\t" << count << " \t total_pre_fix : \t" << tFib2.solidNodeCount << endl;

	tFib2.FreeTrie(tFib2.m_pTrie);
	fout.close();
	fib_number++;
	return 0;
}

int count_ipv6(char *file)
{
	static int fib_num = 1;

	int output = 10;
	int count = 0;
	char pre[129];
	char tmp[129];
	ifstream fin(file);
	if (fin.fail())
	{
		cout << "open failed.." << endl;
		return -1;
	}
	memset(pre,0,sizeof(pre));
	memset(tmp,0,sizeof(tmp));

	int totalnum = 0;

	while (!fin.eof())
	{
		if (strlen(pre) == 0)
		{
			cout << "sssssssssssssss" << endl;
			fin >> pre;
			totalnum++;
			if (fin.eof())
			{
				break;
			}
		}
		
		fin >> tmp;
		totalnum++;
		pre[strlen(pre) - 1] = '1';
		if (strcmp(pre, tmp) == 0)
		{
			if (output>0)
			{

				cout << pre << endl;
				cout << tmp << endl;
				cout << "------------------------" << endl;
				output--;
			}
			count++;
		}
		
		memset(pre, 0, sizeof(pre));
		strcpy(pre, tmp);
		memset(tmp,0,sizeof(tmp));
		



	}

	cout << "IpV6_Fib " << fib_num << ":\t" << count  << "\t totalnum: \t"<<totalnum<<endl;
	ofstream fout("result725",ios::app);
	fout << "IpV6_Fib " << fib_num << ":\t" << count << "\t totalnum: \t" << totalnum << endl;
	fout.close();


		fib_num++;
	return 0;
}

int TrieMerging(CFib *pFib1,char *ribFile2)
{
	char *ribFile2_IP="ribFile2_IP.txt";       

	CFib tFib2= CFib();	

	unsigned int iEntryCount1 = 0;
	unsigned int iEntryCount2 = 0;

	if(bIpFormat)
	{
		iEntryCount2=tFib2.BuildFibFromFile(ribFile2);
	}
	else
	{
		tFib2.ConvertBinToIP(ribFile2,ribFile2_IP);
		iEntryCount2=tFib2.BuildFibFromFile(ribFile2_IP);
	}

	//tFib2.GetNodeCounts();

	FILE *fp=fopen("result.txt","a+");
	printf("FIB %s: solid nodes=\t%d\tleaf nodes=\t%d\tall nodes=\t%d\n",ribFile2,tFib2.solidNodeCount,tFib2.leafNodeCount,tFib2.allNodeCount);
	fprintf(fp,"FIB %s: solid nodes=\t%d\tleaf nodes=\t%d\tall nodes=\t%d\n",ribFile2,tFib2.solidNodeCount,tFib2.leafNodeCount,tFib2.allNodeCount);

	pFib1->outline(&tFib2);  //get the profile of the merged trie

	tFib2.FreeTrie(tFib2.m_pTrie);	

	//pFib1->GetNodeCounts(); // just for MergedTrie

	printf("The merged trie: solid nodes=\t%d\tleaf nodes=\t%d\tall nodes=\t%d\n\n",pFib1->solidNodeCount,pFib1->leafNodeCount,pFib1->allNodeCount);
	fprintf(fp,"The merged trie: solid nodes=\t%d\tleaf nodes=\t%d\tall nodes=\t%d\n\n",pFib1->solidNodeCount,pFib1->leafNodeCount,pFib1->allNodeCount);

	//printf("Fib %s =\t%d\t After trie merging,prefix nodes=\t%u\tallNodeCount=\t%d\n",ribFile2,tFib2.solidNodeCount,pFib1->solidNodeCount,pFib1->allNodeCount);
	//fprintf(fp,"Fib %s =\t%d\t After trie merging,prefix nodes=\t%u\tallNodeCount=\t%d\n",ribFile2,tFib2.solidNodeCount,pFib1->solidNodeCount,pFib1->allNodeCount);

	fclose(fp);

	return 0;
}

void pushToMergedTrie(CFib *pFib1,char *ribFile2,int number)
{
	printf("pushing %s...\n",ribFile2);
	char *ribFile2_IP="ribFile2_IP.txt";
	CFib tFib2= CFib();	

	unsigned int iEntryCount1 = 0;
	unsigned int iEntryCount2 = 0;

	if(bIpFormat)
	{
		iEntryCount2=tFib2.BuildFibFromFile(ribFile2);
	}
	else
	{
		tFib2.ConvertBinToIP(ribFile2,ribFile2_IP);
		iEntryCount2=tFib2.BuildFibFromFile(ribFile2_IP);
	}

	pFib1->pushToOutline(&tFib2,number);  //number starts from 1
	tFib2.FreeTrie(tFib2.m_pTrie);
}

int main(int argc, char * argv[])
{
	

	/*
	char file_str_ip[] = "str_ip.txt";
    trans_strip_to_intip(file_str_ip);
	system("pause");
	return 0;
	*/

	//time_t t_s_main,t_e_mergedtrie;
	//time(&t_s_main);

	CFib MergedFib= CFib();	

	//ShowTipInfo();
	//printf("The size of TrieNode is %d bytes\n",sizeof(TrieNode));

	//system("dir /B /ON rrc*.txt > fileList");
	
	FILE * fList=fopen("ipv6_fileList.txt","r");
	
	while (!feof(fList)) 
	{
		char filename[200];
		memset(filename,0,sizeof(filename));
		fscanf(fList,"%s",filename);
		cout << "filename is: "<<filename << endl;
		if (strlen(filename) > 3)
		{
			//count_ip_cpr(filename);
			count_ipv6(filename);
		}
			//TrieMerging(&MergedFib,filename);

	}
	fclose(fList);

	system("pause");
	return 0;

	int number=0;   //the number of RIBs
	fList=fopen("fileList","r");
	while (!feof(fList)) 
	{
		number++;
		char filename[200];
		memset(filename,0,sizeof(filename));
		fscanf(fList,"%s",filename);

		if(strlen(filename)>3)
			pushToMergedTrie(&MergedFib,filename,number);

	}
	fclose(fList);

	number--;  

	MergedFib.routerNum=number;

	MergedFib.OutputAlltables(MergedFib.m_pTrie);

	//MergedFib.OutputMergedTrie("MergedTrie(only_port).txt",true);
	//MergedFib.OutputMergedTrie("MergedTrie.txt",false);
	//MergedFib.OutputMergedTrie2("MergedTrie_ce.txt", false);
	MergedFib.GetNodeCounts();

	ofstream fout("file_MergedTrie_Info",ios::app); 

	StandardBF *(a_BloomFilter[16]);

	EXSBF *extraBF= new EXSBF();

	printf("The final merged trie:\nsolid nodes=\t%d\tleaf nodes=\t%d\tall nodes=\t%d\n",MergedFib.solidNodeCount,MergedFib.leafNodeCount,MergedFib.allNodeCount);
	
	fout << "The final merged trie:\nsolid nodes=\t" << MergedFib.solidNodeCount << " leaf nodes= " << MergedFib.leafNodeCount << " all nodes= " <<  MergedFib.allNodeCount << endl;

	char *filename1 = "MergedTrie.txt";
	char *filename = "MergedTrie_ce.txt";

	for (int i=0;i<33;i++)
	{
		printf("Level %d    Prefix nodes = %d\n",i,MergedFib.numofDepth[i]);
		//printf("Level numofLeaf %d    Prefix nodes = %d\n", i, MergedFib.numofLeaf[i]);
		fout << "Level " << i << " Prefix nodes_num= " << MergedFib.numofDepth[i] << endl;

		
		if (i<16)
			a_BloomFilter[i] = new StandardBF();

	}

	

	for (int i = 0; i<33; i++)
	{
		printf("Level %d    Yellow nodes = %d\n", i, MergedFib.numofYellow[i]);		
		printf("Level %d    single leaf nodes = %d\n", i, MergedFib.numofLeaf_noce[i]);
		printf("Level %d    Mid nodes = %d\n", i, MergedFib.numofMid[i]);
		printf("Level %d    all leaf nodes = %d\n", i, MergedFib.numofDepth[i]);

	//	fout << "Level " << i << " Yellow nodes_num= " << MergedFib.numofYellow[i] << " Single leaf nodes_num : " << MergedFib.numofLeaf_noce[i] << " Mid nodes_num= " << MergedFib.numofMid[i] << endl;
	}
	/*
	ofstream foutR4("result_ratio_ce.txt", ios::app);
	char strtmp[5];
	foutR4 << "---------------------------------------------" << endl;
	foutR4 << "basic Bf Node" << endl;
	foutR4 << "---------------------------------------------" << endl;
	for (int i = 17; i < 33; i++)
	{
		printf("Level %d    Yellow nodes = %d\n", i, MergedFib.numofYellow[i]);

	//	fout << "Level " << i << " Yellow nodes_num= " << MergedFib.numofYellow[i] << " Mid nodes_num= " << MergedFib.numofMid[i] << endl;
		sprintf(strtmp, "%d", i);
		
		foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
		foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofDepth[i] << endl;
	}
	foutR4 << "---------------------------------------------" << endl;
	foutR4 << "c&e Bf " << endl;
	foutR4 << "---------------------------------------------" << endl;
	for (int i = 17; i < 33; i++)
	{
		sprintf(strtmp, "%d", i);
		
		foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
		if (i != 32)
		{
			if (i == 17)
			{
				foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofDepth[i] + MergedFib.numofYellow[i] << endl;
			}
			else
				foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofLeaf_noce[i] + MergedFib.numofYellow[i] << endl;
		}
		else
		{
			foutR4 << setw(14) << strtmp << "\t" << 0 << endl;
		}
	}
	foutR4 << "---------------------------------------------" << endl;
	foutR4 << "Yellow Node" << endl;
	foutR4 << "---------------------------------------------" << endl;
	for (int i = 17; i < 33; i++)
	{
		sprintf(strtmp, "%d", i);

		
		foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
		foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofYellow[i] << endl;
	}
	foutR4 << "---------------------------------------------" << endl;
	foutR4 << "Mid Node" << endl;	
	foutR4 << "---------------------------------------------" << endl;
	for (int i = 17; i < 33; i++)
	{
		sprintf(strtmp, "%d", i);
			
		foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
		foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofMid[i] << endl;
			
	}
	*/

	//foutR4.close();
	fout.close();
	int bf_num = 16;

	//unsigned char bit_array[65535];	
	//	memset(bit_array, 0, sizeof(bit_array));


	
	for (int i = 8; i <= 8; i++)
	{
		//char temp[30] = "rand_ip_int_my";
		//char temp[30] = "rand_32_100M.tr_22";
		char temp[30] = "rrc00_100M.tr";
		construct_bloom_filter(1, filename1, a_BloomFilter, NULL, MergedFib.numofDepth, BF_NUM);
		bf_query2(temp, &MergedFib, a_BloomFilter, NULL);

		
		construct_bloom_filter(2, filename1, a_BloomFilter, NULL, MergedFib.numofDepth, BF_NUM);
		bf_query2(temp, &MergedFib, a_BloomFilter, NULL);
		
	}
	
	for (int i = 8; i <= 8; i++)
    //for(int i = 16;i>4;i--)
	{
		//char temp[30] = "rand_32_100M.tr_22";
		//char temp[30] = "rand_ip_int_my";
		char temp[30] = "rrc00_100M.tr";
		construct_bloom_filterce(1, filename, a_BloomFilter, NULL, MergedFib, BF_NUM - 1);
		bf_queryce(temp, &MergedFib, a_BloomFilter, NULL);

		
		construct_bloom_filterce(2, filename, a_BloomFilter, NULL, MergedFib, BF_NUM - 1);
		bf_queryce(temp, &MergedFib, a_BloomFilter, NULL);

		
	//	construct_bloom_filter(2, filename, a_BloomFilter, NULL, MergedFib.numofDepth, BF_NUM);
		//	char temp[] = "rrctest1M.txt";

///
	//	bf_query2(temp, &MergedFib, a_BloomFilter, NULL);


	//	strcpy(temp, "rrc00_100M.tr_22");
	//	bf_query2(temp, &MergedFib, a_BloomFilter, NULL);


		//construct_bloom_filter(i, filename, a_BloomFilter, extraBF, MergedFib.numofDepth, BF_NUM);
	//	cout << "extra msg: m,n,k " << extraBF->Get_bf_m() << ends << extraBF->Get_bf_n() << " " << extraBF->Get_bf_k() << endl;

		//char temp1[30] = "rand_32_100M.tr_22";

	

		//bf_query2(temp, &MergedFib, a_BloomFilter, extraBF);
	//	strcpy(temp1, "rrc00_100M.tr_22");
	//	bf_query2(temp1, &MergedFib, a_BloomFilter, extraBF);


		////////////////////////////////////////////////////////

		//construct_bloom_filter2(8, filename, a_BloomFilter, NULL, MergedFib.numofDepth, BF_NUM);
		//	char temp[] = "rrctest1M.txt";

		//char temp[30] = "rand_32_100M.tr_22";
	//	bf_query2(temp, &MergedFib, a_BloomFilter, NULL);
		//	strcpy(temp, "rrc00_100M.tr_22");
		//	bf_query2(temp, &MergedFib, a_BloomFilter, NULL);

      /*
		construct_bloom_filter2(4, filename, a_BloomFilter, extraBF, MergedFib.numofDepth, BF_NUM);
		cout << "extra msg: m,n,k " << extraBF->Get_bf_m() << ends << extraBF->Get_bf_n() << " " << extraBF->Get_bf_k() << endl;

		//char temp1[30] = "rand_32_100M.tr_22";

		cout << "k = " << extraBF->Get_bf_k() << endl;
		cout << "m = " << extraBF->Get_bf_m() << endl;
		cout << "n = " << extraBF->Get_bf_n() << endl;
		//	Sleep(10000);

		bf_query2(temp, &MergedFib, a_BloomFilter, extraBF);


		construct_bloom_filter2(8, filename, a_BloomFilter, extraBF, MergedFib.numofDepth, BF_NUM);
		cout << "extra msg: m,n,k " << extraBF->Get_bf_m() << ends << extraBF->Get_bf_n() << " " << extraBF->Get_bf_k() << endl;

		//char temp1[30] = "rand_32_100M.tr_22";

		cout << "k = " << extraBF->Get_bf_k() << endl;
		cout << "m = " << extraBF->Get_bf_m() << endl;
		cout << "n = " << extraBF->Get_bf_n() << endl;
		//	Sleep(10000);

		bf_query2(temp, &MergedFib, a_BloomFilter, extraBF);


		/////////////////////////////////////////////////////
		//construct_bloom_filter2(16, filename, a_BloomFilter, NULL, MergedFib.numofDepth, BF_NUM);
		//	char temp[] = "rrctest1M.txt";

	//	char temp[30] = "rand_32_100M.tr_22";
	//	bf_query2(temp, &MergedFib, a_BloomFilter, NULL);
		//	strcpy(temp, "rrc00_100M.tr_22");
		//	bf_query2(temp, &MergedFib, a_BloomFilter, NULL);


		construct_bloom_filter2(16, filename, a_BloomFilter, extraBF, MergedFib.numofDepth, BF_NUM);
		cout << "extra msg: m,n,k " << extraBF->Get_bf_m() << ends << extraBF->Get_bf_n() << " " << extraBF->Get_bf_k() << endl;

	//	char temp1[30] = "rand_32_100M.tr_22";

		cout << "k = " << extraBF->Get_bf_k() << endl;
		cout << "m = " << extraBF->Get_bf_m() << endl;
		cout << "n = " << extraBF->Get_bf_n() << endl;
		//	Sleep(10000);

		bf_query2(temp, &MergedFib, a_BloomFilter, extraBF);
		
		for (int i = 5; i!=8,i < 16; i++)
		{

			construct_bloom_filter2(i, filename, a_BloomFilter, extraBF, MergedFib.numofDepth, BF_NUM);
			cout << "extra msg: m,n,k " << extraBF->Get_bf_m() << ends << extraBF->Get_bf_n() << " " << extraBF->Get_bf_k() << endl;

			//	char temp1[30] = "rand_32_100M.tr_22";

			cout << "k = " << extraBF->Get_bf_k() << endl;
			cout << "m = " << extraBF->Get_bf_m() << endl;
			cout << "n = " << extraBF->Get_bf_n() << endl;
			//	Sleep(10000);

			bf_query2(temp, &MergedFib, a_BloomFilter, extraBF);
		
		}*/

	}//for
	cout << "extra msg: m,n,k " << extraBF->Get_bf_m() << ends << extraBF->Get_bf_n() << " " << extraBF->Get_bf_k() << endl;

		

	/*
	for (int i = 8; i <= 8; i++)
	{
		construct_bloom_filter(i, filename, a_BloomFilter, extraBF, MergedFib.numofDepth, BF_NUM);
		char temp[] = "rand_32_10M.tr";
		bf_query(temp, &MergedFib, a_BloomFilter, extraBF);
		strcpy(temp, "rrc00_10M.tr");
		bf_query(temp, &MergedFib, a_BloomFilter, extraBF);

	 }*/
	
   // construct_bloom_filter(filename,bit_array,a_BloomFilter,extraBF, MergedFib.numofDepth,BF_NUM);	
	//bf_query("rrctest1M", a_BloomFilter, bit_array, extraBF);


	printf("\nmission complete!\n");
	
//	cout << "####Time MergedTrie  " << t_e_mergedtrie - t_s_main << " Construct BF: " << t_e_bf - t_e_mergedtrie << " BF_query one file: " << t_e_bf_query - t_e_bf << "total time :" << t_e_bf_query -t_s_main << endl;
	
	for (int i = 0; i<16; i++)
	{		
		delete a_BloomFilter[i];
	}
	system("pause");
	return 0;
}



int _main(int argc, char * argv[])
{
	char  a_str[12][30] = { "rrc00_2013.6.8.8.txt",
		"rrc01_2013.6.8.8.txt",
		"rrc03_2013.6.8.8.txt",
		"rrc04_2013.6.8.8.txt",
		"rrc05_2013.6.8.8.txt",
		"rrc06_2013.6.8.8.txt",
		"rrc07_2013.6.8.8.txt",
		"rrc10_2013.6.8.8.txt",
		"rrc11_2013.6.8.8.txt",
		"rrc12_2013.6.8.8.txt",
		"rrc13_2013.6.8.8.txt",
		"rrc14_2013.6.8.8.txt" };

	/*
	char file_str_ip[] = "str_ip.txt";
	trans_strip_to_intip(file_str_ip);
	system("pause");
	return 0;
	*/



	CFib MergedFib = CFib();		
	ofstream foutR4("result_ratio_ce.txt", ios::app);

	for (int f_num = 1; f_num <= 12; f_num++)
	{
		foutR4 << "##################################" << endl;
		foutR4 << "Rib_num : " <<f_num<< endl;

		for (int ii = 0; ii < f_num; ii++)
		{
			cout << "filename is: " << a_str[ii] << endl;
			TrieMerging(&MergedFib, a_str[ii]);

		}
		
		int number = 0;
		for (int ii = 0; ii < f_num; ii++)
		{		 
			number++;
		    pushToMergedTrie(&MergedFib, a_str[ii], number);
		}


			

			number--;

			MergedFib.routerNum = number;

			MergedFib.OutputAlltables(MergedFib.m_pTrie);

			MergedFib.OutputMergedTrie("MergedTrie(only_port).txt", true);
			MergedFib.OutputMergedTrie("MergedTrie.txt", false);
			//MergedFib.OutputMergedTrie2("MergedTrie.ip", false);
			MergedFib.GetNodeCounts();

			ofstream fout("file_MergedTrie_Info", ios::app);

			StandardBF *(a_BloomFilter[16]);

			EXSBF *extraBF = new EXSBF();

			printf("The final merged trie:\nsolid nodes=\t%d\tleaf nodes=\t%d\tall nodes=\t%d\n", MergedFib.solidNodeCount, MergedFib.leafNodeCount, MergedFib.allNodeCount);

			fout << "The final merged trie:\nsolid nodes=\t" << MergedFib.solidNodeCount << " leaf nodes= " << MergedFib.leafNodeCount << " all nodes= " << MergedFib.allNodeCount << endl;

			char *filename = "MergedTrie.txt";

			for (int i = 0; i<33; i++)
			{
				printf("Level %d    Prefix nodes = %d\n", i, MergedFib.numofDepth[i]);
				//printf("Level numofLeaf %d    Prefix nodes = %d\n", i, MergedFib.numofLeaf[i]);
				fout << "Level " << i << " Prefix nodes_num= " << MergedFib.numofDepth[i] << endl;


				if (i<16)
					a_BloomFilter[i] = new StandardBF();

			}


			for (int i = 0; i<33; i++)
			{
				printf("Level %d    Yellow nodes = %d\n", i, MergedFib.numofYellow[i]);
				printf("Level %d    Mid nodes = %d\n", i, MergedFib.numofMid[i]);
				fout << "Level " << i << " Yellow nodes_num= " << MergedFib.numofYellow[i] << " Mid nodes_num= " << MergedFib.numofMid[i] << endl;
			}

			char strtmp[5];
			foutR4 << "---------------------------------------------" << endl;
			foutR4 << "basic Bf Node" << endl;
			foutR4 << "---------------------------------------------" << endl;
			for (int i = 17; i < 33; i++)
			{
				printf("Level %d    Yellow nodes = %d\n", i, MergedFib.numofYellow[i]);

				fout << "Level " << i << " Yellow nodes_num= " << MergedFib.numofYellow[i] << " Mid nodes_num= " << MergedFib.numofMid[i] << endl;
				sprintf(strtmp, "%d", i);

				foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
				foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofDepth[i] << endl;
			}


		
			foutR4 << "---------------------------------------------" << endl;
			foutR4 << "c&e Bf " << endl;
			foutR4 << "---------------------------------------------" << endl;
			for (int i = 17; i < 33; i++)
			{
				sprintf(strtmp, "%d", i);

				foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
				if (i != 32)
				{
					if (i == 17)
					{
						foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofDepth[i] + MergedFib.numofYellow[i] << endl;
					}
					else
						foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofLeaf_noce[i] + MergedFib.numofYellow[i] << endl;
				}
				else
				{
					foutR4 << setw(14) << strtmp << "\t" << 0 << endl;
				}
			}
			foutR4 << "---------------------------------------------" << endl;
			foutR4 << "Yellow Node" << endl;
			foutR4 << "---------------------------------------------" << endl;
			for (int i = 17; i < 33; i++)
			{
				sprintf(strtmp, "%d", i);


				foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
				foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofYellow[i] << endl;
			}
			foutR4 << "---------------------------------------------" << endl;
			foutR4 << "Mid Node" << endl;
			foutR4 << "---------------------------------------------" << endl;
			for (int i = 17; i < 33; i++)
			{
				sprintf(strtmp, "%d", i);

				foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
				foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofMid[i] << endl;

			}
			foutR4 << "---------------------------------------------" << endl;
			foutR4 << "Single leaf Node" << endl;
			foutR4 << "---------------------------------------------" << endl;
			for (int i = 17; i < 33; i++)
			{
				sprintf(strtmp, "%d", i);


				foutR4 << setw(14) << setiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(7);
				foutR4 << setw(14) << strtmp << "\t" << MergedFib.numofLeaf_noce[i] << endl;
			}



		
		
		



	
	}// for f_num;




	

	foutR4.close();
	
	


	

	printf("\nmission complete!\n");

	//	cout << "####Time MergedTrie  " << t_e_mergedtrie - t_s_main << " Construct BF: " << t_e_bf - t_e_mergedtrie << " BF_query one file: " << t_e_bf_query - t_e_bf << "total time :" << t_e_bf_query -t_s_main << endl;


	system("pause");
	return 0;
}