#ifndef STD_BF
#define STD_BF 

#include "hash_function.h"
//#include <cstring>
#include <iostream>
#include <bitset>
#include <stdlib.h>
#include "define.h"
#define BF_Q 1
using namespace std;



class StandardBF 
{
public:	
	StandardBF(unsigned int m, unsigned int k)
	{
		
		if(k > HashNum){
			cerr << "the # of hash functions cannot exceed " << HashNum << endl;
		}
		bf_m = m;
		bf_k = k;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = new uchar [bf_m/8+1];
		memset(bf_base, 0, bf_m/8+1);

		//uint (* tmp_ptr[HashNum])(const unsigned char * str, uint len) = 
		//{ JSHash, OCaml, OAAT, PJWHash, RSHash,  SDBM, Simple, SML, STL,
		//APHash, BKDR, BOB, DEKHash, DJBHash, FNV32, Hsieh};

		uint (* tmp_ptr[HashNum])(const unsigned char * str, uint len) = 
		{BOB1,BOB2,BOB3,BOB4,BOB5,BOB6,BOB7,BOB8,
		BOB9, BOB10, BOB11, BOB12, BOB13, BOB14, BOB15, BOB16, BOB, CRC32, APHash, SBOX };

		for(int i = 0; i < HashNum; i++){
			bf_hfp[i] = tmp_ptr[i];
		}
	
	}
	StandardBF(){
		bf_m = 0;
		bf_k = 0;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = NULL;
		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB1, BOB2, BOB3, BOB4, BOB5, BOB6, BOB7, BOB8,
		BOB9, BOB10, BOB11, BOB12, BOB13, BOB14, BOB15, BOB16, BOB, CRC32, APHash, SBOX };
		for(int i = 0; i < HashNum; i++){
			bf_hfp[i] = tmp_ptr[i];
		}
	}
	~StandardBF()
	{
		delete [] bf_base;
	}
	void initial(unsigned int m, unsigned int k){
		if(k > HashNum){
			cerr << "the # of hash functions cannot exceed " << HashNum << endl;
		}
		bf_m = m;
		bf_k = k;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = new uchar [bf_m/8+1];
		memset(bf_base, 0, bf_m/8+1);
	}
	void reset(){
		memset(bf_base, 0, bf_m/8+1);
	}
	unsigned int insert(const unsigned char * str, unsigned int len){
		unsigned int value;
		for(uint i = 0; i < bf_k; i++){
			value = bf_hfp[i](str, len) % bf_m;
			bf_base[value/8] |= (128 >> (value%8));
		}
		bf_n++;
		return 1;
	}
	unsigned int query(const unsigned char * str, unsigned int len){
		register unsigned int value;
		for(register uint i = 0; i < bf_k; i++)
		{
			value = bf_hfp[i](str, len) % bf_m;
			if(0 == (bf_base[value/8] & (128 >> (value%8))))return 0;
		}
		return 1;
	}

	void outputOHABF(char * filename)
	{
		FILE *fp=fopen(filename,"w");
		for (uint i=0;i<bf_m/8;i++)fprintf(fp,"%d ",bf_base[i]);
		fclose(fp);

	}

	uint Get_bf_m(){return bf_m;}
	uint Get_bf_k(){return bf_k;}
	uint Get_bf_n(){return bf_n;}

	//uint QuerymemAccNum;
	
private:
	uchar * bf_base; //bloom filter base
	
	unsigned int bf_m; //bloom filter length
	unsigned int bf_k; //hash function numbers;
	unsigned int bf_n; //# of elements inserted
	
	//pointers to hash function
	unsigned int (*bf_hfp[HashNum])(const unsigned char * str, unsigned int len);
};


class QWStandardBF
{
public:
	QWStandardBF(unsigned int m, unsigned int k,unsigned int n,unsigned int q)
	{

		if (k > HashNum){
			cerr << "the # of hash functions cannot exceed " << HashNum << endl;
		}
		bf_m = m;
		bf_k = k;
		bf_n = n;
		bf_q = q;
		bf_samew = 0;
		bf_lessw = 0;
		bf_qw_work = 0;
		//QuerymemAccNum = 0;
		bf_base = new uchar[bf_m / 8 + 1];
		memset(bf_base, 0, bf_m / 8 + 1);
		bf_w = new unsigned int[bf_k];
		memset(bf_w,0,bf_k*sizeof(int));


		//uint (* tmp_ptr[HashNum])(const unsigned char * str, uint len) = 
		//{ JSHash, OCaml, OAAT, PJWHash, RSHash,  SDBM, Simple, SML, STL,
		//APHash, BKDR, BOB, DEKHash, DJBHash, FNV32, Hsieh};

		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB1, BOB2, BOB3, BOB4, BOB5, BOB6, BOB7, BOB8,
		BOB9, BOB10, BOB11, BOB12, BOB13, BOB14, BOB15, BOB16 };

		for (int i = 0; i < HashNum; i++){
			bf_hfp[i] = tmp_ptr[i];
		}

	}
	QWStandardBF(){
		bf_m = 0;
		bf_k = 0;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = NULL;
		bf_q = BF_Q;
		bf_samew = 0;
		bf_lessw = 0;
		bf_qw_work = 0;
		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB1, BOB2, BOB3, BOB4, BOB5, BOB6, BOB7, BOB8,
		BOB9, BOB10, BOB11, BOB12, BOB13, BOB14, BOB15, BOB16 };
		for (int i = 0; i < HashNum; i++){
			bf_hfp[i] = tmp_ptr[i];
		}
	}
	~QWStandardBF()
	{
		delete[] bf_base;
		delete[] bf_w;
		
	}
	void initial(unsigned int m, unsigned int k,unsigned int n,unsigned int q){
		if (k > HashNum){
			cerr << "the # of hash functions cannot exceed " << HashNum << endl;
		}
		bf_m = m;
		bf_k = k;
		bf_n = n;
		bf_q = q;
		
		bf_samew = 0;
		bf_lessw = 0;
		bf_qw_work = 0;
		//QuerymemAccNum = 0;
		bf_base = new uchar[bf_m / 8 + 1];
		memset(bf_base, 0, bf_m / 8 + 1);
		bf_w = new unsigned int[bf_k];
		memset(bf_w, 0, bf_k * sizeof(int));
	}
	void reset(){
		memset(bf_base, 0, bf_m / 8 + 1);
	}
	
	void showwindow()
	{
		cout << "bf_w value: ";
		for (uint i = 0; i < bf_k; i++)
		{
			cout <<"w"<<i<<" "<< bf_w[i]<<" " ;
		}
		cout << endl;
	}

	unsigned int get_qw_count_window(unsigned int function_k,unsigned int start_base,unsigned char start_value)
	{
		unsigned int w_size = 1; //  from the next bit of hashfunction-value, besides, consider the condition that the last bit in bf is "0" 
		unsigned char w_value = start_value;
	    int cur_base = start_base;
		unsigned int q_count = bf_q;
		while ( cur_base >= 0)
		{
			while (w_value)
			{
				if (bf_base[cur_base] & w_value)
				{
					q_count--;
					if (q_count==0)
						break;
					w_size++;
					w_value = w_value >> 1;
				}
				else
				{
					w_size++;
					w_value=w_value >> 1;
				}			
			}	//while w_value

			
			if (q_count==0) //flag ==1, find the first "1"
				break;

			w_value = 128;  //here  reset w_value 0x10000000
			cur_base--;
				//continue;
			
		}//while cur_base
	//	if (function_k==9)
		//cout << "function k " << function_k << "  " << w_size << " bf_w " << bf_w[function_k] << endl;
			
		return w_size;
	
	}

	unsigned int insert(const unsigned char * str, unsigned int len){
		unsigned int value;
		unsigned int w_size=0;
		for (uint i = 0; i < bf_k; i++){
			value = bf_hfp[i](str, len) % bf_m;
			bf_base[value / 8] |= (128 >> (value % 8));
			
			//w_size = get_qw_count_window(i, value/8 , 128 >>((value % 8)+1));
			//bf_w[i] = bf_w[i] > w_size ? bf_w[i] : w_size;

		}
		//bf_n++;
		return 1;
	}

	unsigned int q_w_create(const unsigned char * str, unsigned int len){
		unsigned int value;
		unsigned int w_size = 0;
		for (uint i = 0; i < bf_k; i++){
			value = bf_hfp[i](str, len) % bf_m;
		//	bf_base[value / 8] |= (128 >> (value % 8));

			w_size = get_qw_count_window(i, value / 8, 128 >> ((value % 8) + 1));
			bf_w[i] = bf_w[i] > w_size ? bf_w[i] : w_size;

		}
		
		return 1;
	}

	unsigned int query(const unsigned char * str, unsigned int len){
		register unsigned int value;
		for (register uint i = 0; i < bf_k; i++)
		{
			value = bf_hfp[i](str, len) % bf_m;
			if (0 == (bf_base[value / 8] & (128 >> (value % 8))))
			{
				return 0;
			}
			else 
			{
				if (get_qw_count_window(i, value / 8, 128 >> ((value % 8) + 1)) == bf_w[i])
				{
					bf_samew++;
				}
				else if (get_qw_count_window(i, value / 8, 128 >> ((value % 8) + 1)) < bf_w[i])
				{
					bf_lessw++;
				}
				if (get_qw_count_window(i, value / 8, 128 >> ((value % 8) + 1)) >bf_w[i])
				{
				//	cout << "QW IS Working" << endl;
					bf_qw_work++;
					return 0;
				}
			}
		}
		return 1;
	}

	void outputOHABF(char * filename)
	{
		FILE *fp = fopen(filename, "w");
		for (uint i = 0; i<bf_m / 8; i++)fprintf(fp, "%d ", bf_base[i]);
		fclose(fp);

	}

	void show_bf_bit()
	{
	//	ofstream fout("bf_bit_k3_n100.txt");
		cout << "bf_bit: ";
		char bf_bit[9];
		for (int i = bf_m / 8; i >= 0; i--)
		{
			//cout << bf_base[i] << " int "<<int(bf_base[i]) << endl;
			_itoa(bf_base[i],bf_bit,2);
			cout << bf_bit << ends;
		}
		cout << endl;

	}

	uint Get_bf_m(){ return bf_m; }
	uint Get_bf_k(){ return bf_k; }
	uint Get_bf_n(){ return bf_n; }
	uint Get_bf_q(){ return bf_q; }
	uint Get_bf_qw_work(){ return bf_qw_work; }
	uint Get_bf_samew(){ return bf_samew; }
	uint Get_bf_lessw(){ return bf_lessw; }
	unsigned int * Get_bf_w(){ return bf_w;}
	unsigned char * Get_bf_base(){ return bf_base; }
	unsigned int Get_bf_base_num_of_bit_1()
	{
		unsigned int bit1num = 0;
		unsigned char value = 128;
		for (int i = bf_m / 8; i >= 0; i--)
		{
			while (value)
			{
				if (bf_base[i] & value)
					bit1num++;
				value = value >> 1;

			}
			value = 128;
		}
		return bit1num;
	
	
	}
	void show_query(){ cout <<"qw_work "<<bf_qw_work<<" same w: "<< bf_samew << " less w: " << bf_lessw << endl; }

	
	//uint QuerymemAccNum;

private:
	uchar * bf_base; //bloom filter base
	unsigned int *bf_w;  //window 
	unsigned int bf_q;  // q, number of 1

	unsigned int bf_m; //bloom filter length
	unsigned int bf_k; //hash function numbers;
	unsigned int bf_n; //# of elements inserted

	

	unsigned int bf_samew;
	unsigned int bf_lessw;
	unsigned int bf_qw_work;
	//pointers to hash function
	unsigned int(*bf_hfp[HashNum])(const unsigned char * str, unsigned int len);
};




class EXSBF
{
public:
	EXSBF(unsigned int m, unsigned int k)
	{

		if (k > HashNum){
			cerr << "the # of hash functions cannot exceed " << HashNum << endl;
		}
		bf_m = m;
		bf_k = k;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = new uchar[bf_m / 8 + 1];
		memset(bf_base, 0, bf_m / 8 + 1);

		//uint (* tmp_ptr[HashNum])(const unsigned char * str, uint len) = 
		//{ JSHash, OCaml, OAAT, PJWHash, RSHash,  SDBM, Simple, SML, STL,
		//APHash, BKDR, BOB, DEKHash, DJBHash, FNV32, Hsieh};
		
		/*
		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB16, BOB1, BOB2, BOB3, BOB4, BOB5, BOB6, BOB7, BOB8,
		BOB9, BOB10, BOB11, BOB12, BOB13, BOB14, BOB15 };
		*/
		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB16, BOB15, BOB14, BOB13, BOB12, BOB11, BOB10, BOB9,
		BOB8, BOB7, BOB6, BOB5, BOB4, BOB3, BOB2, BOB1 };
		
		for (int i = 0; i < HashNum; i++){
			bf_hfp[i] = tmp_ptr[i];
		}

	}
	EXSBF(){
		bf_m = 0;
		bf_k = 0;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = NULL;
		/*
		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB16, BOB1, BOB2, BOB3, BOB4, BOB5, BOB6, BOB7, BOB8,
		BOB9, BOB10, BOB11, BOB12, BOB13, BOB14, BOB15 };
		
		*/
		uint(*tmp_ptr[HashNum])(const unsigned char * str, uint len) =
		{ BOB16, BOB15, BOB14, BOB13, BOB12, BOB11, BOB10, BOB9,
		BOB8, BOB7, BOB6, BOB5, BOB4, BOB3, BOB2, BOB1 };
		
		
		for (int i = 0; i < HashNum; i++){
			bf_hfp[i] = tmp_ptr[i];
		}
	}
	~EXSBF()
	{
		delete[] bf_base;
	}
	void initial(unsigned int m, unsigned int k){
		if (k > HashNum){
			cerr << "the # of hash functions cannot exceed " << HashNum << endl;
		}
		bf_m = m;
		bf_k = k;
		bf_n = 0;
		//QuerymemAccNum = 0;
		bf_base = new uchar[bf_m / 8 + 1];
		memset(bf_base, 0, bf_m / 8 + 1);
	}
	void reset(){
		memset(bf_base, 0, bf_m / 8 + 1);
	}
	unsigned int insert(const unsigned char * str, unsigned int len){
		unsigned int value;
		for (uint i = 0; i < bf_k; i++){
			value = bf_hfp[i](str, len) % bf_m;
			bf_base[value / 8] |= (128 >> (value % 8));
		}
		bf_n++;
		return 1;
	}
	unsigned int query(const unsigned char * str, unsigned int len){
		register unsigned int value;
		for (register uint i = 0; i < bf_k; i++)
		{
			value = bf_hfp[i](str, len) % bf_m;
			if (0 == (bf_base[value / 8] & (128 >> (value % 8))))return 0;
		}
		return 1;
	}

	void outputOHABF(char * filename)
	{
		FILE *fp = fopen(filename, "w");
		for (uint i = 0; i<bf_m / 8; i++)fprintf(fp, "%d ", bf_base[i]);
		fclose(fp);

	}

	uint Get_bf_m(){ return bf_m; }
	uint Get_bf_k(){ return bf_k; }
	uint Get_bf_n(){ return bf_n; }

	//uint QuerymemAccNum;

private:
	uchar * bf_base; //bloom filter base

	unsigned int bf_m; //bloom filter length
	unsigned int bf_k; //hash function numbers;
	unsigned int bf_n; //# of elements inserted

	//pointers to hash function
	unsigned int(*bf_hfp[HashNum])(const unsigned char * str, unsigned int len);
};

#endif