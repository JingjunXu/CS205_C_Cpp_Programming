#include <iostream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <time.h>
#include <stdio.h>
#include <string>
#include <fstream>
//获取随机数的计算式 
#define randomInt(a,b)       (rand()% (b - a) + (a) )           //取得[a,b)的随机整数
#define randomDouble(RAND_MAX) (rand() / double(RAND_MAX))      //取得0～1之间的随机双精度浮点数
#define randomFloat(RAND_MAX) (rand() / float(RAND_MAX))        //取得0~1之间的随机单精度浮点数
using namespace std;
class CreateRandomNum{
public:
vector<int> create_rand_int(int a, int b,int num)//取得[a,b)的随机整数组
{
  vector<int> res;
  if(num<=0)
  {
    cout << "The input of the number is wrong" << endl;
    return res;
  }
  for(int i = 1 ; i <= num; i++)
  res.push_back(randomInt(a,b));
  return res;
}
vector<float> create_rand_float(int a, int b,int num)//取得(a,b+1)之间的随机单精度浮点数组
{
  vector<float> res;
  if(num<=0)
  {
    cout << "The input of the number is wrong" << endl;
    return res;
  }
  for(int i=1; i<= num; i++)
  res.push_back((randomFloat(RAND_MAX) + randomInt(a,b)));
  return res;
}
//生成随机双精度浮点数数组与生成随机单精度浮点数数组的方法大同小异
vector<double> create_rand_double(int a, int b,int num)//取得(a,b+1)之间的随机双精度浮点数组
{
  vector<double> res;
  if(num<=0)
  {
    cout << "The input of the number is wrong" << endl;
    return res;
  }
  for(int i = 1 ; i <= num; i++)
  res.push_back((randomDouble(RAND_MAX) + randomInt(a,b)));
  return res;
}
};
class DotProduct{
public:
long long dot_product_int (const vector<int> & num1, const vector<int> & num2){
  if(num1.size()!=num2.size())
  {
    cout << "The size of the vectors are not fit!" << endl;
    return 0;
  }
  long long result=0;
  for(int i=0; i<num1.size();i++)
  result += num1[i] * num2[i];
  return result;
}
float dot_product_float (const vector<float> & num1, const vector<float> & num2){
  if(num1.size()!=num2.size())
  {
    cout << "The size of the vectors are not fit!" << endl;
    return 0.0;
  }
  float result=0.0;
  for(int i=0; i<num1.size();i++)
  result += num1[i] * num2[i];
  return result;
}
long double dot_product_double (const vector<double> & num1, const vector<double> & num2){
  if(num1.size()!=num2.size())
  {
    cout << "The size of the vectors are not fit!" << endl;
    return 0.0;
  }
  long double result= 0.0;
  for(int i=0; i<num1.size();i++)
  result += num1[i] * num2[i];
  return result;
}
};
void SaveDataInt(const vector<int> &num, const string fileName);//fileName的后缀为.bin
void SaveDataFloat(const vector<float> &num, const string fileName);//fileName的后缀为.bin
void SaveDataDouble(const vector<double> &num, const string fileName);//fileName的后缀为.bin
int main()
{
  /*************************实验1*************************/
  // fstream csvFile("Exp1_intTimeCostInC.csv", ios::out);
  // // Check if file is open
  // if (!csvFile.is_open()) {
  //   cout << "Error: Unable to open file." << endl;
  // }
  // else
  // {
  //   csvFile << "intTimeCost" << endl;
  //   CreateRandomNum obj;
  //   vector<int> res_int1, res_int2;
  //   clock_t int_start,int_finish;
  //   long long int_result=0;
  //   DotProduct op;
  //   for(int j=1;j<=10;j++)
  //   {
  //     srand((int)time(0));//利用系统时钟，产生不同的随机数种子
  //     res_int1=obj.create_rand_int(0,100,1000);
  //     res_int2=obj.create_rand_int(0,100,1000);
  //     int_result=0;
  //     //进行点乘并计算耗费的时间
  //     int_start=clock();
  //     int_result=op.dot_product_int(res_int1,res_int2);
  //     int_finish=clock();
  //     csvFile << (double)(int_finish-int_start) / CLOCKS_PER_SEC<< endl;
  //   }
  // }
  // // Close CSV file
  // csvFile.close();

  /*************************实验2*************************/
  // CreateRandomNum obj;
  // vector<int> res_int1, res_int2;
  // vector<float> res_float1, res_float2;
  // vector<double> res_double1, res_double2;
  // //获取随机数组
  // srand((int)time(0));//利用系统时钟，产生不同的随机数种子
  // res_int1=obj.create_rand_int(0, 10000, 10);
  // res_int2=obj.create_rand_int(0, 10000, 10);
  // //数据类型转换
  // for(auto t : res_int1)
  // {
  //   res_float1.push_back(static_cast<float>(t));
  //   res_double1.push_back(static_cast<double>(t));
  // }
  // for(auto t : res_int2)
  // {
  //   res_float2.push_back(static_cast<float>(t));
  //   res_double2.push_back(static_cast<double>(t));
  // }
  // //进行点乘并计算耗费的时间
  // DotProduct op;
  // long long int_result;
  // long double float_result, double_result;
  // clock_t int_start, int_finish, float_start, float_finish, double_start, double_finish;

  // int_start=clock();
  // int_result=op.dot_product_int(res_int1,res_int2);
  // int_finish=clock();

  // float_start=clock();
  // float_result=op.dot_product_float(res_float1,res_float2);
  // float_finish=clock();

  // double_start=clock();
  // double_result=op.dot_product_double(res_double1,res_double2);
  // double_finish=clock();
  // //分别输出结果以及所耗费的时间
  // cout << "The Result of int dot product is "<< int_result<<" Time Cost: " << (double)(int_finish-int_start) / CLOCKS_PER_SEC << " s" << endl;
  // cout << "The Result of float dot product is "<< float_result<<" Time Cost: " << (double)(float_finish-float_start) / CLOCKS_PER_SEC << " s" << endl;
  // cout << "The Result of double dot product is "<< double_result<<" Time Cost: " << (double)(double_finish-double_start) / CLOCKS_PER_SEC << " s" << endl;


  /*************************实验3*************************/
  // CreateRandomNum obj;
  // DotProduct op;
  // vector<int> res_int1,res_int2;
  // vector<float> res_float1,res_float2;
  // vector<double> res_double1,res_double2;
  // long long int_result;
  // long double float_result, double_result;
  // clock_t int_start, int_finish, float_start, float_finish, double_start, double_finish;
  // int_result=0,float_result=0,double_result=0;
  // res_int1=obj.create_rand_int(0,100,1000);
  // res_int2=obj.create_rand_int(0,100,1000);
  // res_float1=obj.create_rand_float(0,100,1000);
  // res_float2=obj.create_rand_float(0,100,1000);
  // res_double1=obj.create_rand_double(0,100,1000);
  // res_double2=obj.create_rand_double(0,100,1000);
  // //将随机生成的数据储存进bin文件，用于传给java使用 （实验5与6）
  // SaveDataInt(res_int1,"intdata1.bin");
  // SaveDataInt(res_int2,"intdata2.bin");
  // SaveDataFloat(res_float1,"floatdata1.bin");
  // SaveDataFloat(res_float2,"floatdata2.bin");
  // SaveDataDouble(res_double1,"doubledata1.bin");
  // SaveDataDouble(res_double2,"doubledata2.bin");
  // //进行点乘并计算耗费的时间
  // int_start=clock();
  // int_result=op.dot_product_int(res_int1,res_int2);
  // int_finish=clock();

  // float_start=clock();
  // float_result=op.dot_product_float(res_float1,res_float2);
  // float_finish=clock();

  // double_start=clock();
  // double_result=op.dot_product_double(res_double1,res_double2);
  // double_finish=clock();
  // //分别输出结果以及所耗费的时间
  // cout << "The Result of int dot product is "<< int_result<<" Time Cost: " << (double)(int_finish-int_start) / CLOCKS_PER_SEC << " s" << endl;
  // cout << "The Result of float dot product is "<< float_result<<" Time Cost: " << (double)(float_finish-float_start) / CLOCKS_PER_SEC << " s" << endl;
  // cout << "The Result of double dot product is "<< double_result<<" Time Cost: " << (double)(double_finish-double_start) / CLOCKS_PER_SEC << " s" << endl;

  
  /*************************实验4与实验5（加-O3编译）*************************/
  // fstream csvFile("Exp4_TimeCostInC.csv", ios::out);//实验5时候改一下名称
  // // Check if file is open
  // if (!csvFile.is_open()) {
  //   cout << "Error: Unable to open file." << endl;
  // }
  // else
  // {
  //   csvFile << "实验次数"<<","<<"int" << ","<<"float"<<","<<"double"<<endl;
  //   CreateRandomNum obj;
  //   DotProduct op;
  //   vector<int> res_int1,res_int2;
  //   vector<float> res_float1,res_float2;
  //   vector<double> res_double1,res_double2;
  //   long long int_result;
  //   long double float_result, double_result;
  //   clock_t int_start, int_finish, float_start, float_finish, double_start, double_finish;
  //   res_int1=obj.create_rand_int(0,100,100000);
  //   res_int2=obj.create_rand_int(0,100,100000);
  //   res_float1=obj.create_rand_float(0,100,100000);
  //   res_float2=obj.create_rand_float(0,100,100000);
  //   res_double1=obj.create_rand_double(0,100,100000);
  //   res_double2=obj.create_rand_double(0,100,100000);
  //   for(int j=1;j<=25;j++)
  //   {
  //     int_result=0,float_result=0,double_result=0;
  //     //进行点乘并计算耗费的时间
  //     int_start=clock();
  //     int_result=op.dot_product_int(res_int1,res_int2);
  //     int_finish=clock();

  //     float_start=clock();
  //     float_result=op.dot_product_float(res_float1,res_float2);
  //     float_finish=clock();

  //     double_start=clock();
  //     double_result=op.dot_product_double(res_double1,res_double2);
  //     double_finish=clock();

  //     csvFile <<j<<","<< (double)(int_finish-int_start) / CLOCKS_PER_SEC<<","<<(double)(float_finish-float_start) / CLOCKS_PER_SEC<<","<<(double)(double_finish-double_start) / CLOCKS_PER_SEC<<endl;
  //   }
  // }
  // // Close CSV file
  // csvFile.close();
  
  /*************************实验6*************************/
  CreateRandomNum obj;
  DotProduct op;
  vector<int> res_int1,res_int2;
  vector<float> res_float1,res_float2;
  vector<double> res_double1,res_double2;
  long long int_result;
  long double float_result, double_result;
  clock_t int_start, int_finish, float_start, float_finish, double_start, double_finish;
  res_int1=obj.create_rand_int(0,100,100000000);
  res_int2=obj.create_rand_int(0,100,100000000);
  res_float1=obj.create_rand_float(0,100,100000000);
  res_float2=obj.create_rand_float(0,100,100000000);
  res_double1=obj.create_rand_double(0,100,100000000);
  res_double2=obj.create_rand_double(0,100,100000000);
  //将随机生成的数据储存进bin文件，用于传给java使用
  SaveDataInt(res_int1,"intdata1.bin");
  SaveDataInt(res_int2,"intdata2.bin");
  SaveDataFloat(res_float1,"floatdata1.bin");
  SaveDataFloat(res_float2,"floatdata2.bin");
  SaveDataDouble(res_double1,"doubledata1.bin");
  SaveDataDouble(res_double2,"doubledata2.bin");
  //进行点乘并计算耗费的时间
  int_start=clock();
  int_result=op.dot_product_int(res_int1,res_int2);
  int_finish=clock();

  float_start=clock();
  float_result=op.dot_product_float(res_float1,res_float2);
  float_finish=clock();

  double_start=clock();
  double_result=op.dot_product_double(res_double1,res_double2);
  double_finish=clock();
  
  //分别输出结果以及所耗费的时间
  cout << "The Result of int dot product is "<< int_result<<" Time Cost: " << (double)(int_finish-int_start) / CLOCKS_PER_SEC << " s" << endl;
  cout << "The Result of float dot product is "<< float_result<<" Time Cost: " << (double)(float_finish-float_start) / CLOCKS_PER_SEC << " s" << endl;
  cout << "The Result of double dot product is "<< double_result<<" Time Cost: " << (double)(double_finish-double_start) / CLOCKS_PER_SEC << " s" << endl;

  return 0;
}
void SaveDataInt(const vector<int> &num, const string fileName)
{
  FILE *fp = fopen( fileName.data() , "wb" );
  if( fp )
  {
    fwrite( &num[0], num.size() * sizeof(int), 1, fp );
    fclose( fp );
    fp = NULL;
  }
  else
  cout << "The file is not available" << endl;
}
void SaveDataFloat(const vector<float> &num, const string fileName)
{
  FILE *fp = fopen( fileName.data() , "wb" );
  if( fp )
  {
    fwrite( &num[0], num.size() * sizeof(float), 1, fp );
    fclose( fp );
    fp = NULL;
  }
  else
  cout << "The file is not available" << endl;
}
void SaveDataDouble(const vector<double> &num, const string fileName)
{
  FILE *fp = fopen( fileName.data() , "wb" );
  if( fp )
  {
    fwrite( &num[0], num.size() * sizeof(double), 1, fp );
    fclose( fp );
    fp = NULL;
  }
  else
  cout << "The file is not available" << endl;
}