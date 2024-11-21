#include <iostream>
#include <sys/time.h>//gettimeofday
#include <fstream>//读写CSV文件
#include <cstring>//使用memset
#include <cblas.h>
#include "mul.hpp"
using namespace std;

#define Col 128
#define Row 128

int main()
{
    /******************循环顺序*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "dgemm_rck" << "," << "dgemm_rkc" << "," <<"dgemm_krc" << "," << "dgemm_kcr" << "," << "dgemm_ckr" << "," << "dgemm_crk" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_rck(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_krc(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_kcr(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_ckr(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_crk(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************分块*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "dgemm_rkc" << "," << "dgemm_rkc_4x4" << "," << "dgemm_rkc_16x16" << "," << "dgemm_rkc_32x32" << "," << "dgemm_rkc_64x64" << "," << "dgemm_rkc_128x128" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_4x4(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_16x16(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_64x64(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_128x128(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************分块加解循环*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "rkc_32x32" << "," << "rkc_32x32_unloop4" << "," << "rkc_32x32_unloop8" << "," << "rkc_32x32_unloop16" << "," << "rkc_32x32_unloop32"<< endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

        // gettimeofday(&start, NULL);
        // dgemm_rkc_32x32(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
        // gettimeofday(&end, NULL);
        // file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
        // memset(C, 0, sizeof(double) * Col * Row);

        // gettimeofday(&start, NULL);
        // dgemm_rkc_32x32_unloop4(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
        // gettimeofday(&end, NULL);
        // file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
        // memset(C, 0, sizeof(double) * Col * Row);

        // gettimeofday(&start, NULL);
        // dgemm_rkc_32x32_unloop8(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
        // gettimeofday(&end, NULL);
        // file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
        // memset(C, 0, sizeof(double) * Col * Row);

        // gettimeofday(&start, NULL);
        // dgemm_rkc_32x32_unloop32(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
        // gettimeofday(&end, NULL);
        // file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
        // memset(C, 0, sizeof(double) * Col * Row);

        // gettimeofday(&start, NULL);
        // dgemm_rkc_32x32_unloop32(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
        // gettimeofday(&end, NULL);
        // file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
        // memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************分块加寄存器变量*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "rkc_32x32" << "," << "rkc_32x32_register" << "," << "rkc_32x32_unloop4_register" << "," << "rkc_32x32_unloop8_register" << "," << "rkc_32x32_unloop16_register" << "," << "rkc_32x32_unloop32_register"<< endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_register(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_unloop4_register(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_unloop8_register(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_unloop16_register(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_unloop32_register(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************矩阵打包*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "rkc_32x32" << "," << "rkc_32x32_packed1" << "," << "rkc_32x32_packed2" << "," << "rkc_32x32_packed3" << "," << "rkc_32x32_packed4" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed1(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed2(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed3(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed4(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************内存对齐*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "no aligned" << "," << "aligned" << "," << "explicit aligned" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";
    //     cout << "here" << endl;
    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed5(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed2(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed6(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************多线程*******************/
    // #ifdef _OPENMP
    //     std::cout << "OpenMP is supported" << std::endl;
    // #else
    //     std::cout << "OpenMP is not supported" << std::endl;
    // #endif
    //  ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "OpenMp" << "," << "normal" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed2_OpenMP(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed2(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    /******************avx2*******************/
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "normal" << "," << "openmp" << "," << "openmp+avx2_v1" << "," << "openmp+avx2_v2" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed2(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

        // gettimeofday(&start, NULL);
        // dgemm_rkc_32x32_packed2_OpenMP(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
        // gettimeofday(&end, NULL);
        // file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
        // memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed_OpenMP_avx2_v1(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc_32x32_packed_OpenMP_avx2_v2(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    //对比openblas与final
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "openblas" << "," << "dgemm_final" << "," << "normal" << endl;
    // struct timeval start, end;

    // double* A = Random(Col, Row);
    // double* B = Random(Col, Row);
    // double* C = new double[Col*Row]();

    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    // dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    // dgemm_rkc(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    // dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    // dgemm_rkc(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    // for(int num = 1; num <=1; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);
        
    //     gettimeofday(&start, NULL);
    //     dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_rkc(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    

    //对比转不转置对速度的影响
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();

    //对比转不转置对速度的影响
    // ofstream file("TimeCost.csv");
    // if(!file.is_open())
    // {
    //     cerr << "File " << __FILE__ << " Line " << __LINE__ << endl;
    //     cerr << "Fail to open the file" << endl;
    //     exit(EXIT_FAILURE);
    // }
    // file << Col << "x" << Row << endl;
    // file << "实验次数" << "," << "normal" << "," << "openmp" << "," << "openmp+avx2_v1" << "," << "openmp+avx2_v2" << endl;
    // struct timeval start, end;
    // for(int num = 1; num <=10; num++)
    // {
    //     double* A = Random(Col, Row);
    //     double* B = Random(Col, Row);
    //     double* C = new double[Col*Row]();
    //     file << num << ",";

    //     gettimeofday(&start, NULL);
    //     dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     gettimeofday(&start, NULL);
    //     dgemm_final(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    //     gettimeofday(&end, NULL);
    //     file << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << ",";
    //     memset(C, 0, sizeof(double) * Col * Row);

    //     file << endl;
    //     delete[] A;
    //     delete[] B;
    //     delete[] C;
    // }
    // file.close();


    /******************检测函数正确性(小规模)*******************/
    // double* A = Random(Row, Col);
    // double* B = Random(Row, Col);
    // double* C = new double[Row*Col]();
    // cout << "blas" << endl;
    // cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    // printInfo(Row, Col, C);
    // memset(C, 0, sizeof(double) * Col * Row);

    // cout << "blas" << endl;
    // cblas_dgemm(CblasRowMajor, CblasConjTrans, CblasConjTrans, Col, Row, Row, 1, B, Col, B, Col, 1, C, Col);
    // printInfo(Row, Col, C);
    // memset(C, 0, sizeof(double) * Col * Row);

    // cout << "my" << endl;
    // dgemm_final(RowMajor, Trans, Trans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    // printInfo(Row, Col, C);
    // memset(C, 0, sizeof(double) * Col * Row);

    // cout << "unloop k" << endl;
    // dgemm_rkc_unloopk(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    // printInfo(Row, Col, C);
    // memset(C, 0, sizeof(double) * Col * Row);

    // cout << "unloop r" << endl;
    // dgemm_rkc_unloopr(RowMajor, NoTrans, NoTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);
    // printInfo(Row, Col, C);
    // memset(C, 0, sizeof(double) * Col * Row);

    /******************检测函数正确性(大规模)*******************/
    double* A = Random(Col, Row);
    double* B = Random(Col, Row);
    double* C = new double[Col*Row]();

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    ofstream file1("openblas_output4.txt");
    for(size_t i = 0; i < Col*Row; i++)
        file1 << C[i] << endl;
    file1.close();

    memset(C, 0, sizeof(double) * Col * Row);

    dgemm_final(RowMajor, Trans, Trans, Col, Row, Row, 1, A, Col, B, Col, 1, C, Col);

    ofstream file2("final_output4.txt");
    for(size_t i = 0; i < Col*Row; i++)
        file2 << C[i] << endl;
    file2.close();
    
    return 0;
}