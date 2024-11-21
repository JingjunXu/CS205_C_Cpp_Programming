#ifndef _MUL_H_
#define _MUL_H_
#include <iostream>
#include <cblas.h>

#ifdef WITH_AVX2
#include <immintrin.h>
#endif 

#ifdef WITH_NEON
#include <arm_neon.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#define left(i, j) A[i * K + j]
#define right(i, j) B[i * N + j]
#define res(i, j) C[i * N + j] 
#define trans(i, j) D[i * K + j] 

#define randomInt(a,b)       (rand()% (b - a) + (a) )           //取得[a,b)的随机整数
#define randomDouble(RAND_MAX) (rand() / double(RAND_MAX))      //取得0～1之间的随机双精度浮点数

typedef enum {RowMajor, ColMajor} Order;
typedef enum {NoTrans, Trans, ConjTrans} TransType;

bool dgemm_rck(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_krc(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_kcr(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_ckr(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_crk(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rck_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_neon(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_4x4(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_16x16(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_64x64(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_128x128(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop4(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop8(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop16(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop32(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop4_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop16_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop4_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop8_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop16_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_unloop32_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed1(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed3(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed4(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed5(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed6(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed2_OpenMP(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed_OpenMP_avx2_v1(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_rkc_32x32_packed_OpenMP_avx2_v2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
bool dgemm_final(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc);
double* Random(size_t row, size_t col);
bool printInfo(const size_t row, const size_t col, const double* data);
#endif