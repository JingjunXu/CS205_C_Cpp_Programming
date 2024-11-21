#include <iostream>
#include <cstring>
#include "mul.hpp"
using namespace std;

bool printInfo(const size_t row, const size_t col, const double* data)
{
    if(row == 0 || col == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The input (" << row << ", " << col << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(data == NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    for(size_t r = 0; r < row; r++)
    {
        for(size_t c = 0; c < col; c++)
            cout << data[r * col + c] << " ";
        cout << endl;
    }
    return true;
}

bool dgemm_rck(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r++)
            for(size_t c = 0; c < N; c++)
                for(size_t k =0; k < K; k++)
                    res(r, c) += left(r, k) * right(k, c);
    }
    return true;
}

bool dgemm_rkc(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r++)
            for(size_t k =0; k < K; k++)
                for(size_t c = 0; c < N; c++)
                    res(r, c) += left(r, k) * right(k, c);
    }
    return true;
}

bool dgemm_krc(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t k =0; k < K; k++)
            for(size_t r = 0; r < M; r++)
                for(size_t c = 0; c < N; c++)
                    res(r, c) += left(r, k) * right(k, c);
    }
    return true;
}

bool dgemm_kcr(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t k =0; k < K; k++)
            for(size_t c = 0; c < N; c++)
                for(size_t r = 0; r < M; r++)
                    res(r, c) += left(r, k) * right(k, c);
    }
    return true;
}

bool dgemm_ckr(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t c = 0; c < N; c++)
            for(size_t k =0; k < K; k++)
                for(size_t r = 0; r < M; r++)
                    res(r, c) += left(r, k) * right(k, c);
    }
    return true;
}

bool dgemm_crk(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t c = 0; c < N; c++)
            for(size_t r = 0; r < M; r++)
                for(size_t k =0; k < K; k++)
                    res(r, c) += left(r, k) * right(k, c);
    }
    return true;
}

bool dgemm_rkc_4x4(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=4)
            for(size_t k =0; k < K; k+=4)
                for(size_t c = 0; c < N; c+=4)
                    for(size_t r0 = 0; r0 < 4; r0++)
                        for(size_t k0 = 0; k0 < 4; k0++)
                        {
                            for(size_t c0 = 0; c0 < 4; c0++)
                            {
                                res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
                                // cout << "("<< r + r0 << ", " << c + c0 << ")";
                                // cout << "("<< res(r + r0, c + c0) << ", " << left(r + r0, k + k0) << "," << right(k + k0, c + c0) << ")";
                                cout << "res("<< r + r0 << ", " << c + c0 << ")+=l(" << r + r0 << "," << k + k0<< ")*r(" << k + k0 << "," <<c + c0 << ")    ";
                            }
                            cout << endl;
                        }
        // cout << res(7,7)<< endl;
    }
    return true;
}

bool dgemm_rkc_16x16(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=16)
            for(size_t k =0; k < K; k+=16)
                for(size_t c = 0; c < N; c+=16)
                for(size_t r0 = 0; r0 < 16; r0++)
                    for(size_t k0 = 0; k0 < 16; k0++)
                        for(size_t c0 = 0; c0 < 16; c0++)
                            res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
    }
    return true;
}

bool dgemm_rkc_32x32_unloop4(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0+=4)
                        {
                            res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
                            res(r + r0, c + c0 + 1) += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                            res(r + r0, c + c0 + 2) += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                            res(r + r0, c + c0 + 3) += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                        }
    }
    return true;
}

bool dgemm_rkc_32x32_unloop8(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0+=8)
                        {
                            res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
                            res(r + r0, c + c0 + 1) += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                            res(r + r0, c + c0 + 2) += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                            res(r + r0, c + c0 + 3) += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                            res(r + r0, c + c0 + 4) += left(r + r0, k + k0) * right(k + k0, c + c0 + 4);
                            res(r + r0, c + c0 + 5) += left(r + r0, k + k0) * right(k + k0, c + c0 + 5);
                            res(r + r0, c + c0 + 6) += left(r + r0, k + k0) * right(k + k0, c + c0 + 6);
                            res(r + r0, c + c0 + 7) += left(r + r0, k + k0) * right(k + k0, c + c0 + 7);
                        }
    }
    return true;
}

bool dgemm_rkc_32x32_unloop16(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0+=16)
                        {
                            res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
                            res(r + r0, c + c0 + 1) += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                            res(r + r0, c + c0 + 2) += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                            res(r + r0, c + c0 + 3) += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                            res(r + r0, c + c0 + 4) += left(r + r0, k + k0) * right(k + k0, c + c0 + 4);
                            res(r + r0, c + c0 + 5) += left(r + r0, k + k0) * right(k + k0, c + c0 + 5);
                            res(r + r0, c + c0 + 6) += left(r + r0, k + k0) * right(k + k0, c + c0 + 6);
                            res(r + r0, c + c0 + 7) += left(r + r0, k + k0) * right(k + k0, c + c0 + 7);
                            res(r + r0, c + c0 + 8) += left(r + r0, k + k0) * right(k + k0, c + c0 + 8);
                            res(r + r0, c + c0 + 9) += left(r + r0, k + k0) * right(k + k0, c + c0 + 9);
                            res(r + r0, c + c0 + 10) += left(r + r0, k + k0) * right(k + k0, c + c0 + 10);
                            res(r + r0, c + c0 + 11) += left(r + r0, k + k0) * right(k + k0, c + c0 + 11);
                            res(r + r0, c + c0 + 12) += left(r + r0, k + k0) * right(k + k0, c + c0 + 12);
                            res(r + r0, c + c0 + 13) += left(r + r0, k + k0) * right(k + k0, c + c0 + 13);
                            res(r + r0, c + c0 + 14) += left(r + r0, k + k0) * right(k + k0, c + c0 + 14);
                            res(r + r0, c + c0 + 15) += left(r + r0, k + k0) * right(k + k0, c + c0 + 15);
                        }
    }
    return true;
}

bool dgemm_rkc_32x32_unloop32(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                    {
                        res(r + r0, c) += left(r + r0, k + k0) * right(k + k0, c);
                        res(r + r0, c + 1) += left(r + r0, k + k0) * right(k + k0, c + 1);
                        res(r + r0, c + 2) += left(r + r0, k + k0) * right(k + k0, c + 2);
                        res(r + r0, c + 3) += left(r + r0, k + k0) * right(k + k0, c + 3);
                        res(r + r0, c + 4) += left(r + r0, k + k0) * right(k + k0, c + 4);
                        res(r + r0, c + 5) += left(r + r0, k + k0) * right(k + k0, c + 5);
                        res(r + r0, c + 6) += left(r + r0, k + k0) * right(k + k0, c + 6);
                        res(r + r0, c + 7) += left(r + r0, k + k0) * right(k + k0, c + 7);
                        res(r + r0, c + 8) += left(r + r0, k + k0) * right(k + k0, c + 8);
                        res(r + r0, c + 9) += left(r + r0, k + k0) * right(k + k0, c + 9);
                        res(r + r0, c + 10) += left(r + r0, k + k0) * right(k + k0, c + 10);
                        res(r + r0, c + 11) += left(r + r0, k + k0) * right(k + k0, c + 11);
                        res(r + r0, c + 12) += left(r + r0, k + k0) * right(k + k0, c + 12);
                        res(r + r0, c + 13) += left(r + r0, k + k0) * right(k + k0, c + 13);
                        res(r + r0, c + 14) += left(r + r0, k + k0) * right(k + k0, c + 14);
                        res(r + r0, c + 15) += left(r + r0, k + k0) * right(k + k0, c + 15);
                        res(r + r0, c + 16) += left(r + r0, k + k0) * right(k + k0, c + 16);
                        res(r + r0, c + 17) += left(r + r0, k + k0) * right(k + k0, c + 17);
                        res(r + r0, c + 18) += left(r + r0, k + k0) * right(k + k0, c + 18);
                        res(r + r0, c + 19) += left(r + r0, k + k0) * right(k + k0, c + 19);
                        res(r + r0, c + 20) += left(r + r0, k + k0) * right(k + k0, c + 20);
                        res(r + r0, c + 21) += left(r + r0, k + k0) * right(k + k0, c + 21);
                        res(r + r0, c + 22) += left(r + r0, k + k0) * right(k + k0, c + 22);
                        res(r + r0, c + 23) += left(r + r0, k + k0) * right(k + k0, c + 23);
                        res(r + r0, c + 24) += left(r + r0, k + k0) * right(k + k0, c + 24);
                        res(r + r0, c + 25) += left(r + r0, k + k0) * right(k + k0, c + 25);
                        res(r + r0, c + 26) += left(r + r0, k + k0) * right(k + k0, c + 26);
                        res(r + r0, c + 27) += left(r + r0, k + k0) * right(k + k0, c + 27);
                        res(r + r0, c + 28) += left(r + r0, k + k0) * right(k + k0, c + 28);
                        res(r + r0, c + 29) += left(r + r0, k + k0) * right(k + k0, c + 29);
                        res(r + r0, c + 30) += left(r + r0, k + k0) * right(k + k0, c + 30);
                        res(r + r0, c + 31) += left(r + r0, k + k0) * right(k + k0, c + 31);
                    }
    }
    return true;
}


bool dgemm_rkc_32x32(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
    }
    return true;
}

bool dgemm_rkc_32x32_unloop4_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0+=4)
                            {
                                register double temp1 = res(r + r0, c + c0);
                                register double temp2 = res(r + r0, c + c0 + 1);
                                register double temp3 = res(r + r0, c + c0 + 2);
                                register double temp4 = res(r + r0, c + c0 + 3);
                                temp1 += left(r + r0, k + k0) * right(k + k0, c + c0);
                                temp2 += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                                temp3 += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                                temp4 += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                                res(r + r0, c + c0) = temp1;
                                res(r + r0, c + c0 + 1) = temp2;
                                res(r + r0, c + c0 + 2) = temp3;
                                res(r + r0, c + c0 + 3) = temp4;
                            }
    }
    return true;
}

bool dgemm_rkc_32x32_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0++)
                            {
                                register double temp = res(r + r0, c + c0);
                                temp += left(r + r0, k + k0) * right(k + k0, c + c0);
                                res(r + r0, c + c0) = temp;
                            }
    }
    return true;
}

bool dgemm_rkc_32x32_unloop8_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0+=8)
                        {
                            register double temp1 = res(r + r0, c + c0);
                            register double temp2 = res(r + r0, c + c0 + 1);
                            register double temp3 = res(r + r0, c + c0 + 2);
                            register double temp4 = res(r + r0, c + c0 + 3);
                            register double temp5 = res(r + r0, c + c0 + 4);
                            register double temp6 = res(r + r0, c + c0 + 5);
                            register double temp7 = res(r + r0, c + c0 + 6);
                            register double temp8 = res(r + r0, c + c0 + 7);
                            temp1 += left(r + r0, k + k0) * right(k + k0, c + c0);
                            temp2 += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                            temp3 += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                            temp4 += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                            temp5 += left(r + r0, k + k0) * right(k + k0, c + c0 + 4);
                            temp6 += left(r + r0, k + k0) * right(k + k0, c + c0 + 5);
                            temp7 += left(r + r0, k + k0) * right(k + k0, c + c0 + 6);
                            temp8 += left(r + r0, k + k0) * right(k + k0, c + c0 + 7);
                            res(r + r0, c + c0) = temp1;
                            res(r + r0, c + c0 + 1) = temp2;
                            res(r + r0, c + c0 + 2) = temp3;
                            res(r + r0, c + c0 + 3) = temp4;
                            res(r + r0, c + c0 + 4) = temp5;
                            res(r + r0, c + c0 + 5) = temp6;
                            res(r + r0, c + c0 + 6) = temp7;
                            res(r + r0, c + c0 + 7) = temp8;
                        }
    }
    return true;
}

bool dgemm_rkc_32x32_unloop16_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0+=8)
                        {
                            register double temp1 = res(r + r0, c + c0);
                            register double temp2 = res(r + r0, c + c0 + 1);
                            register double temp3 = res(r + r0, c + c0 + 2);
                            register double temp4 = res(r + r0, c + c0 + 3);
                            register double temp5 = res(r + r0, c + c0 + 4);
                            register double temp6 = res(r + r0, c + c0 + 5);
                            register double temp7 = res(r + r0, c + c0 + 6);
                            register double temp8 = res(r + r0, c + c0 + 7);
                            register double temp9 = res(r + r0, c + c0 + 8);
                            register double temp10 = res(r + r0, c + c0 + 9);
                            register double temp11 = res(r + r0, c + c0 + 10);
                            register double temp12 = res(r + r0, c + c0 + 11);
                            register double temp13 = res(r + r0, c + c0 + 12);
                            register double temp14 = res(r + r0, c + c0 + 13);
                            register double temp15 = res(r + r0, c + c0 + 14);
                            register double temp16 = res(r + r0, c + c0 + 15);
                            temp1 += left(r + r0, k + k0) * right(k + k0, c + c0);
                            temp2 += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                            temp3 += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                            temp4 += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                            temp5 += left(r + r0, k + k0) * right(k + k0, c + c0 + 4);
                            temp6 += left(r + r0, k + k0) * right(k + k0, c + c0 + 5);
                            temp7 += left(r + r0, k + k0) * right(k + k0, c + c0 + 6);
                            temp8 += left(r + r0, k + k0) * right(k + k0, c + c0 + 7);
                            temp9 += left(r + r0, k + k0) * right(k + k0, c + c0 + 8);
                            temp10 += left(r + r0, k + k0) * right(k + k0, c + c0 + 9);
                            temp11 += left(r + r0, k + k0) * right(k + k0, c + c0 + 10);
                            temp12 += left(r + r0, k + k0) * right(k + k0, c + c0 + 11);
                            temp13 += left(r + r0, k + k0) * right(k + k0, c + c0 + 12);
                            temp14 += left(r + r0, k + k0) * right(k + k0, c + c0 + 13);
                            temp15 += left(r + r0, k + k0) * right(k + k0, c + c0 + 14);
                            temp16 += left(r + r0, k + k0) * right(k + k0, c + c0 + 15);
                            res(r + r0, c + c0) = temp1;
                            res(r + r0, c + c0 + 1) = temp2;
                            res(r + r0, c + c0 + 2) = temp3;
                            res(r + r0, c + c0 + 3) = temp4;
                            res(r + r0, c + c0 + 4) = temp5;
                            res(r + r0, c + c0 + 5) = temp6;
                            res(r + r0, c + c0 + 6) = temp7;
                            res(r + r0, c + c0 + 7) = temp8;
                            res(r + r0, c + c0 + 8) = temp9;
                            res(r + r0, c + c0 + 9) = temp10;
                            res(r + r0, c + c0 + 10) = temp11;
                            res(r + r0, c + c0 + 11) = temp12;
                            res(r + r0, c + c0 + 12) = temp13;
                            res(r + r0, c + c0 + 13) = temp14;
                            res(r + r0, c + c0 + 14) = temp15;
                            res(r + r0, c + c0 + 15) = temp16;
                        }
    }
    return true;
}

bool dgemm_rkc_32x32_unloop32_register(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                for(size_t r0 = 0; r0 < 32; r0++)
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0+=8)
                        {
                            register double temp1 = res(r + r0, c + c0);
                            register double temp2 = res(r + r0, c + c0 + 1);
                            register double temp3 = res(r + r0, c + c0 + 2);
                            register double temp4 = res(r + r0, c + c0 + 3);
                            register double temp5 = res(r + r0, c + c0 + 4);
                            register double temp6 = res(r + r0, c + c0 + 5);
                            register double temp7 = res(r + r0, c + c0 + 6);
                            register double temp8 = res(r + r0, c + c0 + 7);
                            register double temp9 = res(r + r0, c + c0 + 8);
                            register double temp10 = res(r + r0, c + c0 + 9);
                            register double temp11 = res(r + r0, c + c0 + 10);
                            register double temp12 = res(r + r0, c + c0 + 11);
                            register double temp13 = res(r + r0, c + c0 + 12);
                            register double temp14 = res(r + r0, c + c0 + 13);
                            register double temp15 = res(r + r0, c + c0 + 14);
                            register double temp16 = res(r + r0, c + c0 + 15);
                            register double temp17 = res(r + r0, c + c0 + 16);
                            register double temp18 = res(r + r0, c + c0 + 17);
                            register double temp19 = res(r + r0, c + c0 + 18);
                            register double temp20 = res(r + r0, c + c0 + 19);
                            register double temp21 = res(r + r0, c + c0 + 20);
                            register double temp22 = res(r + r0, c + c0 + 21);
                            register double temp23 = res(r + r0, c + c0 + 22);
                            register double temp24 = res(r + r0, c + c0 + 23);
                            register double temp25 = res(r + r0, c + c0 + 24);
                            register double temp26 = res(r + r0, c + c0 + 25);
                            register double temp27 = res(r + r0, c + c0 + 26);
                            register double temp28 = res(r + r0, c + c0 + 27);
                            register double temp29 = res(r + r0, c + c0 + 28);
                            register double temp30 = res(r + r0, c + c0 + 29);
                            register double temp31 = res(r + r0, c + c0 + 30);
                            register double temp32 = res(r + r0, c + c0 + 31);
                            temp1 += left(r + r0, k + k0) * right(k + k0, c + c0);
                            temp2 += left(r + r0, k + k0) * right(k + k0, c + c0 + 1);
                            temp3 += left(r + r0, k + k0) * right(k + k0, c + c0 + 2);
                            temp4 += left(r + r0, k + k0) * right(k + k0, c + c0 + 3);
                            temp5 += left(r + r0, k + k0) * right(k + k0, c + c0 + 4);
                            temp6 += left(r + r0, k + k0) * right(k + k0, c + c0 + 5);
                            temp7 += left(r + r0, k + k0) * right(k + k0, c + c0 + 6);
                            temp8 += left(r + r0, k + k0) * right(k + k0, c + c0 + 7);
                            temp9 += left(r + r0, k + k0) * right(k + k0, c + c0 + 8);
                            temp10 += left(r + r0, k + k0) * right(k + k0, c + c0 + 9);
                            temp11 += left(r + r0, k + k0) * right(k + k0, c + c0 + 10);
                            temp12 += left(r + r0, k + k0) * right(k + k0, c + c0 + 11);
                            temp13 += left(r + r0, k + k0) * right(k + k0, c + c0 + 12);
                            temp14 += left(r + r0, k + k0) * right(k + k0, c + c0 + 13);
                            temp15 += left(r + r0, k + k0) * right(k + k0, c + c0 + 14);
                            temp16 += left(r + r0, k + k0) * right(k + k0, c + c0 + 15);
                            temp17 += left(r + r0, k + k0) * right(k + k0, c + c0 + 16);
                            temp18 += left(r + r0, k + k0) * right(k + k0, c + c0 + 17);
                            temp19 += left(r + r0, k + k0) * right(k + k0, c + c0 + 18);
                            temp20 += left(r + r0, k + k0) * right(k + k0, c + c0 + 19);
                            temp21 += left(r + r0, k + k0) * right(k + k0, c + c0 + 20);
                            temp22 += left(r + r0, k + k0) * right(k + k0, c + c0 + 21);
                            temp23 += left(r + r0, k + k0) * right(k + k0, c + c0 + 22);
                            temp24 += left(r + r0, k + k0) * right(k + k0, c + c0 + 23);
                            temp25 += left(r + r0, k + k0) * right(k + k0, c + c0 + 24);
                            temp26 += left(r + r0, k + k0) * right(k + k0, c + c0 + 25);
                            temp27 += left(r + r0, k + k0) * right(k + k0, c + c0 + 26);
                            temp28 += left(r + r0, k + k0) * right(k + k0, c + c0 + 27);
                            temp29 += left(r + r0, k + k0) * right(k + k0, c + c0 + 28);
                            temp30 += left(r + r0, k + k0) * right(k + k0, c + c0 + 29);
                            temp31 += left(r + r0, k + k0) * right(k + k0, c + c0 + 30);
                            temp32 += left(r + r0, k + k0) * right(k + k0, c + c0 + 31);
                            res(r + r0, c + c0) = temp1;
                            res(r + r0, c + c0 + 1) = temp2;
                            res(r + r0, c + c0 + 2) = temp3;
                            res(r + r0, c + c0 + 3) = temp4;
                            res(r + r0, c + c0 + 4) = temp5;
                            res(r + r0, c + c0 + 5) = temp6;
                            res(r + r0, c + c0 + 6) = temp7;
                            res(r + r0, c + c0 + 7) = temp8;
                            res(r + r0, c + c0 + 8) = temp9;
                            res(r + r0, c + c0 + 9) = temp10;
                            res(r + r0, c + c0 + 10) = temp11;
                            res(r + r0, c + c0 + 11) = temp12;
                            res(r + r0, c + c0 + 12) = temp13;
                            res(r + r0, c + c0 + 13) = temp14;
                            res(r + r0, c + c0 + 14) = temp15;
                            res(r + r0, c + c0 + 15) = temp16;
                            res(r + r0, c + c0 + 16) = temp17;
                            res(r + r0, c + c0 + 17) = temp18;
                            res(r + r0, c + c0 + 18) = temp19;
                            res(r + r0, c + c0 + 19) = temp20;
                            res(r + r0, c + c0 + 20) = temp21;
                            res(r + r0, c + c0 + 21) = temp22;
                            res(r + r0, c + c0 + 22) = temp23;
                            res(r + r0, c + c0 + 23) = temp24;
                            res(r + r0, c + c0 + 24) = temp25;
                            res(r + r0, c + c0 + 25) = temp26;
                            res(r + r0, c + c0 + 26) = temp27;
                            res(r + r0, c + c0 + 27) = temp28;
                            res(r + r0, c + c0 + 28) = temp29;
                            res(r + r0, c + c0 + 29) = temp30;
                            res(r + r0, c + c0 + 30) = temp31;
                            res(r + r0, c + c0 + 31) = temp32;
                        }
    }
    return true;
}


bool dgemm_rkc_64x64(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=64)
            for(size_t k =0; k < K; k+=64)
                for(size_t c = 0; c < N; c+=64)
                    for(size_t r0 = 0; r0 < 64; r0++)
                        for(size_t k0 = 0; k0 < 64; k0++)
                            for(size_t c0 = 0; c0 < 64; c0++)
                                res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
    }
    return true;
}

bool dgemm_rkc_128x128(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        for(size_t r = 0; r < M; r+=128)
            for(size_t k =0; k < K; k+=128)
                for(size_t c = 0; c < N; c+=128)
                    for(size_t r0 = 0; r0 < 128; r0++)
                        for(size_t k0 = 0; k0 < 128; k0++)
                            for(size_t c0 = 0; c0 < 128; c0++)
                                res(r + r0, c + c0) += left(r + r0, k + k0) * right(k + k0, c + c0);
    }
    return true;
}

bool dgemm_rkc_32x32_unloop4_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
#ifdef WITH_AVX2
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        //将right转置
        // double* D = new double[N*K];
        // for(size_t k =0; k < K; k++)
        //     for(size_t c = 0; c < N; c++)
        //         trans(c, k) = right(k, c);
        __m256d left_data;
        __m256d right_data1, right_data2, right_data3, right_data4;
        __m256d res_data1, res_data2, res_data3, res_data4;
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 8; k0++)
                            for(size_t c0 = 0; c0 < 32; c0+=4)
                            {
                                size_t s = (c / 32) * 32 + c0;
                                double*temp1=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp2=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp3=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp4=(double*)(aligned_alloc(256,4*sizeof(double)));
                                res_data1 =  _mm256_setzero_pd();
                                res_data2 =  _mm256_setzero_pd();
                                res_data3 =  _mm256_setzero_pd();
                                res_data4 =  _mm256_setzero_pd();
                                left_data = _mm256_loadu_pd(A + ((r / 32) * 32 + r0) * K + k + k0 * 4);
                                right_data1 = _mm256_loadu_pd(B + (s) * K + k0 * 4);
                                right_data2 = _mm256_loadu_pd(B + (s + 1) * K + k0 * 4);
                                right_data3 = _mm256_loadu_pd(B + (s + 2) * K + k0 * 4);
                                right_data4 = _mm256_loadu_pd(B + (s + 3) * K + k0 * 4);
                                res_data1 = _mm256_add_pd(res_data1, _mm256_mul_pd(left_data, right_data1));
                                res_data2 = _mm256_add_pd(res_data2, _mm256_mul_pd(left_data, right_data2));
                                res_data3 = _mm256_add_pd(res_data3, _mm256_mul_pd(left_data, right_data3));
                                res_data4 = _mm256_add_pd(res_data4, _mm256_mul_pd(left_data, right_data4));
                                _mm256_storeu_pd(temp1, res_data1);
                                _mm256_storeu_pd(temp2, res_data2);
                                _mm256_storeu_pd(temp3, res_data3);
                                _mm256_storeu_pd(temp4, res_data4);
                                size_t t1 = c + c0;
                                size_t t2 = r + r0;
                                for(size_t t = 0; t < 4; t++)
                                {
                                    res(t2, t1) += temp1[t];
                                    res(t2, t1 + 1) += temp2[t];
                                    res(t2, t1 + 2) += temp3[t];
                                    res(t2, t1 + 3) += temp4[t];
                                }
                                delete[] temp1;
                                delete[] temp2;
                                delete[] temp3;
                                delete[] temp4;
                            }
    }
    return true;
#else
    cout << "AVX2 is not supported" << endl;
    return false;
#endif
}


bool dgemm_rkc_32x32_unloop16_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
#ifdef WITH_AVX2
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        //将right转置
        // double* D = new double[N*K];
        // for(size_t k =0; k < K; k++)
        //     for(size_t c = 0; c < N; c++)
        //         trans(c, k) = right(k, c);
        __m256d left_data;
        __m256d right_data1, right_data2, right_data3, right_data4, right_data5, right_data6, right_data7, right_data8, 
            right_data9, right_data10, right_data11, right_data12, right_data13, right_data14, right_data15, right_data16;
        __m256d res_data1, res_data2, res_data3, res_data4, res_data5, res_data6, res_data7, res_data8, 
            res_data9, res_data10, res_data11, res_data12, res_data13, res_data14, res_data15, res_data16;
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 8; k0++)
                            for(size_t c0 = 0; c0 < 32; c0+=16)
                            {
                                size_t s = (c / 32) * 32 + c0;
                                double*temp1=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp2=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp3=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp4=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp5=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp6=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp7=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp8=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp9=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp10=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp11=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp12=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp13=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp14=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp15=(double*)(aligned_alloc(256,4*sizeof(double)));
                                double*temp16=(double*)(aligned_alloc(256,4*sizeof(double)));
                                res_data1 =  _mm256_setzero_pd();
                                res_data2 =  _mm256_setzero_pd();
                                res_data3 =  _mm256_setzero_pd();
                                res_data4 =  _mm256_setzero_pd();
                                res_data5 =  _mm256_setzero_pd();
                                res_data6 =  _mm256_setzero_pd();
                                res_data7 =  _mm256_setzero_pd();
                                res_data8 =  _mm256_setzero_pd();
                                res_data9 =  _mm256_setzero_pd();
                                res_data10 =  _mm256_setzero_pd();
                                res_data11 =  _mm256_setzero_pd();
                                res_data12 =  _mm256_setzero_pd();
                                res_data13 =  _mm256_setzero_pd();
                                res_data14 =  _mm256_setzero_pd();
                                res_data15 =  _mm256_setzero_pd();
                                res_data16 =  _mm256_setzero_pd();
                                left_data = _mm256_loadu_pd(A + ((r / 32) * 32 + r0) * K + k + k0 * 4);
                                right_data1 = _mm256_loadu_pd(B + (s) * K + k0 * 4);
                                right_data2 = _mm256_loadu_pd(B + (s + 1) * K + k0 * 4);
                                right_data3 = _mm256_loadu_pd(B + (s + 2) * K + k0 * 4);
                                right_data4 = _mm256_loadu_pd(B + (s + 3) * K + k0 * 4);
                                right_data5 = _mm256_loadu_pd(B + (s + 4) * K + k0 * 4);
                                right_data6 = _mm256_loadu_pd(B + (s + 5) * K + k0 * 4);
                                right_data7 = _mm256_loadu_pd(B + (s + 6) * K + k0 * 4);
                                right_data8 = _mm256_loadu_pd(B + (s + 7) * K + k0 * 4);
                                right_data9 = _mm256_loadu_pd(B + (s + 8) * K + k0 * 4);
                                right_data10 = _mm256_loadu_pd(B + (s + 9) * K + k0 * 4);
                                right_data11 = _mm256_loadu_pd(B + (s + 10) * K + k0 * 4);
                                right_data12 = _mm256_loadu_pd(B + (s + 11) * K + k0 * 4);
                                right_data13 = _mm256_loadu_pd(B + (s + 12) * K + k0 * 4);
                                right_data14 = _mm256_loadu_pd(B + (s + 13) * K + k0 * 4);
                                right_data15 = _mm256_loadu_pd(B + (s + 14) * K + k0 * 4);
                                right_data16 = _mm256_loadu_pd(B + (s + 15) * K + k0 * 4);
                                res_data1 = _mm256_add_pd(res_data1, _mm256_mul_pd(left_data, right_data1));
                                res_data2 = _mm256_add_pd(res_data2, _mm256_mul_pd(left_data, right_data2));
                                res_data3 = _mm256_add_pd(res_data3, _mm256_mul_pd(left_data, right_data3));
                                res_data4 = _mm256_add_pd(res_data4, _mm256_mul_pd(left_data, right_data4));
                                res_data5 = _mm256_add_pd(res_data5, _mm256_mul_pd(left_data, right_data5));
                                res_data6 = _mm256_add_pd(res_data6, _mm256_mul_pd(left_data, right_data6));
                                res_data7 = _mm256_add_pd(res_data7, _mm256_mul_pd(left_data, right_data7));
                                res_data8 = _mm256_add_pd(res_data8, _mm256_mul_pd(left_data, right_data8));
                                res_data9 = _mm256_add_pd(res_data9, _mm256_mul_pd(left_data, right_data9));
                                res_data10 = _mm256_add_pd(res_data10, _mm256_mul_pd(left_data, right_data10));
                                res_data11 = _mm256_add_pd(res_data11, _mm256_mul_pd(left_data, right_data11));
                                res_data12 = _mm256_add_pd(res_data12, _mm256_mul_pd(left_data, right_data12));
                                res_data13 = _mm256_add_pd(res_data13, _mm256_mul_pd(left_data, right_data13));
                                res_data14 = _mm256_add_pd(res_data14, _mm256_mul_pd(left_data, right_data14));
                                res_data15 = _mm256_add_pd(res_data15, _mm256_mul_pd(left_data, right_data15));
                                res_data16 = _mm256_add_pd(res_data16, _mm256_mul_pd(left_data, right_data16));
                                _mm256_storeu_pd(temp1, res_data1);
                                _mm256_storeu_pd(temp2, res_data2);
                                _mm256_storeu_pd(temp3, res_data3);
                                _mm256_storeu_pd(temp4, res_data4);
                                _mm256_storeu_pd(temp5, res_data5);
                                _mm256_storeu_pd(temp6, res_data6);
                                _mm256_storeu_pd(temp7, res_data7);
                                _mm256_storeu_pd(temp8, res_data8);
                                _mm256_storeu_pd(temp9, res_data9);
                                _mm256_storeu_pd(temp10, res_data10);
                                _mm256_storeu_pd(temp11, res_data11);
                                _mm256_storeu_pd(temp12, res_data12);
                                _mm256_storeu_pd(temp13, res_data13);
                                _mm256_storeu_pd(temp14, res_data14);
                                _mm256_storeu_pd(temp15, res_data15);
                                _mm256_storeu_pd(temp16, res_data16);
                                size_t t1 = c + c0;
                                size_t t2 = r + r0;
                                for(size_t t = 0; t < 4; t++)
                                {
                                    res(t2, t1) += temp1[t];
                                    res(t2, t1 + 1) += temp2[t];
                                    res(t2, t1 + 2) += temp3[t];
                                    res(t2, t1 + 3) += temp4[t];
                                    res(t2, t1 + 4) += temp5[t];
                                    res(t2, t1 + 5) += temp6[t];
                                    res(t2, t1 + 6) += temp7[t];
                                    res(t2, t1 + 7) += temp8[t];
                                    res(t2, t1 + 8) += temp9[t];
                                    res(t2, t1 + 9) += temp10[t];
                                    res(t2, t1 + 10) += temp11[t];
                                    res(t2, t1 + 11) += temp12[t];
                                    res(t2, t1 + 12) += temp13[t];
                                    res(t2, t1 + 13) += temp14[t];
                                    res(t2, t1 + 14) += temp15[t];
                                    res(t2, t2 + 15) += temp16[t];
                                }
                                delete[] temp1;
                                delete[] temp2;
                                delete[] temp3;
                                delete[] temp4;
                                delete[] temp5;
                                delete[] temp6;
                                delete[] temp7;
                                delete[] temp8;
                                delete[] temp9;
                                delete[] temp10;
                                delete[] temp11;
                                delete[] temp12;
                                delete[] temp13;
                                delete[] temp14;
                                delete[] temp15;
                                delete[] temp16;
                            }
    }
    return true;
#else
    cout << "AVX2 is not supported" << endl;
    return false;
#endif
}

bool dgemm_rkc_avx2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
#ifdef WITH_AVX2
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        //将right转置
        double* D = new double[N*K];
        for(size_t k =0; k < K; k++)
            for(size_t c = 0; c < N; c++)
                trans(c, k) = right(k, c);
        __m256d left_data, right_data1, right_data2, right_data3, right_data4, res_data1, res_data2, res_data3, res_data4;
        for(size_t r = 0; r < M; r++)
        {
            for(size_t k =0; k < (K / 4); k++)
            {
                for(size_t c = 0; c < N; c+=4)
                {
                    double*temp1=(double*)(aligned_alloc(256,4*sizeof(double)));
                    double*temp2=(double*)(aligned_alloc(256,4*sizeof(double)));
                    double*temp3=(double*)(aligned_alloc(256,4*sizeof(double)));
                    double*temp4=(double*)(aligned_alloc(256,4*sizeof(double)));
                    res_data1 =  _mm256_setzero_pd();
                    res_data2 =  _mm256_setzero_pd();
                    res_data3 =  _mm256_setzero_pd();
                    res_data4 =  _mm256_setzero_pd();
                    left_data = _mm256_loadu_pd(A + r * K + k * 4);
                    right_data1 = _mm256_loadu_pd(D + c * K + k * 4);
                    right_data2 = _mm256_loadu_pd(D + (c + 1) * K + k * 4);
                    right_data3 = _mm256_loadu_pd(D + (c + 2) * K + k * 4);
                    right_data4 = _mm256_loadu_pd(D + (c + 3) * K + k * 4);
                    res_data1 = _mm256_add_pd(res_data1, _mm256_mul_pd(left_data, right_data1));
                    res_data2 = _mm256_add_pd(res_data2, _mm256_mul_pd(left_data, right_data2));
                    res_data3 = _mm256_add_pd(res_data3, _mm256_mul_pd(left_data, right_data3));
                    res_data4 = _mm256_add_pd(res_data4, _mm256_mul_pd(left_data, right_data4));
                    _mm256_storeu_pd(temp1, res_data1);
                    _mm256_storeu_pd(temp2, res_data2);
                    _mm256_storeu_pd(temp3, res_data3);
                    _mm256_storeu_pd(temp4, res_data4);
                    for(size_t t = 0; t < 4; t++)
                    {
                        res(r, c) += temp1[t];
                        res(r, c + 1) += temp2[t];
                        res(r, c + 2) += temp3[t];
                        res(r, c + 3) += temp4[t];
                    }
                    delete[] temp1;
                    delete[] temp2;
                    delete[] temp3;
                    delete[] temp4;
                }
            }
        }    
    }
    return true;
#else
    cout << "AVX2 is not supported" << endl;
    return false;
#endif
}

bool dgemm_rkc_32x32_packed1(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            register double temp = left(r + r0, k + k0);
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res(r + r0, c + c0) += temp * right_packed[k0 * 32 + c0];
                        }
                }
        delete[] right_packed;
    }
    return true;
}

bool dgemm_rkc_32x32_packed2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            register double temp1 = left_packed[r0 * 32 + k0];
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res(r + r0, c + c0) += temp1 * right_packed[k0 * 32 + c0];
                        }
                    }
                }
        delete[] right_packed;
        delete[] left_packed;
    }
    return true;
}


bool dgemm_rkc_32x32_packed2_OpenMP(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        #pragma omp parallel for schedule(dynamic)
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
            {
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            register double temp1 = left_packed[r0 * 32 + k0];
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res(r + r0, c + c0) += temp1 * right_packed[k0 * 32 + c0];
                        }
                    }
                }
            }
        delete[] right_packed;
        delete[] left_packed;
    }
    return true;
}

bool dgemm_rkc_32x32_packed_OpenMP_avx2_v1(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        #pragma omp parallel for schedule(dynamic)
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
            {
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    __m256d right_data, left_data, res_data;
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            left_data = _mm256_broadcast_sd(&left_packed[r0 * 32 + k0]);
                            double* temp = (double*)aligned_alloc(32, 4 * sizeof(double));
                            for(size_t c0 = 0; c0 < 32; c0+=4)
                            {
                                res_data = _mm256_setzero_pd();
                                right_data = _mm256_loadu_pd(&right_packed[k0 * 32 + c0]);
                                res_data = _mm256_add_pd(res_data, _mm256_mul_pd(left_data, right_data));
                                _mm256_storeu_pd(temp, res_data);
                                for(int i = 0; i < 4; i++)
                                    res(r + r0, c + c0 + i) += temp[i]; 
                            }
                            delete[] temp;
                        }
                    }
                }
            }
        delete[] right_packed;
        delete[] left_packed;
    }
    return true;
}

bool dgemm_rkc_32x32_packed_OpenMP_avx2_v2(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        #pragma omp parallel for schedule(dynamic)
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
            {
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    __m256d right_data1, left_data, res_data1, right_data2, res_data2, right_data3, res_data3, right_data4, res_data4;
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            left_data = _mm256_broadcast_sd(&left_packed[r0 * 32 + k0]);
                            double* temp1 = (double*)aligned_alloc(32, 4 * sizeof(double));
                            double* temp2 = (double*)aligned_alloc(32, 4 * sizeof(double));
                            double* temp3 = (double*)aligned_alloc(32, 4 * sizeof(double));
                            double* temp4 = (double*)aligned_alloc(32, 4 * sizeof(double));
                            for(size_t c0 = 0; c0 < 32; c0+=16)
                            {
                                res_data1 = _mm256_setzero_pd();
                                res_data2 = _mm256_setzero_pd();
                                res_data3 = _mm256_setzero_pd();
                                res_data4 = _mm256_setzero_pd();
                                right_data1 = _mm256_loadu_pd(&right_packed[k0 * 32 + c0]);
                                right_data2 = _mm256_loadu_pd(&right_packed[k0 * 32 + c0 + 4]);
                                right_data3 = _mm256_loadu_pd(&right_packed[k0 * 32 + c0 + 8]);
                                right_data4 = _mm256_loadu_pd(&right_packed[k0 * 32 + c0 + 12]);
                                res_data1 = _mm256_add_pd(res_data1, _mm256_mul_pd(left_data, right_data1));
                                res_data2 = _mm256_add_pd(res_data2, _mm256_mul_pd(left_data, right_data2));
                                res_data3 = _mm256_add_pd(res_data3, _mm256_mul_pd(left_data, right_data3));
                                res_data4 = _mm256_add_pd(res_data4, _mm256_mul_pd(left_data, right_data4));
                                _mm256_storeu_pd(temp1, res_data1);
                                _mm256_storeu_pd(temp2, res_data2);
                                _mm256_storeu_pd(temp3, res_data3);
                                _mm256_storeu_pd(temp4, res_data4);
                                for(int i = 0; i < 4; i++)
                                {
                                    res(r + r0, c + c0 + i) += temp1[i]; 
                                    res(r + r0, c + c0 + i + 4) += temp2[i]; 
                                    res(r + r0, c + c0 + i + 8) += temp3[i]; 
                                    res(r + r0, c + c0 + i + 12) += temp4[i]; 
                                }
                            }
                            delete[] temp1;
                            delete[] temp2;
                            delete[] temp3;
                            delete[] temp4;
                        }
                    }
                }
            }
        delete[] right_packed;
        delete[] left_packed;
    }
    return true;
}


bool dgemm_rkc_32x32_packed5(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = new double[32*32];
        double* left_packed = new double[32*32];
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            register double temp1 = left_packed[r0 * 32 + k0];
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res(r + r0, c + c0) += temp1 * right_packed[k0 * 32 + c0];
                        }
                    }
                }
        delete[] right_packed;
        delete[] left_packed;
    }
    return true;
}

bool dgemm_rkc_32x32_packed6(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(lda != K || ldb != N)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The parameters are not fit ( " << lda << ", " << K << ") (" << ldb << ", " << N << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        __attribute__((aligned(1024))) double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        __attribute__((aligned(1024))) double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            register double temp1 = left_packed[r0 * 32 + k0];
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res(r + r0, c + c0) += temp1 * right_packed[k0 * 32 + c0];
                        }
                    }
                }
        delete[] right_packed;
        delete[] left_packed;
    }
    return true;
}

bool dgemm_rkc_32x32_packed3(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0 || (M != N) || (M != K))
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* res_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        memset(res_packed, 0, 32 * 32 * sizeof(double));
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            register double temp1 = left_packed[r0 * 32 + k0];
                            for(size_t c0 = 0; c0 < 32; c0++)
                                res_packed[r0 * 32 + c0] += temp1 * right_packed[k0 * 32 + c0];
                        }
                    }
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            res(r + r0, c + c0) += res_packed[r0 * 32 + c0];
                }  
        delete[] right_packed;
        delete[] left_packed;
        delete[] res_packed;
    }
    return true;
}

bool dgemm_rkc_32x32_packed4(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0 || (M != N) || (M != K))
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(order == RowMajor && TransA == NoTrans && TransB == NoTrans)
    {
        double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        double* res_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
        memset(res_packed, 0, 32 * 32 * sizeof(double));
        for(size_t r = 0; r < M; r+=32)
            for(size_t k =0; k < K; k+=32)
                for(size_t c = 0; c < N; c+=32)
                {
                    for(size_t k0 = 0; k0 < 32; k0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t k0 = 0; k0 < 32; k0++)
                            left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                    for(size_t r0 = 0; r0 < 32; r0++)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                        {
                            double* pright = &right_packed[k0 * 32];
                            double* pres = &res_packed[r0 * 32];
                            register double temp1 = left_packed[r0 * 32 + k0];
                            for(size_t c0 = 0; c0 < 32; c0++)
                                *(pres++) += temp1 * *(pright++);
                        }
                    }
                    double* pres = res_packed;
                    for(size_t r0 = 0; r0 < 32; r0++)
                        for(size_t c0 = 0; c0 < 32; c0++)
                            res(r + r0, c + c0) += *(pres++);
                }  
        delete[] right_packed;
        delete[] left_packed;
        delete[] res_packed;
    }
    return true;
}

double* Random(size_t row, size_t col)//随机生成0到255的double数组
{
    if(row == 0 || col == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The input (" << row << ", " << col << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    } 
    double* data = new double[row*col];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    srand((unsigned)time(NULL));
    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0; i < row*col; i++)
        *(data + i) = randomDouble(RAND_MAX) + randomInt(0, 10);
    return data;
}

bool dgemm_final(Order order, TransType TransA, TransType TransB, const size_t M, const size_t N, const size_t K, const double alpha, const double* A, const size_t lda, 
    const double* B, const size_t ldb, const double beta, double* C, const size_t ldc)
{
    if(M == 0 || N == 0 || K == 0 || lda == 0 || ldb == 0 || (M != N) || (M != K))
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the matrix are wrong ( " << M << ", " << N << ", " << K << ", " << lda << ", " << ldb << ")" << endl;
        exit(EXIT_FAILURE);
    }
    if(A == NULL || B == NULL || C == NULL)
    {
        cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "NULL pointer!" << endl;
        exit(EXIT_FAILURE);
    }
    if(M < 512)
    {
        if(TransA == NoTrans && TransB == NoTrans)
        {
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r++)
                for(size_t k =0; k < K; k++)
                    for(size_t c = 0; c < N; c++)
                        res(r, c) = beta * res(r, c) + alpha * left(r, k) * right(k, c);
        }
        else if((TransA == Trans || TransB ==  ConjTrans) && TransB == NoTrans)
        {
            #pragma omp parallel for schedule(dynamic)
            for(size_t k =0; k < K; k++)
                for(size_t r = 0; r < M; r++)
                    for(size_t c = 0; c < N; c++)
                        res(r, c) = beta * res(r, c) + alpha * left(r, k) * right(k, c);
        }
        else if(TransA == NoTrans && (TransB == Trans || TransB == ConjTrans))
        {
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r++)
                for(size_t c = 0; c < N; c++)
                    for(size_t k =0; k < K; k++)
                        res(r, c) = beta * res(r, c) + alpha * left(r, k) * right(k, c);
        }
        else if((TransA == Trans || TransB ==  ConjTrans) && (TransB == Trans || TransB == ConjTrans))
        {
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r++)
                for(size_t k =0; k < K; k++)
                    for(size_t c = 0; c < N; c++)
                        res(r, c) = beta * res(r, c) + alpha * left(r, k) * right(k, c);
        }
        else
        {
            cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
            cerr << "Wrong Matrix TransType input" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        if(TransA == NoTrans && TransB == NoTrans)
        {
            double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r+=32)
                for(size_t k =0; k < K; k+=32)
                {
                    for(size_t c = 0; c < N; c+=32)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0++)
                                right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                            for(size_t k0 = 0; k0 < 32; k0++)
                                left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                        {
                            for(size_t k0 = 0; k0 < 32; k0++)
                            {
                                register double temp1 = left_packed[r0 * 32 + k0];
                                for(size_t c0 = 0; c0 < 32; c0++)
                                    res(r + r0, c + c0) = alpha * (temp1 * right_packed[k0 * 32 + c0]) + beta * res(r + r0, c + c0);
                            }
                        }
                    }
                }
            delete[] right_packed;
            delete[] left_packed;
        }
        else if((TransA == Trans || TransB ==  ConjTrans) && TransB == NoTrans)
        {
            double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r+=32)
                for(size_t k =0; k < K; k+=32)
                {
                    for(size_t c = 0; c < N; c+=32)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0++)
                                right_packed[k0 * 32 + c0] = right(k + k0, c + c0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                            for(size_t k0 = 0; k0 < 32; k0++)
                                left_packed[r0 * 32 + k0] = left(k + k0, r + r0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                        {
                            for(size_t k0 = 0; k0 < 32; k0++)
                            {
                                register double temp1 = left_packed[r0 * 32 + k0];
                                for(size_t c0 = 0; c0 < 32; c0++)
                                    res(r + r0, c + c0) = alpha * (temp1 * right_packed[k0 * 32 + c0]) + beta * res(r + r0, c + c0);
                            }
                        }
                    }
                }
            delete[] right_packed;
            delete[] left_packed;
        }
        else if(TransA == NoTrans && (TransB == Trans || TransB == ConjTrans))
        {
            double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r+=32)
                for(size_t k =0; k < K; k+=32)
                {
                    for(size_t c = 0; c < N; c+=32)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0++)
                                right_packed[k0 * 32 + c0] = right(c + c0, k + k0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                            for(size_t k0 = 0; k0 < 32; k0++)
                                left_packed[r0 * 32 + k0] = left(r + r0, k + k0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                        {
                            for(size_t k0 = 0; k0 < 32; k0++)
                            {
                                register double temp1 = left_packed[r0 * 32 + k0];
                                for(size_t c0 = 0; c0 < 32; c0++)
                                    res(r + r0, c + c0) = alpha * (temp1 * right_packed[k0 * 32 + c0]) + beta * res(r + r0, c + c0);
                            }
                        }
                    }
                }
            delete[] right_packed;
            delete[] left_packed;
        }
        else if((TransA == Trans || TransB ==  ConjTrans) && (TransB == Trans || TransB == ConjTrans))
        {
            double* right_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            double* left_packed = (double*)aligned_alloc(1024, 32 * 32 * sizeof(double));
            #pragma omp parallel for schedule(dynamic)
            for(size_t r = 0; r < M; r+=32)
                for(size_t k =0; k < K; k+=32)
                {
                    for(size_t c = 0; c < N; c+=32)
                    {
                        for(size_t k0 = 0; k0 < 32; k0++)
                            for(size_t c0 = 0; c0 < 32; c0++)
                                right_packed[k0 * 32 + c0] = right(c + c0, k + k0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                            for(size_t k0 = 0; k0 < 32; k0++)
                                left_packed[r0 * 32 + k0] = left(k + k0, r + r0);
                        for(size_t r0 = 0; r0 < 32; r0++)
                        {
                            for(size_t k0 = 0; k0 < 32; k0++)
                            {
                                register double temp1 = left_packed[r0 * 32 + k0];
                                for(size_t c0 = 0; c0 < 32; c0++)
                                    res(r + r0, c + c0) = alpha * (temp1 * right_packed[k0 * 32 + c0]) + beta * res(r + r0, c + c0);
                            }
                        }
                    }
                }
            delete[] right_packed;
            delete[] left_packed;
        }
        else
        {
            cerr << "File " << __FILE__ << " Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
            cerr << "Wrong Matrix TransType input" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return true;
}
