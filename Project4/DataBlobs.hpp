#ifndef _DATABLOBS_H_
#define _DATABLOBS_H_
#include <iostream>
#include <cstring>//利用memcpy
#include <opencv2/opencv.hpp>//读图片
#include <cblas.h>//矩阵乘法

#ifdef WITH_AVX2
#include <immintrin.h>
#endif 

// #ifdef WITH_NEON
// #include <arm_neon.h>
// #endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace cv;

template<typename T>
class DataBlobs
{
private:                                                                                                       
    size_t row;
    size_t col;
    size_t channel;
    T* data;
    size_t* ref_count;
    inline static size_t DataBlobs_num = 0;
public:
    //构造器
    //默认构造器
    DataBlobs();

    //完成半初始化的构造器
    DataBlobs(const size_t row, const size_t col, const size_t channel);

    //完成初始化的构造器
    DataBlobs(const size_t row, const size_t col, const size_t channel, T* data, size_t* ref_count);

    //复制构造器
    DataBlobs(const DataBlobs& blob);

    //以图片为数据的构造器
    DataBlobs(const string filename);
    //将rgb三通道的数据单独获取
    bool getRed(DataBlobs& r) const;
    bool getGreen(DataBlobs& g) const;
    bool getBlue(DataBlobs& b) const;

    //以输入rgb三个通道的数据为基础的构造器
    static bool CombineRGB(const DataBlobs& r, const DataBlobs& g, const DataBlobs& b, DataBlobs& rgb);

    //将数据以图片的形式输出
    bool output(const string filename) const;
    static bool output(const DataBlobs& blob, const string filename);

    //卷积操作
    //输入卷积对象，输入卷积核，输入输出对象
    static bool padding(const DataBlobs& input, const size_t ker_r, const size_t ker_c, DataBlobs& pad);
    static bool cnn(const DataBlobs& pad, const DataBlobs& kernel, DataBlobs& output);
    static bool cnn_unloop(const DataBlobs& pad, const DataBlobs& kernel, DataBlobs& output);
    static bool cnn_avx2_ker3x3(const DataBlobs& pad, const DataBlobs& kernel, DataBlobs& output);
    static bool cnn_neon_ker3x3(const DataBlobs& pad, const DataBlobs& kernel, DataBlobs& output);
    static bool rearrange_img2Col(const DataBlobs& input, const size_t ker_r, const size_t ker_c, DataBlobs& rearrange);
    static bool cnn_img2Col(const DataBlobs& rearrange, const size_t input_row, const size_t input_col, const DataBlobs& kernel, DataBlobs& output);
    static bool cnn_neon(const DataBlobs& pad, const DataBlobs& kernel, DataBlobs& output);

    //析构函数
    ~DataBlobs();

    //重载运算符
    //重载赋值运算符
    DataBlobs& operator=(const DataBlobs& blob);

    //重载加法运算符
    DataBlobs operator+(T num) const;
    DataBlobs operator+(const DataBlobs& blob) const;
    DataBlobs& operator+=(T num);
    DataBlobs& operator+=(const DataBlobs& blob);
    //友元函数
    template<typename Tp>
    friend DataBlobs operator+(Tp num, const DataBlobs<Tp>& blob){
        return blob + num;
    }

    //重载减法运算符
    DataBlobs operator-(T num) const;
    DataBlobs operator-(const DataBlobs& blob) const;
    DataBlobs& operator-=(T num);
    DataBlobs& operator-=(const DataBlobs& blob);
    //友元函数
    template<typename Tp>
    friend DataBlobs operator-(Tp num, const DataBlobs<Tp>& blob){
        return blob - num;
    }

    //重载除法运算
    DataBlobs operator/(T num) const;
    DataBlobs& operator/=(T num);

    //重载乘法运算
    DataBlobs operator*(T num) const;
    DataBlobs operator*(const DataBlobs& blob) const;
    DataBlobs& operator*=(T num);
    DataBlobs& operator*=(const DataBlobs& blob);
    //友元函数
    template<typename Tp>
    friend DataBlobs operator*(Tp num,const DataBlobs<Tp>& blob){
        return blob * num;
    }

    //直接访问数据块的数据
    T& operator()(size_t r, size_t c, size_t ch);

    //重载自增运算符
    //前缀
    DataBlobs& operator++();
    //后缀
    DataBlobs operator++(int);
    
    //重载自减运算符
    //前缀
    DataBlobs& operator--();
    //后缀
    DataBlobs operator--(int);

    //重载比较运算符
    //相等
    bool operator==(const DataBlobs& blob);
    //不等
    bool operator!=(const DataBlobs& blob);

    //重载<<运算符
    template<typename Tp>
    friend std::ostream& operator<<(std::ostream& , DataBlobs<Tp>& blob);

    //辅助函数
    //返回数据块的总数
    static size_t getDataBlobsNum();
    //返回数据块的行数
    const size_t getRow(){
        return this->row;
    }
    //返回数据块的列数
    const size_t getCol(){
        return this->col;
    }
    //返回数据块的通道数
    const size_t getChannel(){
        return this->channel;
    }
    //返回数据块的数据
    T* getData(){
        return this->data;
    }
    //返回数据块被引用的次数
    size_t* getRefCount(){
        return this->ref_count;
    }
    //返回这个矩阵是否为空
    bool empty(){
        if(this->row == 0||this->col == 0||this->channel == 0||this->data == NULL||
            this->ref_count == NULL)
            return true;
        else
            return false;
    }
};

template<typename T>
size_t DataBlobs<T>::getDataBlobsNum(){
    return DataBlobs_num;
}

template<typename T>
bool DataBlobs<T>::padding(const DataBlobs& input, const size_t ker_r, const size_t ker_c, DataBlobs& pad){
    if(input.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(input.row == 0 || input.col == 0 || input.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << input.row << ", " << input.col << ", " << input.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(ker_r == 0 || ker_c == 0 || ker_r / 2 == 0 || ker_c / 2 ==0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << ker_r << ", " << ker_c << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    size_t input_row = input.row; 
    size_t input_col = input.col;
    size_t pad_row = input_row + (ker_r - 1);
    size_t pad_col = input_col + (ker_c - 1);
    T* pad_data = new T[pad_row*pad_col]();//初始化为0
    if(pad_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for padding data" << endl;
        exit(EXIT_FAILURE);
    }
    size_t row_bias = (ker_r - 1) / 2;
    size_t col_bias = (ker_c - 1) / 2;
    for(int i = row_bias; i < input_row + row_bias; i++){
        memcpy(pad_data + i*pad_col + col_bias, input.data + (i - col_bias)*input_col, input_col*sizeof(T));
    }
    pad.row = pad_row;
    pad.col = pad_col;
    pad.channel = 1;
    if(pad.ref_count==NULL){
        pad.data = pad_data;
        pad.ref_count = new size_t(1);
    }
    else if(*(pad.ref_count)==1){
        delete[] pad.data;
        pad.data = pad_data;
    }
    else{
        pad.data = pad_data;
        --*pad.ref_count;
        pad.ref_count = new size_t(1);
    }
    return true;
}

template<typename T>
bool DataBlobs<T>::cnn(const DataBlobs<T>& pad, const DataBlobs<T>& kernel, DataBlobs<T>& output){
    if(pad.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(pad.row == 0 || pad.col == 0 || pad.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << pad.row << ", " << pad.col << ", " << pad.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the kernel (2nd) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.row / 2 == 0 || kernel.col / 2 ==0 || kernel.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << kernel.row << ", " << kernel.col << ", " << kernel.channel << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    size_t pad_row = pad.row; 
    size_t pad_col = pad.col;
    size_t kernel_row = kernel.row;
    size_t kernel_col = kernel.col;
    size_t input_row = pad_row - (kernel_row - 1);
    size_t input_col = pad_col - (kernel_col - 1);
    T* output_data = new T[input_row*input_col];
    if(output_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for padding data" << endl;
        exit(EXIT_FAILURE);
    }
    //卷积操作
    for(size_t t = 0; t < input_row*input_col; t++)
    {
        size_t j = t / input_col;
        size_t i = t % input_col;
        for(size_t k =0; k<kernel_col * kernel_row;k++)
        {
            *(output_data+ t) += (*(pad.data + (j * pad_col) + i + ((k / kernel_col)*pad_col) + (k % kernel_col))) * (*(kernel.data + k));
        }
    }
    //输出结果
    output.row = input_row;
    output.col = input_col;
    output.channel = 1;
    if(output.ref_count==NULL){
        output.data = output_data;
        output.ref_count = new size_t(1);
    }
    else if(*(output.ref_count)==1){
        delete[] output.data;
        output.data = output_data;
    }
    else{
        output.data = output_data;
        --*output.ref_count;
        output.ref_count = new size_t(1);
    }
    return true;
}

template<typename T>
bool DataBlobs<T>::cnn_unloop(const DataBlobs<T>& pad, const DataBlobs<T>& kernel, DataBlobs<T>& output){
    if(pad.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(pad.row == 0 || pad.col == 0 || pad.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << pad.row << ", " << pad.col << ", " << pad.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the kernel (2nd) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.row / 2 == 0 || kernel.col / 2 ==0 || kernel.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << kernel.row << ", " << kernel.col << ", " << kernel.channel << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    size_t pad_row = pad.row; 
    size_t pad_col = pad.col;
    size_t kernel_row = kernel.row;
    size_t kernel_col = kernel.col;
    size_t input_row = pad_row - (kernel_row - 1);
    size_t input_col = pad_col - (kernel_col - 1);
    T* kernel_data = kernel.data;
    T* pad_data = pad.data;
    T* output_data = new T[input_row*input_col];
    if(output_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for padding data" << endl;
        exit(EXIT_FAILURE);
    }
    //卷积操作
    for(size_t i = 0; i < input_row; i++)
    {
        size_t stride = i*input_col;
        size_t j = 0;
        for(; j < (input_col/8)*8; j+=8)
        {
            size_t out_index = stride + j;
            for(size_t m = 0; m<kernel_row;m++)
            {
                for(size_t n = 0,k = 0; n < kernel_col; n++,k++)
                {
                    size_t pad_index = (i * pad_col) + j + (m*pad_col) + n;
                    *(output_data+out_index) += (*(pad_data + pad_index)) * (*(kernel_data+k));
                    *(output_data+out_index+1) += (*(pad_data + pad_index + 1)) * (*(kernel_data+k));
                    *(output_data+out_index+2) += (*(pad_data + pad_index + 2)) * (*(kernel_data+k));
                    *(output_data+out_index+3) += (*(pad_data + pad_index + 3)) * (*(kernel_data+k));
                    *(output_data+out_index+4) += (*(pad_data + pad_index + 4)) * (*(kernel_data+k));
                    *(output_data+out_index+5) += (*(pad_data + pad_index + 5)) * (*(kernel_data+k));
                    *(output_data+out_index+6) += (*(pad_data + pad_index + 6)) * (*(kernel_data+k));
                    *(output_data+out_index+7) += (*(pad_data + pad_index + 7)) * (*(kernel_data+k));
                }
            }
        }
        if((input_col % 8) != 0)
        {
            for(int w=j;w<(input_col % 8)+j;w++)
            {
                for(int m = 0; m<kernel_row;m++)
                {
                    for(int n = 0,k = 0; n < kernel_col; n++,k++)
                    {
                        *(output_data+(stride)+w) += (*(pad_data + (i * pad_col) + w + m*pad_col + n)) * (*(kernel_data+k));
                    }
                }
            }
        }
    }
    //输出结果
    output.row = input_row;
    output.col = input_col;
    output.channel = 1;
    if(output.ref_count==NULL){
        output.data = output_data;
        output.ref_count = new size_t(1);
    }
    else if(*(output.ref_count)==1){
        delete[] output.data;
        output.data = output_data;
    }
    else{
        output.data = output_data;
        --*output.ref_count;
        output.ref_count = new size_t(1);
    }
    return true;
}

template<typename T>
bool DataBlobs<T>::cnn_avx2_ker3x3(const DataBlobs<T>& pad, const DataBlobs<T>& kernel, DataBlobs<T>& output){
    if(pad.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(pad.row == 0 || pad.col == 0 || pad.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << pad.row << ", " << pad.col << ", " << pad.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the kernel (2nd) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.row != 3|| kernel.col != 3 || kernel.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << kernel.row << ", " << kernel.col << ", " << kernel.channel << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    size_t pad_row = pad.row; 
    size_t pad_col = pad.col;
    size_t kernel_row = kernel.row;
    size_t kernel_col = kernel.col;
    size_t input_row = pad_row - (kernel_row - 1);
    size_t input_col = pad_col - (kernel_col - 1);
    size_t num = 32 / sizeof(T);
    T* kernel_data = kernel.data;
    T* input_data = pad.data;
    T* output_data = new T[input_row*input_col];
    if(output_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
#ifdef WITH_AVX2
    __m256 im11,im12,im13,im21,im22,im23,im31,im32,im33;
    __m256 ker11,ker12,ker13,ker21,ker22,ker23,ker31,ker32,ker33;
    __m256 res;
    T* Ker11=static_cast<T*>(aligned_alloc(256,32));
    T* Ker12=static_cast<T*>(aligned_alloc(256,32));
    T* Ker13=static_cast<T*>(aligned_alloc(256,32));
    T* Ker21=static_cast<T*>(aligned_alloc(256,32));
    T* Ker22=static_cast<T*>(aligned_alloc(256,32));
    T* Ker23=static_cast<T*>(aligned_alloc(256,32));
    T* Ker31=static_cast<T*>(aligned_alloc(256,32));
    T* Ker32=static_cast<T*>(aligned_alloc(256,32));
    T* Ker33=static_cast<T*>(aligned_alloc(256,32));
    for(int i=0;i<num;i++)
    {
        *(Ker11+i)=*(kernel_data);
        *(Ker12+i)=*(kernel_data+1);
        *(Ker13+i)=*(kernel_data+2);
        *(Ker21+i)=*(kernel_data+3);
        *(Ker22+i)=*(kernel_data+4);
        *(Ker23+i)=*(kernel_data+5);
        *(Ker31+i)=*(kernel_data+6);
        *(Ker32+i)=*(kernel_data+7);
        *(Ker33+i)=*(kernel_data+8);
    }
    ker11 = _mm256_loadu_ps(Ker11);
    ker12 = _mm256_loadu_ps(Ker12);
    ker13 = _mm256_loadu_ps(Ker13);
    ker21 = _mm256_loadu_ps(Ker21);
    ker22 = _mm256_loadu_ps(Ker22);
    ker23 = _mm256_loadu_ps(Ker23);
    ker31 = _mm256_loadu_ps(Ker31);
    ker32 = _mm256_loadu_ps(Ker32);
    ker33 = _mm256_loadu_ps(Ker33);
    T*temp=static_cast<T*>(aligned_alloc(256,32));
    for(size_t i=0; i<input_row;i++)
    {
        size_t j=0;
        for(; j<(input_col/num)*num; j+=num)
        {
            res =  _mm256_setzero_ps();
            im11 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (0*pad_col) + 0);
            im12 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (0*pad_col) + 1);
            im13 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (0*pad_col) + 2);
            im21 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (1*pad_col) + 0);
            im22 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (1*pad_col) + 1);
            im23 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (1*pad_col) + 2);
            im31 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (2*pad_col) + 0);
            im32 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (2*pad_col) + 1);
            im33 = _mm256_loadu_ps(input_data + (i * pad_col) + j + (2*pad_col) + 2);
            res = _mm256_add_ps(res, _mm256_mul_ps(im11, ker11));
            res = _mm256_add_ps(res, _mm256_mul_ps(im12, ker12));
            res = _mm256_add_ps(res, _mm256_mul_ps(im13, ker13));
            res = _mm256_add_ps(res, _mm256_mul_ps(im21, ker21));
            res = _mm256_add_ps(res, _mm256_mul_ps(im22, ker22));
            res = _mm256_add_ps(res, _mm256_mul_ps(im23, ker23));
            res = _mm256_add_ps(res, _mm256_mul_ps(im31, ker31));
            res = _mm256_add_ps(res, _mm256_mul_ps(im32, ker32));
            res = _mm256_add_ps(res, _mm256_mul_ps(im33, ker33));
            _mm256_storeu_ps(temp, res);
            for(size_t k=0;k<num;k++)
                *(output_data+(i*input_col)+j+k) = *(temp+k);
        }
        if((input_col % num) != 0)
        {
            for(size_t w=j;w<(input_col % num)+j;w++)
            {
                for(size_t k =0; k<kernel_col*kernel_row;k++)
                {
                    *(output_data+(i*input_col)+w) += (*(input_data + (i * pad_col) + w + (k / kernel_col)*pad_col + (k % kernel_col))) * (*(kernel_data+k));
                }
            }
        }
    }
    output.row = input_row;
    output.col = input_col;
    output.channel = 1;
    if(output.ref_count==NULL){
        output.data = output_data;
        output.ref_count = new size_t(1);
    }
    else if(*(output.ref_count)==1){
        delete[] output.data;
        output.data = output_data;
    }
    else{
        output.data = output_data;
        --*output.ref_count;
        output.ref_count = new size_t(1);
    }
    return true;
#else
    cerr << "AVX2 is not supported" << endl;
    return false;
#endif
}

template<typename T>
bool DataBlobs<T>::cnn_neon_ker3x3(const DataBlobs<T>& pad, const DataBlobs<T>& kernel, DataBlobs<T>& output){
    if(pad.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(pad.row == 0 || pad.col == 0 || pad.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << pad.row << ", " << pad.col << ", " << pad.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the kernel (2nd) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.row != 3|| kernel.col != 3 || kernel.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << kernel.row << ", " << kernel.col << ", " << kernel.channel << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    size_t pad_row = pad.row; 
    size_t pad_col = pad.col;
    size_t kernel_row = kernel.row;
    size_t kernel_col = kernel.col;
    size_t input_row = pad_row - (kernel_row - 1);
    size_t input_col = pad_col - (kernel_col - 1);
    size_t num = 16 / sizeof(T);
    T* kernel_data = kernel.data;
    T* input_data = pad.data;
    T* output_data = new T[input_row*input_col];
    if(output_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << " Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
#ifdef WITH_NEON
    float32x4_t im11,im12,im13,im21,im22,im23,im31,im32,im33;
    float32x4_t ker11,ker12,ker13,ker21,ker22,ker23,ker31,ker32,ker33;
    float32x4_t res;
    T* Ker11=static_cast<T*>(aligned_alloc(128,16));
    T* Ker12=static_cast<T*>(aligned_alloc(128,16));
    T* Ker13=static_cast<T*>(aligned_alloc(128,16));
    T* Ker21=static_cast<T*>(aligned_alloc(128,16));
    T* Ker22=static_cast<T*>(aligned_alloc(128,16));
    T* Ker23=static_cast<T*>(aligned_alloc(128,16));
    T* Ker31=static_cast<T*>(aligned_alloc(128,16));
    T* Ker32=static_cast<T*>(aligned_alloc(128,16));
    T* Ker33=static_cast<T*>(aligned_alloc(128,16));
    for(int i=0;i<num;i++)
    {
        *(Ker11+i)=*(kernel_data);
        *(Ker12+i)=*(kernel_data+1);
        *(Ker13+i)=*(kernel_data+2);
        *(Ker21+i)=*(kernel_data+3);
        *(Ker22+i)=*(kernel_data+4);
        *(Ker23+i)=*(kernel_data+5);
        *(Ker31+i)=*(kernel_data+6);
        *(Ker32+i)=*(kernel_data+7);
        *(Ker33+i)=*(kernel_data+8);
    }
    ker11 = vld1q_f32(Ker11);
    ker12 = vld1q_f32(Ker12);
    ker13 = vld1q_f32(Ker13);
    ker21 = vld1q_f32(Ker21);
    ker22 = vld1q_f32(Ker22);
    ker23 = vld1q_f32(Ker23);
    ker31 = vld1q_f32(Ker31);
    ker32 = vld1q_f32(Ker32);
    ker33 = vld1q_f32(Ker33);
    T*temp=static_cast<T*>(aligned_alloc(128,16));
    for(size_t i=0; i<input_row;i++)
    {
        size_t j=0;
        for(; j<(input_col/num)*num; j+=num)
        {
            res =  vdupq_n_f32(0);
            im11 = vld1q_f32(input_data + (i * pad_col) + j + (0*pad_col) + 0);
            im12 = vld1q_f32(input_data + (i * pad_col) + j + (0*pad_col) + 1);
            im13 = vld1q_f32(input_data + (i * pad_col) + j + (0*pad_col) + 2);
            im21 = vld1q_f32(input_data + (i * pad_col) + j + (1*pad_col) + 0);
            im22 = vld1q_f32(input_data + (i * pad_col) + j + (1*pad_col) + 1);
            im23 = vld1q_f32(input_data + (i * pad_col) + j + (1*pad_col) + 2);
            im31 = vld1q_f32(input_data + (i * pad_col) + j + (2*pad_col) + 0);
            im32 = vld1q_f32(input_data + (i * pad_col) + j + (2*pad_col) + 1);
            im33 = vld1q_f32(input_data + (i * pad_col) + j + (2*pad_col) + 2);
            res = vaddq_f32(res, vmulq_f32(im11, ker11));
            res = vaddq_f32(res, vmulq_f32(im12, ker12));
            res = vaddq_f32(res, vmulq_f32(im13, ker13));
            res = vaddq_f32(res, vmulq_f32(im21, ker21));
            res = vaddq_f32(res, vmulq_f32(im22, ker22));
            res = vaddq_f32(res, vmulq_f32(im23, ker23));
            res = vaddq_f32(res, vmulq_f32(im31, ker31));
            res = vaddq_f32(res, vmulq_f32(im32, ker32));
            res = vaddq_f32(res, vmulq_f32(im33, ker33));
            vst1q_f32(temp, res);
            for(size_t k=0;k<num;k++)
                *(output_data+(i*input_col)+j+k) = *(temp+k);
        }
        if((input_col % num) != 0)
        {
            for(size_t w=j;w<(input_col % num)+j;w++)
            {
                for(size_t k =0; k<kernel_col*kernel_row;k++)
                {
                    *(output_data+(i*input_col)+w) += (*(input_data + (i * pad_col) + w + (k / kernel_col)*pad_col + (k % kernel_col))) * (*(kernel_data+k));
                }
            }
        }
    }
    output.row = input_row;
    output.col = input_col;
    output.channel = 1;
    if(output.ref_count==NULL){
        output.data = output_data;
        output.ref_count = new size_t(1);
    }
    else if(*(output.ref_count)==1){
        delete[] output.data;
        output.data = output_data;
    }
    else{
        output.data = output_data;
        --*output.ref_count;
        output.ref_count = new size_t(1);
    }
    return true;
#else
    cerr << "NEON is not supported" << endl;
    return 0.0;
#endif
}

template<typename T>
bool DataBlobs<T>::rearrange_img2Col(const DataBlobs& input, const size_t ker_r, const size_t ker_c, DataBlobs& rearrange){
    if(input.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(input.row == 0 || input.col == 0 || input.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << input.row << ", " << input.col << ", " << input.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(ker_r / 2 == 0 || ker_c / 2 ==0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << ker_r << ", " << ker_c << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    size_t input_row = input.row; 
    size_t input_col = input.col;
    size_t pad_row = input_row + (ker_r - 1);
    size_t pad_col = input_col + (ker_c - 1);
    size_t re_row = pad_row * pad_col;
    size_t re_col = ker_r * ker_c;
    T* pad_data = new T[pad_row*pad_col]();//初始化为0
    if(pad_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for padding data" << endl;
        exit(EXIT_FAILURE);
    }
    T* re_data = new T[re_row*re_col];
    if(re_data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for padding data" << endl;
        delete[] pad_data;
        exit(EXIT_FAILURE);
    }
    size_t row_bias = (ker_r - 1) / 2;
    size_t col_bias = (ker_c - 1) / 2;
    //padding
    for(int i = row_bias; i < input_row + row_bias; i++){
        memcpy(pad_data + i*pad_col + col_bias, input.data + (i - col_bias)*input_col, input_col*sizeof(T));
    }
    //rearrange data
    size_t t = 0;
    for(size_t i = 0; i < input_row; i++){
        for(size_t j = 0; j < input_col; j++){
            for(size_t m = 0; m < ker_r; m++){
                for(size_t n = 0; n < ker_c; n++, t++){
                    re_data[t] = *(pad_data + i*pad_col + j + m*ker_c + n);
                }
            }
        }
    }
    rearrange.row = re_row;
    rearrange.col = re_col;
    rearrange.channel = 1;
    if(rearrange.ref_count==NULL){
        rearrange.data = re_data;
        rearrange.ref_count = new size_t(1);
    }
    else if(*(rearrange.ref_count)==1){
        delete[] rearrange.data;
        rearrange.data = re_data;
    }
    else{
        rearrange.data = re_data;
        --*rearrange.ref_count;
        rearrange.ref_count = new size_t(1);
    }
    delete[] pad_data;
    return true;
}

template<typename T>
bool DataBlobs<T>::cnn_img2Col(const DataBlobs& rearrange, const size_t input_row, const size_t input_col, const DataBlobs& kernel, DataBlobs& output){
    if(rearrange.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input (1st) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(rearrange.row == 0 || rearrange.col == 0 || rearrange.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the input (1st) parameter is (" << rearrange.row << ", " << rearrange.col << ", " << rearrange.channel << ") invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the kernel (2nd) parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(kernel.row / 2 == 0 || kernel.col / 2 ==0 || kernel.channel != 1){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the kernel (2nd) parameter is (" << kernel.row << ", " << kernel.col << ", " << kernel.channel << ") invalid (The size must be odd)" << endl;
        exit(EXIT_FAILURE);
    }
    DataBlobs<T> ker((kernel.row)*(kernel.col),1,1,kernel.data,kernel.ref_count);
    DataBlobs<T> res;
    res = rearrange * ker;
    output.row = input_row;
    output.col = input_col;
    output.channel = 1;
    if(output.ref_count==NULL){
        output.data = res.data;
        output.ref_count = res.ref_count;
        ++*(output.ref_count);
    }
    else if(*(output.ref_count)==1){
        delete[] output.data;
        output.data = res.data;
        delete[] output.ref_count;
        output.ref_count = res.ref_count;
        ++*(output.ref_count);
    }
    else{
        output.data = res.data;
        --*output.ref_count;
        output.ref_count = res.ref_count;
        ++*(output.ref_count);
    }
    return true;
}

//默认构造器
template<typename T>
DataBlobs<T>::DataBlobs(){
    // cout << "Constructor() is invoked" << endl;
    this->row = 0;
    this->col = 0;
    this->channel = 0;
    this->data = NULL;
    this->ref_count = NULL;
    DataBlobs_num++;
}

//半初始化构造器
template<typename T>
DataBlobs<T>::DataBlobs(const size_t row, const size_t col, const size_t channel){
    // cout << "Constructor(row,col,channel) is invoked" << endl;
    if(row == 0 || col == 0 || channel == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(row,col,channel)" << endl;
        cerr << "The input ( "<< row << ", "<< col << ", " << channel << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[row*col*channel];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(row,col,channel)" << endl;
        cerr << "Fail to allocate memory for the data" << endl;
        exit(EXIT_FAILURE);
    }
    memset(data, 0, sizeof(T)*row*col*channel);
    this->row = row;
    this->col = col;
    this->channel = channel;
    this->data = data;
    this->ref_count = new size_t(1);
    DataBlobs_num++;
}

//完成初始化的构造器
template<typename T>
DataBlobs<T>::DataBlobs(const size_t row, const size_t col, const size_t channel, T* data, size_t* ref_count){
    // cout << "Constructor(row,col,channel,data,ref_count) is invoked" << endl;
    if(row == 0 || col == 0 || channel == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(row,col,channel,data,ref_count)" << endl;
        cerr << "The input ( "<< row << ", "<< col << ", " << channel << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(row,col,channel,data,ref_count)" << endl;
        cerr << "The 4th parameter has no valid data" << endl;
        exit(EXIT_FAILURE);
    }
    if(ref_count==NULL)
        this->ref_count = new size_t(1);
    else
    {
        this->ref_count = ref_count;
        ++*this->ref_count;
    }
    this->row = row;
    this->col = col;
    this->channel = channel;
    this->data = data;
    DataBlobs_num++;
}

//复制构造器
template<typename T>
DataBlobs<T>::DataBlobs(const DataBlobs<T>& blob){
    // cout << "Copy Constructor(blob) is invoked" << endl;
    this->row = blob.row;
    this->col = blob.col;
    this->channel = blob.channel;
    this->data = blob.data;
    if((blob.ref_count)==NULL)
        this->ref_count = NULL;
    else
    {
        this->ref_count = blob.ref_count;
        ++*this->ref_count;
    }
    DataBlobs_num++;
}

//以图片为数据的构造器
template<typename T>
DataBlobs<T>::DataBlobs(const string filename){//mian里面调用的时候要手动进行显示的实例化
    // cout << "Constructor(filename) is invoked" << endl;
    Mat image = imread(filename);
    if(image.empty()){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(filename)" << endl;
        cerr << "The parameter is invalid or Fail to load the image" << endl;
        exit(EXIT_FAILURE);
    }
    if((image.channels() != 1) && (image.channels() != 3)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(filename)" << endl;
        cerr << "The channel of the image is" << image.channels() << ", which is invalid (only valid for 1 or 3)" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(image.rows)*(image.cols)*(image.channels())];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "(filename)" << endl;
        cerr << "Fail to allocate memory for the data" << endl;
        exit(EXIT_FAILURE);
    }
    this->row = image.rows;
    this->col = image.cols;
    this->channel = image.channels();
    size_t row_stride = (this->col)*(this->channel);
    if(this->channel == 1)
        for(int i = 0; i < this->row; i++)
            for(int j = 0; j  < this->col; j++)
                data[i*row_stride + j] =(T)image.at<uchar>(i,j);
    else if(this->channel == 3)//将BGR通道转换储存为RGB通道
        for(int i = 0; i < this->row; i++)
            for(int j = 0; j < this-> col; j++)
            {
                data[i*row_stride + j*3] = static_cast<T>(image.at<Vec3b>(i,j)[2]);
                data[i*row_stride + j*3 + 1] = static_cast<T>(image.at<Vec3b>(i,j)[1]);
                data[i*row_stride + j*3 + 2] = static_cast<T>(image.at<Vec3b>(i,j)[0]);
            }
    this->data = data;
    this->ref_count = new size_t(1);
    DataBlobs_num++;
}

//将rgb三通道单独获取
template<typename T>
bool DataBlobs<T>::getRed(DataBlobs<T>& r) const{
    if(this->channel!=3){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The channel of the Datablobs is " << this->channel <<"(only valid for channel 3)" << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((this->row == 0) || (this->col == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the Datablobs is (" << this->row << ", " << this->col << "), which is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    int row_stride = (this->col)*(this->channel);
    for(int i = 0; i < this->row; i++)
        for(int j = 0; j < this-> col; j++)
            data[i*this->col + j] = *(this->data + i*row_stride + j*3);
    if((r.ref_count!=NULL)&&(--*r.ref_count)==0)
    {
        delete[] r.data;
        delete[] r.ref_count;
    }
    r.row = this->row;
    r.col = this->col;
    r.channel = 1;
    r.data = data;
    r.ref_count = new size_t(1);
    return true;
}

template<typename T>
bool DataBlobs<T>::getGreen(DataBlobs<T>& g) const{
    if(this->channel!=3){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The channel of the Datablobs is " << this->channel <<"(only valid for channel 3)" << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((this->row == 0) || (this->col == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the Datablobs is (" << this->row << ", " << this->col << "), which is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    int row_stride = (this->col)*(this->channel);
    for(int i = 0; i < this->row; i++)
        for(int j = 0; j < this-> col; j++)
            data[i*this->col + j] = *(this->data + i*row_stride + j*3 + 1);
    if((g.ref_count!=NULL)&&(--*g.ref_count)==0)
    {
        delete[] g.data;
        delete[] g.ref_count;
    }
    g.row = this->row;
    g.col = this->col;
    g.channel = 1;
    g.data = data;
    g.ref_count = new size_t(1);
    return true;
}

template<typename T>
bool DataBlobs<T>::getBlue(DataBlobs& b) const{
    if(this->channel!=3){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The channel of the Datablobs is " << this->channel <<"(only valid for channel 3)" << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((this->row == 0) || (this->col == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the Datablobs is (" << this->row << ", " << this->col << "), which is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    int row_stride = (this->col)*(this->channel);
    for(int i = 0; i < this->row; i++)
        for(int j = 0; j < this-> col; j++)
            data[i*this->col + j] = *(this->data + i*row_stride + j*3 + 2);
    if((b.ref_count!=NULL)&&(--*b.ref_count)==0)
    {
        delete[] b.data;
        delete[] b.ref_count;
    }
    b.row = this->row;
    b.col = this->col;
    b.channel = 1;
    b.data = data;
    b.ref_count = new size_t(1);
    return true;
}

//合并rgb三通道的数据
template<typename T>
bool DataBlobs<T>::CombineRGB(const DataBlobs<T>& r, const DataBlobs<T>& g, const DataBlobs<T>& b,DataBlobs<T>& rgb){
    if((r.row != g.row)||(g.row != b.row)||(r.col != g.col)||(g.col != b.col)||
        (r.row == 0)||(r.col == 0)||(g.row == 0)||(g.col == 0)||(b.row == 0)||(b.col == 0))
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The sizes of the r ("<< r.row << ", " << r.col <<"), g ("<< g.row << ", " << g.col <<"), b ("<< b.row << ", " << b.col <<") Datablobs are invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((r.channel != 1)||(g.channel != 1)||(b.channel != 1)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The channels of the r (" << r.channel << "), g (" << g.channel << "), b (" << b.channel << ") Datablobs are invalid (they must be 1)" << endl;
        exit(EXIT_FAILURE);
    }
    if((r.data == NULL)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the 1st parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((g.data == NULL)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the 2nd parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((b.data == NULL)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the 3nd parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(r.row)*(r.col)*3];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    if((rgb.ref_count!=NULL)&&(--*rgb.ref_count)==0)
    {
        delete[] b.data;
        delete[] b.ref_count;
    }
    rgb.row = r.row;
    rgb.col = r.col;
    rgb.channel = 3;
    size_t row_stride = (rgb.col) * 3;
    size_t row_stride_ = (rgb.col);
    for(int i = 0; i < rgb.row; i++)
        for(int j = 0; j < rgb.col; j++){
            data[i*row_stride + j*3] = *(r.data + i*row_stride_ + j);
            data[i*row_stride + j*3 + 1] = *(g.data + i*row_stride_ + j);
            data[i*row_stride + j*3 + 2] = *(b.data + i*row_stride_ + j);
        }
    rgb.data = data;
    rgb.ref_count = new size_t(1);
    return true;
}

//将数据以图片的形式输出
template<typename T>
bool DataBlobs<T>::output(const string filename) const{
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((this->channel != 1) && (this->channel != 3)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The channel of the image is "<< this->channel <<", which is invalid(must be 1 or 3)" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned char* data = new unsigned char[(this->row)*(this->col)*(this->channel)];
    if(data == NULL){
        cerr << "In output() Function" << endl;
        cerr << "Fail to allocate memory" << endl;
        exit(EXIT_FAILURE);
    }
    size_t row_stride = (this->col)*(this->channel);
    if(this->channel == 3){
        for(int i = 0; i < this->row; i++)//将RGB转化成BGR
            for(int j = 0; j < this-> col; j++)
            {
                data[i*row_stride + j*3] = static_cast<unsigned char>(*(this->data + i*row_stride + j*3 + 2));
                data[i*row_stride + j*3 + 1] = static_cast<unsigned char>(*(this->data + i*row_stride + j*3 + 1));
                data[i*row_stride + j*3 + 2] = static_cast<unsigned char>(*(this->data + i*row_stride + j*3));
            }
        Mat image(this->row, this->col, CV_8UC3, data);
        imwrite(filename, image);
    }
    else if(this->channel == 1){
        for(int i = 0; i < (this->row)*(this->col); i++)
            data[i] = static_cast<unsigned char>(*(this->data + i));
        Mat image(this->row, this->col, CV_8UC1, data);
        imwrite(filename, image);
    }
    delete[] data;
    return true;
}

template<typename T>
bool DataBlobs<T>::output(const DataBlobs& blob, const string filename){
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if((blob.channel != 1) && (blob.channel != 3)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The channel of the image is "<< blob.channel <<", which is invalid(must be 1 or 3)" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned char* data = new unsigned char[(blob.row)*(blob.col)*(blob.channel)];
    if(data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    size_t row_stride = (blob.col)*(blob.channel);
    if(blob.channel == 3){
        for(int i = 0; i < blob.row; i++)//将RGB转化成BGR
            for(int j = 0; j < blob.col; j++)
            {
                data[i*row_stride + j*3] = static_cast<unsigned char>(*(blob.data + i*row_stride + j*3 + 2));
                data[i*row_stride + j*3 + 1] = static_cast<unsigned char>(*(blob.data + i*row_stride + j*3 + 1));
                data[i*row_stride + j*3 + 2] = static_cast<unsigned char>(*(blob.data + i*row_stride + j*3));
            }
        Mat image(blob.row, blob.col, CV_8UC3, data);
        imwrite(filename, image);
    }
    else if(blob.channel == 1){
        for(int i = 0; i < (blob.row)*(blob.col); i++)
            data[i] = static_cast<unsigned char>(*(blob.data + i));
        Mat image(blob.row, blob.col, CV_8UC1, data);
        imwrite(filename, image);
    }
    delete[] data;
    return true;
}

//析构函数
template<typename T>
DataBlobs<T>::~DataBlobs(){
    // cout << "Destructor() is invoked" << endl;
    DataBlobs_num--;
    if(ref_count!=NULL&&(--*ref_count)==0)
    {
        delete[] this->data;
        this->data = NULL;
        delete[] this->ref_count;
        this->ref_count = NULL;
    }
}

//运算符重载
//重载赋值运算符
template<typename T>
DataBlobs<T>& DataBlobs<T>::operator=(const DataBlobs<T>& blob){
    if(this == &blob)
        return *this;
    this->row = blob.row;
    this->col = blob.col;
    this->channel = blob.channel;
    if(this->ref_count == NULL){
        this->data = blob.data;
        this->ref_count = blob.ref_count;
    }
    else if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = blob.data;
        delete[] this->ref_count;
        this->ref_count = blob.ref_count;
    }
    else{
        this->data = blob.data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = blob.ref_count;
    }
    if(this->ref_count != NULL)
        ++*this->ref_count;
    return *this;
}

//重载加法运算符
template<typename T>
DataBlobs<T> DataBlobs<T>::operator+(T num) const{//需要创建中间量来传递数值，但是中间量只起到传递数值的作用，它的析构函数不会被自动调用
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] + num;
    DataBlobs B(this->row, this->col, this->channel, data, NULL);
    return B;
}

template<typename T>
DataBlobs<T> DataBlobs<T>::operator+(const DataBlobs<T>& blob) const{
    if((this->row != blob.row)||(this->col != blob.col)||(this->channel != blob.channel)||
        (this->row == 0)||(this->col == 0)||(this->channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the two Datablobs (" << this->row << ", " << this->col <<", " << this->channel <<") ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") should be the same." << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] + *(blob.data + i);
    DataBlobs B(this->row, this->col, this->channel, data, NULL);
    return B;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator+=(T num){
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] + num;
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator+=(const DataBlobs<T>& blob){
    if((this->row != blob.row)||(this->col != blob.col)||(this->channel != blob.channel)||
        (this->row == 0)||(this->col == 0)||(this->channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the two Datablobs (" << this->row << ", " << this->col <<", " << this->channel <<") ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") should be the same." << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] + *(blob.data+i);
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}

//重载减法运算符
template<typename T>
DataBlobs<T> DataBlobs<T>::operator-(T num) const{
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] - num;
    DataBlobs B(this->row, this->col, this->channel, data, NULL);
    return B;
}

template<typename T>
DataBlobs<T> DataBlobs<T>::operator-(const DataBlobs<T>& blob) const{
    if((this->row != blob.row)||(this->col != blob.col)||(this->channel != blob.channel)||
        (this->row == 0)||(this->col == 0)||(this->channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the two Datablobs (" << this->row << ", " << this->col <<", " << this->channel <<") ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") should be the same." << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] - *(blob.data + i);
    DataBlobs B(this->row, this->col, this->channel, data, NULL);
    return B;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator-=(T num){
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] - num;
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator-=(const DataBlobs<T>& blob){
    if((this->row != blob.row)||(this->col != blob.col)||(this->channel != blob.channel)||
        (this->row == 0)||(this->col == 0)||(this->channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the two Datablobs (" << this->row << ", " << this->col <<", " << this->channel <<") ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") should be the same." << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] - *(blob.data+i);
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}

//重载除法运算符
template<typename T>
DataBlobs<T> DataBlobs<T>::operator/(T num) const{
    if(num == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The input parameter cannot be 0" << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] / num;
    DataBlobs B(this->row, this->col, this->channel, data, NULL);
    return B;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator/=(T num){
    if(num == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The input parameter cannot be 0" << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] / num;
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}

//重载乘法运算符
template<typename T>
DataBlobs<T> DataBlobs<T>::operator*(T num) const{
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] * num;
    DataBlobs B(this->row, this->col, this->channel, data, NULL);
    return B;
}

template<typename T>
DataBlobs<T> DataBlobs<T>::operator*(const DataBlobs& blob) const{//使用openblas进行计算
    if((this->col != blob.row)||(this->row == 0)||(this->col == 0)||(blob.col == 0)||(this->channel != 1)||(blob.channel != 1)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the two Datablobs (" << this->row << ", " << this->col <<", " << this->channel <<") ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") should be the same." << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(blob.col)];
    float* left = new float[(this->col)*(this->row)];
    float* right = new float[(blob.row)*(blob.col)];
    float* res = new float[(this->row)*(blob.col)];
    if((data == NULL)||(left == NULL)||(right == NULL)||(res == NULL)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->col)*(this->row); i++)
        left[i] = static_cast<float>(*(this->data + i));
    for(int i = 0; i < (blob.row)*(blob.col); i++)
        right[i] = static_cast<float>(*(blob.data + i));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->row, blob.col, this->col, 1.0f, left, this->col, right, blob.col, 0, res, blob.col);
    for(int i = 0; i < (this->row)*(blob.col); i++)
        data[i] = static_cast<T>(res[i]);
    DataBlobs B(this->row, blob.col, 1, data, NULL);
    delete[] left;
    delete[] right;
    delete[] res;
    return B;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator*=(T num){
    if(this->data==NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(this->col)*(this->channel)];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->row)*(this->col)*(this->channel); i++)
        data[i] = this->data[i] * num;
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}

template<typename T>
DataBlobs<T>& DataBlobs<T>::operator*=(const DataBlobs& blob){
    if((this->col != blob.row)||(this->row == 0)||(this->col == 0)||(blob.col == 0)||(this->channel != 1)||(blob.channel != 1)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the two Datablobs (" << this->row << ", " << this->col <<", " << this->channel <<") ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") should be the same." << endl;
        exit(EXIT_FAILURE);
    }
    if(this->data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the input parameter is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    T* data = new T[(this->row)*(blob.col)];
    float* left = new float[(this->col)*(this->row)];
    float* right = new float[(blob.row)*(blob.col)];
    float* res = new float[(this->row)*(blob.col)];
    if((data == NULL)||(left == NULL)||(right == NULL)||(res == NULL)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (this->col)*(this->row); i++)
        left[i] = static_cast<float>(*(this->data + i));
    for(int i = 0; i < (blob.row)*(blob.col); i++)
        right[i] = static_cast<float>(*(blob.data + i));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->row, blob.col, this->col, 1.0f, left, this->col, right, blob.col, 0, res, blob.col);
    for(int i = 0; i < (this->row)*(blob.col); i++)
        data[i] = static_cast<T>(res[i]);
    if(*(this->ref_count)==1){
        delete[] this->data;
        this->data = data;
    }
    else{
        this->data = data;
        --*this->ref_count;
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    delete[] left;
    delete[] right;
    delete[] res;
    return *this;
}

//直接访问数据块数据
template<typename T>
T& DataBlobs<T>::operator()(size_t r, size_t c, size_t ch){
    if((r >= this->row) || (c >= this->col) || (ch > this->channel)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The index is out of range" << endl;
        exit(EXIT_FAILURE);
    }
    return this->data[r * this->col * this->channel + c * this->channel + (ch-1)];
}

//重载自增运算符
//前缀
template<typename T>
DataBlobs<T>& DataBlobs<T>::operator++(){
    if(this->data == NULL || this->ref_count == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((this->row == 0)||(this->col == 0)||(this->channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the Datablobs ( "<< row << ", "<< col << ", " << channel << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    size_t num = (this->row)*(this->col)*(this->channel);
    if(*(this->ref_count) == 1){
        for(int i = 0; i < num; i++)
            *(this->data + i) = *(this->data + i) + 1;
    }
    else{
        T* data = new T[num];
        if(data == NULL){
            cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
            cerr << "Fail to allocate memory for the data" << endl;
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < num; i++)
            data[i] = *(this->data + i) + 1;
        this->data = data;
        --*(this->ref_count);
        this->ref_count == NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}
//后缀
template<typename T>
DataBlobs<T> DataBlobs<T>::operator++(int){
    DataBlobs old = *this;
    operator++();
    return old;
}
//重载自减运算符
//前缀
template<typename T>
DataBlobs<T>& DataBlobs<T>::operator--(){
    if(this->data == NULL || this->ref_count == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((this->row == 0)||(this->col == 0)||(this->channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the Datablobs ( "<< row << ", "<< col << ", " << channel << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    size_t num = (this->row)*(this->col)*(this->channel);
    if(*(this->ref_count) == 1){
        for(int i = 0; i < num; i++)
            *(this->data + i) = *(this->data + i) - 1;
    }
    else{
        T* data = new T[num];
        if(data == NULL){
            cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
            cerr << "Fail to allocate memory for the data" << endl;
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < num; i++)
            data[i] = *(this->data + i) - 1;
        this->data = data;
        --*(this->ref_count);
        this->ref_count = NULL;
        this->ref_count = new size_t(1);
    }
    return *this;
}
//后缀
template<typename T>
DataBlobs<T> DataBlobs<T>::operator--(int){
    DataBlobs old = *this;
    operator--();
    return old;
}

//重载比较运算符
//相等
template<typename T>
bool DataBlobs<T>::operator==(const DataBlobs& blob){
    if(this->row != blob.row)
        return false;
    if(this->col != blob.col)
        return false;
    if(this->channel != blob.channel)
        return false;
    size_t num = (this->row)*(this->col)*(this->channel);
    for(int i=0; i<num; i++)
        if(*(this->data + i) != *(blob.data + i))
            return false;
    return true;
}

template<typename T>
bool DataBlobs<T>::operator!=(const DataBlobs& blob){
    if(*this == blob)
        return false;
    else
        return true;
}

template<typename T>
std::ostream& operator<<(std::ostream & os, DataBlobs<T>& blob){
    if(blob.data == NULL){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The data of the Datablobs is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if((blob.row == 0)||(blob.col == 0)||(blob.channel == 0)){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The size of the Datablobs ("<< blob.row << ", " << blob.col << ", " << blob.channel <<") is invalid." << endl;
        exit(EXIT_FAILURE);
    }
    os << "row: " << blob.row << " col: " << blob.col << " channel: " << blob.channel << endl;
    os << "ref_count: " << *blob.ref_count << " "<< blob.ref_count << endl;
    for(size_t c = 1; c <= blob.channel; c++){
        for(size_t i = 0; i < blob.row; i++){
            for(size_t j = 0; j < blob.col; j++)
                os << blob(i, j, c) << " ";
            os << endl;
        }
        os << endl;
    }
    return os;
}
#endif