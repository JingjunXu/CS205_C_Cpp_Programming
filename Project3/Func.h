#pragma once
#include <stdbool.h>//使用bool数据类型增强可读性
typedef struct {
    size_t width;
    size_t height;
    size_t channels;
    float* data1;
    float* data2;
    float* data3;
} Image;

typedef struct {
    size_t width;
    size_t height;
    size_t channels;
    float* data1;
    float* data2;
    float* data3;
} Kernel;

//传入数据
bool input_jpeg_image(const char* filename, Image* image);
bool output_jpeg_image(const char* filename, Image* image);
//初始化卷积核
bool kernel_init(Kernel* kernel, const size_t height, const size_t width, const size_t ch, float* data1, float* data2, float* data3);
//填充0使得使得输出的图像与输入的图像尺寸一样
float* padding(const float* input_data, const size_t input_height, const size_t input_width, const size_t kernel_height, const size_t kernel_weight);
//cnn函数
//四重循环
bool cnnFunction_v1(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);//input的应该是经过padding的
//三重循环
bool cnnFunction_v2(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnnFunction_v3(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
//二重循环
bool cnnFunction_v4(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
//解循环
bool cnn_unloop(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnn_omp_unloop(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
//SIMD:512
bool cnn_avx2_kernel3x3_512(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
//SIMD:256
bool cnn_avx2_kernel1x1(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnn_avx2_omp_kernel1x1(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnn_avx2_kernel3x3(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnn_avx2_omp_kernel3x3(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnn_avx2_kernel5x5(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
bool cnn_avx2_omp_kernel5x5(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
//使用openblas
bool cnn_openblas(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data);
/********************总的函数************************/
bool Convolution(const Image* input_image, const Kernel* kernel, Image* output_image);
//辅助函数
bool printinfo(const float* data, const size_t height, const size_t width, const char* name);
float* Random(const size_t height, const size_t width);
bool test(bool (*cnn)(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data),const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data,
    const size_t kernel_height,const size_t kernel_width,FILE* file,const char* name);
bool Free_Image(Image* image);