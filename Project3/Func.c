#include <stdio.h>
#include <stdlib.h>
#include <string.h>//使用memcpy函数来提高padding的效率
#include <jpeglib.h>//读取jpg
#include <stdbool.h>//使用bool数据类型增强可读性
#include <sys/time.h>//gettimeofday
#include <time.h>//获取随机数种子
#include <cblas.h>
#include "Func.h"

#ifdef WITH_AVX2
#include <immintrin.h>
#endif 

#ifdef _OPENMP
#include <omp.h>
#endif

bool input_jpeg_image(const char* filename, Image* image)
{
    if(image==NULL) {
        printf("Error Image(structure) input!");
        return false;
    }
    //打开JPEG文件，创建JPEG解码器对象与JPEG错误对象
    FILE *infile;
    if ((infile = fopen(filename, "rb")) == NULL) {
        printf("Error opening JPEG file!");
        return false;
    }
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    // 设置基本的参数
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);
    //开始解压
    jpeg_start_decompress(&cinfo);

    //将数据传入结构体image
    int row_stride = cinfo.output_width * cinfo.output_components;
    unsigned char *buffer = (unsigned char*)malloc(row_stride*sizeof(unsigned char));
    if(cinfo.output_components==1)
    {
        cinfo.out_color_space = JCS_GRAYSCALE; //设置输出颜色空间为灰度图像
        cinfo.quantize_colors = FALSE;
        float *Data = (float*)aligned_alloc(256,row_stride*cinfo.output_height*sizeof(float));
        while (cinfo.output_scanline < cinfo.output_height) {
            int i = cinfo.output_scanline;
            jpeg_read_scanlines(&cinfo, &buffer, 1);
            for(int j=0;j<row_stride;j++)
                Data[i*row_stride+j]=(float)buffer[j];
        }
        image->data1 = Data;
        image->data2 = NULL;
        image->data3 = NULL;
    }
    else if(cinfo.output_components==3)
    {
        cinfo.out_color_space = JCS_RGB; //设置输出颜色空间为RGB
        // 分别存储每个通道的像素值
        float* red_channel = (float*)aligned_alloc(256,cinfo.output_width * cinfo.output_height * sizeof(float));
        float* green_channel = (float*)aligned_alloc(256,cinfo.output_width * cinfo.output_height * sizeof(float));
        float* blue_channel = (float*)aligned_alloc(256,cinfo.output_width * cinfo.output_height * sizeof(float));
        while (cinfo.output_scanline < cinfo.output_height)
        {
            jpeg_read_scanlines(&cinfo, &buffer, 1);
            int offset = (cinfo.output_scanline - 1) * cinfo.output_width;
            for (int i = 0, j = 0; i < row_stride; i += 3, j++) {
                red_channel[offset + j] =(float) buffer[i];
                green_channel[offset + j] =(float) buffer[i + 1];
                blue_channel[offset + j] =(float) buffer[i + 2];
            }
        }
        image->data1 = red_channel;
        image->data2 = green_channel;
        image->data3 = blue_channel;
    }
    else
    {
        printf("The channels of the image must be 1 or 3");
        return false;
    }
    image->width = cinfo.output_width;
    image->height = cinfo.output_height;
    image->channels = cinfo.output_components;

    //完成解压进程，并释放所有相关资源
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    free(buffer);

    return true;
}

bool output_jpeg_image(const char* filename, Image* image)
{
    if(image==NULL) {
        printf("Error Image(structure) input!");
        return false;
    }
    //创建JPEG压缩器对象,JPEG错误对象与打开文件
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    FILE* outfile;
    if ((outfile = fopen(filename, "wb")) == NULL) {
        printf("Error opening JPEG file!");
        return false;
    }

    //将输出流与压缩对象关联起来
    jpeg_stdio_dest(&cinfo, outfile);

    //设置基本参数
    cinfo.image_width = image->width;
    cinfo.image_height = image->height;
    cinfo.input_components = image->channels;
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 90, TRUE);
    //启动压缩进程
    jpeg_start_compress(&cinfo, TRUE);

    //将image数据传出到压缩对象中
    JSAMPROW row_pointer[1];
    int row_stride = cinfo.image_width * cinfo.input_components;
    unsigned char *buffer = (char*)malloc(row_stride*sizeof(unsigned char));
    if(image->channels==1)
    {
        cinfo.in_color_space = JCS_GRAYSCALE;
        while (cinfo.next_scanline < cinfo.image_height) {
            int i = cinfo.next_scanline;
            for(int j=0;j<row_stride;j++)
                buffer[j]=(unsigned char)image->data1[i*row_stride+j];
            row_pointer[0] = &buffer[0];
            jpeg_write_scanlines(&cinfo, row_pointer, 1);
        }
    }
    else if(image->channels==3)
    {
        cinfo.in_color_space = JCS_RGB;
        while (cinfo.next_scanline < cinfo.image_height) {
            int offset = cinfo.next_scanline * cinfo.image_width;
            for (int i = 0, j = 0; i < row_stride; i += 3, j++) {
                buffer[i] = (unsigned char)image->data1[offset + j];
                buffer[i + 1] =(unsigned char)image->data2[offset + j];
                buffer[i + 2] =(unsigned char)image->data3[offset + j];
            }
            row_pointer[0] = &buffer[0];
            jpeg_write_scanlines(&cinfo, row_pointer, 1);
        }
    }
    else
    {
        printf("The channels of the image must be 1 or 3");
        return false;
    }
    //完成压缩进程，并释放所有相关资源
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
    
    return true;
}

float* padding(const float* input_data, const size_t input_height, const size_t input_width, const size_t kernel_height, const size_t kernel_width)
{
    if(input_data==NULL)
    {
        printf("NULL pointer!");
        return FALSE;
    }
    size_t pad_width = input_width + (kernel_width - 1);
    size_t pad_height = input_height + (kernel_height - 1);
    float* pad_data=(float*)calloc(pad_width*pad_height,sizeof(float));
    if(pad_data==NULL)
    {
        printf("Error Calloc for the padding data!");
        return NULL;
    }
    for(int i = ((kernel_width - 1) / 2); i < input_height + ((kernel_width - 1) / 2); i++){
        memcpy(pad_data + i * pad_width  + ((kernel_width - 1) / 2) , input_data + (i - ((kernel_width - 1) / 2)) * input_width , input_width  * sizeof(float));
    }
    return pad_data;
}

bool kernel_init(Kernel* kernel, const size_t height, const size_t width, const size_t ch, float* Data1, float* Data2, float* Data3)//如果channel等于3的时候,data2与data3为NULL
{
    if(ch==1&&(Data1!=NULL))
    {
        kernel->data1 = Data1;
        kernel->data2 = NULL;
        kernel->data3 = NULL;
    }
    else if(ch==3&&(Data1!=NULL)&&(Data2!=NULL)&&(Data3!=NULL))
    {
        kernel->data1 = Data1;
        kernel->data2 = Data2;
        kernel->data3 = Data3;
    }
    else
    {
        printf("Error Input!");
        return false;
    }
    kernel->width = width;
    kernel->height = height;
    kernel->channels = ch;
    return true;
}

bool cnnFunction_v1(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)//input_data是padding之后的，input_height与input_width是input_image的size,只能通道1或者增加传参
{
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
    for(int j = 0; j < input_height; j++)
    {
        for(int i = 0; i < input_width;i++)
        {
            for(int m = 0; m<kernel_height;m++)
            {
                for(int n = 0,k = 0; n < kernel_width; n++,k++)
                {
                    *(output_data+(j*input_width)+i) += (*(input_data + (j * pad_width) + i + (m*pad_width) + n)) * (*(kernel_data+k));
                }
            }
        }
    }
    return true;
}

bool cnnFunction_v2(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
    for(int j = 0; j < input_height; j++)
    {
        for(int i = 0; i < input_width;i++)
        {
            for(int k =0; k<kernel_height * kernel_width;k++)
            {
                *(output_data+(j * input_width)+i) += (*(input_data + (j * pad_width) + i + ((k / kernel_width)*pad_width) + (k % kernel_width))) * (*(kernel_data+k));
            }
        }
    }
    return true;
}

bool cnnFunction_v3(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
    for(int t = 0; t < input_height*input_width; t++)
    {
        int j = t / input_width;
        int i = t % input_width;
        for(int m = 0; m<kernel_height;m++)
        {
            for(int n = 0,k = 0; n < kernel_width; n++,k++)
            {
                *(output_data+(j*input_width)+i) += (*(input_data + (j * pad_width) + i + (m*pad_width) + n)) * (*(kernel_data+k));
            }
        }
    }
    return true;
}

bool cnnFunction_v4(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)//input_data是padding之后的，input_height与input_width是input_image的size,只能通道1或者增加传参
{
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
    for(int t = 0; t < input_height*input_width; t++)
    {
        int j = t / input_width;
        int i = t % input_width;
        for(int k =0; k<kernel_height * kernel_width;k++)
        {
            *(output_data+(j * input_width)+i) += (*(input_data + (j * pad_width) + i + ((k / kernel_width)*pad_width) + (k % kernel_width))) * (*(kernel_data+k));
        }
    }
    return true;
}

bool cnn_unloop(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
    for(int j = 0; j < input_height; j++)
    {
        int stride = j*input_width;
        int i = 0;
        for(; i < (input_width/8)*8;i+=8)
        {
            for(int m = 0; m<kernel_height;m++)
            {
                for(int n = 0,k = 0; n < kernel_width; n++,k++)
                {
                    *(output_data+(stride)+i) += (*(input_data + (j * pad_width) + i + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+1)) += (*(input_data + (j * pad_width) + (i+1) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+2)) += (*(input_data + (j * pad_width) + (i+2) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+3)) += (*(input_data + (j * pad_width) + (i+3) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+4)) += (*(input_data + (j * pad_width) + (i+4) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+5)) += (*(input_data + (j * pad_width) + (i+5) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+6)) += (*(input_data + (j * pad_width) + (i+6) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+7)) += (*(input_data + (j * pad_width) + (i+7) + (m*pad_width) + n)) * (*(kernel_data+k));
                }
            }
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int m = 0; m<kernel_height;m++)
                {
                    for(int n = 0,k = 0; n < kernel_width; n++,k++)
                    {
                        *(output_data+(stride)+w) += (*(input_data + (j * pad_width) + w + m*pad_width + n)) * (*(kernel_data+k));
                    }
                }
            }
        }
    }
    return true;
}

bool cnn_omp_unloop(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
    #pragma omp parallel for
    for(int j = 0; j < input_height; j++)
    {
        int stride = j*input_width;
        int i = 0;
        for(; i < (input_width/8)*8;i+=8)
        {
            for(int m = 0; m<kernel_height;m++)
            {
                for(int n = 0,k = 0; n < kernel_width; n++,k++)
                {
                    *(output_data+(stride)+i) += (*(input_data + (j * pad_width) + i + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+1)) += (*(input_data + (j * pad_width) + (i+1) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+2)) += (*(input_data + (j * pad_width) + (i+2) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+3)) += (*(input_data + (j * pad_width) + (i+3) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+4)) += (*(input_data + (j * pad_width) + (i+4) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+5)) += (*(input_data + (j * pad_width) + (i+5) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+6)) += (*(input_data + (j * pad_width) + (i+6) + (m*pad_width) + n)) * (*(kernel_data+k));
                    *(output_data+(stride)+(i+7)) += (*(input_data + (j * pad_width) + (i+7) + (m*pad_width) + n)) * (*(kernel_data+k));
                }
            }
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int m = 0; m<kernel_height;m++)
                {
                    for(int n = 0,k = 0; n < kernel_width; n++,k++)
                    {
                        *(output_data+(stride)+w) += (*(input_data + (j * pad_width) + w + m*pad_width + n)) * (*(kernel_data+k));
                    }
                }
            }
        }
    }
    return true;
}


bool cnn_avx2_kernel1x1(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(kernel_height!=1||kernel_width!=1)
    {
        printf("The size of the kernel must be 1x1!");
        return false;
    }
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
#ifdef WITH_AVX2
    __m256 im, ker, res;
    float*Ker=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(int i=0;i<8;i++)
    {
        *(Ker+i)=*(kernel_data);
    }
    ker = _mm256_loadu_ps(Ker);
    float*temp=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(size_t j=0; j<input_height;j++)
    {
        size_t i=0;
        for(; i<(input_width/8)*8; i+=8)
        {
            res =  _mm256_setzero_ps();
            im = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width));
            res = _mm256_add_ps(res, _mm256_mul_ps(im, ker));
            _mm256_storeu_ps(temp, res);
            for(size_t k=0;k<8;k++)
                *(output_data+(j*input_width)+i+k) = *(temp+k);
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int m = 0; m<kernel_height;m++)
                {
                    for(int n = 0,k = 0; n < kernel_width; n++,k++)
                    {
                        *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + m*pad_width + n)) * (*(kernel_data+k));
                    }
                }
            }
        }
    }
    return TRUE;
#else
    printf("AVX2 is not supported");
    return FALSE;
#endif
}

bool cnn_avx2_omp_kernel1x1(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(kernel_height!=1||kernel_width!=1)
    {
        printf("The size of the kernel must be 1x1!");
        return FALSE;
    }
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return FALSE;
    }
    size_t pad_width=input_width+(kernel_width-1);
#ifdef WITH_AVX2
    __m256 im, ker, res;
    float*Ker=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(int i=0;i<8;i++)
    {
        *(Ker+i)=*(kernel_data);
    }
    ker = _mm256_loadu_ps(Ker);
    float*temp=(float*)(aligned_alloc(256,8*sizeof(float)));
    #pragma omp parallel for
    for(size_t j=0; j<input_height;j++)
    {
        size_t i=0;
        for(; i<(input_width/8)*8; i+=8)
        {
            res =  _mm256_setzero_ps();
            im = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width));
            res = _mm256_add_ps(res, _mm256_mul_ps(im, ker));
            _mm256_storeu_ps(temp, res);
            for(size_t k=0;k<8;k++)
                *(output_data+(j*input_width)+i+k) = *(temp+k);
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int m = 0; m<kernel_height;m++)
                {
                    for(int n = 0,k = 0; n < kernel_width; n++,k++)
                    {
                        *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + m*pad_width + n)) * (*(kernel_data+k));
                    }
                }
            }
        }
    }
    return TRUE;
#else
    printf("AVX2 is not supported");
    return FALSE;
#endif
}

bool cnn_avx2_kernel3x3(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(kernel_height!=3||kernel_width!=3)
    {
        printf("The size of the kernel must be 3x3!");
        return false;
    }
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return false;
    }
    size_t pad_width=input_width+(kernel_width-1);
#ifdef WITH_AVX2
    __m256 im11,im12,im13,im21,im22,im23,im31,im32,im33;
    __m256 ker11,ker12,ker13,ker21,ker22,ker23,ker31,ker32,ker33;
    __m256 res;
    float*Ker11=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker12=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker13=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker21=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker22=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker23=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker31=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker32=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker33=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(int i=0;i<8;i++)
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
    float*temp=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(size_t j=0; j<input_height;j++)
    {
        size_t i=0;
        for(; i<(input_width/8)*8; i+=8)
        {
            res =  _mm256_setzero_ps();
            im11 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 0);
            im12 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 1);
            im13 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 2);
            im21 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 0);
            im22 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 1);
            im23 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 2);
            im31 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 0);
            im32 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 1);
            im33 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 2);
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
            for(size_t k=0;k<8;k++)
                *(output_data+(j*input_width)+i+k) = *(temp+k);
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int k =0; k<kernel_height*kernel_width;k++)
                {
                    *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + (k / kernel_width)*pad_width + (k % kernel_width))) * (*(kernel_data+k));
                }
            }
        }
    }
    return TRUE;
#else
    printf("AVX2 is not supported");
    return FALSE;
#endif
}

// bool cnn_avx2_kernel3x3_512(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
//     const size_t kernel_width, float* output_data)
// {
//     if(kernel_height!=3||kernel_width!=3)
//     {
//         printf("The size of the kernel must be 3x3!");
//         return false;
//     }
//     if(input_width<16)
//     {
//         printf("The width fo the data must bigger than 8!");
//         return false;
//     }
//     if(input_data==NULL||kernel_data==NULL||output_data==NULL)
//     {
//         printf("NULL Pointer!");
//         return false;
//     }
//     size_t pad_width=input_width+(kernel_width-1);
// #ifdef WITH_AVX2
//     __m512 im11,im12,im13,im21,im22,im23,im31,im32,im33;
//     __m512 ker11,ker12,ker13,ker21,ker22,ker23,ker31,ker32,ker33;
//     __m512 res;
//     float*Ker11=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker12=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker13=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker21=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker22=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker23=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker31=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker32=(float*)(aligned_alloc(512,16*sizeof(float)));
//     float*Ker33=(float*)(aligned_alloc(512,16*sizeof(float)));
//     for(int i=0;i<16;i++)
//     {
//         *(Ker11+i)=*(kernel_data);
//         *(Ker12+i)=*(kernel_data+1);
//         *(Ker13+i)=*(kernel_data+2);
//         *(Ker21+i)=*(kernel_data+3);
//         *(Ker22+i)=*(kernel_data+4);
//         *(Ker23+i)=*(kernel_data+5);
//         *(Ker31+i)=*(kernel_data+6);
//         *(Ker32+i)=*(kernel_data+7);
//         *(Ker33+i)=*(kernel_data+8);
//     }
//     ker11 = _mm512_loadu_ps(Ker11);
//     ker12 = _mm512_loadu_ps(Ker12);
//     ker13 = _mm512_loadu_ps(Ker13);
//     ker21 = _mm512_loadu_ps(Ker21);
//     ker22 = _mm512_loadu_ps(Ker22);
//     ker23 = _mm512_loadu_ps(Ker23);
//     ker31 = _mm512_loadu_ps(Ker31);
//     ker32 = _mm512_loadu_ps(Ker32);
//     ker33 = _mm512_loadu_ps(Ker33);
//     float*temp=(float*)(aligned_alloc(512,16*sizeof(float)));
//     for(size_t j=0; j<input_height;j++)
//     {
//         size_t i=0;
//         for(; i<(input_width/16)*16; i+=16)
//         {
//             res =  _mm512_setzero_ps();
//             im11 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 0);
//             im12 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 1);
//             im13 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 2);
//             im21 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 0);
//             im22 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 1);
//             im23 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 2);
//             im31 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 0);
//             im32 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 1);
//             im33 = _mm512_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 2);
//             res = _mm512_add_ps(res, _mm512_mul_ps(im11, ker11));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im12, ker12));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im13, ker13));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im21, ker21));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im22, ker22));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im23, ker23));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im31, ker31));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im32, ker32));
//             res = _mm512_add_ps(res, _mm512_mul_ps(im33, ker33));
//             _mm512_storeu_ps(temp, res);
//             for(size_t k=0;k<16;k++)
//                 *(output_data+(j*input_width)+i+k) = *(temp+k);
//         }
//         if((input_width % 16) != 0)
//         {
//             for(int w=i;w<(input_width % 16)+i;w++)
//             {
//                 for(int k =0; k<kernel_height*kernel_width;k++)
//                 {
//                     *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + (k / kernel_width)*pad_width + (k % kernel_width))) * (*(kernel_data+k));
//                 }
//             }
//         }
//     }
//     return TRUE;
// #else
//     printf("AVX2 is not supported");
//     return FALSE;
// #endif
// }

bool cnn_avx2_omp_kernel3x3(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(kernel_height!=3||kernel_width!=3)
    {
        printf("The size of the kernel must be 3x3!");
        return FALSE;
    }
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return FALSE;
    }
    size_t pad_width=input_width+(kernel_width-1);
#ifdef WITH_AVX2
    __m256 im11,im12,im13,im21,im22,im23,im31,im32,im33;
    __m256 ker11,ker12,ker13,ker21,ker22,ker23,ker31,ker32,ker33;
    __m256 res;
    float*Ker11=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker12=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker13=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker21=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker22=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker23=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker31=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker32=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker33=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(int i=0;i<8;i++)
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
    float*temp=(float*)(aligned_alloc(256,8*sizeof(float)));
    #pragma omp parallel for
    for(size_t j=0; j<input_height;j++)
    {
        size_t i=0;
        for(; i<(input_width/8)*8; i+=8)
        {
            res =  _mm256_setzero_ps();
            im11 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 0);
            im12 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 1);
            im13 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 2);
            im21 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 0);
            im22 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 1);
            im23 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 2);
            im31 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 0);
            im32 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 1);
            im33 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 2);
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
            for(size_t k=0;k<8;k++)
                *(output_data+(j*input_width)+i+k) = *(temp+k);
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int k =0; k<kernel_height*kernel_width;k++)
                {
                    *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + (k / kernel_width)*pad_width + (k % kernel_width))) * (*(kernel_data+k));
                }
            }
        }
    }
    return TRUE;
#else
    printf("AVX2 is not supported\n");
    return FALSE;
#endif
}

bool cnn_avx2_kernel5x5(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)//加速的是单此卷积,由于内存的不连续，感觉很难下手
{
    if(kernel_height!=5||kernel_width!=5)
    {
        printf("The size of the kernel must be 5x5");
        return FALSE;
    }
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return FALSE;
    }
    size_t pad_width=input_width+(kernel_width-1);
#ifdef WITH_AVX2
    __m256 im11,im12,im13,im14,im15,im21,im22,im23,im24,im25,im31,im32,im33,im34,im35,im41,im42,im43,im44,im45,im51,im52,im53,im54,im55;
    __m256 ker11,ker12,ker13,ker14,ker15,ker21,ker22,ker23,ker24,ker25,ker31,ker32,ker33,ker34,ker35,ker41,ker42,ker43,ker44,ker45,ker51,ker52,ker53,ker54,ker55;
    __m256 res;
    float*Ker11=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker12=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker13=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker14=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker15=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker21=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker22=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker23=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker24=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker25=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker31=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker32=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker33=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker34=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker35=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker41=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker42=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker43=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker44=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker45=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker51=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker52=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker53=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker54=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker55=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(int i=0;i<8;i++)
    {
        *(Ker11+i)=*(kernel_data);
        *(Ker12+i)=*(kernel_data+1);
        *(Ker13+i)=*(kernel_data+2);
        *(Ker14+i)=*(kernel_data+3);
        *(Ker15+i)=*(kernel_data+4);
        *(Ker21+i)=*(kernel_data+5);
        *(Ker22+i)=*(kernel_data+6);
        *(Ker23+i)=*(kernel_data+7);
        *(Ker24+i)=*(kernel_data+8);
        *(Ker25+i)=*(kernel_data+9);
        *(Ker31+i)=*(kernel_data+10);
        *(Ker32+i)=*(kernel_data+11);
        *(Ker33+i)=*(kernel_data+12);
        *(Ker34+i)=*(kernel_data+13);
        *(Ker35+i)=*(kernel_data+14);
        *(Ker41+i)=*(kernel_data+15);
        *(Ker42+i)=*(kernel_data+16);
        *(Ker43+i)=*(kernel_data+17);
        *(Ker44+i)=*(kernel_data+18);
        *(Ker45+i)=*(kernel_data+19);
        *(Ker51+i)=*(kernel_data+20);
        *(Ker52+i)=*(kernel_data+21);
        *(Ker53+i)=*(kernel_data+22);
        *(Ker54+i)=*(kernel_data+23);
        *(Ker55+i)=*(kernel_data+24);
    }
    ker11 = _mm256_loadu_ps(Ker11);
    ker12 = _mm256_loadu_ps(Ker12);
    ker13 = _mm256_loadu_ps(Ker13);
    ker14 = _mm256_loadu_ps(Ker14);
    ker15 = _mm256_loadu_ps(Ker15);
    ker21 = _mm256_loadu_ps(Ker21);
    ker22 = _mm256_loadu_ps(Ker22);
    ker23 = _mm256_loadu_ps(Ker23);
    ker24 = _mm256_loadu_ps(Ker24);
    ker25 = _mm256_loadu_ps(Ker25);
    ker31 = _mm256_loadu_ps(Ker31);
    ker32 = _mm256_loadu_ps(Ker32);
    ker33 = _mm256_loadu_ps(Ker33);
    ker34 = _mm256_loadu_ps(Ker34);
    ker35 = _mm256_loadu_ps(Ker35);
    ker41 = _mm256_loadu_ps(Ker41);
    ker42 = _mm256_loadu_ps(Ker42);
    ker43 = _mm256_loadu_ps(Ker43);
    ker44 = _mm256_loadu_ps(Ker44);
    ker45 = _mm256_loadu_ps(Ker45);
    ker51 = _mm256_loadu_ps(Ker51);
    ker52 = _mm256_loadu_ps(Ker52);
    ker53 = _mm256_loadu_ps(Ker53);
    ker54 = _mm256_loadu_ps(Ker54);
    ker55 = _mm256_loadu_ps(Ker55);
    float*temp=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(size_t j=0; j<input_height;j++)
    {
        size_t i=0;
        for(; i<(input_width/8)*8; i+=8)
        {
            //获取连续内存
            res =  _mm256_setzero_ps();
            im11 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 0);
            im12 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 1);
            im13 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 2);
            im14 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 3);
            im15 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 4);
            im21 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 0);
            im22 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 1);
            im23 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 2);
            im24 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 3);
            im25 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 4);
            im31 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 0);
            im32 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 1);
            im33 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 2);
            im34 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 3);
            im35 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 4);
            im41 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 0);
            im42 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 1);
            im43 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 2);
            im44 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 3);
            im45 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 4);
            im51 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 0);
            im52 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 1);
            im53 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 2);
            im54 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 3);
            im55 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 4);
            res = _mm256_add_ps(res, _mm256_mul_ps(im11, ker11));
            res = _mm256_add_ps(res, _mm256_mul_ps(im12, ker12));
            res = _mm256_add_ps(res, _mm256_mul_ps(im13, ker13));
            res = _mm256_add_ps(res, _mm256_mul_ps(im14, ker14));
            res = _mm256_add_ps(res, _mm256_mul_ps(im15, ker15));
            res = _mm256_add_ps(res, _mm256_mul_ps(im21, ker21));
            res = _mm256_add_ps(res, _mm256_mul_ps(im22, ker22));
            res = _mm256_add_ps(res, _mm256_mul_ps(im23, ker23));
            res = _mm256_add_ps(res, _mm256_mul_ps(im24, ker24));
            res = _mm256_add_ps(res, _mm256_mul_ps(im25, ker25));
            res = _mm256_add_ps(res, _mm256_mul_ps(im31, ker31));
            res = _mm256_add_ps(res, _mm256_mul_ps(im32, ker32));
            res = _mm256_add_ps(res, _mm256_mul_ps(im33, ker33));
            res = _mm256_add_ps(res, _mm256_mul_ps(im34, ker34));
            res = _mm256_add_ps(res, _mm256_mul_ps(im35, ker35));
            res = _mm256_add_ps(res, _mm256_mul_ps(im41, ker41));
            res = _mm256_add_ps(res, _mm256_mul_ps(im42, ker42));
            res = _mm256_add_ps(res, _mm256_mul_ps(im43, ker43));
            res = _mm256_add_ps(res, _mm256_mul_ps(im44, ker44));
            res = _mm256_add_ps(res, _mm256_mul_ps(im45, ker45));
            res = _mm256_add_ps(res, _mm256_mul_ps(im51, ker51));
            res = _mm256_add_ps(res, _mm256_mul_ps(im52, ker52));
            res = _mm256_add_ps(res, _mm256_mul_ps(im53, ker53));
            res = _mm256_add_ps(res, _mm256_mul_ps(im54, ker54));
            res = _mm256_add_ps(res, _mm256_mul_ps(im55, ker55));
            _mm256_storeu_ps(temp, res);
            for(size_t k=0;k<8;k++)
                *(output_data+(j*input_width)+i+k) = *(temp+k);
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int k =0; k<kernel_height*kernel_width;k++)
                {
                    *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + (k / kernel_width)*pad_width + (k % kernel_width))) * (*(kernel_data+k));
                }
            }
        }
    }
    return TRUE;
#else
    printf("AVX2 is not supported\n");
    return FALSE;
#endif
}

bool cnn_avx2_omp_kernel5x5(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(kernel_height!=5||kernel_width!=5)
    {
        printf("The size of the kernel must be 5x5");
        return FALSE;
    }
    if(input_width<8)
    {
        printf("The width fo the data must bigger than 8!");
        return false;
    }
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return FALSE;
    }
    size_t pad_width=input_width+(kernel_width-1);
#ifdef WITH_AVX2
    __m256 im11,im12,im13,im14,im15,im21,im22,im23,im24,im25,im31,im32,im33,im34,im35,im41,im42,im43,im44,im45,im51,im52,im53,im54,im55;
    __m256 ker11,ker12,ker13,ker14,ker15,ker21,ker22,ker23,ker24,ker25,ker31,ker32,ker33,ker34,ker35,ker41,ker42,ker43,ker44,ker45,ker51,ker52,ker53,ker54,ker55;
    __m256 res;
    float*Ker11=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker12=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker13=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker14=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker15=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker21=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker22=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker23=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker24=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker25=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker31=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker32=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker33=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker34=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker35=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker41=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker42=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker43=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker44=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker45=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker51=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker52=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker53=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker54=(float*)(aligned_alloc(256,8*sizeof(float)));
    float*Ker55=(float*)(aligned_alloc(256,8*sizeof(float)));
    for(int i=0;i<8;i++)
    {
        *(Ker11+i)=*(kernel_data);
        *(Ker12+i)=*(kernel_data+1);
        *(Ker13+i)=*(kernel_data+2);
        *(Ker14+i)=*(kernel_data+3);
        *(Ker15+i)=*(kernel_data+4);
        *(Ker21+i)=*(kernel_data+5);
        *(Ker22+i)=*(kernel_data+6);
        *(Ker23+i)=*(kernel_data+7);
        *(Ker24+i)=*(kernel_data+8);
        *(Ker25+i)=*(kernel_data+9);
        *(Ker31+i)=*(kernel_data+10);
        *(Ker32+i)=*(kernel_data+11);
        *(Ker33+i)=*(kernel_data+12);
        *(Ker34+i)=*(kernel_data+13);
        *(Ker35+i)=*(kernel_data+14);
        *(Ker41+i)=*(kernel_data+15);
        *(Ker42+i)=*(kernel_data+16);
        *(Ker43+i)=*(kernel_data+17);
        *(Ker44+i)=*(kernel_data+18);
        *(Ker45+i)=*(kernel_data+19);
        *(Ker51+i)=*(kernel_data+20);
        *(Ker52+i)=*(kernel_data+21);
        *(Ker53+i)=*(kernel_data+22);
        *(Ker54+i)=*(kernel_data+23);
        *(Ker55+i)=*(kernel_data+24);
    }
    ker11 = _mm256_loadu_ps(Ker11);
    ker12 = _mm256_loadu_ps(Ker12);
    ker13 = _mm256_loadu_ps(Ker13);
    ker14 = _mm256_loadu_ps(Ker14);
    ker15 = _mm256_loadu_ps(Ker15);
    ker21 = _mm256_loadu_ps(Ker21);
    ker22 = _mm256_loadu_ps(Ker22);
    ker23 = _mm256_loadu_ps(Ker23);
    ker24 = _mm256_loadu_ps(Ker24);
    ker25 = _mm256_loadu_ps(Ker25);
    ker31 = _mm256_loadu_ps(Ker31);
    ker32 = _mm256_loadu_ps(Ker32);
    ker33 = _mm256_loadu_ps(Ker33);
    ker34 = _mm256_loadu_ps(Ker34);
    ker35 = _mm256_loadu_ps(Ker35);
    ker41 = _mm256_loadu_ps(Ker41);
    ker42 = _mm256_loadu_ps(Ker42);
    ker43 = _mm256_loadu_ps(Ker43);
    ker44 = _mm256_loadu_ps(Ker44);
    ker45 = _mm256_loadu_ps(Ker45);
    ker51 = _mm256_loadu_ps(Ker51);
    ker52 = _mm256_loadu_ps(Ker52);
    ker53 = _mm256_loadu_ps(Ker53);
    ker54 = _mm256_loadu_ps(Ker54);
    ker55 = _mm256_loadu_ps(Ker55);
    float*temp=(float*)(aligned_alloc(256,8*sizeof(float)));
    #pragma omp parallel for
    for(size_t j=0; j<input_height;j++)
    {
        size_t i=0;
        for(; i<(input_width/8)*8; i+=8)
        {
            res =  _mm256_setzero_ps();
            im11 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 0);
            im12 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 1);
            im13 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 2);
            im14 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 3);
            im15 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (0*pad_width) + 4);
            im21 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 0);
            im22 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 1);
            im23 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 2);
            im24 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 3);
            im25 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (1*pad_width) + 4);
            im31 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 0);
            im32 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 1);
            im33 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 2);
            im34 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 3);
            im35 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (2*pad_width) + 4);
            im41 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 0);
            im42 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 1);
            im43 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 2);
            im44 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 3);
            im45 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (3*pad_width) + 4);
            im51 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 0);
            im52 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 1);
            im53 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 2);
            im54 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 3);
            im55 = _mm256_loadu_ps(input_data + (j * pad_width) + i + (4*pad_width) + 4);
            res = _mm256_add_ps(res, _mm256_mul_ps(im11, ker11));
            res = _mm256_add_ps(res, _mm256_mul_ps(im12, ker12));
            res = _mm256_add_ps(res, _mm256_mul_ps(im13, ker13));
            res = _mm256_add_ps(res, _mm256_mul_ps(im14, ker14));
            res = _mm256_add_ps(res, _mm256_mul_ps(im15, ker15));
            res = _mm256_add_ps(res, _mm256_mul_ps(im21, ker21));
            res = _mm256_add_ps(res, _mm256_mul_ps(im22, ker22));
            res = _mm256_add_ps(res, _mm256_mul_ps(im23, ker23));
            res = _mm256_add_ps(res, _mm256_mul_ps(im24, ker24));
            res = _mm256_add_ps(res, _mm256_mul_ps(im25, ker25));
            res = _mm256_add_ps(res, _mm256_mul_ps(im31, ker31));
            res = _mm256_add_ps(res, _mm256_mul_ps(im32, ker32));
            res = _mm256_add_ps(res, _mm256_mul_ps(im33, ker33));
            res = _mm256_add_ps(res, _mm256_mul_ps(im34, ker34));
            res = _mm256_add_ps(res, _mm256_mul_ps(im35, ker35));
            res = _mm256_add_ps(res, _mm256_mul_ps(im41, ker41));
            res = _mm256_add_ps(res, _mm256_mul_ps(im42, ker42));
            res = _mm256_add_ps(res, _mm256_mul_ps(im43, ker43));
            res = _mm256_add_ps(res, _mm256_mul_ps(im44, ker44));
            res = _mm256_add_ps(res, _mm256_mul_ps(im45, ker45));
            res = _mm256_add_ps(res, _mm256_mul_ps(im51, ker51));
            res = _mm256_add_ps(res, _mm256_mul_ps(im52, ker52));
            res = _mm256_add_ps(res, _mm256_mul_ps(im53, ker53));
            res = _mm256_add_ps(res, _mm256_mul_ps(im54, ker54));
            res = _mm256_add_ps(res, _mm256_mul_ps(im55, ker55));
            _mm256_storeu_ps(temp, res);
            for(size_t k=0;k<8;k++)
                *(output_data+(j*input_width)+i+k) = *(temp+k);
        }
        if((input_width % 8) != 0)
        {
            for(int w=i;w<(input_width % 8)+i;w++)
            {
                for(int k =0; k<kernel_height*kernel_width;k++)
                {
                    *(output_data+(j*input_width)+w) += (*(input_data + (j * pad_width) + w + (k / kernel_width)*pad_width + (k % kernel_width))) * (*(kernel_data+k));
                }
            }
        }
    }
    return TRUE;
#else
    printf("AVX2 is not supported");
    return FALSE;
#endif
}

bool cnn_openblas(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data)
{
    if(input_data==NULL||kernel_data==NULL||output_data==NULL)
    {
        printf("NULL Pointer!");
        return FALSE;
    }
    float* data = (float*)malloc(sizeof(float)*kernel_height*kernel_width*input_height*input_width);
    if(data == NULL)
    {
        printf("Fail to allocate memory in cnn_openblas\n");
        return false;
    }
    size_t pad_width = input_width + (kernel_width - 1);
    int t = 0;
    for(int i = 0; i < input_height; i++){
        for(int j = 0; j < input_width; j++){
            for(int m = 0; m < kernel_height; m++){
                for(int n = 0; n < kernel_width; n++, t++){
                    data[t] = *(input_data + i*pad_width + j + m*kernel_width + n);
                }
            }
        }
    }
    struct timeval start, end;
    gettimeofday(&start, NULL);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, input_height*input_width, 1, kernel_height*kernel_width, 1.0f, data, kernel_height*kernel_width, kernel_data, 1, 0, output_data, 1);
    gettimeofday(&end, NULL);
    printf("%lf\n", (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0);
    // cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->row, blob.col, this->col, 1.0f, left, this->col, right, blob.col, 0, res, blob.col);
    free(data);
    return true;
}

bool Convolution(const Image* input_image, const Kernel* kernel, Image* output_image)
{
    //Safe Detection
    if(input_image==NULL||kernel==NULL||output_image==NULL)
    {
        printf("NULL pointer!\n");
        return FALSE;
    }
    if(input_image->channels!=kernel->channels)
    {
        printf("The channels of the image is not fit with the kernel!\n");
        return FALSE;
    }
    if(input_image->channels!=3&&input_image->channels!=1)
    {
        printf("The channels should be 1 or 3.\n");
        return FALSE;
    }

    if(input_image->channels==1)
    {
        //padding
        float* pad_data = padding(input_image->data1,input_image->height,input_image->width,kernel->height,kernel->width);

        //cnn
        float* data1 = (float*)malloc(sizeof(float)*input_image->height*input_image->width);
        if(kernel->height==1&&kernel->width==1)
        {
            cnn_avx2_kernel1x1(pad_data,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            output_image->data1 = data1;
        }
        else if(kernel->height==3&&kernel->width==3)
        {
            cnn_openblas(pad_data,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            // cnn_avx2_kernel3x3(pad_data,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            output_image->data1 = data1;
        }
        else if(kernel->height==5&&kernel->width==5)
        {
            cnn_avx2_kernel5x5(pad_data,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            output_image->data1 = data1;
        }
        output_image->data2 = NULL;
        output_image->data3 = NULL;
        free(pad_data);
    }
    else if(input_image->channels==3)
    {
        //padding
        float* pad_data1 = padding(input_image->data1,input_image->height,input_image->width,kernel->height,kernel->width);
        float* pad_data2 = padding(input_image->data2,input_image->height,input_image->width,kernel->height,kernel->width);
        float* pad_data3 = padding(input_image->data3,input_image->height,input_image->width,kernel->height,kernel->width);

        //cnn
        float* data1 = (float*)malloc(sizeof(float)*input_image->height*input_image->width);
        float* data2 = (float*)malloc(sizeof(float)*input_image->height*input_image->width);
        float* data3 = (float*)malloc(sizeof(float)*input_image->height*input_image->width);
        if(kernel->height==1&&kernel->width==1)
        {
            cnn_avx2_kernel1x1(pad_data1,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            cnn_avx2_kernel1x1(pad_data2,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data2);
            cnn_avx2_kernel1x1(pad_data3,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data3);
            output_image->data1 = data1;
            output_image->data2 = data2;
            output_image->data3 = data3;
        }
        else if(kernel->height==3&&kernel->width==3)
        {
            cnn_avx2_kernel3x3(pad_data1,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            cnn_avx2_kernel3x3(pad_data2,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data2);
            cnn_avx2_kernel3x3(pad_data3,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data3);
            output_image->data1 = data1;
            output_image->data2 = data2;
            output_image->data3 = data3;
        }
        else if(kernel->height==5&&kernel->width==5)
        {
            cnn_avx2_kernel5x5(pad_data1,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data1);
            cnn_avx2_kernel5x5(pad_data2,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data2);
            cnn_avx2_kernel5x5(pad_data3,input_image->height,input_image->width,kernel->data1,kernel->height,kernel->width,data3);
            output_image->data1 = data1;
            output_image->data2 = data2;
            output_image->data3 = data3;
        }
        else
        {
            printf("The size of the kernel must be 1x1, 3x3 or 5x5");
        }
        free(pad_data1);
        free(pad_data2);
        free(pad_data3);
    }
    output_image->width = input_image->width;
    output_image->height = input_image->height;
    output_image->channels = input_image->channels;
    return true;
}

float* Random(const size_t height, const size_t width)
{
    float* data = (float*)malloc(sizeof(float)*height*width);
    if(data==NULL)
    {
        printf("Malloc Error");
        return NULL;
    }
    srand((unsigned)time(NULL));
    for(int i=0;i<height*width;i++)
        *(data + i) = ((rand() / (float)RAND_MAX) + (rand()% 254 ));
    return data;
}

bool printinfo(const float* data, const size_t height, const size_t width, const char* name)
{
    if(data==NULL)
    {
        printf("NULL Pointer!");
        return FALSE;
    }
    printf("%s:\n",name);
    for(int j=0; j<height;j++)
    {
        for(int i=0; i<width; i++)
            printf("%.2f ",*(data + j*width + i));
        printf("\n");
    }
    printf("\n");
    return true;
}

bool test(bool (*cnn)(const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data, const size_t kernel_height, 
    const size_t kernel_width, float* output_data),const float* input_data, const size_t input_height, const size_t input_width, const float* kernel_data,
    const size_t kernel_height,const size_t kernel_width,FILE* file,const char* name)
{
    
    bool flag = false;
    float*pad_data = padding(input_data,input_height,input_width,kernel_height,kernel_width);
    if(pad_data==NULL)
    {
        printf("Padding Error");
        return false;
    }
    float* cnn_data = (float*)calloc(input_height*input_width,sizeof(float));
    if(cnn_data==NULL)
    {
        printf("cnn_data calloc error!");
        return false;
    }
    struct timeval start, end;
    int ret;
    ret = gettimeofday(&start, NULL);
    if (ret == -1)
    {
        printf("Error: gettimeofday()\n");
        return ret;
    }
    flag = cnn(pad_data, input_height, input_width, kernel_data, kernel_height, kernel_width, cnn_data);
    if(!flag)
    {
        printf("cnn Error!");
        return false;
    }
    ret = gettimeofday(&end, NULL);
    if (ret == -1)
    {
        printf("Error: gettimeofday()\n");
        return ret;
    }
    fprintf(file,"%lf,",(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0);
    // printinfo(cnn_data,input_height,input_width,name);
    // printf("Time Cost : %lf s\n",(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0);
    // printf("\n\n");
    free(cnn_data);
    free(pad_data);
}

bool Free_Image(Image* image)
{
    if(image==NULL)
    {
        printf("NULL pointer!");
        return false;
    }
    if(image->data1!=NULL)
    {
        free(image->data1);
        image->data1=NULL;
    }
    if(image->data2!=NULL)
    {
        free(image->data2);
        image->data2=NULL;
    }
    if(image->data3!=NULL)
    {
        free(image->data3);
        image->data3=NULL;
    }
    return true;
}