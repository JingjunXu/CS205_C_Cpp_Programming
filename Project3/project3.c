#include <stdio.h>
#include <stdlib.h>
#include <string.h>//使用memcpy函数来提高padding的效率
#include <jpeglib.h>//读取jpg
#include <stdbool.h>//使用bool数据类型增强可读性
#include <sys/time.h>//gettimeofday
#include <time.h>//获取随机数种子
#include "Func.h"

#ifdef WITH_AVX2
#include <immintrin.h>
#endif 

#ifdef WITH_NEON
#include <arm_neon.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


int main()
{
    //Convolution operation to an image
    // bool flag = false;
    // Image input_image;
    // flag=input_jpeg_image("gray.jpg",&input_image);
    // if(!flag)
    // {
    //     printf("input error\n");
    //     return -1;
    // }
    // else
    //     printf("input success!\n");
    // Kernel kernel;
    // float kernel_data[9] = {-1, 0 ,1, -2 ,0, 2, -1, 0, 1};
    // flag=kernel_init(&kernel,3,3,input_image.channels,kernel_data,kernel_data,kernel_data);
    // if(!flag)
    // {
    //     printf("kernel initialization error\n");
    //     return -1;
    // }
    // else
    //     printf("kernel initialization success!\n");
    // Image output_image;
    // flag = Convolution(&input_image,&kernel,&output_image);
    // if(!flag)
    // {
    //     printf("cnn error\n");
    //     return -1;
    // }
    // else
    //     printf("cnn success!\n");
    // flag = output_jpeg_image("rgboutput.jpg",&output_image);
    // if(!flag)
    // {
    //     printf("output error");
    //     return -1;
    // }
    // else
    //     printf("output success!\n");
    // flag = Free_Image(&input_image);
    // if(!flag)
    // {
    //     printf("free image error\n");
    //     return -1;
    // }
    // else
    //     printf("free image success!\n");

    //experiment of efficiency
    FILE *fp;
    fp = fopen("TimeCost.csv", "a");
    if(fp == NULL) {
        printf("Failed to open file\n");
        return 1;
    }
    for(int num = 1; num<=10; num++)
    {
        //获取随机数
        float* image_data = (float*) malloc(128*128*sizeof(float));
        image_data = Random(128,128);
        float* kernel_data = (float*)malloc(3*3*sizeof(float));
        kernel_data = Random(3,3);
        fprintf(fp,"%d,",num);
        test(cnnFunction_v1,image_data,128,128,kernel_data,3,3,fp,"cnn_v1");
        // test(cnn_unloop,image_data,4096,4096,kernel_data,3,3,fp,"cnn_avx2_kernel3x3");
        // test(cnn_omp_unloop,image_data,4096,4096,kernel_data,3,3,fp,"cnn_avx2_kernel3x3");
        test(cnn_avx2_kernel3x3,image_data,128,128,kernel_data,3,3,fp,"cnn_avx2_kernel3x3");
        // test(cnn_avx2_omp_kernel3x3,image_data,4096,4096,kernel_data,3,3,fp,"cnn_avx2_omp_kernel3x3");
        test(cnn_openblas,image_data,128,128,kernel_data,3,3,fp,"cnn_openblas");
        // test(cnn_avx2_kernel3x3_512,image_data,4096,4096,kernel_data,3,3,fp,"cnn_avx2_kernel3x3_128");
        fprintf(fp,"\n");
        free(image_data);
        free(kernel_data);
    }
    fclose(fp);
    return 0;
}
