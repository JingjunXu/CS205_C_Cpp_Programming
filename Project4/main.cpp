#include<iostream>
#include <opencv2/opencv.hpp>//读写图片，做卷积
#include <sys/time.h>//gettimeofday
#include <fstream>//读写CSV文件
#include "DataBlobs.hpp"
using namespace std;
using namespace cv;

/**************由于进行的实验过多，这里没有展示所有的实验的代码，大部分已经放在报告里面**************/

#define randomInt(a,b)       (rand()% (b - a) + (a) )           //取得[a,b)的随机整数
#define randomDouble(RAND_MAX) (rand() / double(RAND_MAX))      //取得0～1之间的随机双精度浮点数
#define randomFloat(RAND_MAX) (rand() / float(RAND_MAX))        //取得0~1之间的随机单精度浮点数

float* Random(size_t row, size_t col);
int main()
{
    //Convolutions
    // float* kernel_data = new float[9]{-1, 0 ,1, -2 ,0, 2, -1, 0, 1};
    // DataBlobs kernel(3,3,1,kernel_data,NULL);
    // DataBlobs<float> input("gray.jpg");
    // DataBlobs<float> r1,g1,b1,r2,g2,b2,r3,g3,b3,Output;
    // input.getRed(r1);
    // input.getGreen(g1);
    // input.getBlue(b1);
    // DataBlobs<float>::padding(r1,3,3,r2);
    // DataBlobs<float>::padding(g1,3,3,g2);
    // DataBlobs<float>::padding(b1,3,3,b2);
    // DataBlobs<float>::cnn(r2,kernel,r3);
    // DataBlobs<float>::cnn(g2,kernel,g3);
    // DataBlobs<float>::cnn(b2,kernel,b3);
    // DataBlobs<float>::cnn_unloop(r2,kernel,r3);
    // DataBlobs<float>::cnn_unloop(g2,kernel,g3);
    // DataBlobs<float>::cnn_unloop(b2,kernel,b3);
    // DataBlobs<float>::cnn_avx2_ker3x3(r2,kernel,r3);
    // DataBlobs<float>::cnn_avx2_ker3x3(g2,kernel,g3);
    // DataBlobs<float>::cnn_avx2_ker3x3(b2,kernel,b3);
    // DataBlobs<float>::rearrange_img2Col(r1,3,3,r2);
    // DataBlobs<float>::rearrange_img2Col(g1,3,3,g2);
    // DataBlobs<float>::rearrange_img2Col(b1,3,3,b2);
    // DataBlobs<float>::cnn_img2Col(r2,input.getRow(),input.getCol(),kernel,r3);
    // DataBlobs<float>::cnn_img2Col(g2,input.getRow(),input.getCol(),kernel,g3);
    // DataBlobs<float>::cnn_img2Col(b2,input.getRow(),input.getCol(),kernel,b3);
    // DataBlobs<float>::CombineRGB(r3,g3,b3,Output);
    // Output.output("cnn_img2Col_output.jpg");

    //Convolution using opencv
    // cv::Mat Kernel = (cv::Mat_<float>(3, 3) << -1, 0, 1,
    //                                         -2, 0, 2,
    //                                         -1, 0, 1);
    // cv::Mat image = cv::imread("gray.jpg");
    // cv::Mat result;
    // cv::filter2D(image, result, -1, Kernel);
    // imwrite("opencv_output.jpg", result);

    //x86 arm comparison
    //matrix 乘法
    // ofstream outFile;
    // outFile.open("TimeCost_x86.csv", std::ios::out | std::ios::trunc);
    // for(int num = 1; num<=10; num++)
    // {
    //     float* data1 = Random(4096,4096);
    //     float* data2 = Random(4096,4096);
    //     DataBlobs d1(4096,4096,1,data1,NULL);
    //     DataBlobs d2(4096,4096,1,data2,NULL);
    //     outFile << num << ",";
    //     struct timeval start, end;
    //     gettimeofday(&start, NULL);
    //     d1 * d2;
    //     gettimeofday(&end, NULL);
    //     outFile << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << endl;
    // }
    // outFile.close();

    //cnn_simd_avx2
    // ofstream outFile;
    // outFile.open("TimeCost_x86.csv", std::ios::out | std::ios::trunc);
    // float* kernel_data = new float[9]{-1, 0 ,1, -2 ,0, 2, -1, 0, 1};
    // DataBlobs kernel(3,3,1,kernel_data,NULL);
    // DataBlobs<float> res;
    // for(int num = 1; num<=10; num++)
    // {
    //     float* data1 = Random(128,128);
    //     DataBlobs d1(128,128,1,data1,NULL);
    //     outFile << num << ",";
    //     struct timeval start, end;
    //     gettimeofday(&start, NULL);
    //     DataBlobs<float>::cnn_avx2_ker3x3(d1,kernel,res);
    //     gettimeofday(&end, NULL);
    //     outFile << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << endl;
    // }
    // outFile.close();

    //cnn_simd_arm
    // ofstream outFile;
    // outFile.open("TimeCost_arm.csv", std::ios::out | std::ios::trunc);
    // float* kernel_data = new float[9]{-1, 0 ,1, -2 ,0, 2, -1, 0, 1};
    // DataBlobs kernel(3,3,1,kernel_data,NULL);
    // DataBlobs<float> res;
    // for(int num = 1; num<=10; num++)
    // {
    //     float* data1 = Random(128,128);
    //     DataBlobs d1(128,128,1,data1,NULL);
    //     outFile << num << ",";
    //     struct timeval start, end;
    //     gettimeofday(&start, NULL);
    //     DataBlobs<float>::cnn_neon_ker3x3(d1,kernel,res);
    //     gettimeofday(&end, NULL);
    //     outFile << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0 << endl;
    // }
    // outFile.close();

    return 0;
}

float* Random(size_t row, size_t col){
    if(row == 0 || col == 0){
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "The input (" << row << ", " << col << ") is invalid" << endl;
        exit(EXIT_FAILURE);
    }
    float* data = new float[row*col];
    if(data==NULL)
    {
        cerr << "File " << __FILE__ << "Line " << __LINE__ << "Func " << __FUNCTION__ << "()" << endl;
        cerr << "Fail to allocate memory for data" << endl;
        exit(EXIT_FAILURE);
    }
    srand((unsigned)time(NULL));
    for(int i=0;i<row*col;i++)
        *(data + i) = ((rand() / (float)RAND_MAX) + (rand()% 254 ));
    return data;
}