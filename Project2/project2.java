import java.util.ArrayList;
import java.io.FileWriter;
import java.io.IOException;
public class project2 {
    public static void PrintVectorStatistic (ArrayList num)
    {
        for(int i=0;i<num.size();i++)
        {
            System.out.print(num.get(i));
            System.out.print(" ");
        }
        System.out.println();
    }
    public static void main(String[] args) 
    {
        /*************************实验1*************************/
        // //由于java循环多次之后会加速，故要一次次做实验
        // ArrayList<Integer> res_int1=new ArrayList<Integer>(), res_int2=new ArrayList<Integer>();
        // CreateRandomNum obj = new CreateRandomNum();
        // res_int1=obj.create_rand_int(0, 100, 10000000);
        // res_int2=obj.create_rand_int(0, 100, 10000000);
        // long int_result,int_start, int_finish;
        // int_start= System.nanoTime();
        // int_result=DotProduct.DotProductInteger(res_int1,res_int2);
        // int_finish= System.nanoTime();
        // System.out.println("The Result of int dot product is " + int_result + " Time Cost: " + (double)(int_finish-int_start)/1000000000 + " s");

        /*************************实验2*************************/
        // ArrayList<Integer> res_int1=new ArrayList<Integer>(), res_int2=new ArrayList<Integer>();
        // ArrayList<Float> res_float1=new ArrayList<Float>(), res_float2=new ArrayList<Float>();
        // ArrayList<Double> res_double1=new ArrayList<Double>(), res_double2=new ArrayList<Double>();
        // CreateRandomNum obj = new CreateRandomNum();
        // res_int1=obj.create_rand_int(0, 100, 1000);
        // res_int2=obj.create_rand_int(0, 100, 1000);
        // //数据类型转换
        // for(int i=0;i<res_int1.size();i++)
        // {
        //     res_float1.add((float)(res_int1.get(i)));
        //     res_double1.add((double)(res_int1.get(i)));
        // }
        // for(int i=0;i<res_int2.size();i++)
        // {
        //     res_float2.add((float)(res_int2.get(i)));
        //     res_double2.add((double)(res_int2.get(i)));
        // }
        // long int_result;
        // double float_result, double_result;
        // long int_start, int_finish, float_start=0, float_finish=0, double_start=0, double_finish=0;
        // //进行点乘并计算耗费的时间
        // int_start= System.nanoTime();
        // int_result=DotProduct.DotProductInteger(res_int1,res_int2);
        // int_finish= System.nanoTime();

        // float_start= System.nanoTime();
        // float_result=DotProduct.DotProductFloat(res_float1,res_float2);
        // float_finish= System.nanoTime();

        // double_start= System.nanoTime();
        // double_result=DotProduct.DotProductDouble(res_double1,res_double2);
        // double_finish= System.nanoTime();
        // //分别输出结果以及所耗费的时间
        // System.out.println("The Result of int dot product is " + int_result + " Time Cost: " + (double)(int_finish-int_start)/1000000000 + " s");
        // System.out.println("The Result of float dot product is " + float_result + " Time Cost: " + (double)(float_finish-float_start)/1000000000 + " s");
        // System.out.println("The Result of double dot product is " + double_result + " Time Cost: " + (double)(double_finish-double_start)/1000000000 + " s");

        /*************************实验3*************************/
        // ArrayList<Integer> res_int1=new ArrayList<Integer>(), res_int2=new ArrayList<Integer>();
        // ArrayList<Float> res_float1=new ArrayList<Float>(), res_float2=new ArrayList<Float>();
        // ArrayList<Double> res_double1=new ArrayList<Double>(), res_double2=new ArrayList<Double>();
        // ReadFileFromC.fileReadInt(res_int1, "intdata1.bin");
        // ReadFileFromC.fileReadInt(res_int2, "intdata2.bin");
        // ReadFileFromC.fileReadFloat(res_float1, "floatdata1.bin");
        // ReadFileFromC.fileReadFloat(res_float2, "floatdata2.bin");
        // ReadFileFromC.fileReadDouble(res_double1, "doubledata1.bin");
        // ReadFileFromC.fileReadDouble(res_double2, "doubledata2.bin");
        // long int_result;
        // double float_result, double_result;
        // long int_start, int_finish, float_start=0, float_finish=0, double_start=0, double_finish=0;
        // //进行点乘并计算耗费的时间
        // int_start= System.nanoTime();
        // int_result=DotProduct.DotProductInteger(res_int1,res_int2);
        // int_finish= System.nanoTime();

        // float_start= System.nanoTime();
        // float_result=DotProduct.DotProductFloat(res_float1,res_float2);
        // float_finish= System.nanoTime();

        // double_start= System.nanoTime();
        // double_result=DotProduct.DotProductDouble(res_double1,res_double2);
        // double_finish= System.nanoTime();
        // //分别输出结果以及所耗费的时间
        // System.out.println("The Result of int dot product is " + int_result + " Time Cost: " + (double)(int_finish-int_start)/1000000000 + " s");
        // System.out.println("The Result of float dot product is " + float_result + " Time Cost: " + (double)(float_finish-float_start)/1000000000 + " s");
        // System.out.println("The Result of double dot product is " + double_result + " Time Cost: " + (double)(double_finish-double_start)/1000000000 + " s");

        /*************************实验4*************************/
        // String csvFile = "Exp4_TimeCostInJ.csv";
        // String csvSplitBy = ",";
        // try (FileWriter fw = new FileWriter(csvFile)) {
        //     ArrayList<Integer> res_int1=new ArrayList<Integer>(), res_int2=new ArrayList<Integer>();
        //     ArrayList<Float> res_float1=new ArrayList<Float>(), res_float2=new ArrayList<Float>();
        //     ArrayList<Double> res_double1=new ArrayList<Double>(), res_double2=new ArrayList<Double>();
        //     CreateRandomNum obj = new CreateRandomNum();
        //     long int_result;
        //     double float_result, double_result;
        //     long int_start, int_finish, float_start, float_finish, double_start, double_finish;
        //     fw.append("实验次数").append(csvSplitBy).append("int").append(csvSplitBy).append("float").append(csvSplitBy).append("double").append("\n");
        //     for(int i=1;i<=25;i++)
        //     {
        //         res_int1=obj.create_rand_int(0, 100, 100000);
        //         res_int2=obj.create_rand_int(0, 100, 100000);
        //         res_float1=obj.create_rand_float(0, 100, 100000);
        //         res_float2=obj.create_rand_float(0, 100, 100000);
        //         res_double1=obj.create_rand_double(0, 100, 100000);
        //         res_double2=obj.create_rand_double(0, 100, 100000);
        //         //进行点乘并计算耗费的时间
        //         int_start= System.nanoTime();
        //         int_result=DotProduct.DotProductInteger(res_int1,res_int2);
        //         int_finish= System.nanoTime();

        //         float_start= System.nanoTime();
        //         float_result=DotProduct.DotProductFloat(res_float1,res_float2);
        //         float_finish= System.nanoTime();

        //         double_start= System.nanoTime();
        //         double_result=DotProduct.DotProductDouble(res_double1,res_double2);
        //         double_finish= System.nanoTime();
        //         fw.append(String.valueOf(i)).append(csvSplitBy);
        //         fw.append(String.valueOf((double)(int_finish-int_start)/1000000000.0)).append(csvSplitBy);
        //         fw.append(String.valueOf((double)(float_finish-float_start)/1000000000.0)).append(csvSplitBy);
        //         fw.append(String.valueOf((double)(double_finish-double_start)/1000000000.0)).append(csvSplitBy).append("\n");
        //     }
        //     fw.close();
        // } catch (IOException e) {
        //     e.printStackTrace();
        // }

        /*************************实验6*************************/
        // ArrayList<Integer> res_int1=new ArrayList<Integer>(), res_int2=new ArrayList<Integer>();
        // ArrayList<Float> res_float1=new ArrayList<Float>(), res_float2=new ArrayList<Float>();
        // ArrayList<Double> res_double1=new ArrayList<Double>(), res_double2=new ArrayList<Double>();
        // ReadFileFromC.fileReadInt(res_int1, "intdata1.bin");
        // ReadFileFromC.fileReadInt(res_int2, "intdata2.bin");
        // ReadFileFromC.fileReadFloat(res_float1, "floatdata1.bin");
        // ReadFileFromC.fileReadFloat(res_float2, "floatdata2.bin");
        // ReadFileFromC.fileReadDouble(res_double1, "doubledata1.bin");
        // ReadFileFromC.fileReadDouble(res_double2, "doubledata2.bin");
        // long int_result;
        // double float_result, double_result;
        // long int_start, int_finish, float_start=0, float_finish=0, double_start=0, double_finish=0;
        // //进行点乘并计算耗费的时间
        // int_start= System.nanoTime();
        // int_result=DotProduct.DotProductInteger(res_int1,res_int2);
        // int_finish= System.nanoTime();

        // float_start= System.nanoTime();
        // float_result=DotProduct.DotProductFloat(res_float1,res_float2);
        // float_finish= System.nanoTime();

        // double_start= System.nanoTime();
        // double_result=DotProduct.DotProductDouble(res_double1,res_double2);
        // double_finish= System.nanoTime();
        // //分别输出结果以及所耗费的时间
        // System.out.println("The Result of int dot product is " + int_result + " Time Cost: " + (double)(int_finish-int_start)/1000000000 + " s");
        // System.out.println("The Result of float dot product is " + float_result + " Time Cost: " + (double)(float_finish-float_start)/1000000000 + " s");
        // System.out.println("The Result of double dot product is " + double_result + " Time Cost: " + (double)(double_finish-double_start)/1000000000 + " s");
    }
}