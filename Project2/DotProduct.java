import java.util.ArrayList;
public class DotProduct {
    public static long DotProductInteger(ArrayList<Integer> num1, ArrayList<Integer> num2)
    {
        if(num1.size()!=num2.size())
        {
            System.out.println("The size of the vectors are not fit!");
            return 0;
        }
        long result=0;
        for(int i=0;i<num1.size();i++)
        result += num1.get(i) * num2.get(i);
        return result;
    }
    public static long DotProductLong(ArrayList<Long> num1, ArrayList<Long> num2)
    {
        if(num1.size()!=num2.size())
        {
            System.out.println("The size of the vectors are not fit!");
            return 0;
        }
        long result=0;
        for(int i=0;i<num1.size();i++)
        result += num1.get(i) * num2.get(i);
        return result;
    }
    public static float DotProductFloat(ArrayList<Float> num1,ArrayList<Float> num2)
    {
        if(num1.size()!=num2.size())
        {
            System.out.println("The size of the vectors are not fit!");
            return 0;
        }
        float result= 0;
        for(int i=0;i<num1.size();i++)
        result += num1.get(i) * num2.get(i);
        return result;
    }
    public static double DotProductDouble(ArrayList<Double> num1, ArrayList<Double> num2)
    {
        if(num1.size()!=num2.size())
        {
            System.out.println("The size of the vectors are not fit!");
            return 0;
        }
        double result=0;
        for(int i=0;i<num1.size();i++)
        result += num1.get(i)*num2.get(i);
        return result;
    }
}
