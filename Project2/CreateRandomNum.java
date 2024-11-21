import java.util.Random;
import java.util.ArrayList;
public class CreateRandomNum {
    public static ArrayList<Integer> create_rand_int(int a,int b,int num)//取得[a,b)的随机整数组
    {
        if(num<=0)
        {
            System.out.println("Wrong input of the number");
            return null;
        }
        ArrayList<Integer> res = new ArrayList<Integer>();
        Random temp = new Random();
        for(int i=1;i<=num;i++)
        res.add(temp.nextInt(b-a)+a);
        return res;
    }
    public static ArrayList<Float> create_rand_float(int a, int b, int num)//取得(a,b+1)之间的随机单精度浮点数组
    {
        if(num<=0)
        {
            System.out.println("Wrong input of the number");
            return null;
        }
        ArrayList<Float> res= new ArrayList<>();
        Random temp= new Random();
        for(int i=0;i<num;i++)
        res.add(temp.nextFloat()+temp.nextInt(b-a)+a);
        return res;
    }
    public static ArrayList<Double> create_rand_double(int a,int b,int num)//取得(a,b+1)之间的随机双精度浮点数组
    {
        if(num<=0)
        {
            System.out.println("Wrong input of the number");
            return null;
        }
        ArrayList<Double> res =new ArrayList<>();
        Random temp =new Random();
        for(int i=1; i<=num; i++)
        res.add(temp.nextDouble()+temp.nextInt(b-a)+a);
        return res;
    }
}
