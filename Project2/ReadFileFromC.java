import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.MappedByteBuffer;
import java.util.ArrayList;
public class ReadFileFromC {
    public static void fileReadInt (ArrayList<Integer> num, String filename)//fileName的后缀为.bin
    {
        File file= new File(filename);
        if(!file.exists())
        {
            System.out.println("file not exist");
            return;
        }
        try
        {
            // 打开文件
            FileInputStream fis = new FileInputStream(filename);
            FileChannel fileChannel = fis.getChannel();
            // 将文件映射到内存
            MappedByteBuffer buffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0, fileChannel.size());
            buffer.order(ByteOrder.LITTLE_ENDIAN); // 设置缓冲区的字节顺序为小端字节序
            // 读取数据
            while (buffer.hasRemaining()) {
                num.add( buffer.getInt());
            }
            // 关闭文件
            fileChannel.close();
            fis.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
    }
    public static void fileReadFloat (ArrayList<Float> num, String filename)//fileName的后缀为.bin
    {
        File file= new File(filename);
        if(!file.exists())
        {
            System.out.println("file not exist");
            return;
        }
        try
        {
            // 打开文件
            FileInputStream fis = new FileInputStream(filename);
            FileChannel fileChannel = fis.getChannel();

            // 将文件映射到内存
            MappedByteBuffer buffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0, fileChannel.size());
            buffer.order(ByteOrder.LITTLE_ENDIAN); // 设置缓冲区的字节顺序为小端字节序

            // 读取数据
            while (buffer.hasRemaining()) {
                num.add( buffer.getFloat());
            }
            
            // 关闭文件
            fileChannel.close();
            fis.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
    }
    public static void fileReadDouble (ArrayList<Double> num, String filename)//fileName的后缀为.bin
    {
        File file= new File(filename);
        if(!file.exists())
        {
            System.out.println("file not exist");
            return;
        }
        try
        {
            // 打开文件
            FileInputStream fis = new FileInputStream(filename);
            FileChannel fileChannel = fis.getChannel();

            // 将文件映射到内存
            MappedByteBuffer buffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0, fileChannel.size());
            buffer.order(ByteOrder.LITTLE_ENDIAN); // 设置缓冲区的字节顺序为小端字节序

            // 读取数据
            while (buffer.hasRemaining()) {
                num.add(buffer.getDouble());
            }
            
            // 关闭文件
            fileChannel.close();
            fis.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
    }
}