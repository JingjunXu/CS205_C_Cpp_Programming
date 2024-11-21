#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define divideprecision 20
#define outputprecision 4
typedef struct bignum//Store decimal digits
{
    int Sign;
    int Dot_position;//Relative to the value in Data
    int Length;
    int* Data;
}Bignum;
typedef struct Scentificnum//Storing numbers from Scientific notation
{
    int Sign;
    double Data1;//Data before E or e
    double Data2;//Data after E or e
}SNnum;
int TestIfNum (const char* str);//Check whether the input of the number is correct
int TestIfEEEE (const char* str);//Check whether the expression form of digital machine is Scientific notation
void Initialize(Bignum* pnum,char* str);//Initialize decimal digits
int TestSign(const char* num);//Detect the sign of the number
int Dot_pos(const char* num);//Detect the position of the decimal point in the number
int CalcLen(char* str,int Dot_pos);//Calculate the length of data
int* CpyData(char* str,int Len);//Store the numbers into an int array
Bignum BigNumAdd(Bignum* pnum1,Bignum* pnum2);//Add large numbers
Bignum BigNumSubtract(Bignum* pnum1,Bignum* pnum2);//Subtract large numbers
Bignum BigNumMultiple(Bignum* pnum1,Bignum* pnum2);//Multiply large numbers
Bignum BigNumDivide(Bignum* pnum1,Bignum* pnum2);//Divide large numbers
void Printout(Bignum result);//Output decimal digits
void SNInitialize(SNnum* pnum,char* str);//Number that initializes Scientific notation
int TestEePosition(const char*str);//Detect the position of E or e of the number in Scientific notation
double CpyData1(char*str,int index);//Store the number before E or e
double CpyData2(char*str,int index,int Length);//Store the number after E or e
SNnum SNnumMultiple(SNnum* pnum1,SNnum* pnum2);//Multiplying Scientific notation Numbers
SNnum SNnumDivide(SNnum* pnum1,SNnum* pnum2);//Divide numbers by Scientific notation
void SNprintout(SNnum result);//Output the number of scientific and technological methods
int main(int argc,char** argv)
{
    if(TestIfNum(argv[1])==0||TestIfNum(argv[3])==0)//Check whether the input number is compliant
    printf("The input cannot be interpreted as numbers!");
    else if(TestIfEEEE(argv[1])==1||TestIfEEEE(argv[3])==1)//Check whether the input number is scientific and technological method
    {
        SNnum num1,num2,result;
        SNnum* pnum1=&num1,*pnum2=&num2;
        SNInitialize(pnum1,argv[1]);
        SNInitialize(pnum2,argv[3]);
        if(strcmp(argv[2],"*")==0)//Perform multiplication operation of Scientific notation
        {
            result=SNnumMultiple(pnum1,pnum2);
            printf("%s * %s = ",argv[1],argv[3]);
            SNprintout(result);
        }
        else if(strcmp(argv[2],"/")==0)//Perform division operation of Scientific notation
        {
            if(pnum2->Data1==0)//Determine whether the divisor is zero
            printf("A number cannot be divied by zero.");
            else
            {
                result=SNnumDivide(pnum1,pnum2);
                printf("%s * %s = ",argv[1],argv[3]);
                SNprintout(result);
            }
        }
        else if(strcmp(argv[2],"+")==0||strcmp(argv[2],"-")==0)
        printf("Recommend to enter numbers in decimal to compelete the operation");
        else
        printf("The input cannot be interpreted as an operator!");
    }
    else
    {
        Bignum num1,num2,result;
        Bignum* pnum1=&num1,*pnum2=&num2;
        Initialize(pnum1,argv[1]);
        Initialize(pnum2,argv[3]);
        if(strcmp(argv[2],"+")==0){
            if(pnum1->Sign==1&&pnum2->Sign==1){
                result=BigNumAdd(pnum1,pnum2);
            }
            else if(pnum1->Sign==-1&&pnum2->Sign==-1){
                result=BigNumAdd(pnum1,pnum2);
                result.Sign=-1;
            }
            else if(pnum1->Sign==-1&&pnum2->Sign==1){
                result=BigNumSubtract(pnum2,pnum1);
            }
            else{
                result=BigNumSubtract(pnum1,pnum2);
            }
            printf("%s + %s = ",argv[1],argv[3]);
            Printout(result);
        }
        else if(strcmp(argv[2],"-")==0){
            if(pnum1->Sign==1&&pnum2->Sign==1){
                result=BigNumSubtract(pnum1,pnum2);
            }
            else if(pnum1->Sign==-1&&pnum2->Sign==-1){
                result=BigNumAdd(pnum1,pnum2);
                if(result.Sign==-1)
                result.Sign=1;
                else
                result.Sign=-1;
            }
            else if(pnum1->Sign==-1&&pnum2->Sign==1){
                result=BigNumAdd(pnum2,pnum1);
                result.Sign=-1;
            }
            else{
                result=BigNumAdd(pnum1,pnum2);
            }
            printf("%s - %s = ",argv[1],argv[3]);
            Printout(result);
        }
        else if(strcmp(argv[2],"*")==0){
            result=BigNumMultiple(pnum1,pnum2);
            printf("%s * %s = ",argv[1],argv[3]);
            Printout(result);
        }
        else if(strcmp(argv[2],"/")==0){
            int sign=0;
            for(int i=0;i<pnum2->Length;i++)
            {
                if(pnum2->Data[i]!=0)
                {
                    sign=1;
                    break;
                }
            }
            if(sign==1){
                result=BigNumDivide(pnum1,pnum2);
                printf("%s / %s = ",argv[1],argv[3]);
                Printout(result);
            }
            else
            printf("A number cannot be divied by zero.");
        }
        else
        printf("The input cannot be interpreted as an operator!");
    }
    return 0;
}
int TestIfNum (const char* str)
{
    int sign=1,flag=0;
    if(strspn(str,"0123456789.+-Ee")!=strlen(str))//Check whether there are any undesirable characters
    return sign=0;
    char E_e[2]={'\0'};
    if(strstr(str,"e"))//Check whether it is Scientific notation (flag==0 means yes)
    strcpy(E_e,"e");
    else if(strstr(str,"E"))
    strcpy(E_e,"E");
    else//Detect whether the number entered in decimal is correct
    {
        int k=0,cnt=0;
        while(str[k]!='\0')
        {
            //The presence of plus and minus signs affects the position of the decimal point
            int signOfPM=0;
            if(str[k]=='+'||str[k]=='-')//Check whether the plus/minus sign is correct
            {
                signOfPM=1;
                if(k!=0)//Detection of the position and number of plus and minus signs
                return sign=0;
            }
            if(str[k]=='.')//Check whether the decimal point is correct
            {
                if(signOfPM==0)//Check whether the decimal point position is correct
                {
                    if(k==0)
                    return sign=0;
                }
                else if(signOfPM=1)
                {
                    if(k==1)
                    return sign=0;
                }
                cnt++;
                if(cnt>1)//Check whether the decimal point position is correct
                return sign=0;
            }
            k++;
        }
        flag=1;
    }
    if(flag==0)
    {
        int jmax=strspn(str,"0123456789+-."),cnt=0;
        if(jmax==0)
        return sign=0;
        for(int j=0;j<jmax;j++)//Check whether the string before e or E is correct
        {
            int signOfPM=0;
            if(str[j]=='+'||str[j]=='-')//Check whether the plus/minus sign is correct
            {
                signOfPM=1;
                if(j!=0)//Detection of the position and number of plus and minus signs
                return sign=0;
            }
            if(str[j]=='.')//Check whether the decimal point is correct
            {
                if(signOfPM==0)//Check whether the decimal point position is correct
                {
                    if(j==0)
                    return sign=0;
                }
                else if(signOfPM=1)
                {
                    if(j==1)
                    return sign=0;
                }
                cnt++;
                if(cnt>1)//Check whether the number is correct
                return sign=0;
            }
        }
        char *strr=(char*)calloc(100,sizeof(char));//Get the string after e or E
        strcpy(strr,strstr(str,E_e));
        if(strlen(strr) > 1)//Check whether there is a number after e or E
        {
            int i=0,cnt=0;
            while(strr[i] != '\0')
            {
                int signOfPM=0;
                if(strr[i]=='+'||strr[i]=='-')//Check whether the plus/minus sign is correct
                {
                    signOfPM=1;
                    if(i!=1)//Detection of the position and number of plus and minus signs
                    return sign=0;
                }
                if(strr[i]=='.')
                {
                    if(signOfPM==0)//Check whether the decimal point position is correct
                    {
                        if(i==0)
                        return sign=0;
                    }
                    else if(signOfPM=1)
                    {
                        if(i==1)
                        return sign=0;
                    }    
                    cnt++;
                    if(cnt>1)//Check whether the number of decimal points is correct
                    return sign=0;
                }
                if(i>0 && (strr[i]=='e'||strr[i]=='E'))//Check whether the number of e or E is correct
                return sign=0;
                i++;
            }
        }
        else
        return sign=0;
    }
    return sign;
}
int TestIfEEEE (const char* str)
{
    int i=0,sign=0;
    while(str[i]!='\0')
    {
        if(str[i]=='E'||str[i]=='e')
        return sign=1;
        i++;
    }
    return sign;
}
void Initialize(Bignum* pnum,char* str)
{
    pnum->Sign=TestSign(str);
    pnum->Dot_position=Dot_pos(str);
    pnum->Length=CalcLen(str,pnum->Dot_position);
    pnum->Data=CpyData(str,pnum->Length);
}
int TestSign(const char* str)
{
    int sign=1;
    if(str[0]=='-')
    return sign=-1;
    return sign;
}
int Dot_pos(const char* str)
{
    int dot_pos=0;
    while(str[dot_pos]!='\0')
    {
        if(str[dot_pos]=='.')
        break;
        dot_pos++;
    }
    if(str[0]=='+'||str[0]=='-')
    dot_pos--;
    return dot_pos;
}
int CalcLen(char* str,int Dot_pos)
{
    int Len=strlen(str);
    if(str[0]=='+'||str[0]=='-')
    Len--;
    if(Len!=Dot_pos)
    Len--;
    return Len;
}
int* CpyData(char* str,int Len)
{
    int* Data=(int*)calloc(Len,sizeof(int));
    for(int i=0,j=0;i<Len;j++)
    {
        if(str[j]=='+'||str[j]=='-'||str[j]=='.')
        continue;
        else
        {
            Data[i]=str[j]-'0';
            i++;
        }
    }
    return Data;
}
Bignum BigNumAdd(Bignum* pnum1,Bignum* pnum2)//By default, two positive numbers are added
{
    Bignum result;
    int *num1,*num2,*tempp,*sum;
    int Len=0;
    /*Align the decimals of the two numbers, store them in reverse order in num1 and num2,
      create sum to store the results, and initially update the position of the decimal point of the results*/
    if((pnum1->Length-pnum1->Dot_position)>=(pnum2->Length-pnum2->Dot_position))
    {
        num1=(int*)calloc(pnum1->Length,sizeof(int));
        num2=(int*)calloc(pnum2->Length+(pnum1->Length-pnum1->Dot_position)-(pnum2->Length-pnum2->Dot_position),sizeof(int));
        for(int i=pnum1->Length-1,j=0;i>=0;i--,j++)
        num1[j]=pnum1->Data[i];
        for(int i=pnum2->Length-1,j=(pnum1->Length-pnum1->Dot_position)-(pnum2->Length-pnum2->Dot_position);i>=0;i--,j++)
        num2[j]=pnum2->Data[i];
        if(pnum1->Dot_position>=pnum2->Dot_position)
        {
            Len=pnum1->Length;
            tempp=(int*)calloc(Len,sizeof(int));
            for(int i=0;i<pnum2->Length+(pnum1->Length-pnum1->Dot_position)-(pnum2->Length-pnum2->Dot_position);i++)
            tempp[i]=num2[i];
            num2=tempp;
            result.Dot_position=pnum1->Dot_position;
        }
        else
        {
            Len=pnum2->Length+(pnum1->Length-pnum1->Dot_position)-(pnum2->Length-pnum2->Dot_position);
            tempp=(int*)calloc(Len,sizeof(int));
            for(int i=0;i<pnum1->Length;i++)
            tempp[i]=num1[i];
            num1=tempp;
            result.Dot_position=pnum2->Dot_position;
        }
        sum=(int*)calloc(Len,sizeof(int));
    }
    else
    {
        num1=(int*)calloc(pnum1->Length+(pnum2->Length-pnum2->Dot_position)-(pnum1->Length-pnum1->Dot_position),sizeof(int));
        num2=(int*)calloc(pnum2->Length,sizeof(int));
        for(int i=pnum1->Length-1,j=(pnum2->Length-pnum2->Dot_position)-(pnum1->Length-pnum1->Dot_position);i>=0;i--,j++)
        num1[j]=pnum1->Data[i];
        for(int i=pnum2->Length-1,j=0;i>=0;i--,j++)
        num2[j]=pnum2->Data[i];
        if(pnum1->Dot_position>=pnum2->Dot_position)
        {
            Len=pnum1->Length+(pnum2->Length-pnum2->Dot_position)-(pnum1->Length-pnum1->Dot_position);
            tempp=(int*)calloc(Len,sizeof(int));
            for(int i=0;i<pnum2->Length;i++)
            tempp[i]=num2[i];
            num2=tempp;
            result.Dot_position=pnum1->Dot_position;
        }
        else
        {
            Len=pnum2->Length;
            tempp=(int*)calloc(Len,sizeof(int));
            for(int i=0;i<pnum1->Length+(pnum2->Length-pnum2->Dot_position)-(pnum1->Length-pnum1->Dot_position);i++)
            tempp[i]=num1[i];
            num1=tempp;
            result.Dot_position=pnum2->Dot_position;
        }
        sum=(int*)calloc(Len,sizeof(int));
    }
    int temp=0,carry=0;
    for(int i=0;i<Len;i++)//Add two numbers and complete carry
    {
        temp=num1[i]+num2[i]+carry;
        sum[i]=temp%10;
        carry=temp/10;
    }
    int* data;
    if(carry>0)//If the number of results increases, store the results in a larger space
    {
        result.Length=Len+1;//Update result data length
        result.Dot_position++;//Update the position of the decimal point of the result
        data=(int*)calloc(Len+1,sizeof(int));//Store data in positive order
        data[0]=carry;
        for(int i=1,j=Len-1;j>=0;i++,j--)
        data[i]=sum[j];
    }
    else
    {
        result.Length=Len;//Update result data length
        data=(int*)calloc(Len,sizeof(int));//Update the position of the decimal point of the result
        for(int i=0,j=Len-1;j>=0;i++,j--)//Store data in positive order
        data[i]=sum[j];
    }
    result.Data=data;
    result.Sign=1;
    return result;
}
Bignum BigNumSubtract(Bignum* pnum1,Bignum* pnum2)//Subtract two positive numbers by default
{
    Bignum result;
    Bignum *big,*small;//Distinguish the size of two numbers and judge the result symbol
    if(pnum1->Dot_position>pnum2->Dot_position)
    {
        result.Sign=1;
        big=pnum1;
        small=pnum2;
    }
    else if(pnum1->Dot_position<pnum2->Dot_position)
    {
        result.Sign=-1;
        big=pnum2;
        small=pnum1;
    }
    else
    {
        int sign=0,i=0;
        while(pnum1->Data[i]!='\0'||pnum2->Data[i]!='\0')
        {
            if(pnum1->Data[i]>pnum2->Data[i])
            {
                result.Sign=-1;
                big=pnum1;
                small=pnum2;
                sign=1;
                break;
            }
            else if(pnum1->Data[i]<pnum2->Data[i])
            {
                result.Sign=-1;
                big=pnum2;
                small=pnum1;
                sign=1;
                break;
            }
            i++;
        }
        if(sign==0)
        {
            if(pnum1->Length>pnum2->Length)
            {
                result.Sign=1;
                big=pnum1;
                small=pnum2;
            }
            else if(pnum1->Length<pnum2->Length)
            {
                result.Sign=-1;
                big=pnum2;
                small=pnum1;
            }
            else
            {
                int num[1]={0};
                result.Data=num;
                result.Sign=1;
                result.Length=1;
                return result;
            }
        }
    }
    /*Align the two numbers with decimals, store the large numbers in reverse order in minuend 
    and the small numbers in reverse order in subtrahend, so as to subtract the small numbers from the large numbers*/
    int *minuend,*subtrahend;
    int small_len=0,temp=0,borrow=0,big_len=0;
    if((big->Length-big->Dot_position)<(small->Length-small->Dot_position))
    {
        minuend=(int*)calloc((small->Length-small->Dot_position)+big->Dot_position,sizeof(int));//123.56
        for(int i=(small->Length-small->Dot_position)-(big->Length-big->Dot_position),j=big->Length-1;j>=0;j--,i++)
        minuend[i]=big->Data[j];
        subtrahend=(int*)calloc(small->Length,sizeof(int));
        for(int i=0,j=small->Length-1;j>=0;j--,i++)
        subtrahend[i]=small->Data[j];
        small_len=small->Length;
        big_len=(small->Length-small->Dot_position)+big->Dot_position;
    }
    else
    {
        minuend=(int*)calloc(big->Length,sizeof(int));
        for(int i=0,j=big->Length-1;j>=0;j--,i++)
        minuend[i]=big->Data[j];
        subtrahend=(int*)calloc((big->Length-big->Dot_position)+small->Dot_position,sizeof(int));
        for(int i=(big->Length-big->Dot_position)-(small->Length-small->Dot_position),j=small->Length-1;j>=0;j--,i++)
        subtrahend[i]=small->Data[j];
        small_len=(big->Length-big->Dot_position)+small->Dot_position;
        big_len=big->Length;
    }
    //Subtract decimals from large numbers, complete the borrowing, and store the results in minuend
    for(int i=0;i<small_len;i++)
    {
        temp=minuend[i]-subtrahend[i]-borrow;
        if(temp<0)
        {
            minuend[i]=temp+10;
            borrow=1;
        }
        else
        {
            minuend[i]=temp;
            borrow=0;
        }
    }
    if(borrow==1)
    minuend[small_len]-=borrow;
    //Optimize the result output and delete the redundant zero
    int findex=0,bindex=big_len-1,Length=0,dot_position=0;//Find the index in reverse order
    while(minuend[findex]==0)
    findex++;
    while(minuend[bindex]==0)
    bindex--;
    //Derive the positive order index from the reverse order index
    int Findex=big_len-1-bindex,Bindex=big_len-1-findex;
    if(Findex>=big->Dot_position)//Determine the position of the decimal point
    {
        dot_position=1;
        bindex++;
    }
    else
    dot_position=big->Dot_position-Findex;
    Length=bindex-findex+1;//Length of confirmation result
    int* data=(int*)calloc(Length,sizeof(int));//Store the optimized data in positive order
    for(int i=bindex,j=0;i>=findex;i--,j++)
    data[j]=minuend[i];
    result.Data=data;
    result.Length=Length;
    result.Dot_position=dot_position;
    return result;
}
Bignum BigNumMultiple(Bignum* pnum1,Bignum* pnum2)
{
    Bignum result;
    int* num1=(int*)calloc(pnum1->Length,sizeof(int));//Store data in reverse order
    int* num2=(int*)calloc(pnum2->Length,sizeof(int));
    for(int i=0,j=pnum1->Length-1;i<pnum1->Length;i++,j--)
    num1[i]=pnum1->Data[j];
    for(int i=0,j=pnum2->Length-1;i<pnum2->Length;i++,j--)
    num2[i]=pnum2->Data[j];
    int* data=(int*)calloc(pnum1->Length+pnum2->Length,sizeof(int));//Create result array
    for(int i=0;i<pnum1->Length;i++)//Multiply
    {
        for(int j=0;j<pnum2->Length;j++)
        data[i+j]+=num1[i]*num2[j];
    }
    for(int i=0;i<pnum1->Length+pnum2->Length-1;i++)//Carry forward
    {
        if(data[i]>9)
        data[i+1]+=data[i]/10;
        data[i] = data[i] % 10;
    }
    int index=0;//Optimize the result output and delete the redundant zero
    for(index=pnum1->Length+pnum2->Length-1;index>=0;index--)
    {
        if(data[index]==0)
        continue;
        else
        break;
    }
    int* final_data=(int*)calloc(index+1,sizeof(int));//Store the optimized results in positive order
    for(int i=0,j=index;i<=index;i++,j--)
    final_data[i]=data[j];
    result.Data=final_data;
    result.Sign=pnum1->Sign*pnum2->Sign;//Update symbols of data
    result.Length=index+1;//Update data length
    //Update decimal point position
    result.Dot_position=index+1-(pnum1->Length-pnum1->Dot_position)-(pnum2->Length-pnum2->Dot_position);
    return result;
}
Bignum BigNumDivide(Bignum* pnum1,Bignum* pnum2)
{
    Bignum result;
    int* divisor=(int*)calloc(pnum2->Length,sizeof(int));//Default divisor length is less than 20
    for(int i=0;i<pnum2->Length;i++)
    divisor[i]=pnum2->Data[i];
    //Store data of divisor and dividend in positive order
    int* dividend=(int*)calloc(divideprecision-1+pnum1->Length,sizeof(int));
    for(int i=0;i<pnum1->Length;i++)
    dividend[i]=pnum1->Data[i];
    int* data=(int*)calloc(divideprecision-1+pnum1->Length,sizeof(int));
    int cnt=0,borrow=0,left=0,sign=0,k=0,t=0,w=0,i=0,j=0;
    //Convert division into subtraction for operation
    for(i=pnum2->Length-1,w=0;i<divideprecision-1+pnum1->Length;i++,w++)
    {
        cnt=0,sign=0;       
        while(1)
        {
            if(left==0)
            {
                for(k=w,j=0;k<pnum2->Length+w;k++,j++)
                {
                    if(dividend[k]>divisor[j])
                    break;
                    else if(dividend[k]<divisor[j])
                    {
                        left=dividend[w];
                        sign=1;
                        break;
                    }
                }
            }
            if(sign==1)
            break;
            borrow=0;
            for(j=pnum2->Length-1,t=pnum2->Length-1+w;j>=0;j--,t--)
            {
                int temp=dividend[t]-divisor[j]-borrow;
                if(temp<0)
                {
                    dividend[t]=temp+10;
                    borrow=1;
                }
                else
                {
                    dividend[t]=temp;
                    borrow=0;
                }
            }
            cnt++;
            left-=borrow;
        }
        data[i]=cnt;
    }
    result.Dot_position=pnum1->Dot_position+pnum2->Length-pnum2->Dot_position;
    int index=0;
    for(index=0;index<result.Dot_position;index++)
    {
        if(data[index]!=0)
        break;
    }
    if(index==result.Dot_position)
    index--;
    int* final_data=(int*)calloc(divideprecision,sizeof(int));
    for(i=index,j=0;i<index+divideprecision;i++,j++)
    final_data[j]=data[i];
    result.Data=final_data;
    result.Length=divideprecision;
    result.Sign=pnum1->Sign*pnum2->Sign;
    result.Dot_position-=index;
    return result;
}
void Printout(Bignum result)
{
    if(result.Sign==-1)//Symbol of output result
    printf("-");
    /*The position of the last non-zero number is counted because the integer is output in decimal
    form due to the precision in large number division*/
    int index=result.Length-1;
    while(result.Data[index]==0&&index>=result.Dot_position)
    index--;
    if(index>=result.Dot_position)//Judge if it is a decimal
    {
        int *data=(int*)calloc(result.Dot_position+outputprecision,sizeof(int)),*dataa=NULL;
        int i=0;
        while(i<result.Length&&i<result.Dot_position+outputprecision)//Store data in positive order
        {
            data[i]=result.Data[i];
            i++;
        }
        //When the number after the decimal point of the result is greater than the output digit of the decimal point
        if(result.Dot_position+outputprecision<result.Length)
        {
            if(result.Data[result.Dot_position+outputprecision]>=5)//Rounding
            {
                data[result.Dot_position+outputprecision-1]++;
                int carry=0,temp=0;
                for(i=result.Dot_position+outputprecision-1;i>=0;i--)//Carry forward
                {
                    temp=data[i]+carry;
                    data[i]=temp%10;
                    carry=temp/10;
                }
                if(carry!=0)//If the number of results increases, store the results in a larger space
                {
                    dataa=(int*)calloc(result.Dot_position+outputprecision+1,sizeof(int));
                    dataa[0]=carry;
                    for(int i=1;i<result.Dot_position+outputprecision+1;i++)
                    dataa[i]=data[i-1];
                    result.Dot_position++;
                    data=dataa;
                }
            }
        }
        for(i=0;i<result.Dot_position;i++)//Output Decimal
        printf("%d",data[i]);
        printf(".");
        for(i=result.Dot_position;i<result.Dot_position+outputprecision;i++)
        printf("%d",data[i]);
    }
    else
    {
        for(int i=0;i<result.Dot_position;i++)//Output integer
        printf("%d",result.Data[i]);
    }
}
void SNInitialize(SNnum* pnum,char* str)
{
    int Ee_position=0,Length=0;
    Ee_position=TestEePosition(str);
    Length=strlen(str);
    pnum->Sign=TestSign(str);
    pnum->Data1=CpyData1(str,Ee_position);
    pnum->Data2=CpyData2(str,Ee_position,Length);
}
int TestEePosition(const char*str)
{
    int index=1;
    if(strstr(str,"e"))
    index=strspn(str,"0123456789.+-");
    else if(strstr(str,"E"))
    index=strspn(str,"0123456789.+-");
    return index;
}
double CpyData1(char*str,int index)
{
    char*str1;
    double result=0;
    str1=(char*)calloc(index,sizeof(char));
    for(int i=0;i<index;i++)
    str1[i]=str[i];
    result=atof(str1);
    return result;
}
double CpyData2(char*str,int index,int Length)
{
    char* str2;
    double result=0;
    if(Length==index)
    return result;
    str2=(char*)calloc(Length-index-1,sizeof(char));
    for(int i=index+1,j=0;i<Length;i++,j++)
    str2[j]=str[i];
    result=atof(str2);
    return result;
}
SNnum SNnumMultiple(SNnum* pnum1,SNnum* pnum2)
{
    SNnum result;
    result.Sign=pnum1->Sign*pnum2->Sign;
    result.Data1=pnum1->Data1*pnum2->Data1;
    result.Data2=pnum1->Data2+pnum2->Data2;
    while(result.Data1>=10)
    {
        result.Data1/=10;
        result.Data2++;
    }
    return result;
}
SNnum SNnumDivide(SNnum* pnum1,SNnum* pnum2)
{
    SNnum result;
    result.Sign=pnum1->Sign*pnum2->Sign;
    result.Data1=pnum1->Data1/pnum2->Data1;
    result.Data2=pnum1->Data2-pnum2->Data2;
    while(result.Data1<1)
    {
        result.Data1*=10;
        result.Data2--;
    }
    return result;
}
void SNprintout(SNnum result)
{
    if(result.Data1==0)//Determine whether the number is zero
    printf("0");//If it is zero, no output e or E is required
    else//Output base
    {
        printf("%.2lf",result.Data1);
        if(result.Data2!=0)//Output index
        printf("e%.2lf",result.Data2);
    }
}