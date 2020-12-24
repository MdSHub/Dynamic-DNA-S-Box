import numpy as np
from math import *
from os import urandom as _urandom




# Code to generate random index values of S-box array 
def Randombits(k):
    numbytes = (k + 7) // 8                      
    x = int.from_bytes(_urandom(numbytes), 'big')
    return x >> (numbytes * 8 - k)
def SboxindexGenerator(a, b):
    stop=b+1
    start=a
    n = stop - start
    k = n.bit_length()
    r=Randombits(k) 
    while r >= n:
        r = Randombits(k)
    return start + r


#   Code of Non-Linearity given by Musheer Sir,in which we pass 1-Dimensional   S-box array of size 256 

    
def fwt(f):  
    import math
    wf = []
    for x in f:
        if x == 0:
            wf.append(1)
        else:  
            wf.append(-1)
    order = len(f)  # order = 2^n
    n = int(math.log(order, 2))
    size = int(math.floor(order / 2))
    while size >= 1:
        left = 0
        while left < order - 1:
            for p in range(int(size)):
                right = left + int(size)
                a = wf[left]
                b = wf[right]
                wf[left] = a + b
                wf[right] = a - b
                left = left + 1
            left = right + 1
        size = int(math.floor(size / 2))
    
    return wf


def bf_nonlinearity(f, n):
    import math
    fw = fwt(f)
    #find modulus of each element in Walsh transform
    for i in range(len(fw)):
        fw[i] = abs(fw[i])
    # nonlinearity from the Walsh transform
    x = ((2 ** (n - 1)) - (max(fw) / 2))
    # print"\tNL of function is",
    # print x
    return x
 
def binary1(num, length):
    binary_string_list = list(format(num, '0{}b'.format(length)))
    #print("num,binary String List",num,binary_string_list)
    return [int(digit) for digit in binary_string_list]
     
# function of NL which is called by our Main function on 1 dimensional S-box array of size 256 
def non_linearity(S):     
    import math
    order = len(S)
    n = int(math.log(order, 2))
    nl_array = []  # nl_array[] stores calculated NL for each function yi
    for bitno in range(n):
        f = []
        for index in range(order):  #for each element in Sbox
            binary_value = binary1(S[index], n)
            bit = binary_value[bitno]
            f.append(bit)
        
        bfnl = bf_nonlinearity(f, n)
        nl_array.append(bfnl)
    return sum(nl_array)/len(nl_array)







def check(Sbox,k):
    p=Sbox[k]
    k=k+1
    while(k<256):
        if(Sbox[k] == p):
            return False
        k=k+1
    return True 

def SboxGenerator(M,lenJ):
    
    
    Sbox=np.zeros(256 ,dtype=int)
    Sbox[0]=(M ^ lenJ) % 256
    for i in range(256):
        value=(M ^ (lenJ+i)) % 256  
        # F=i*M*lenJ % 255 + 1
        # F=(F*Z*i) % 255 + 1
        index = SboxindexGenerator(0,255)
        Sbox[index]=value
        for k in range(256):      
            p=Sbox[k]
            if(check(Sbox,k)):
                Sbox[k]=p
            else:
                Sbox[k]=Sbox[k]+1 
        
    return Sbox


# for subtraction between adjacent characters of DNA string, input is of 32 bit and output is of 16 bit
def sub1(dnas,subtraction):
    s1=""
    i=1
    while(i < len(dnas)):
        a=dnas[i-1]
        b=dnas[i]
        if(a=='A'):
            a=0
        elif(a=='T'):
            a=1    
        elif(a=='C'):
            a=2
        else:
            a=3
        if(b=='A'):
            b=0
        elif(b=='T'):
            b=1    
        elif(b=='C'):
            b=2
        else:
            b=3        
        i=i+2
        s1+=subtraction[a][b]
    return s1


# for subtraction between  DNA string s1 and s2 of 16  bit each , and  output is of 16 bit
def sub2(s1,s2,subtraction):
    s3=""
    i=0

    while(i < len(s1)):
        a=s1[i]
        b=s2[i]
        if(a=='A'):
            a=0
        elif(a=='T'):
            a=1    
        elif(a=='C'):
            a=2
        else:
            a=3
        if(b=='A'):
            b=0
        elif(b=='T'):
            b=1    
        elif(b=='C'):
            b=2
        else:
            b=3        
        i=i+1
        s3+=subtraction[b][a]
    return s3    


# for addition between adjacent characters of DNA string, input is of 32 bit and output is of 16 bit
def add(dnas,addition):
    s1=""
    i=0
    while(i < len(dnas)):
        a=dnas[i]
        b=dnas[i+1]
        if(a=='A'):
            a=0
        elif(a=='T'):
            a=1    
        elif(a=='C'):
            a=2
        else:
            a=3
        if(b=='A'):
            b=0
        elif(b=='T'):
            b=1    
        elif(b=='C'):
            b=2
        else:
            b=3        
        i=i+2
        s1+=addition[a][b]
    return s1



# for XOR operation  between  DNA string s1 and s2 of 16  bit each , and  output is of 16 bit
def xor(s1,s2,xorTable):

    s3=""
    i=0
    l1=len(s1)
    l2=len(s2)
    if(l1 > l2):
        l=l2
    else:
        l=l1    

    while(i < l):
        a=s1[i]
        b=s2[i]
        if(a=='C'):
            a=0
        elif(a=='G'):
            a=1    
        elif(a=='A'):
            a=2
        else:
            a=3
        if(b=='C'):
            b=0
        elif(b=='G'):
            b=1    
        elif(b=='A'):
            b=2
        else:
            b=3        
        i=i+1
        s3+=xorTable[b][a]
    return s3 



# Convert binary form of string to DNA string
def dna(binarys,rules):
    a=""
    i=1
    j=0
    r=""

    while(i < len(binarys)):
        r=binarys[i-1]
        r+=binarys[i]
        p=int(r,2)
        a+=rules[p][j%8]
        i=i+2
        j=j+1
    return a


def binary(s1):
    a=""
    for x in s1:
        if(x!=0):
            a+=bin(x).replace("0b","0")
        else:
            a+=bin(x).replace("0b","00")
                
    return a

def missing(arr):
    i=0
    miss=np.zeros(257, dtype=int)
    while(i < 256):
        miss[arr[i]]+=1
        i+=1
    
    a=[]
    for k in range(257):
        if (miss[k] == 0):
            a.append(k)
    return a

def bijective(Sbox):

    arr=np.array(Sbox)
    arr=np.sort(arr)
    b=missing(arr)
    k=0
    b1=1
    while(k < 256):
        p=Sbox[k]
        q=0
        index=275
        while(q < 256):
            if(q == k):
                q+=1   
            elif(Sbox[q]== p):
                index=q   
                break
            else:
                q+=1

        if(index < 256):
            Sbox[index]=b[b1]
            b1+=1 
            
        k+=1
    return Sbox    



def Main():
    rules=np.array([['A','A','C','C','G','G','T','T'],['C','G','A','T','A','T','C','G'],
                   ['G','C','T','A','T','A','G','C'], ['T','T','G','G','C','C','A','A']])
    
    addition=np.array([['T','G','A','C'],['G','C','T','A'],
                      ['A','T','C','G'],['C','A','G','T']])
                      
    subtraction=np.array([['C','G','A','T'],['A','C','T','G'],
                         ['G','T','C','A'],['T','A','G','C']])
    xorTable=np.array([['C','G','A','T'],['G','C','T','A'],
                         ['A','T','C','G'],['T','A','G','C']])

    answer=np.zeros(1000,dtype=float)
    
        
    """ String generation with the help of chaos method
        x(n+1)=c*x(n) % 1

        where c=random integer between 1 to 256
        and x(0)=0.12345 """
    
    # Genrate first string s1
    
    
    
    x1=np.zeros(64,dtype=float)
    y1=np.zeros(64,dtype=int)
    x1[0]=0.12345
    s1=np.zeros(8,dtype=int)
    k=1
    c=SboxindexGenerator(1,256)
    if(c%2==0):
        c=c+1
    
    
    l=0
    while(k<64):
        
        x1[k]=(c*x1[k-1])%1
        if(x1[k] >= 0.5):
            y1[k]=1
        else:
            y1[k]=0    
        k+=1
        if(k %8 == 0):
            p=k
            n=""
            i=k-8
            while(i<k):
                
                n+=str(y1[i])
                i+=1
            s1[l]=int(n,2)
            l=l+1



    # Generate Second string s2        
    x2=np.zeros(64,dtype=float)
    y2=np.zeros(64,dtype=int)
    s2=np.zeros(8,dtype=int)
    k=1
    x2[0]=0.12345
    c=SboxindexGenerator(1,256)
    if(c%2==0):
        c=c+1
        
    
    l=0
    while(k<64):
        
        x2[k]=(c*x2[k-1])%1
        if(x2[k] >= 0.5):
            y2[k]=1
        else:
            y2[k]=0    
        k+=1
        if(k %8 == 0):
            p=k
            n=""
            i=k-8
            while(i<k):
                
                n+=str(y2[i])
                i+=1
            s2[l]=int(n,2)
            l=l+1

    
    bin1=binary(s1)
    bin2=binary(s2)
    
    dnas1=dna(bin1,rules)
    
    dnas2=dna(bin2,rules)
    
    #dnaSub1=add(dnas1,addition)    
    dnaSub1=sub1(dnas1,subtraction)
    
    #dnaSub2=add(dnas2,addition)    
    dnaSub2=sub1(dnas2,subtraction)
    
    dna3=xor(dnaSub1,dnaSub2,xorTable)    
    #dna3=sub2(dnaSub1,dnaSub2,subtraction)
    
    

    k=1
    root=dna3[0]



    while(k < len(dna3)):
        if(k%2==1):
            root=ord(root) & ord(dna3[k])
        else:
            root=ord(root) ^ ord(dna3[k])    
        k+=1
        root=chr(root)
        
    
    
    
    root=ord(root)
    

    Run1=1
    i=0
    s=""
    while(i<8):
        s+=str(s1[i])
        i+=1
    
    i=0
    while(i < len(s)):
        if(s[i] == '0'):
            i+=1
            continue
        Run1=Run1*(int(s[i]))
        
        i+=1
    z=""
    i=0   
    while(i<8):
        z+=str(s2[i])
        i+=1
    
    i=0
    Run2=1
    while(i < len(z)):
        if(z[i]=='0'):
            i+=1
            continue
        Run2=Run2*(int(z[i]))
        i+=1
    
    
    Run=Run1+Run2
    id=0
    max=0.000
    SboxSeed=root

    

    while(id <= 3):
        
        Sbox=SboxGenerator(SboxSeed,len(s)+len(z))
        Sbox=bijective(Sbox)
        a=non_linearity(Sbox)
        
        if(a>max):
            max=a
            max1=SboxSeed
            new=Sbox
            
        id+=1
        if(id==1):
            SboxSeed=Run
        elif(id==2):
            SboxSeed=Run+root
        else:
            SboxSeed=Run*root        
        
        
    print("S1 =",s1)
    print("S2 =",s2)
    # print("\n")    
    # print("The Generated Sbox ")
    # print("----------------------------------------------------------")    
    # print(new.reshape(16,16))

    # print("\n")
    # print("NL = : ",max)
    # print("\n")    
    
    return new     
        
                
if __name__=="__main__":
    Main()