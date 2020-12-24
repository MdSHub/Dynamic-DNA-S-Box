# Dynamic-DNA-S-Box
# dynamic key Dependent 8x8 Substitution box.
The Design is  inspired by the Biological DNA techniques. The new DNA-base S-Box is used in order to reduce the calculations process in addition to make use of the unknown randomness of DNA bases in creating the S-box.
S-box are the only non linear component of the Block Ciphers. A lot of security of block ciphers depend on it.
Chaotic map is a method used to generate random sequence of floating point numbers.
X(n+1) = (C*X(n)) mod 1


Algorithm..
The initial value of x(0)=0.12345 and C=137(a random value between 1 to 256)
Input: //We generate two random strings with the help of chaotic map and consider them as our input to generate the Sbox. 
Output: new S-Box of size 8x8 .  
      Begin: 
Convert string1 to DNA coding. 
Convert string2 to DNA coding.
By subtraction table find new half string from string1.
By subtraction table find new half string from string2.
Implement xor operation between results of 3 and 4 to find new string3.
      r=‘A’ or ‘C’ or ‘T’ or ‘G’
For k =1 to  string3.length -1
     Begin
     if( k is odd ) then R=ascii( R ) Bitwise AND ascii(string3(k))
     else R=ascii( R ) bitwise OR ascii( string3(k) )
     End for k
7.  For i=0 to string1.length -1
     Begin 
               Run1=Run1 x  string1( i )
     End for i
For p=0 to string2.length -1
     Begin
               Run2=Run2 x string2( p )
     End for p
     9.    Run = Run1 + Run2
       SboxSeed = R
       Param=length( string1+ string2 )
  For id=0 to 3 
       Begin
           Sbox[ ] = SboxGenerator( SboxSeed,param)
           //Sbox is a linear array of size 256
           Sbox[ ] = bijective( Sbox)
           // bijective  is a function which remove inconsistencies 
           // of Sbox 
       End for id
11. Sbox[ ] = optimization( Sbox )
     //To optimize the non Linearity of Sbox by rotating rows        
    //and columns of Sbox with the help of chaotic map.   





     

     



















