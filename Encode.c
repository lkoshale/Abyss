#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>



//order d c b a

typedef struct Pair_{
    unsigned long a;
    unsigned long b;
    unsigned long c;
    unsigned long d;
} Pair;



Pair genPair(char* kmer){
    Pair ans;
    ans.a =0;ans.b=0;
    ans.c=0;ans.d=0;
    
    int count = 0;
    int index = strlen(kmer) -1;
    unsigned long shiftop = 1; 
    for( ;count<32 && index>=0;index--){
        switch(kmer[index]){
            case 'A': ans.d+=shiftop*0;    
                    break;
            case 'T': ans.d+=shiftop*1;
                    break;
            case 'G': ans.d+=shiftop*2;
                    break;
            case 'C': ans.d+=shiftop*3;
                    break;
        }
        shiftop=shiftop<<2;     //multiply by 4
      //  printf("%lu \n",shiftop);
        count++;
    }

    shiftop = 1;
    count = 0;
    for( ;count<32 && index>=0;index--){
        switch(kmer[index]){
            case 'A': ans.c+=shiftop*0;    
                    break;
            case 'T': ans.c+=shiftop*1;
                    break;
            case 'G': ans.c+=shiftop*2;
                    break;
            case 'C': ans.c+=shiftop*3;
                    break;
        }
        shiftop=shiftop<<2;       //multiply by 4
        count++;
    }

    shiftop = 1;
    count = 0;
    for( ;count<32 && index>=0;index--){
        switch(kmer[index]){
            case 'A': ans.b+=shiftop*0;    
                    break;
            case 'T': ans.b+=shiftop*1;
                    break;
            case 'G': ans.b+=shiftop*2;
                    break;
            case 'C': ans.b+=shiftop*3;
                    break;
        }
        shiftop=shiftop<<2;       //multiply by 4
        count++;
    }

    shiftop = 1;
    count = 0;
    for( ;count<32 && index>=0;index--){
        switch(kmer[index]){
            case 'A': ans.a+=shiftop*0;    
                    break;
            case 'T': ans.a+=shiftop*1;
                    break;
            case 'G': ans.a+=shiftop*2;
                    break;
            case 'C': ans.a+=shiftop*3;
                    break;
        }
        shiftop=shiftop<<2;      //multiply by 4
        count++;
    }

    return ans;
}

char getBase(int val){
    char ch = '$';
    switch(val){
        case 0: ch= 'A';
                break;
        case 1: ch= 'T';
                break;
        case 2: ch= 'G';
                break;
        case 3: ch= 'C';
                break;
    }

    return ch;
}

char* getDNAstring(Pair p,int sa,int sb,int sc,int sd){
    char* str = (char*)malloc(sizeof(char)*150);
    sprintf(str," ");
    unsigned long k = 3;
    for(int i=0;i<64 && i<sd;i+=2){
        int d = (p.d>>i)& k;
        sprintf(str,"%c%s",getBase(d),str);
    }

    for(int i=0;i<64&& i<sc;i+=2){
        int d = (p.c>>i)& k;
        sprintf(str,"%c%s",getBase(d),str);
    }

    for(int i=0;i<64 && i<sb;i+=2){
        int d = (p.b>>i)& k;
        sprintf(str,"%c%s",getBase(d),str);
    }

    for(int i=0;i<64 && i <sa;i+=2){
        int d = (p.a>>i)& k;
        sprintf(str,"%c%s",getBase(d),str);
    }


    return str;
}

//mod power of two as bitwise and 
// num % 4 =   num & 000011  //get last two bit

int main(){

    char buff[100]="AATC";
    Pair p = genPair(buff);

    printf("%lu %lu %lu %lu\n",p.a,p.b,p.c,p.d);

    p.c = 1;
    int k = strlen(buff)*2 - 2;
    unsigned long end = 3;
    unsigned long start = end<<62;
    //extract last val from privous and add at begining
    // same as below expr
    unsigned long a1 = (p.d>>2) | ( (p.c & end)<<k) ;
    printf("%lu %lu %lu\n",start,end,a1);
    
   // Pair t; t.d=65; t.c=0;t.b=0;t.a=0;
  //  printf("%s\n",getDNAstring(p,0,0,0,4));


    return 0;
}