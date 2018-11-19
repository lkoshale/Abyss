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
	
	unsigned long e;
	int count = 0,i=0;
	int index = strlen(kmer) -1;
	//printf("%d\n",index);
	unsigned long shiftop = 1; 

	//converting rightmost part of kmer into long and storing in d
	for(i=0;i<4 && index>=0 ;i++){
		shiftop =1;
		count=0;
		e=0;
		for(;count<32 && index>=0;index--){
			switch(kmer[index]){
				case 'A': e+=shiftop*0;    
						break;
				case 'T': e+=shiftop*1;
						break;
				case 'G': e+=shiftop*2;
						break;
				case 'C': e+=shiftop*3;
						break;
			}
			shiftop=shiftop<<2;     //multiply by 4
			//  printf("%lu \n",shiftop);
			count++;
		}

		switch(i){
			case 0: ans.d = e;break;
			case 1: ans.c = e;break;
			case 2: ans.b = e;break;
			case 3: ans.a = e;break;
		}
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

    char buff[100]="CT";
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

     unsigned long lim = 1;
     unsigned long i = 3;
     unsigned long d1 = (p.d<<2)| i;
     printf("%lu\n",d1);   
     d1 = d1 & (  (lim<< (strlen(buff)*2))-1);
     printf("%lu\n",d1);
   // Pair t; t.d=65; t.c=0;t.b=0;t.a=0;
  //  printf("%s\n",getDNAstring(p,0,0,0,4));


    return 0;
}