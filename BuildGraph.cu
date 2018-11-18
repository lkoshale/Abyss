#include <stdio.h>
#include <cuda.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <stdbool.h>

// map gnome seq to graph
// read k-mer from file


/*
    mapping
    A -> 0;
    T -> 1;
    G -> 2;
    C -> 3;

*/

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
    




__global__ void buildGraph(unsigned long* Node, unsigned int* Edge, unsigned int N,int SA,int SB,int SC,int SD){
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(idx < N){
        unsigned int k = 4*idx;
        unsigned long a = Node[k]; 
        unsigned long b = Node[k+1]
        unsigned long c = Node[k+2];
        unsigned long d = Node[k+3];

        unsigned int e = 8*idx;

        unsigned long end = 3;
        unsigned long start = end<<62;

        //at begining
        for(unsigned long i=0;i<4;i++){
            //add last base of C to first of D
            //atleast size of 1 base 2 bits
            unsigned long a1,b1,c1,d1;
            a1=b1=c1=d1=0;

            if(SD>1 && SC>1)
                d1 = (d>>2) | ( (c & end)<<(SD-2));
            if(SC>1 && SB>1)
                c1= (c>>2) | ( (b & end)<<(SC-2));
            if(SB>1 && SA>1)
                b1 = (b>>2) | ( (a & end)<<(SB-2))
            
            //Add here diff bases
            if(SA>1)
                a1 = (a>>2) | (i<<(SA-2));
            else if(SB>1)
                b1 = (b>>2) | (i<<(SB-2));
            else if(SC>1)
                c1 = (c>>2) | (i<<(SC-2));
            else if(SD>1)
                d1 = (d>>2) | (i<<(SD-2));

            //search for a1,b1,c1,d1

        }

        //at end
        for(unsigned long i=0;i<4;i++){
            unsigned long a1,b1,c1,d1;
            a1=b1=c1=d1=0;
            unsigned long lim = 1;
            
            if(SA>1)
                a1 = (a<<2) | ( (b & start)>>62 );
            if(SB>1)
                b1 = (b<<2) | ( (c & start)>>62 );
            if(SC>1)
                c1 = (c<<2) | ( (d & start)>>62 );
            if(SD>1)
                d1 = (d<<2)| i ;
            
            //removing additional bits if not full limit
            if(SA>1 && SA<64)
                a1 = a1 & ((lim<<SA)-1);
            if(SA>1 && SA<64)
                b1 = b1 & ((lim<<SB)-1);
            if(SA>1 && SA<64)
                c1 = c1 & ((lim<<SC)-1);
            if(SA>1 && SA<64)
                d1 = d1 & ((lim<<SD)-1);

            //search now
        }
    } 
}
 

// each k-mer starts at 1,8,16,... mapping  array for edges

//takes comnd line arg K (k-mer)
int main(int argc,char* argv){

    int K,N;
    if(argc < 3){
        printf("Argument K missing !!\n");
        exit(0);
    }


    K = atoi(argv[1]); 
    unsigned int N = 100;       //taking it as given
    /********************
    *   Data structure
    *   Node representing a sequence => 4 64bit num
    *   adjacency                    => 8 32bit num
    *
    *********************/
    

    File* fptr = fopen("data.txt","r");

    // Array containing Nodes each 4 digit represents one seq
    unsigned long* Node = (unsigned long*)malloc(sizeof(unsigned long)*N*4);

    // Array contain edges for each node as 8 per Node
    unsigned int*  Edge = (unsigned int*) malloc(sizeof(unsigned int)*N*8);

 
    char buffer[100];
    unsigned int Nindex = 0;
    while(fscanf(fptr,"%s\n",buffer)!=NULL){
        if(strlen(buffer)==K){
            Pair p = genPair(buffer);
            Node[Nindex] = p.a;
            Node[Nindex+1]=p.b;
            Node[Nindex+2]=p.c;
            Node[Nindex+3]=p.d;
            Nindex+=4;
        }else{
            // handel out of order strings
        }

    }
 
    //****for gpu*******

    unsigned long* DNode;
    unsigned long DNsize = sizeof(unsigned long)*N*4;
    
    unsigned int* DEdge;
    unsigned long DEsize = sizeof(unsigned int)*N*8;
    
    cudaMalloc(&DNode,DNsize);
    cudaMalloc(&DEdge,DEsize);

    //copy node array to gpu
    cudaMemcpy(DNode,Node,DNsize,cudaMemcpyHostToDevice);

    //set all val to -1 
    cudaMemset(DEdge,-1,DEsize);

    //invoke kernel

    int numThreads = 1024;
    int numBlocks = (N+numThreads-1)/numThreads;



    //end kernel

    cudaMemcpy(h_edge,d_edge,size*8*sizeof(int),cudaMemcpyDeviceToHost);

    //node contains valid node and its adjacency info in 
    // edge array .


    //TODO
    // Traverse the graph and build the assembly
    


}
