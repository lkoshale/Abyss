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

int getIndex(char* kmer){
    int ans = 0;
    for(int i=0;i<strlen(kmer);i++){

    }

    return ans;
}


// each k-mer starts at 1,8,16,... mapping  array for edges

//takes comnd line arg K (k-mer)
int main(int argc,char* argv){

    int K;
    if(argc < 2){
        printf("Argument K missing !!\n");
        exit(0);
    }

    K = atoi(argv[1]);
    int size = (int)pow(4,K);       //all possible k-mers 

    //array that stores if the node is valid or not i.e present or not
    bool* h_node = (bool*)malloc(sizeof(bool)*size);

    File* fptr = fopen("data.txt","r");

    char buffer[100];
    while(fscanf(fptr,"%s\n",buffer)!=NULL){
        if(strlen(buffer)==K){
            int idx = getIndex(buffer);
            h_node[idx] = true;
        }else{

        }
    }

    // output of gpu kernel stored in it
    //array of int* each 8 consecutive block belongs to a node
    int* h_edge = (int*)malloc(sizeof(int)*size*8);
    
    //****for gpu*******

    int* d_edge;
    cudaMalloc(&d_edge,size*8*sizeof(int));
    bool* d_node;
    cudaMalloc(&d_node,size*sizeof(int));

    //copy node array to gpu
    cudaMemcpy(d_node,h_node,size*sizeof(int),cudaMemcpyHostToDevice);

    //set all val to -1 
    cudaMemset(d_edge,-1,size*8*sizeof(int));

    //invoke kernel

    //end kernel

    cudaMemcpy(h_edge,d_edge,size*8*sizeof(int),cudaMemcpyDeviceToHost);

    //node contains valid node and its adjacency info in 
    // edge array .


    //TODO
    // Traverse the graph and build the assembly
    



    
    
    



}