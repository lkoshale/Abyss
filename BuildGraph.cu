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

// each k-mer starts at 1,8,16,... mapping  array for edges

//node array if k-mer is valid or not i.e represent a node of graph if 1

int main(int argc,char* argv){

    int K;
    if(argc < 2){
        printf("Argument K missing !!\n");
        exit(0);
    }

    K = atoi(argv[1]);
    int size = (int)pow(4,K);

    //array of node*
    int* d_edge = (int*)malloc(sizeof(int)*size*8);

    bool* d_node = (bool*)malloc(sizeof(bool)*size);

    File* fptr = fopen("data.txt","r");

    char buffer[100];
    while(fscanf(fptr,"%s\n",buffer)!=NULL){
        if(strlen(buffer)==K){
            int idx = getIndex(buffer);
            d_node[idx] = true;
        }else{

        }
    }

    // input as scanned in arrays
    



}