#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <cuda.h>

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <stack>

using namespace std;
// map gnome seq to graph
// read k-mer from file

//GPU KERNEL here
__global__ void buildGraph(unsigned long* Node,unsigned int* Edge,unsigned int N,int SA,int SB,int SC,int SD,unsigned int* stt,unsigned int* enn){
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(idx < N){
        unsigned int k = 4*idx;
        unsigned long a = Node[k]; 
        unsigned long b = Node[k+1];
        unsigned long c = Node[k+2];
        unsigned long d = Node[k+3];

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
                b1 = (b>>2) | ( (a & end)<<(SB-2));
            
            //Add here diff bases
            if(SA>1)
                a1 = (a>>2) | (i<<(SA-2));
            else if(SB>1)
                b1 = (b>>2) | (i<<(SB-2));
            else if(SC>1)
                c1 = (c>>2) | (i<<(SC-2));
            else if(SD>1)
                d1 = (d>>2) | (i<<(SD-2));

            int m,l,r;
            l=0;
            r=(int)N-1;
		    while ( l <= r) 
		    { 
		        m = l + (r-l)/2; 
		  	
		        // Check if x is present at mid 
		        if (Node[4*m] == a1  && Node[4*m+1] == b1  && Node[4*m+2] == c1  && Node[4*m+3] == d1 ){
		        	Edge[8*idx+i]=m; 
		        	break;
		        }
		        else if(Node[4*m] >a1) r = m-1;
		        else if(Node[4*m] <a1) l = m +1;
		        else{
		        	if(Node[4*m+1] >b1) r=m-1;
		        	else if(Node[4*m+1] <b1) l = m +1;
		        	else{
		        		if(Node[4*m+2] >c1) r=m-1;
		        		else if(Node[4*m+2] <c1) l = m +1;
		        		else{
		        			if(Node[4*m+3] >d1) r=m-1;
		        			else l = m +1;
							//all possible case is already taken care.
		        		}
		        	}

		        }
		        
		    } 
        }
       // if(idx==8)printf("%lu\n",d);
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
            if(SB>1 && SA<64)
                a1 = a1 & ((lim<<SA)-1);
            if(SC>1 && SB<64)
                b1 = b1 & ((lim<<SB)-1);
            if(SD>1 && SC<64)
                c1 = c1 & ((lim<<SC)-1);
            if(SD<64)
                d1 = d1 & ((lim<<SD)-1);

            //search now
            
            int m,l,r;
            l=0;
            r=(int)N-1;
           // if(idx==8) printf("idx:%u\t = a1:%lu\tb1:%lu\tc1:%lu\td1:%lu\n",idx+1 ,a1,b1,c1,d1);
		    while ( l <= r) 
		    { 
		        m = l + (r-l)/2; 
		       // if(idx==8)printf(" %d ",m);
		  
		        // Check if x is present at mid 
		        if (Node[4*m] == a1  && Node[4*m+1] == b1  && Node[4*m+2] == c1  && Node[4*m+3] == d1 ){
		        	Edge[8*idx+4+i]=m; 
		        	
		        
		        //	if(idx==8) printf("found\n");
		        	break;
		        }
		        else if(Node[4*m] >a1) r = m-1;
		        else if(Node[4*m] <a1) l = m +1;
		        else{
		        	if(Node[4*m+1] >b1) r=m-1;
		        	else if(Node[4*m+1] <b1) l = m +1;
		        	else{
		        		if(Node[4*m+2] >c1) r=m-1;
		        		else if(Node[4*m+2] <c1) l = m +1;
		        		else{
		        			if(Node[4*m+3] >d1) r=m-1;
		        			else l = m +1;
							//all possible case is already taken care.
		        		}
		        	}
		        }
		  		

		    }
		    
        }

        if(Edge[idx*8]==UINT_MAX && Edge[idx*8+1]==UINT_MAX && Edge[idx*8+2]==UINT_MAX && Edge[idx*8+3]==UINT_MAX) *stt = idx;
        if(Edge[idx*8+4]==UINT_MAX && Edge[idx*8+5]==UINT_MAX && Edge[idx*8+6]==UINT_MAX && Edge[idx*8+7]==UINT_MAX) *enn = idx;

    } 
}



/*
    mapping
    A -> 0;
    C -> 1;
    G -> 2;
    T -> 3;

*/

// we are representing kmer as as a pair of 4 unsigned long
// at max there can be 128 base pairs in kmer 
typedef struct Pair_{      
    unsigned long a;
    unsigned long b;
    unsigned long c;
    unsigned long d;
} Pair;

class Vertex {
public:
    string contig;
    vector<Vertex*> Edges;

    Vertex(){
        contig="";
    }

    void addContig(string st1){
        contig+=st1[st1.size()-1];
    }

    void addEdge(Vertex* v){
        Edges.push_back(v);
    }
};



//global vars
unordered_map<unsigned int,Vertex*>Graph;
set<unsigned int> Removed;
int SA,SB,SC,SD;

void add_all_edges(unsigned int* Edge,Vertex* v,int i,int N);
void rec_merge(unsigned int* Edge,unsigned long* Node,int N,Vertex* v,int id,int KS);
string getDNA(unsigned long* Node,int id,int N);


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
				case 'C': e+=shiftop*1;
						break;
				case 'G': e+=shiftop*2;
						break;
				case 'T': e+=shiftop*3;
						break;
			}
			shiftop=shiftop<<2;     //multiply by 4
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
    
    int A=0,B=0,C=0,D=0,l;
    int l=2*K;

    for(int i=0;i<4 && l>0 ;i++){
    	int m=0;
    	if(l<=64) {m=l;l=0;}
    	else{ 
    		m=64;
    		l=l-64;
    	}
    	switch(i){
    		case 0: D=m;break;
    		case 1: C=m;break;
    		case 2: B=m;break;
    		case 3: A=m;break;
    	}	
    	
    }

    printf("%d %d %d %d\n",A,B,C,D );
	SA=A;SB=B;SC=C;SD=D;



    File* fptr = fopen("ReadsData.txt","r");

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

    unsigned int* start = (unsigned int*)malloc(sizeof(unsigned int));
	unsigned int* end = (unsigned int*)malloc(sizeof(unsigned int));

 
    //****for gpu*******
    unsigned long* DNode;
    unsigned long DNsize = sizeof(unsigned long)*N*4;
    
    unsigned int* DEdge;
    unsigned long DEsize = sizeof(unsigned int)*N*8;
    
    unsigned int* Dstart;
    unsigned int* Dend;

    cudaMalloc(&DNode,DNsize);
    cudaMalloc(&DEdge,DEsize);
    cudaMalloc(&Dstart,sizeof(unsigned int));
    cudaMalloc(&Dend,sizeof(unsigned int));

    //copy node array to gpu
    cudaMemcpy(DNode,Node,DNsize,cudaMemcpyHostToDevice);
    cudaMemcpy(Dstart,start,sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(Dend,end,sizeof(unsigned int),cudaMemcpyHostToDevice);


    //set all val to -1 
    cudaMemset(DEdge,UINT_MAX,DEsize);

    //invoke kernel

    int numThreads = 1024;
    int numBlocks = (N+numThreads-1)/numThreads;


    buildGraph<<<numBlocks,numThreads>>>(DNode,DEdge,N,SA,SB,SC,SD,Dstart,Dend);

    //end kernel
    cudaMemcpy(DEdge,Edge,DEsize,cudaMemcpyDeviceToHost);

    //node contains valid node and its adjacency info in 
    // edge array .


    //TODO
    // Traverse the graph and build the assembly
    
	for(unsigned int i=0;i<N;i++){
		
		set<unsigned int>::iterator sit;
		sit=Removed.find(i);
		if(sit!=Removed.end())
			continue;

		Vertex* v;
		unordered_map<unsigned int,Vertex*>::iterator it;
		it = Graph.find(i);
		if(it==Graph.end()){
			v = new Vertex;
			Graph.insert(pair<unsigned int,Vertex*>(i,v));
		}
		else{
			v = it->second;
		}

		//only single outgoing edge
		
		rec_merge(Edge,Node,N,v,i,K);
		
	}

	for(auto sit=Removed.begin();sit!=Removed.end();sit++){
		unordered_map<unsigned int,Vertex*>::iterator it;
		it=Graph.find(*sit);
		if(it!=Graph.end()){
			Graph.erase(it);
		}
	}


		//print graph
	cout<<"Graph\n";
	for(auto it=Graph.begin();it!=Graph.end();it++){
		cout<<it->second->contig<<" :";
		for(int j=0;j<it->second->Edges.size();j++){
			cout<<it->second->Edges[j]->contig<<" , ";
		}
		cout<<"\n";
	}	



	unordered_map<unsigned int,Vertex*>::iterator it2;
	it2=Graph.find(*start);
	if(it2==Graph.end()){
		cout<<"Error\n";
		exit(0);
	}
	Vertex* st = it2->second;
	Vertex* en = NULL;
	for(auto it=Graph.begin();it!=Graph.end();it++){
		if(it->second->Edges.size()==0 && it->second->contig.size()>0)
			en=it->second;
	}	

	if(en==NULL){
		cout<<"Error\n";
		exit(0);
	}
	 

	if(Graph.size()==1){
		for(auto it=Graph.begin();it!=Graph.end();it++){
			cout<<it->second->contig<<"\n";
			return 0;
		}	
	}

	//find the path
	//make dummy edge from end to start
	en->addEdge(st);

	//if more than one vertex
	cout<<"start tarversing\n";
	vector<Vertex*>path;
	stack<Vertex*> curr_path;
	Vertex* curr_v= st;
	curr_path.push(st);

	while(!curr_path.empty()){
		if(curr_v->Edges.size()>0){
			curr_path.push(curr_v);
			Vertex* next_v = curr_v->Edges[curr_v->Edges.size()-1];
			curr_v->Edges.pop_back();

			curr_v=next_v;
		}
		else{
			path.push_back(curr_v);
			curr_v =  curr_path.top();
			curr_path.pop();
		}
	}

	cout<<"path\n";
	for(int i=path.size()-1;i>=0;i--){
		cout<<path[i]->contig<<"\n";
	}

	return 0;



}

void rec_merge(unsigned int* Edge,unsigned long* Node,int N,Vertex* v,int id,int KS){
	
	set<unsigned int>::iterator sit;
	sit=Removed.find(id);
	if(sit!=Removed.end()){
		unsigned int parent;
		for(int j=0;j<4;j++){
			if(Edge[8*id+j]!=UINT_MAX){
				parent=Edge[8*id+j];
				break;
			}
		}
		unordered_map<unsigned int,Vertex*>::iterator itr;
		itr=Graph.find(parent);
		if(itr!=Graph.end()){
			//add all nodes to vertex
			for(int k=0;k<(itr->second)->Edges.size();k++){
				v->addEdge(itr->second->Edges[k]);
			}

			v->contig+=itr->second->contig.substr(KS);

			Graph.erase(itr);
			
			return;
		}

	}
	
	int out=0;
	unsigned int oid;
	for(int j=0;j<4;j++){
		if(Edge[8*id+4+j]!=UINT_MAX && out==0){
			oid=Edge[8*id+4+j];
			out++;
		}
	}


	bool flag_merge = false;

	if(out==1){
		int in=0;
		for(int j=0;j<4;j++){
			if(Edge[8*oid+j]!=UINT_MAX)
				in++;
		}

		if(in==1){
			// cout<<v->contig<<"-";
			//recurse
			if(v->contig.size()==0)
				v->contig = getDNA(Node,id,N);
			else
				v->addContig(getDNA(Node,id,N));
			
			// cout<<v->contig<<"\n";
			rec_merge(Edge,Node,N,v,oid,KS);
			Removed.insert(oid);
			flag_merge=true;
		}

	}

	if(flag_merge==false){	
		// cout<<v->contig<<"-";

		if(v->contig.size()==0)
			v->contig = getDNA(Node,id,N);
		else
			v->addContig(getDNA(Node,id,N));
		
	   // cout<<v->contig<<"\n";

		add_all_edges(Edge,v,id,N);	
	}	

}

void add_all_edges(unsigned int* Edge,Vertex* v,int i,int N){
	if(i>=N)
		return;
	
	for(int j=0;j<4;j++){
		if(Edge[8*i+4+j]!=UINT_MAX){
			Vertex* temp;
			unordered_map<unsigned int,Vertex*>::iterator itr;
			itr = Graph.find(Edge[8*i+4+j]);
			if(itr==Graph.end()){
				temp = new Vertex;
				Graph.insert(pair<unsigned int,Vertex*>(Edge[8*i+4+j],temp));
				
			}
			else{
				temp = itr->second;
				
			}	
			v->addEdge(temp);
		}
	}
}

string getDNA(unsigned long* Node,int id,int N){
	string dna="";
	unsigned long end = 3;
	unsigned long a,b,c,d;
	a=Node[4*id];
	b=Node[4*id+1];
	c=Node[4*id+2];
	d=Node[4*id+3];	
	for(int i=0;i< SD;i+=2){
		unsigned long temp = (d>>i) & end;
		char ch='$';
		switch(temp){
			case 0:ch='A';break;
			case 1:ch='C';break;
			case 2:ch='G';break;
			case 3:ch='T';break;
		}
		dna=ch+dna;
	}

	for(int i=0;i< SC;i+=2){
		unsigned long temp = (c>>i) & end;
		char ch='$';
		switch(temp){
			case 0:ch='A';break;
			case 1:ch='C';break;
			case 2:ch='G';break;
			case 3:ch='T';break;
		}
		dna=ch+dna;
	}

	for(int i=0;i< SB;i+=2){
		unsigned long temp = (b>>i) & end;
		char ch='$';
		switch(temp){
			case 0:ch='A';break;
			case 1:ch='C';break;
			case 2:ch='G';break;
			case 3:ch='T';break;
		}
		dna=ch+dna;
	}

	for(int i=0;i< SA;i+=2){
		unsigned long temp = (a>>i) & end;
		char ch='$';
		switch(temp){
			case 0:ch='A';break;
			case 1:ch='C';break;
			case 2:ch='G';break;
			case 3:ch='T';break;
		}
		dna=ch+dna;
	}

	return dna;
}
