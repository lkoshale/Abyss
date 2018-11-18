#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include<cstdlib>
#include<cstdio>

#define K 10


using namespace std;

int main(){
    
    string genome="";
    set<string>Read;

    string line;
    getline(cin,line);
    //ignore first line
    while(getline(cin,line)){
        genome+=line;
    }

    for(int i=0;i+K<genome.size();i++){
        string temp = genome.substr(i,K);
        Read.insert(temp);
    }

    for(auto it=Read.begin();it!=Read.end();it++){
        cout<<*it<<"\n";
    }

    return 0;   
}