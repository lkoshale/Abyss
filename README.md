# Abyss


# TODO LIST
    - Get data for input i.e some 1000 k-mer with known genome.
    - Understand and find code for assembling DNA from de-bruijin graph.
    - Map k-mer to memory and build de-bruijin graph on GPU.
    - Get DNA sequence using graph.

# To genrate Reads
 - to set K, change the #define in <b>GenDataFasta.cpp<b>
 - compile g++ GenDataFasta.cpp -std=c++11
 - Run: ./a.out < inputfile >ouputfile
