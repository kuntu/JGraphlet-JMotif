# JMotif
## Introduction
A Java program that
1. Provides a fast algorithm to compute motif census in large scale network
2. Given a graph, provides a analytical model to compute average motif frequencies in random graphs (faster than simulation methods) with the same 
  * in/out-degree-pair sequence 
  * number of nodes and edges 
  * number of mutual dyad(reciprocal edge), asymetric dyad (asymetric edge) and null dyad (null edge). 
3. Provides methods to generate random graphs with the same 
* in/out-degree-pair sequence, 
* (reciprocal, incoming, outgoing)-degree-tripplet sequence, 
* number of nodes and edges 
* number of mutual dyad(reciprocal edge), asymetric dyad (asymetric edge) and null dyad (null edge).
4. Provides simulation methods to compute motif census in random graphs with the same 
* an in/out-degree-pair sequence, 
* (reciprocal, incoming, outgoing)-degree-tripplet sequence, 
* number of nodes and edges 
* number of mutual dyad(reciprocal edge), asymetric dyad (asymetric edge) and null dyad (null edge).

The 16 triad are named using M-A-N labeling (refer to the end of this readme file and triadLabeling.png), the census is in the order of [003 012 102 021D 021U 021C 111D 111U 030T 030C 201 120D 120U 120 C 210 300]
## Running example
### 1. Compile files of motif census for static network and temporal networks
For example, in command window under folder bin/

javac -d ./ -sourcepath ../ ../StaticGraphExperiment.java

javac -d ./ -sourcepath ../ ../TemporalGraphExperiment.java

### 2. Prepare configuration file for batch processing of data files
1. create a file to indicate whether to compute motif census in a network or to compute average motif frequencies in random graph. Refer to files configs/commands/ for more details
2. create a file to list all the networks to be computed. Refer to files in configs/dataSources/ for more details
3. prepare the data file format according to instruction in configs/dataSources/
- make sure the paths are correct.
### 3. obtain motif census (in network or random graphs)/generate random graphs/ using configuration files

Example: in command window under folder bin/

java StaticGraphExperiment ../configs/dataSources/staticGraphList.cfg ../configs/commands/triadStaticGraph.cfg

java TemporalGraphExperiment ../configs/dataSources/dynamicGraphList.cfg ../configs/commands/triadDynamicGraph.cfg

## M-A-N Labeling of Triad
A M-A-N label consists
of three digits and an optional character. The first digit represents the number of
mutual connections between two nodes. The second digit is the number of
asymmetric connections. The last digit denotes the number of null dyads. The
optional character represents the directions of edges: "D" for down; "U" for up;
"T" for transition and "C" for cyclic. See triadLabeling.png Image for details
