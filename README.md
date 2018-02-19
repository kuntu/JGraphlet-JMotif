# JMotif
## Introduction
1. Provide a fast algorithm to compute motif census in large scale network
2. Provide a analytical model to compute average motif frequencies in random graphs give an in/out-degree-pair sequence. 
- The 16 triad are named using M-A-N labeling (refer to the end of this readme file and triadLabeling.png), the census is in the order of [003 012 102 021D 021U 021C 111D 111U 030T 030C 201 120D 120U 120 C 210 300]
## Running example
### Compile files of motif census for static network and temporal networks
For example, in command window under folder bin/

javac -d ./ -sourcepath ../ ../TemporalGraphExperiment.java

javac -d ./ -sourcepath ../ ../StaticGraphExperiment.java
### Prepare configuration file for batch processing of data files
1. create a file to indicate whether to compute motif census in a network or to compute average motif frequencies in random graph. Refer to files configs/commands/ for more details
2. create a file to list all the networks to be computed. Refer to files in configs/dataSources/ for more details
3. prepare the data file format according to instruction in configs/dataSources/
- make sure the paths are correct.
### obtain motif census using configuration files

Example: in command window under folder bin/

- java TemporalGraphExperiment ../configs/dataSources/dynamicGraphList.cfg ../configs/commands/triadDynamicGraph.cfg

- java TemporalGraphExperiment ../configs/dataSources/staticGraphList.cfg ../configs/commands/triadStaticGraph.cfg

### M-A-N Labeling of Triad
A M-A-N label consists
of three digits and an optional character. The first digit represents the number of
mutual connections between two nodes. The second digit is the number of
asymmetric connections. The last digit denotes the number of null dyads. The
optional character represents the directions of edges: "D" for down; "U" for up;
"T" for transition and "C" for cyclic. See triadLabeling.png Image for details
