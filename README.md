# JMotif
## Introduction
1. Provide a fast algorithm to compute motif census in large scale network
2. Provide a analytical model to compute average motif frequencies in random graphs give an in/out-degree-pair sequence. 
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

for example, in command window under folder bin/

- java TemporalGraphExperiment ../configs/dataSources/dynamicGraphList.cfg ../configs/commands/triadDynamicGraph.cfg

- java TemporalGraphExperiment ../configs/dataSources/dynamicGraphList.cfg ../configs/commands/triadDynamicGraph.cfg
