# MSREDUCE
## Description of the files :
 * **HCD-DS-2.tar.gz** : HCS-DS-2 data set.
 * **msreduce.java** : Java implementation of MS-REDUCE method.
  
for verification of the results we made use of crux toolkit, that can be found here : http://cruxtoolkit.sourceforge.net/ 
# Usage
a basic java installation is required to build and run, build using:

```javac msreduce.java```

run on a test dataset using:

```
java msreduce <input folder> <xx> <output folder> 

xx = percentage data to be retained e.g. 10, 20, 30. 

input folder = any of the unzipped test data set directories e.g. HCD-DS-2

output folder = path to output directory
```
