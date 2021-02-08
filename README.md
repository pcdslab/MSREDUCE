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
# Related Article
Muaaz Gul Awan, Fahad Saeed, MS-REDUCE: an ultrafast technique for reduction of big mass spectrometry data for high-throughput processing, Bioinformatics, Volume 32, Issue 10, 15 May 2016, Pages 1518–1526, https://doi.org/10.1093/bioinformatics/btw023
