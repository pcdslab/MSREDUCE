
*Bug Fixes*
In the most recent commit some bug fixes were performed in the quantization method.
The fixed version gives improved performance for the UPS2 data set (discussed int he published reprot).
Minor improvements in the performance over other datasets were also observed.

to compile :
javac msreduce.java


to execute:
java msreduce <input folder> <xx> <output folder>
xx = percentage data to be retained e.g. 10, 20, 30.