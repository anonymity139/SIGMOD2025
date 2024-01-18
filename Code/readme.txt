Readme (Interactive Learning for Diverse Top-k Set)
=========================
This package contains all source codes for 
a. Algorithm TDIA 
	1. It only works for the special case of IDT.
	2. The code is in folder TDIA
b. Algorithm HDIA 
	1. It works for the general case of IDT. 
	2. The code is in folder TDIA
c. Algorithm RH 
	1. It is an adapted existing algorithm.
	2. The code is in folder RH.
d. Algorithm ActiveRanking 
	1. It is an adapted existing algorithm.
	2. The code is in folder ActiveRanking.
e. Algorithm Preference-Learning 
	1. It is an adapted existing algorithm.
	2. The code is in folder PreferenceLearning.
f. Algorithm UH-Simplex 
	1. It is an adapted existing algorithm.
	2. The code is in folder Simplex.
g. Algorithm singlePass 
	1. It is an adapted existing algorithm.
	2. The code is in folder singlePass.

Make sure there is a folder called "input/", a folder called "output/" and a file called 
"config.txt" under the working directory.
They will be used for storing the input/output files, some intermediate results and the 
input parameters.

Usage Step
==========
a. Compilation
	mkdir build
	cd build
	cmake ..
	make

	You will need to install the GLPK package (for solving LPs) at first.
	See GLPK webpage <http://www.gnu.org/software/glpk/glpk.html>.
	Then update the path in CMakeLists.txt
		set(INC_DIR /usr/local/Cellar/glpk/5.0/include)
		set(LINK_DIR /usr/local/Cellar/glpk/5.0/lib)
	Update path "/usr/local/Cellar/glpk/5.0" to the path you install the GLPK package
	
b. Execution
	./run

c. Config
	The config file contains the input parameters (whose format will be described in Appendix A).

c. Input
	The input file contains the dataset (whose format will be described in Appendix B).
	
d. Output
	The output will be shown on the console (whose format will be described in Appendix C).

Example
=======
Sample input (input/Anti2d100k.txt) are provided. The dataset is described by two scoring attributes. 
Try: ./run



Appendix A. Format of Config File
------------------------------------
The format is: AlgName DatasetName k e[1] e[2] ... e[d]
AlgName - the name of the algorithm
DatasetName - the name of the dataset
k - parameter k indicating the size of the output
u[1] - the first attribute value of the user's utility vector
u[2] - the second attribute value of the user's utility vector
...
u[d] - the d-th attribute value of the user's utility vector
For example, you might see
-----------------------
TDIA  Anti2d100k.txt  10  0.5 0.5 
-----------------------


Appendix B. Format of Input File
------------------------------------
The format of the first line is: n d
n - the dataset size, integer
d - the number of scoring attributes in the dataset, integer
The format of the following n lines is
-----------------------------------------------------------------
<scoring attribute 1> <scoring attribute 2> ... <scoring attribute d> <sensitive attribute> 
-----------------------------------------------------------------
Each line corresponds to a point.


Appendix C. Format of Console Output
-----------------------------------------------------------------
The format of the output is
-----------------------------------------------------------------------------------
|      Algorithm | # of Questions |           Time |       Accuracy | Point #ID |
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
|   Ground Truth |              - |              - |              - |    PtID-0 |
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
|        AlgName |        Q-count | Execution Time |       Accuracy |    PtID-1 |
-----------------------------------------------------------------------------------
where AlgName is the name of the algorithm,
Q-count is the number of questions asked by the algorithm,
Execution Time is the execution time of the algorithm,
PtID-0 is the tuples' ID of the user's favorite point (ground truth),
PtID-1 is the tuples' ID of the user's favorite point returned by the algorithm,
Accuracy is the utility similarity of tuples between PtID-0 and PtID-0. 

For example, you might see:
-----------------------------------------------------------------------------------
|      Algorithm | # of Questions |           Time |       Accuracy | Point #ID |
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
|   Ground Truth |              - |              - |              - |     75145 |
|              - |              - |              - |              - |      9607 |
|              - |              - |              - |              - |     71679 |
|              - |              - |              - |              - |     40073 |
|              - |              - |              - |              - |     66976 |
|              - |              - |              - |              - |     61970 |
|              - |              - |              - |              - |     26812 |
|              - |              - |              - |              - |     35040 |
|              - |              - |              - |              - |     79884 |
|              - |              - |              - |              - |     11880 |
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
|           TDIA |              7 |       1.019308 |       1.000000 |     75145 |
|              - |              - |              - |              - |     66976 |
|              - |              - |              - |              - |     71679 |
|              - |              - |              - |              - |      9607 |
|              - |              - |              - |              - |     61970 |
|              - |              - |              - |              - |     40073 |
|              - |              - |              - |              - |     26812 |
|              - |              - |              - |              - |     35040 |
|              - |              - |              - |              - |     11880 |
|              - |              - |              - |              - |     79884 |
-----------------------------------------------------------------------------------
