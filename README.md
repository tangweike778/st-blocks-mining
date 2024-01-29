# st-blocks-mining

# RPA PA RA
The startup function for RA, PA, and RPA is "Converge.main". It requires seven arguments to be input as follows: "InputPath OutputPath delimeter k pro1 pro2 pro3"  
InputPath:File path for input data  
OutputPath:Path for result output  
delimeter:delimeter  
k:Number of spatiotemporal blocks I hope to find  
pro1,pro2,pro3:The proportion of density, time, and space in score function f. The proportion of pro1 in f is pro1/(pro1 + pro2 + pro3).  
We have provided a test data in data folder. you can run the test data by entering the parameters "./data/test.txt ./output/test 10 1 1 1" in main function".  

# EA
The startup function for EA is "Converge.main". It requires seven arguments to be input as follows: "InputPath OutputPath delimeter k h pro1 pro2 pro3"  
InputPath:File path for input data  
OutputPath:Path for result output  
delimeter:delimeter  
k:Number of spatiotemporal blocks I hope to find  
h:How many intervals are the spatiotemporal range divided into.  
pro1,pro2,pro3:The proportion of density, time, and space in score function f. The proportion of pro1 in f is pro1/(pro1 + pro2 + pro3).  
We have provided a test data in data folder. you can run the test data by entering the parameters "./data/test.txt ./output/test 10 100 1 1 1" in main function.  
