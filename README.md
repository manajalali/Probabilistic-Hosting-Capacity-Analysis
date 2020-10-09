
** Reference paper:

	S.  Taheri,  M. Jalali, V.  Kekatos,  and  L.  Tong,  “Fast Probabilistic Hosting Capacity Analysisfor Active Distribution Systems,” 
	IEEE   Trans.   Smart   Grid.,   2020,   (submitted).   [Online].   Available:https://arxiv.org/abs/2002.01980 

Notes: 

1- The Master code will solve the OPF using MP-OPF Algorithm of the paper. This sciprt uses the Pre_load() functions to generate the needed
   matrices and used the Theta_maker() funtion to generate the library of thetas. Other functions used in these three files are explained
   in the respective mfiles.

2- You can choose between the IEEE 123-bus system and the 8500-bus system by selecting the correct preload function (lines 12 and 13). 
   The IEEE 123-bus preload file requires the feeder data in table format in digital numbers. The IEEE 8500-bus file reads the excel 
   files downloaded from the IEEE website directly. Please note that you can use the 123-bus file to generate the feeder data for feeders 
   that have numerical bus IDs and the 8500-bus file to generate feeder matrices and relabeling feeders with string-type bus IDs.



