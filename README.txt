Source code to reproduce bioRxiv manuscript "Tracking the outbreak. An optimized delimiting survey strategy 
for Xylella fastidiosa in Alicante, Spain" by E. Lázaro et al., 2020.

For questions, comments or remarks about the code please contact E. Lázaro (lazaro_ele@gva.es).

Code has been structured in three main folders, 1) "code" and 2) "data". 

Folder "code" contains all the R code to reproduce all the outputs (tables and figures) 
showed in the manuscript and in the supplementary material file. It is organised by means 
of the following scripts which must be executed in the order displayed:

- a) three_strategy_methods.R --> Contains all the R code to implement the three-phase delimiting strategy.

- b) two_strategy_methods.R --> Contains all the R code to implement the two-phase delimiting strategy.

- c) incidence.R --> Contains all the R code to implement the incidence modeling of the whole database and also
for the datasubsets created to chech sampling intensity.

- d) article.R --> Contains all the R code to reproduce article outputs. It requires the objects 
created by the execution of the three previous scripts (three_strategy_methods.R, 
two_strategy_methods.R and incidence.R) which are prepared to be stored in the folder "results".

- e) supplementary.R --> Contains all the R code to reproduce supplementary material outputs. It is based on
the objects created by the execution of the three previous scripts (three_strategy_methods.R, 
two_strategy_methods.R and incidence.R) which are prepared to be stored in the folder "results".

- f) Bdiclcpomodel_stack.R --> Contains a function to carry out model selection based on WAIC criterion.

Folder "data" contains all the datasets and spatial objects used to implement article methodology. 
It is organised in different subfolders according to the nature of the objects. All the "data" files are 
called in the different "code" scripts and some of them are created by the execution of them.

IMPOTANT: Previously to execute "code" scripts it is recommendable to create an additional folder
named "results" as te same level of "code" and "data" folders. This folder "results" can be organised as
follows to store correctly output objects.

	|		     |models(subsubfolfer)	
	|incidence(subfolder)|
	|		     |data_subsets_models(subsubfolder)	
results |three_phase(subfolder)
(Folder)|
	|two_phase(subdolder)
	|	
	|

Once the "results" folder had been created, outputs are ready to be stored as "code" scripts are executed.
Note that figures, are prepared to be stored in the default path  originated by the working 
directory allocation, you can create and additional subfolder named "img" within the "results" folder, 
change the storing pathfile in the corresponding figure code lines and store them all together.

