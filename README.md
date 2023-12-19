              
                  ___  ____   _ _____ _____ _   _______ ___________ _____ _____ _____ ______ 
                  |  \/  | | | /  ___|  ___| \ | |  _  \_   _| ___ \  _  /  ___|  ___|| ___ \
                  | .  . | | | \ `--.| |__ |  \| | | | | | | | |_/ / | | \ `--.| |__  | |_/ /
                  | |\/| | | | |`--. \  __|| . ` | | | | | | |  __/| | | |`--. \  __| |    / 
                  | |  | | |_| /\__/ / |___| |\  | |/ / _| |_| |   \ \_/ /\__/ / |____| |\ \ 
                  \_|  |_/\___/\____/\____/\_| \_/___/  \___/\_|    \___/\____/\____(_)_| \_|
              
                                                                                           
README

Author: Gustaw Eriksson

Date: 19-12-2022

Version: 1.0

Contact: gustaw.eriksson@ki.se

# Introduction
Muscle, endometrium and adipose in R (MUSENDIPOSE.R) is a downstream data analysis pipeline based on Seurat developed to analysis CellRanger output for a on-going research study applying 10x Genomics based single-nuclei (sn) RNA-seq on endometrium tissue. The pipeline has been developed to handle datasets of large numbers of cells from multiple study groups. Depending on the size of the dataset, high performance cluster environment has to be used to run the pipeline. 

The pipeline is still under development and two parts of it is available in this repository. Each part and script is described below:

### Project directory
Create a project directory and copy the scripts to it.

### Input directory
Create a data directory in the project directory as such: MUSENDIPOSER/Data/. Within the /Data directory, create a directory named /CellRanger_Count.  

In Data/CellRanger_Count, create a folder for each sample (Data/CellRanger_Count/Sample_X) containing Cell Ranger output for each sample and the required "filtered_feature_bc_matrix" folder which is generated after running Cell Ranger.    Create the directory by running:

```
mkdir MUSENDIPOSER  
mkdir MUSENDIPOSER/Data  
mkdir MUSENDIPOSER/Data/CellRanger_Count  
cp -r CellRanger_sample_X MUSENDIPOSER/Data/CellRanger_Count
```
Or by running Create_directory.sh script from your selected directory:

```
bash Create_directory.sh
```
  
The directory structure should be:

/MUSENDIPOSER  
/MUSENDIPOSER/Data  
/MUSENDIPOSER/Data/CellRanger_Count

### Output directory
The scripts in the pipeline will automatically generate a /Output directory in the project directory. Within it, output folders for each script and analysis will be outputted.

### Running the scripts
Now you are ready to run the scripts. Run the scripts in the numbered order and make sure to manually check the output between each scripts.

# Acknowledgement

The computations were enabled by resources provided by the Swedish National Infrastructure for Computing (SNIC) at UPPMAX partially funded by the Swedish Research Council through grant agreement no. 2018-05973. 
