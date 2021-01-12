# easyFlow
easyFlow : An easy to do unsupervised analysis for mass and flow cytometry data.
From fcs data to statistics !

## How To easyflow

### Creation of a Virtual Machine (VM)

The complete easyFlow pipeline was tested on Linux Ubuntu operating system only (16.04 LTS 64bit). 

To create a VM running on Ubuntu, please follow the link :

https://github.com/VallinP/HowTo_VirtualMachine


### Install easyFlow package

devtools::install_github("VallinP/easyFlow")


### easyFlow analysis

Before to start, you have to perform a quick analysis of your data to :
- Compensate your data
- Select events of interrest (ie, live cells, CD45+ ....) 
- Save the workspace in the same folder of your fcs.


Now, you are ready to start an easyFlow analysis.


1- Initialization

  library(easyflow)
  easyFlow_initialization()


2- Setup parameters


3- Analysis


4- Finalizing analysis


5- Interpretation

