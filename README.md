# easyFlow
easyFlow : An easy to do unsupervised analysis for mass and flow cytometry data.
From fcs data to statistics !

## How To easyflow

### Installation 

#### Creation of an Virtual machine 

The complete easyFlow pipeline was tested only on Linux Ubuntu operating system (16.04 LTS 64bit). 

To create a virtual machine running on Ubuntu, please follow these steps :

1/ Download and install a virtual machine creator (ie VMware >16, or Virtualbox)

2/ Download Ubuntu 16.04 LTS 64bit at the offical website (https://releases.ubuntu.com/16.04/).

3/ Create a new virtual machine with at least 16Gb of allocated RAM, 4 cores CPU, 150Go of allocated hard drive disk space.

4/ Start the VM (virtual machine)

5/ Create an admin account

6/ If necessary, change regional parameters (system settings >text entry >Add >France)

7/ Add console to your shortcuts

8/ Start console

9/ Upgrade the OS by entering the following command: 

    sudo apt full-upgrade

10/Restart the VM


#### Installation of the prerequisite

1/ Installing the Default JRE/JDK
The easiest option for installing Java is using the version packaged with Ubuntu. Specifically, this will install OpenJDK 8, the latest and recommended version.

First, update the package index.

    sudo apt-get update

You can install the JDK with the following command:

    sudo apt-get install default-jdk

2/ Setting the JAVA_HOME Environment Variable

Many programs, such as Java servers, use the JAVA_HOME environment variable to determine the Java installation location. To set this environment variable, we will first need to find out where Java is installed. You can do this by executing the same command as in the previous section:

    sudo update-alternatives --config java

Copy the path from your preferred installation and then open /etc/environment using nano or your favorite text editor.

    sudo nano /etc/environment

create a new line and paste : JAVA_HOME="PATH"
ie : JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"

Save and exit the file, and reload it.

    source /etc/environment

You can now test whether the environment variable has been set by executing the following command:

    echo $JAVA_HOME

This will return the path you just set.

3/ Restart the VM


#### Linux libraries installation

    sudo apt install libssl-dev curl
    sudo apt-get install libcurl4-openssl-dev 
    sudo apt-get install libcurl4-gnutls-dev
    sudo apt-get install libmagick++-dev
    sudo apt-get install libcairo2-dev
    sudo apt-get install unixodbc-dev
    sudo apt-get install libpq-dev
    sudo apt-get install libxml2-dev
    sudo apt-get install libv8-3.14-dev
    sudo apt-get install libudunits2-dev
    sudo apt-get install libpoppler-cpp-dev
    sudo apt-get install libnetcdf-dev
    sudo apt-get install libgsl-dev
    sudo apt-get install libgdal-dev
    sudo apt-get install libgeos-dev
    sudo apt-get install jags
    sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev 
    sudo apt-get install libgmp3-dev
    sudo apt-get install libfftw3-dev
    sudo apt-get install libcr-dev mpich mpich-doc
    
    sudo apt-get install aptitude
    sudo add-apt-repository ppa:marutter/c2d4u3.5

    sudo apt full-upgrade
    
    sudo apt-get upgrade


#### Rstudio and Rpackages installation

##### Install and update Linux libraries

    sudo apt-get install build-essential

    sudo apt-get install fort77

    sudo apt-get install xorg-dev

    sudo apt-get install liblzma-dev  libblas-dev gfortran

    sudo apt-get install gcc-multilib

    sudo apt-get install gobjc++

    sudo apt-get install aptitude

    sudo aptitude install libreadline-dev

    sudo aptitude install libcurl4-openssl-dev

    sudo apt-get install default-jdk

    sudo apt-get install texlive-latex-base

    sudo apt-get install libcairo2-dev 

#### Install r-base for ubuntu

    sudo apt-get install r-base

Update any R libraries installed via APT.

    sudo apt full-upgrade

    sudo apt-get upgrade

#### Install newest version of R from source

    wget https://cran.r-project.org/src/base/R-3/R-3.4.0.tar.gz

    ./configure --prefix=/home/"your profile"/R/R-3.4.0 --with-x=yes --enable-R-shlib=yes --with-cairo=yes

    make

NEWS.pdf file is missing and will make installation crash.

    touch doc/NEWS.pdf

    make install

Do not forget to update your PATH

    export PATH=~/R/R-3.4.0/bin:$PATH

    export RSTUDIO_WHICH_R=~/R/R-3.4.0/bin/R

Install libjpeg62

    sudo apt-get install libjpeg62

for some reason it prompted me to do 'sudo apt-get -f install' after. I did and it worked...


#### R packages Installation

##### Upadate librairies and start R

    sudo apt update

    sudo R

##### Load R packages list 

    packages <- read.csv("~/Documents/packages.csv")[,1]

#################
# CRAN packages #
#################

##### List missing packages

    new.packages <- packages[!(as.vector(packages) %in% as.vector(installed.packages()[,"Package"]))]

    length(new.packages)

##### Install missing packages 

    if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

#########################
# Packages Bioconductor #
#########################

##### Install Bioconductor

    source("https://bioconductor.org/biocLite.R")

    biocLite()                  ## R version 3.0 or later

##### List missing packages
'skip new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

    new.packages <- packages[!(as.vector(packages) %in% as.vector(installed.packages()[,"Package"]))]

    length(new.packages)

##### Install missing packages
'skip if(length(new.packages)) biocLite(new.packages, dependencies = TRUE)

    for (i in 1:length(new.packages)) {

      biocLite(new.packages[i], dependencies = TRUE)

    }

###################
# Packages Github #
###################

##### List missing packages

    new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

##### Install missing packages
'skip if(length(new.packages)) biocLite(new.packages, dependencies = TRUE)

    library("devtools")

    source("https://bioconductor.org/biocLite.R")

    biocLite("flowCore", dependencies = TRUE)

    biocLite("impute", dependencies = TRUE)

    library("devtools")

    install_github("ramnathv/rCharts")

    install_github("nolanlab/Rclusterpp")

    install_github("tchitchek-lab/SPADEVizR")

    install_github("kevinushey/data.table.extras")

    install_github("mul118/shinyGridster")

    install_github("kevinushey/Kmisc")

    install_github("vqv/ggbiplot")

    install_github('daattali/shinyjs')

# Quit R software

    cd ./Documents/Scripts/RchyOptimyx/

    sudo su -c "R -e \"source('./install_RchyOptimyx.R')\""




