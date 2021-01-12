# easyFlow
easyFlow : An easy to do unsupervised analysis for mass and flow cytometry data.
From fcs data to statistics !

## How To easyflow

### Creation of a Virtual Machine (VM)

The complete easyFlow pipeline was tested only on Linux Ubuntu operating system (16.04 LTS 64bit). 

To create a VM running on Ubuntu, please follow these steps :

1/ Download and install a VM creator (ie VMware >16, or Virtualbox)

2/ Download Ubuntu 16.04 LTS 64bit at the offical website (https://releases.ubuntu.com/16.04/).

3/ Create a new virtual machine with at least 16Gb of allocated RAM, 4 cores CPU, 150Go of allocated hard drive disk space.

4/ Start the VM

5/ Create an admin account

6/ If necessary, change regional parameters (system settings >text entry >Add >France)

7/ Add terminal to your shortcuts

8/ Start a terminal

9/ Upgrade the OS by entering the following command: 

    sudo apt full-upgrade

10/Restart the VM


### Java installation

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


### Linux libraries installation

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


### R software and R packages installation

#### Install and update Linux libraries

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

#### CRAN packages

##### Install missing packages

    new.packages <- packages[!(as.vector(packages) %in% as.vector(installed.packages()[,"Package"]))]
    length(new.packages)
    if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

#########################

#### Packages Bioconductor 

##### Install Bioconductor

    source("https://bioconductor.org/biocLite.R")
    biocLite()                  ## R version 3.0 or later

##### Install missing packages

    new.packages <- packages[!(as.vector(packages) %in% as.vector(installed.packages()[,"Package"]))]
    length(new.packages)
    for (i in 1:length(new.packages)) {
      biocLite(new.packages[i], dependencies = TRUE)
    }

###################

#### Packages Github 

##### Install missing packages

    new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
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

Quit R software and start a terminal

    cd ./Documents/Scripts/RchyOptimyx/
    sudo su -c "R -e \"source('./install_RchyOptimyx.R')\""


### Install Nginx

    sudo apt-get update
    sudo apt-get -y install nginx

### Install RStudio 

Let’s install some pre-requisites:

    sudo apt-get -y install gdebi-core

Download Rstudio and install it.

    wget https://download1.rstudio.org/desktop/xenial/amd64/rstudio-1.2.1335-amd64.deb
    sudo gdebi rstudio-1.2.1335-amd64.deb

For rJava and ReporteRs installation you also need to do this:

    sudo apt-get install libxml2-dev
    sudo R CMD javareconf

Install fonts as well.

    sudo apt-get install t1-xfree86-nonfree ttf-xfree86-nonfree ttf-xfree86-nonfree-syriac xfonts-75dpi xfonts-100dpi

Install RStudio Server

Allow Firewall port 8787 :

    sudo ufw allow 8787

Download the latest RStudio Server — consult RStudio Downloads page to get the URL for the latest version. Then install the file you downloaded. These next two lines are using the latest version as of writing this post.

    wget https://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb 
    sudo gdebi rstudio-server-1.1.463-amd64.deb 

Done! By default, Rstudio hostname -I uses port 8787, so to access RStudio go to 
http://123.456.1.2:8787 and you should be greeted with an RStudio login page. (If you forgot what your droplet’s IP is, you can find out by running hostname -I) 

Verify Rstudio installation

    sudo rstudio-server verify-installation

### Install Shiny Server

You can safely skip this step if you don’t use Shiny and aren’t interested in being able to host Shiny apps yourself. But don’t forget that Shiny Server can also be used to host Rmarkdown files, not just shiny apps. This means that even if you don’t develop shiny apps you might still have a use for Shiny Server if you want to host interactive Rmarkdown documents.

To install Shiny Server, first install the shiny package:

    sudo su - -c "R -e \"install.packages('shiny', repos='http://cran.rstudio.com/')\""

(Note again that we’re installing shiny in a way that will make it available to all users, as I explained above).

Just like when we installed RStudio, again we need to get the URL of the latest Shiny Server from the Shiny Server downloads page, download the file, and then install it. These are the two commands using the version that is most up-to-date right now:

    wget https://download3.rstudio.org/ubuntu-12.04/x86_64/shiny-server-1.5.6.875-amd64.deb
    sudo gdebi shiny-server-1.5.6.875-amd64.deb

Shiny Server is now installed and running. Assuming there were no problems, if you go to http://123.456.1.2:3838/ you should see Shiny Server’s default homepage, which includes some instructions and two Shiny apps. If you see an error on the bottom Shiny app, it’s probably because you don’t have the rmarkdown R package installed (the instructions on the default Shiny Server page mention this). After installing rmarkdown in R, the bottom Shiny app should work as well. Don’t forget to install rmarkdown so that it will be available to all users as described above. I suggest you read through the instructions page at http://123.456.1.2:3838/.

### Install nodejs

Install & Update script
    
    Install curl
    sudo apt-get install libcurl4-openssl-dev
    sudo apt install curl

To install or update nvm, you can use the install script using curl:

    curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.34.0/install.sh | bash

The script clones the nvm repository to ~/.nvm and adds the source line to your profile (~/.bash_profile, ~/.zshrc, ~/.profile, or ~/.bashrc).

Note: If the environment variable $XDG_CONFIG_HOME is present, it will place the nvm files there.
    
    export NVM_DIR="${XDG_CONFIG_HOME/:-$HOME/.}nvm"
    [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh" # This loads nvm

Note: You can add --no-use to the end of the above script (...nvm.sh --no-use) to postpone using nvm until you manually use it.

You can customize the install source, directory, profile, and version using the NVM_SOURCE, NVM_DIR, PROFILE, and NODE_VERSION variables. Eg: curl ... | 
NVM_DIR="path/to/nvm". Ensure that the NVM_DIR does not contain a trailing slash.

NB. The installer can use git, curl, or wget to download nvm, whatever is available.

Note: On Linux, after running the install script, if you get nvm: command not found or see no feedback from your terminal after you type:

    command -v nvm

simply close your current terminal, open a new terminal, and try verifying again.

To install a specific version of node:
    
    nvm install 10.15.3 



