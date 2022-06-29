# How to install FEniCS?

NOTE: in May 2022, a newer version of FEniCS - FEniCSx version 0.4 - was released. However, the work in this repository was developed using the lastest stable release of legacy FEniCS version 2019.1.0, which was released in April 2019.

A convenient way to install FEniCS across operating platforms is through Anaconda. [Anaconda](https://www.anaconda.com/) is a good cross-platform computing environment to perform Python/R data science and machine learning applications with thousands of open-source packages and libraries. 

Alternative ways to install FEniCS could be through Ubuntu, Docker, or building from source. Details could be found in this [link](https://fenicsproject.org/download/archive/).

## 1. To install Anaconda on Linux (or Ubuntu subsystem on Windows),
### Download and install
```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
```
### Initiate conda with different shells

Right after the installation of Anaconda, you will be asked to initiate conda with the command
```
conda init
```
Or, you can use 
```
conda init --all
```
to initiate conda for all the shells (bash, tcsh, fish, xonsh, zsh, powershell, etc.. NOTE, the command will modify .bashrc and/or .tcshrc). 

### Alternatives
You may want to initiate conda with the following lines
```
source $Anaconda_root_path"/etc/profile.d/conda.csh"
```
for tcsh or, 
```
source $Anaconda_root_path"etc/profile.d/conda.sh"
```
for bash shell.

### Disk usage
A. The installation of Anaconda may need a minimum of 4.5GB. <br />

## 2. To install FEniCS via Anaconda
After the installation of Anaconda, use the following command to create a new conda environment called 'fenicsproject'
```
conda create -n fenicsproject -c conda-forge fenics
source activate fenicsproject
```
## 3. To install FEniCS via Ubuntu
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
```
