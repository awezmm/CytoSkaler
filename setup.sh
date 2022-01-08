#!/bin/bash

DIR="distribution/"
if [ -d "$DIR" ]; then
   echo "Error! A 'distribution' folder already exists in this directory. Either delete the 'distribution' folder or setup this program in a different directory."
   exit 1
fi


#curl https://people.cs.uchicago.edu/~aamohsin/distribution.tar.gz > distribution.tar.gz && tar xvzf distribution.tar.gz && rm distribution.tar.gz && cd distribution && mv * ../ && cd .. && rm -rf distribution
#cat DistributionParts/* > distribution.tar.gz && rm -rf DistributionParts && tar xvzf distribution.tar.gz && rm distribution.tar.gz && cd distribution && mv * ../ && cd .. && rm -rf distribution
cat DistributionParts/* > distribution.tar.gz && rm -rf DistributionParts && tar xvzf distribution.tar.gz && rm distribution.tar.gz && cp distribution/* . && rm -rf distribution


echo " " >> ~/.bashrc
echo "# MATLAB 2019A Runtime Path Append" >> ~/.bashrc
echo "export DYLD_LIBRARY_PATH=/Applications/MATLAB/MATLAB_Runtime/v96/runtime/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/sys/os/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/bin/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/extern/bin/maci64" >> ~/.bashrc
echo " " >> ~/.bashrc


echo " " >> ~/.cshrc
echo "# MATLAB 2019A Runtime Path Append" >> ~/.cshrc
echo "setenv DYLD_LIBRARY_PATH" >> ~/.cshrc
echo "Applications/MATLAB/MATLAB_Runtime/v96/runtime/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/sys/os/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/bin/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/extern/bin/maci64" >> ~/.cshrc
echo " " >> ~/.cshrc


echo " " >> ~/.bash_profile
echo "# MATLAB 2019A Runtime Path Append" >> ~/.bash_profile
echo "export DYLD_LIBRARY_PATH=/Applications/MATLAB/MATLAB_Runtime/v96/runtime/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/sys/os/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/bin/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/extern/bin/maci64" >> ~/.bash_profile
echo " " >> ~/.bash_profile

echo " " >> ~/.zshenv
echo "# MATLAB 2019A Runtime Path Append" >> ~/.zshenv
echo "export DYLD_LIBRARY_PATH=/Applications/MATLAB/MATLAB_Runtime/v96/runtime/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/sys/os/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/bin/maci64:/Applications/MATLAB/MATLAB_Runtime/v96/extern/bin/maci64" >> ~/.zshenv
echo " " >> ~/.zshenv


pwd > ~/.CytoSkalerInfo

