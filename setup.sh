#!/bin/bash

curl https://people.cs.uchicago.edu/~aamohsin/distribution.tar.gz > distribution.tar.gz

tar xvzf distribution.tar.gz

rm distribution.tar.gz

cd distribution

mv * ../

cd ..

rm -rf distribution

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


rm -f ~/.CytoSkalerInfo
pwd >> ~/.CytoSkalerInfo

