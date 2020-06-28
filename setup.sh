#!/bin/bash

curl https://people.cs.uchicago.edu/~aamohsin/distribution.tar.gz > distribution.tar.gz

tar xvzf distribution.tar.gz

rm distribution.tar.gz

cd distribution

mv * ../

cd ..

rm distribution



