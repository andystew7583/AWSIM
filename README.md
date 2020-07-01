# AWSIM
AndreW Stewart Isopycnal Model

Because vanity has no limits.

This package includes:
1. The core model code, written in C (see model_code/).
2. A library of useful Matlab files for configuring and analyzing experiments (see matlab_common/).
3. Some example experiment configuration scripts (see examples/).

To get stared open one of the example experiment configurations using Matlab ('examples/ACC' is probably the simplest, and will be assumed for the purposes of this example). Create an experiment using setparams_XXX.m as follows:

>> setparams_ACC('.','test_ACC')

This will create an experiment folder 'test_ACC' within the examples/ACC/ directory. You can create the experiment elsewhere by changing the first argument of setparams_ACC. Now use a unix terminal to navigate to your experiment folder and compile using.

>> sh Make.sh

Th
