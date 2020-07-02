# AWSIM
AndreW Stewart Isopycnal Model

Because vanity has no limits.

This package includes:
1. The core model code, written in C (see model_code/).
2. A library of useful Matlab files for configuring and analyzing experiments (see matlab_common/).
3. Some example experiment configuration scripts (see examples/).

The example experiments serve as an introduction to working with the model. To get started, open one of the example experiment configurations using Matlab. The 'examples/ACC' is probably the simplest, and will be assumed for the purposes of this example. Create an experiment using setparams_ACC.m as follows:

>> setparams_ACC('.','test_ACC')

This will create an experiment folder 'test_ACC' within the examples/ACC/ directory. You can create the experiment elsewhere by changing the first argument of setparams_ACC. Now use a unix terminal to navigate to your experiment folder and compile using:

>> sh Make.sh

This script compiles the code using gcc (you will need to have gcc installed for this to work, or you will need to modify this script for your C compiler). Then you can run the model using

>> sh Run.sh

As the model runs, you can visualize the model fields, for example using the anim.m script provided in the matlab_common/ folder:

>> anim('.','test_ACC','z',1,-1,-1)

(See anim.m for information on input parameters.) This will play a movie of the relative vorticity in the upper isopycnal layer of this test case. The anim.m script also serves as an example of how to analyze the model output. 
