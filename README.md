# Fixative modelling & EPG fitting
This respository contains the code associated with the manuscript titled:\
*A method to remove the influence of fixative concentration on post-mortem T2 maps using a Kinetic Tensor model*

Code is as follows:

*FixativeModelling.m*\
Code to perform KI and KT modelling

*EPG* folder\
Folder contains code for step one and step two of the EPG fitting as described in the manuscript. Step one estimates T2 and B1 from multi-echo TSE data. Step two estimates T2 only, with regularisation on T2 to prevent spurious results in regions of low B1 and low SNR. For further details please see the Supporting Information of the manuscript.

Fitting is initialised using the wrapper scripts.

The EPG fitting is based on code written by Matthias Weigel, and the depiction and discussion of Extended Phase Graphs in the following publication:

Weigel M. J Magn Reson Imaging 2015; 41: 266-295. DOI: 10.1002/jmri.24619\
*Extended Phase Graphs: Dephasing, RF Pulses, and Echoes - Pure and Simple*

The original EPG code which this fitting is based on can be obtained by contacting Matthias at epg@matthias-weigel.net
                                                                                                                           
