# Fixative modelling & EPG fitting
This respository contains the code associated with the publication:\
*Tendler, Benjamin C., et al. "A method to remove the influence of fixative concentration on post-mortem T2 maps using a Kinetic Tensor model." bioRxiv (2020). https://doi.org/10.1101/2020.09.16.299784*


Code is as follows:

**FixativeModelling.m**\
Code to perform KI and KT modelling

**EPG folder**\
Folder contains code for step one and step two of the EPG fitting as described in the publication. Step one estimates T<sub>2</sub> and B<sub>1</sub> from multi-echo TSE data. Step two estimates T<sub>2</sub> only, with regularisation on T<sub>2</sub> to prevent spurious results in regions of low B<sub>1</sub> and low SNR. For further details please see the Supporting Information.

EPG fitting is initialised using the wrapper scripts *EPG_wrapper_step1.m* and *EPG_wrapper_step2.m*

The EPG fitting is based on the "cp_cpmg_epg_domain_fplus_fminus" MATLAB function written by Matthias Weigel as part of his EPG software, associated with the following publication:

Weigel M. J Magn Reson Imaging 2015; 41: 266-295. DOI: 10.1002/jmri.24619\
*Extended Phase Graphs: Dephasing, RF Pulses, and Echoes - Pure and Simple*

EPG software written by Matthias can be obtained by contacting Matthias at epg@matthias-weigel.net
