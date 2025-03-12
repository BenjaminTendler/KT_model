# Fixative modelling & EPG fitting
This respository contains the code associated with the publication:\
*Tendler, Benjamin C., et al. "A method to remove the influence of fixative concentration on postmortem T2 maps using a kinetic tensor model." Human Brain Mapping (2021).  https://doi.org/10.1002/hbm.25661*

Any questions please contact benjamin.tendler@ndcn.ox.ac.uk

Details are as follows:

**FixativeModelling.m**\
Code to perform KI and KT modelling as described in the publication.

**EPG folder**\
Folder contains code for step one and step two of the EPG fitting as described in the publication. Step one estimates T<sub>2</sub> and B<sub>1</sub> from multi-echo TSE data. Step two estimates T<sub>2</sub> only, with regularisation on T<sub>2</sub> to prevent spurious results in regions of low B<sub>1</sub> and low SNR. For further details please see the Supporting Information.

EPG fitting is initialised using the wrapper scripts *EPG_wrapper_step1.m* and *EPG_wrapper_step2.m*

The EPG fitting is based on the "*cp_cpmg_epg_domain_fplus_fminus.m*" MATLAB function written by Matthias Weigel, available at https://github.com/matthias-weigel/EPG.

