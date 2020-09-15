function EPG_wrapper_step2(out_path,TSE_path,TE_path,mask_path,B1_path,med_path)
%%
% This code estimates T2 from the TSE signal via Extended Phase Graphs (EPG),
% with regularisation on T2 using the median T2 estimate across the sample from step 1

%Details of this fitting approach are described in the manuscript:
%A method to remove the influence of fixative concentration on post-mortem T2 maps using a Kinetic Tensor model

% The fitting code is based on the depiction and discussion of Extended
% Phase Graphs in the following publication:
% 
% Weigel M. J Magn Reson Imaging 2015; 41: 266-295. DOI: 10.1002/jmri.24619
% "Extended Phase Graphs: Dephasing, RF Pulses, and Echoes - Pure and Simple"  

% The EPG fitting function is based on the "cp_cpmg_epg_domain_fplus_fminus" 
% MATLAB function written by Matthias Weigel as part of his EPG software, 
% which can be obtained by contacting Matthias at epg@matthias-weigel.net 

% Wrapper script & modifications to perform fitting written by Benjamin Tendler 
% Contact benjamin.tendler@ndcn.ox.ac.uk

%%
%Output nifti files correspond to:
% T2map         Estimated T2 map
% T2map_SE      Estimated T2 standard error map
% T2_sig_recon  Reconstructed TSE signal
% T2_sig_res    TSE Residual (experimental - reconstructed)
% R2map         R2 map (1/T2map)
% fval          MATLAB lsqnonlin squared 2-norm residual at solution
% exitflag      MATLAB lsqnonlin exitflag

%Inputs are as follows:

% out_path      Output path where estimates will be saved
% TSE_path      Path to 4D TSE data, where each image corresponds to a different TE
% TE_path       Path to text file containing echo times in ms. NB the EPG model will used the difference between the first two echoes as the effective echo spacing
% mask_path     Path to mask of data (region where fitting will be performed)
% B1_path       Path to B1 map (from step 1)
% med_path      Path to text file containing median T2 value over dataset

%%
%Read in data
%TSE, TEs, mask, B1 & median
TSE=niftiread(TSE_path);fclose all;
TEs=importdata(TE_path)*1000;fclose all;
mask=niftiread(mask_path);fclose all;
B1=niftiread(B1_path);fclose all;
brain_med=importdata(med_path);fclose all;
%%
%Define output arrays
%T2 map, T2 error, Residual, Exit flag, reconstructed TSE signal
T2=zeros(size(TSE,1),size(TSE,2),size(TSE,3));
T2_SE=zeros(size(TSE,1),size(TSE,2),size(TSE,3));
fval=zeros(size(TSE,1),size(TSE,2),size(TSE,3));
exitflag=zeros(size(TSE,1),size(TSE,2),size(TSE,3));
TSE_recon=zeros(size(TSE));
%%
%Define fitting options
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off');
%%
%Define regularisation lambda
lambda=2;
%%
%Perform fitting (voxelwise)
for m=1:size(TSE,3)
    for l=1:size(TSE,2)
        for k=1:1:size(TSE,1)
            %Only fit voxels which are contained within mask
            if mask(k,l,m)~=0
                %Extract voxel data to be fit
                data=double(squeeze(TSE(k,l,m,:)));
                %Normalise data (to avoid fitting for S0)
                data_norm=data./(sum(data(:).^2)).^0.5; data_norm(isnan(data_norm))=0; data_norm(data_norm==inf)=0; data_norm(data_norm==-inf)=0;           
                %Define fitting function
                f=@(x)EPG_fitting_T2_reg(x,[data_norm;lambda*brain_med],double(B1(k,l,m)),TEs,1,(sum(data(:).^2)).^0.5,lambda);
                %Perform fitting
                [fit_out,fval(k,l,m),residual,exitflag(k,l,m),~,~,J]=lsqnonlin(f,50,0,inf,options);
                %Reconstruct TSE signal from fit parameters - multiply by normalised signal to be of same order as original data
                TSE_recon(k,l,m,:)=EPG_fitting_T2_reg(fit_out,data_norm,double(B1(k,l,m)),TEs,0,1,0).*(sum(data(:).^2)).^0.5;
                %Pass T2 estimate to output array
                T2(k,l,m)=fit_out(1);
                %Estimate T2 standard error
                [~,SE]=nlparci_SE(fit_out,residual,'jacobian',J);
                %Pass error estimate to array
                T2_SE(k,l,m)=SE(1);
            end
        end
    end
end

%%
%Generate R2
R2=1./T2;R2(isnan(R2))=0;R2(isinf(R2))=0;
%Correct nans
T2(isnan(T2))=0;
T2_SE(isnan(T2_SE))=0;
TSE_recon(isnan(TSE_recon))=0;
%%
%Write files
niftiwrite(T2,[out_path,'_T2map'],'Compressed',true)
niftiwrite(T2_SE,[out_path,'_T2map_SE'],'Compressed',true)
niftiwrite(TSE_recon,[out_path,'_T2_sig_recon'],'Compressed',true)
niftiwrite(double(TSE)-TSE_recon,[out_path,'_T2_sig_res'],'Compressed',true)
niftiwrite(R2,[out_path,'_R2map'],'Compressed',true)
niftiwrite(fval,[out_path,'_fval'],'Compressed',true)
niftiwrite(exitflag,[out_path,'_exitflag'],'Compressed',true)

