function FixativeModelling(out_path,L1_path,L2_path,L3_path,V1_path,V2_path,V3_path,mask_path,tau,iterations,model,flow,conc_frac,str_end,seed_path)
%%
% This code simulates fixative influx/outflux in tissue using a Kinetic
% Tensor (KT) or Kinetic Isotropic (KI) model

%Details of these models can be found in the manuscript titled:
%A method to remove the influence of fixative concentration on post-mortem T2 maps using a Kinetic Tensor model

%Written by Benjamin Tendler - contact benjamin.tendler@ndcn.ox.ac.uk

%%
% Outputs correspond to the resulting concentration maps

% Inputs are as follows:

% out_path      - folder where data will be output
% L1/2/3_path   - Paths to L1/2/3 (eigenvalue) estimates of diffusion tensor
% V1/2/3        - Paths to V1/2/3 (eigenvector) estimates of diffusion tensor
% mask_path     - Path to tissue mask
% tau           - Duration of each iteration (s)
% iterations    - Number of iterations
% model         - 1 is KI model, 2 is KT model
% flow          - 1 models outflow, 2 models inflow 
% conc_frac     - For outflow modelling, output is instead a map of the duration (in s) for each voxel to reach a concentration equal to conc_frac (between 0 and 1) - this is not constrained by the number of iterations. Setting this value to 0 ignores this option 
% str_end       - Append the output file names with this at the end
% seed_path     - When modelling inflow, seed mask defines is where the pool of fixative exists - if set to 0, will be the inverse of the tissue mask
%%
%Load in data
%Read in diffusion data
L1=niftiread(L1_path); 
L2=niftiread(L2_path); 
L3=niftiread(L3_path); 
V1=niftiread(V1_path); 
V2=niftiread(V2_path); 
V3=niftiread(V3_path); 
%Read in mask
mask=niftiread(mask_path); 
%Get mask header info
mask_info=niftiinfo(mask_path);
%Get seed path - if seed path equal to 0, fixative pool will be inverse of the tissue mask
if isequal(seed_path,0) || isequal(seed_path,'0')
    seed=0;
else
    seed=single(niftiread(seed_path));
end
%%
%Get x,y and z voxel dimensions
dim(1)=mask_info.PixelDimensions(1);
dim(2)=mask_info.PixelDimensions(2);
dim(3)=mask_info.PixelDimensions(3);
%Define concentration condition as off (by default)
conc_cond=0;
%%
%Define inflow or outflow
if flow==1
    concentration=mask;
    out_str_flow='_outflow';
elseif flow==2
    if seed == 0
        concentration=ones(size(mask))-mask;
    else
        concentration=seed;
    end
    out_str_flow='_inflow';
    %Define condition if estimating the time it takes for each voxel to reach a concentration
    if conc_frac~=0
        out_str_flow=['_inflow_time_to_concentration_days'];
        iterations=1;
        conc_cond=1;
    end
end
%%
%Define individual components of the diffusion tensor - model 1 is KI,
%Model 2 is KT
if model==1
    mean_diffusivity=(L1+L2+L3)/3;
    mean_diffusivity=mean(mean_diffusivity(mean_diffusivity~=0));
    D{1,1}=mean_diffusivity*mask;
    D{1,2}=mask*0;
    D{1,3}=mask*0;
    D{2,1}=mask*0;
    D{2,2}=mean_diffusivity*mask;
    D{2,3}=mask*0;
    D{3,1}=mask*0;
    D{3,2}=mask*0;
    D{3,3}=mean_diffusivity*mask;
    out_str='KI';
    clear mean_diffusivity L1 L2 L3
elseif model==2
    D{1,1}=(V1(:,:,:,1).*L1.*V1(:,:,:,1)+V2(:,:,:,1).*L2.*V2(:,:,:,1)+V3(:,:,:,1).*L3.*V3(:,:,:,1)).*mask;
    D{1,2}=(V1(:,:,:,1).*L1.*V1(:,:,:,2)+V2(:,:,:,1).*L2.*V2(:,:,:,2)+V3(:,:,:,1).*L3.*V3(:,:,:,2)).*mask;
    D{1,3}=(V1(:,:,:,1).*L1.*V1(:,:,:,3)+V2(:,:,:,1).*L2.*V2(:,:,:,3)+V3(:,:,:,1).*L3.*V3(:,:,:,3)).*mask;
    D{2,1}=(V1(:,:,:,2).*L1.*V1(:,:,:,1)+V2(:,:,:,2).*L2.*V2(:,:,:,1)+V3(:,:,:,2).*L3.*V3(:,:,:,1)).*mask;
    D{2,2}=(V1(:,:,:,2).*L1.*V1(:,:,:,2)+V2(:,:,:,2).*L2.*V2(:,:,:,2)+V3(:,:,:,2).*L3.*V3(:,:,:,2)).*mask;
    D{2,3}=(V1(:,:,:,2).*L1.*V1(:,:,:,3)+V2(:,:,:,2).*L2.*V2(:,:,:,3)+V3(:,:,:,2).*L3.*V3(:,:,:,3)).*mask;
    D{3,1}=(V1(:,:,:,3).*L1.*V1(:,:,:,1)+V2(:,:,:,3).*L2.*V2(:,:,:,1)+V3(:,:,:,3).*L3.*V3(:,:,:,1)).*mask;
    D{3,2}=(V1(:,:,:,3).*L1.*V1(:,:,:,2)+V2(:,:,:,3).*L2.*V2(:,:,:,2)+V3(:,:,:,3).*L3.*V3(:,:,:,2)).*mask;
    D{3,3}=(V1(:,:,:,3).*L1.*V1(:,:,:,3)+V2(:,:,:,3).*L2.*V2(:,:,:,3)+V3(:,:,:,3).*L3.*V3(:,:,:,3)).*mask;
    out_str='KT';
    clear V1 V2 V3 L1 L2 L3
end
%%
%Calculate D~1,2,3 
for k=1:3
    for l=1:3
        D_shift{k,l}=(circshift(D{k,l},1,l)-circshift(D{k,l},-1,l))/(2*dim(l));
    end
end
D_tilde{1}=D_shift{1,1}+D_shift{1,2}+D_shift{1,3};
D_tilde{2}=D_shift{2,1}+D_shift{2,2}+D_shift{2,3};
D_tilde{3}=D_shift{3,1}+D_shift{3,2}+D_shift{3,3};
clear D_shift
%%
%Perform KI/KT modelling
%Set up output for duration of time for voxels to reach concentration defined by conc_frac
concentration_duration=zeros(size(concentration));
iter=1;
no_vox=sum(mask(:));
while no_vox~=0
    for k=1:iterations
        %Define shifted components of the concentration by one voxel along x,y and z (here 2,2,2 corresponds to
        %the original image, 2,2,3 corresponds to one voxel shifted in positive z, 2,1,2 corresponds to one voxel shifted in negative y etc)
        for l=-1:1
            for m=-1:1
                for n=-1:1
                    concentration_shift{l+2,m+2,n+2}=circshift(circshift(circshift(concentration,l,1),m,2),n,3);
                end
            end
        end
        %Calculate individual components of concentration (second derivative)
        d2concentration{1,1}=(concentration_shift{1,2,2}+concentration_shift{3,2,2}-2*concentration)/(dim(1).^2);
        d2concentration{2,2}=(concentration_shift{2,1,2}+concentration_shift{2,3,2}-2*concentration)/(dim(2).^2);
        d2concentration{3,3}=(concentration_shift{2,2,1}+concentration_shift{2,2,3}-2*concentration)/(dim(3).^2);
        d2concentration{1,2}=((concentration_shift{1,1,2}+concentration_shift{3,3,2})-(concentration_shift{3,1,2}+concentration_shift{1,3,2}))/(4*dim(1)*dim(2));
        d2concentration{1,3}=((concentration_shift{1,2,1}+concentration_shift{3,2,3})-(concentration_shift{3,2,1}+concentration_shift{1,2,3}))/(4*dim(1)*dim(3));
        d2concentration{2,3}=((concentration_shift{2,1,1}+concentration_shift{2,3,3})-(concentration_shift{2,3,1}+concentration_shift{2,1,3}))/(4*dim(2)*dim(3));
        %Calculate first derivative components
        dconcentration{1}=(concentration_shift{3,2,2}-concentration_shift{1,2,2})/(2*dim(1));
        dconcentration{2}=(concentration_shift{2,3,2}-concentration_shift{2,1,2})/(2*dim(2));
        dconcentration{3}=(concentration_shift{2,2,3}-concentration_shift{2,2,1})/(2*dim(3));
        %Form concentration
        concentration=(D{1,1}.*d2concentration{1,1}+D{2,2}.*d2concentration{2,2}+D{3,3}.*d2concentration{3,3}+2*D{1,2}.*d2concentration{1,2}+2*D{1,3}.*d2concentration{1,3}+2*D{2,3}.*d2concentration{2,3}+D_tilde{1}.*dconcentration{1}+D_tilde{2}.*dconcentration{2}+D_tilde{3}.*dconcentration{3})*tau+concentration;
        %Set concentration to zero outside brain
        if flow==1
            concentration=concentration.*mask;
        elseif flow==2
            concentration(mask==0)=0;
            concentration(seed==1)=1;
        end
        %Remove errors
        concentration(concentration<0)=0;
        concentration(concentration>1)=1;
        %Output concentration distribution per day
        if mod(k,round(86400/tau))-1==0 && k~=1
            niftiwrite(concentration,[out_path,'/concentration_',out_str,out_str_flow,str_end,'_day_',num2str(round(k/round(86400/tau)))],'Compressed',true)
        end
        if conc_cond ==0

            ['iteration ', num2str(iter), ' of ', num2str(iterations)]
        else
            concentration_duration(concentration_duration ==0 & concentration >= conc_frac & mask==1)=tau*iter/(24*60*60);
            no_vox=sum(mask(:))-sum(mask(concentration_duration>0),'all');
            ['iteration ', num2str(iter),': ',num2str(no_vox),' voxels remaining with concentration less than ', num2str(conc_frac)]
            %Prevent getting stuck on a few erroneous voxels by breaking loop if no change in concentration for 1000 iterations
            no_vox_arr(iter)=no_vox;
            if iter>1000
                if no_vox == no_vox_arr(iter-1000)
                    concentration_duration(concentration_duration ==0 & mask==1)=inf;
                    ['unable to reach condition for ',num2str(no_vox),' voxels: set to time of infinity']
                    no_vox=0;
                end
            end
        end
        iter=iter+1;
    end
    if conc_frac==0 | flow==1
        no_vox=0;
    end
end
clear D d2concentration dconcentration concentration_shift D_tilde
%%
%Output maps
if conc_frac==0 | flow==1
    niftiwrite(concentration,[out_path,'/concentration_',out_str,out_str_flow,str_end],'Compressed',true)
else
    niftiwrite(concentration_duration,[out_path,'/concentration_',out_str,out_str_flow,str_end],'Compressed',true)
end
