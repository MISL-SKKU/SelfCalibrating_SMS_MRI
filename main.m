%% Initialize
clear all;close all;clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load SMS data in 3D k-space
% load('smsdata.mat')

% *********************
%   Data description:
% *********************
%    raw: A 4D (size n1 x n2 x nc x vs) array of measured 3D k-space data
%           to be reconstructed (with Extended controlled aliasing).
%    nc : number of receiver coils
%    n1, n2, n3: number of encodings along ky, kx, and kz, respectively
%    ** n3 represents the number of slice experiencing controlled aliasing.
%    ** n3 = mb (mb factor) x ppa (in-plane acceleration)


[n1,n2,nc,n3,~]=size(raw);


% Set basic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [**]This part should be modified to the characteristics of the data.
mb = 3;
ppa = 2;
caipi = 3; 
% vs = mb*ppa; % for one-step reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Without extended controlled aliasing, we apply it by upsampling in k-space.
% if (n3 == mb)&&(ppa>1)
%     tmp = zeros(n1,n2,nc,n3*ppa);
%     tmp(:,:,:,1:ppa:end) = raw;
%     raw=tmp;
% end


% extract sampling mask
mask_raw = (squeeze(raw(:,100,1,:))~=0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Decomposition & Extended Self-Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The measured SMS data is decomposed into two sets of data for 
%  1) SMS imaging           --- mbaraw    (n1 x n2 x nc)
%  2) SMS calibration       --- sms_calib (s1 x s2 x nc x n3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extended Sampling Pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find self-calibrating region
cal_idx = find(sum(mask_raw,2)>1);

% Extract shift pattern to match CAIPI condition between SMS imaging and
% SMS calibration
tmp = mask_raw(1:n3,1:ppa:end);
mask = repmat(tmp,ceil(n1/n3),1);
% mask(n1,:)=0; mask = mask(1:n1,:); % if size of mask != mask_raw

% Create mask for only sms data without self-calibrating data
mask_vs = zeros(n1,n3);
mask_vs(:,1:ppa:end)=mask;
for p=2:ppa
    mask_vs(p:end,p:ppa:end)=mask(1:end-(p-1),:);
end
% Display sampling pattern corresponding to CAIPI shift pattern
% figure,imshow(mask_vs(1:12,:),[],'initialMagnification',300)

% Create mask for sms data with self-calibrating data
mask_vs_sc = mask_vs;
mask_vs_sc(cal_idx,1:ppa:end)=1;
% figure,imshow(mask_vs_sc,[],'initialMagnification',300)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract only SMS data without self-calibrating data
raw_vs = zeros(size(raw));
for c=1:mb*ppa
    raw_vs(mask_vs(:,c)~=0,:,:,c)=raw(mask_vs(:,c)~=0,:,:,c);
end
% Generate SMS dataset
mbraw = sum(raw_vs,4);
% Display SMS image
% figure,imshow(sos(ifft2c(mbraw),3),[])


% Extract SMS calibration data in the self-calibrating region
sccal_3d = raw(cal_idx,:,:,:);
sms_calib = ifft(sccal_3d,[],4);

% Without self-calibration, this part can be modified for external calib.
% sms_calib = external_calib;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extended Self-Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply linear gradient to k-space of each slice for FOV shifts 
%    to match CAIPI condition between SMS imaging and SMS calibration
sms_calib_shft=shift_by_mask(sms_calib,mask_vs(cal_idx,:),1);
% Display calibration image
% figure,imshow(sos(ifft2c(zpad(sms_calib_shft(:,:,:,2),n1,n2,nc)),3),[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One-Step Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set reconstrction parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [**]This part should be adjusted to the characteristics of the data.
params.mb=mb;params.caipi=caipi; params.ppa=ppa;
params.kSize=[7,7]; params.nul=350; params.r=90; params.iter=40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct complementary null space using calibration data.
[Vs,Nsc]=constructionSPSL(sms_calib_shft,params,0);

% SMS reconstruction
res = gnrHSL_fast(mbraw,Vs,Nsc,raw,mask_vs,params);
% Apply linear gradient to k-space of each slice for FOV shifts
res_hsl_1s = shift_by_mask(res,mask_vs,0);

% Display result image
for si=1:mb
figure(si),imshow(sos(ifft2c(zpad(res_hsl_1s(:,:,:,si),n1,n2,nc)),3),[])
end
