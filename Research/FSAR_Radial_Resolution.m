%% *****************************************************
%Function: Analyzing resolution w.r.t axial distance.
%Author:Omkar Pradhan
%Date:1/30/2017
%Inherited from: SARAmbigFunc_Approx_v2b (Local filename)
%Local filenam: SARRes_v1
%% *****************************************************

% relPth = 'C:\Users\Capmunk';%change this depending upon computer used
relPth = 'D:';%change this depending upon computer used
% relPth = 'C:\Users\aamir207';

gglDrvPth = 'Google Drive';
FigSvPth = 'VALKYRIE\Imgs_Fgrs\Figures\Sims';
MatSvPth = 'VALKYRIE\MIT_062015_DataPrcss\Rdr_Sim';

%% ***********************Radar specs*********
c0 = 3*1e8;
f = 1e9;
w = 2*pi*f;
v = c0/sqrt(3.17);

rho_o = 0; %distance to target in meters
z_o = (25:50:525); % z (vertical) distance to target from radar position d_o
% z_o = (25); % z (vertical) distance to target from radar position d_o
Delrho = (-60:0.5e-2:60);% radial difference between target point and imaged point

StrtDpth = 2; Stp = (11:5:26);
    tmp = linspace(2,11,604);
    DpthIncr = diff(tmp(1:2));

resInd = zeros(length(z_o),length(Stp));
res = zeros(length(z_o),length(Stp));
for ss = 1:length(Stp)
    StpDpth = Stp(ss);
    
    % DpthIncr=1.5*1e-2; % depth increment
    d = StrtDpth:DpthIncr:StpDpth+1; % radar positions
    N = length(d);
    M= 25; % every sub aperture has M depths;
    P = floor(N/M);
    a1= N*DpthIncr;   b1= 1.0; c1= 1; %for StpDpth = 20 m
    a2= 8.842;   b2=0.8813; c2=1.003;%StpDpth = 10 m
    
    ChiEst1 = zeros(length(Delrho),length(z_o));
    ChiEst2 = ChiEst1;
    ChiExct = ChiEst1;
    for ii= 1:length(z_o)
        X1 = c1-a1/z_o(ii)^b1;%use curve fit toolbox
        X2 = c2-a2/z_o(ii)^b2;%use curve fit toolbox
        
        [phn,sincArg,ph1,ph2] = AmbFuncArgs(w,v,z_o(ii),rho_o,Delrho,DpthIncr,d,d,M);
        ChiExct(:,ii) = sum(exp(1i*(2*w/v)*phn),2);%exact ambiguity function
        
        ChiEst1(:,ii) = N*sinc((w/v)*sincArg*N./X1/pi);
        [~,resInd(ii,ss)]=min(abs(((abs(ChiEst1(:,ii))/max(abs(ChiEst1(:,ii)))).^2-0.5)));
        ChiEst2(:,ii) = N*sinc(ph1*P./X2/pi).*sinc(ph2*M/pi);
        ii;
    end
    res(:,ss) = abs(Delrho(1,resInd(:,ss))).'*2;

    ss
end
% indz=1;figure;plot(Delrho-rho_o,(abs(ChiEst1(:,indz))/max(abs(ChiEst1(:,indz)))).^2);
% legend('Est1-single sinc');title(sprintf('rho_o=%d, z_o=%d, N=%d',rho_o,z_o(indz),N));
% 
% indz=length(z_o);figure;plot(Delrho-rho_o,([abs(ChiExct(:,indz)),abs(ChiEst1(:,indz)),abs(ChiEst2(:,indz))]));
% legend('Exact','Est1-single sinc','Est2-2 sincs');title(sprintf('rho_o=%d, z_o=%d, N=%d',rho_o,z_o(indz),N));

figure;plot(z_o,res,'--*');

