%%**************************************************************************************************************************************
%Function: Numerical computation of a FSAR ambiguity function (AF) in a 2D (radial-axial) grid.
%Author:Omkar Pradhan
%Date: 6/5/2016
%Local filename: SARThryAA_3D_v3c
%Description:
%AA - represents TxRx antenna ID (i.e. monostatic in this case) - 
% - Compute ambiguity function from pure analytical form of pulse modulation and exponential term (see pg.9 in res.notes#2)
% - The 2D grid is centered at zTrgt, rhoTrgt with appropriate distance windows in the 2 axes. this is done to investigate AF at very
distant points without blowing up the array size.
% - error can be added to the depth increments
%%**************************************************************************************************************************************

toSave = 1; %% Flag to save(1) or not(0) the resulting figures

MtNm1 = 'ChiSq';%Name: Ambig.Func. (AF) TxRx antenna (AA or BB or AB etc.)
RhoZNm = 'RhoZ';
SvInd = 1;%use this to index .mat files in case repeating same calculations
%_mmddyy

% relPth = 'C:\Users\Capmunk';%change this depending upon computer used
relPth = 'D:';%change this depending upon computer used
% relPth = 'C:\Users\aamir207';

gglDrvPth = 'Google Drive';
FigSvPth = 'VALKYRIE\Imgs_Fgrs\Figures\Sims';
MatSvPth = 'VALKYRIE\MIT_062015_DataPrcss\Processed_Mat_01082017_1';
%% Environment parameters
%Fif MHz signal sampled at Fs MHz

c0 = 3e8;%in m/s
Ep0 = 8.8542e-12;
EpIcerel = 3.17;
Eprel = EpIcerel;%last relative permittivity represents the scatterer.

%% target parameters
%Unit reflector position
xTrgt = (0);  zTrgt = -25;   yTrgt = 0;% in meters
rhoTrgt = sqrt(xTrgt^2+yTrgt^2);
Err_d = 5*1e-2; %Depth measurement error
Ne = 300;
% zStcks = 4;%  stack volumes along -ve z axis (X and Y limits remain the same)

%% Radar TxRx parameters

Fc = 870e6;%tx(carrier) frequency in hertz ( not the downconverted freq received by ADC)
Flo = 820e6;%local osc. frequency
Fif = Fc-Flo;
Fs = 500e6;%sampling frequency
Ts = 1/Fs;%sampling time period
Tprf = 3e-6;% Time period of pulse repetition frequency
Ton = 150e-9;% On time of pulse

sAB = 20 * 1e-2; % in meters (diameter of cryobot)- distance  bet. antennas A & B (180 deg apart)

%% derived parameters from TxRx parameters

w=2*pi*Fif;
Ns=round(Tprf/Ts);% # of samples for entire transmit-receive period
PhPerSmp = 2*pi/(Fs/Fif);%phase in radians / sample (Fs/Fif is # of samples for 2pi rad)

% t=linspace(0,2*Ton,2*Ton*Fs);%time vector for sine wave construction
t=(-Ton:Ts:Ton);%time vector for triangular wave construction

tNs = length(t);%# of samples in constructed signal for use later

lmbdaEff = c0/sqrt(real(Eprel))/Fc;%Wavelength in media
k = 2*pi/lmbdaEff;% wave number
% FarFldDist = 2*lmbdaEff^2/4;

%% SAR parameters
LstDpth = 600;
dpths = 1:LstDpth;

StrtDpth = 2;% in meters
StpDpth = 11;% in meters

tmp = linspace(StrtDpth,StpDpth,604);
difDpth = diff(tmp(1:2));
rdrDpths = -linspace(StrtDpth,StrtDpth+LstDpth*difDpth,LstDpth)';%radar depth in meters
rdrDpths(Ne:end) = rdrDpths(Ne:end) + Err_d;

%% AF reconstruction parameters- X,Y,Z volume grid
rhoCell = 5;
dim=2;
rhoWndw = rhoCell*3;%window around target (in meters) to be imaged.
zWndw = 10;
nPnts = 200*3;%# of pixels in each axis
nZ=3;
rhoStrt = 0; %yStrt = 0;
rhoEnd = 0; %yEnd = 0;
rhoMid=0;%yMid=0;
zMid = zTrgt;

% maxPnts = 100;
% pxlSz = 1e-2;%Prefered pixel size
rho = single(linspace(rhoStrt,rhoStrt+rhoWndw,nPnts)); %Range bins
% y = single(linspace(yMid-wndw,yMid+wndw,nPnts)); %Range bins
pxlSiz = abs(diff(rho(1:2)));
z = single(linspace(zMid+nZ*zWndw,zMid-nZ*zWndw,nZ*nPnts)); %Range bins
% z = single(linspace(zMid-nZ*wndw,zMid+nZ*wndw,nZ*nPnts)); %Range bins
% difZ = diff(rho(1:2));
% z = -single(abs(rdrDpths(end))+5:difZ:(abs(rdrDpths(end))+5+nPnts*nZ*difZ));
% z(end) = [];
phi = single(linspace(0,90,length(rho))*pi/180);

% rdrDpths = -rdrDpths;
[Rho,Zk] = meshgrid(rho,z);
Rho = reshape(single(Rho),nPnts^dim*nZ,1);
% Yj = reshape(single(Yj),nPnts^dim*nZ,1);
Zk = reshape(single(Zk),nPnts^dim*nZ,1);
% Phi = reshape(single(phi),nPnts^dim*nZ,1);
%% 2. Analytic weighting windows 
AngPnts = (0:-1:-180);
hmngWndw = hamming(numel(AngPnts));

%% Sinusoidal signal - reflected sig. and refernce sig.
vNs = 0:Ns;%vector of smaple points
rA = single(sqrt(rhoTrgt.^2+(zTrgt-rdrDpths).^2));%distance of target from antenna A
triWv = [tripuls(t,2*Ton),zeros(1,Ns-tNs)];%unit magnitude triangular function
%AA
delNs = single(2*rA/(c0/sqrt(real(Eprel)))/Ts);
intrpNs = interp1(vNs',vNs',delNs,'nearest');%interpolate Lag values to find the integral values closest to delNs
% [~,iNs] = ismember(intrpNs,vNs');
% sigDel=exp(1i*(2*pi*Fif*repmat(t,length(pdAA),1)+repmat(pdAA,1,length(t))));%delayed reference signal
% SigDel = [sigDel,zeros(length(pdAA),Ns-size(sigDel,2))];%delayed ref. sig. padded to same size as DechoAA
PulsAutoCorr = zeros(LstDpth,Ns);
for dd = 1:LstDpth
% shift the traingular wave to center at delay equivalent to target
% rountrip distance.
PulsAutoCorr(dd,:) = circshift(triWv,intrpNs(dd)-round(tNs/2),2);
end


%% Coherent aggregation.
%Reference signal has to be correct phase delayed in space and time. Then
%coherently add matched filtered output from range bins for each depth. On
%a time/sample scale matched filter output is the same for any pair of
%echo and ref. sig. irrespective of starting phase and is shifted only
%according to time delay (rountrip wave travel time to target)

%pre-allocate for speed
%D contains echos stacked in a vector with each row corresp. to an R(x,y,z) and changes for every dd
% Dtmp = sparse(zeros(nPnts^(dim)*zStcks,Ns));
%D contains echos stacked in a vector with each row corresp. to an R(x,y,z) and changes for every dd
% D = single(ones(nPnts^(dim)*zStcks,tNs)+1i*ones(nPnts^(dim)*zStcks,tNs));
%reference signal for matched filter impulse response
% sigRef = single(ones(nPnts^(dim)*zStcks,tNs)+1i*ones(nPnts^(dim)*zStcks,tNs));
%matched filter ouput is stored in this. (recursively added at each depth)
mfWt = (zeros(nPnts^(dim)*nZ,1));
mf = mfWt;
pdAA = single(zeros(nPnts^(dim)*nZ,1));
delNsAA = pdAA;
intrpNs = uint16(pdAA);
tstrt = tic;
tic
parfor dd = 1:LstDpth
%     for aa = 1:nStcks
    %find phase delays for reference signal
    riA = sqrt(Rho.^2+(Zk-rdrDpths(dpths(dd))).^2);
    [phi,~] = cart2pol(Rho,Zk-rdrDpths(dpths(dd)));
    wt = interp1(AngPnts,hmngWndw,phi*180/pi,'spline');% weighted sum with hamming window

    
    %     pdAA = single(mod(2*k*riA,2*pi)); % phase diff. in rdians
    pdAA = 2*k*riA; % phase diff. in rdians
    
    %find sample delay for time of travel
    delNsAA = single(2*riA/(c0/sqrt(real(Eprel)))/Ts);
    %interpolate  to pick out the waveform beginning at correct
    %lag value. 1 is added since indexing begins from 1
    intrpNs = interp1(vNs',vNs',delNsAA,'nearest')+1;
%     [~,iNs] = ismember(intrpNs,vNs');
    
    mf = mf+PulsAutoCorr(dd,intrpNs)'.*exp(1i*mod(pdAA-2*k*rA(dd),2*pi));

%     mfWt = mfWt+PulsAutoCorr(dd,intrpNs)'.*exp(1i*mod(pdAA-2*k*rA(dd),2*pi)).*wt;
    
    dd
end
toc

%% replicate on a constant z plane to show azimuthal ambiguity
% MFWt = reshape(mfWt,nPnts*nZ,nPnts)/LstDpth;
MF = reshape(mf,nPnts*nZ,nPnts)/LstDpth;

[~,iZ] = min(abs(z-(zTrgt)));
[~,iRho] = min(abs(rho-rhoTrgt));

% MFaz = repmat(MFWt(iZ,:)',1,length(phi));%each row is for a specific phi value
%interpolat to x-y
[PHI,RHO] = meshgrid(phi,rho);
[X,Y]= pol2cart(PHI,RHO);


%% Plot setup


[Rho,Z] = meshgrid(rho,z);
% Zslc = single(zTrgt);
% Xslc = single(rhoTrgt);
% Yslc = single(yTrgt);
RhoZ = [rho(1),z(1);rho(end),z(end);length(rho),length(z)]; %for storing



ChiFig1=figure('Color',[1,1,1],'name','Ambiguity Function');

figure(ChiFig1);
% plot(rho(iRho),z(iZ),'Marker','+','LineWidth',20,'MarkerSize',30,...
%     'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
% l=legend(sprintf('%d dB',round(20*log(abs(MF(iZ,iRho))))),'Location','NorthEast');
% l.FontName = 'Arial';
% l.FontSize = 20;
% l.Box = 'off';
% hold on

% mesh(Rho,Z-zTrgt,20*log10(abs(MF)./max(abs(MF(:)))));colormap jet;shading interp
mesh(Rho,Z-zTrgt,(abs(MF)./max(abs(MF(:)))).^2);colormap jet;shading interp
xlabel(sprintf('\\rho (m)'),'fontname','Times New Roman','fontsize',14);
ylabel('z (m)','fontname','Times New Roman','fontsize',14);
% zlabel('|\\chi|^2','fontname','Times New Roman','fontsize',14);
set(gca,'XLim',[rho(1),rho(end)],'XTick',(rho(1):1:rho(end)));
set(gca,'YLim',[z(end)-zTrgt,z(1)-zTrgt],'YTick',(z(end)-zTrgt:5:z(1)-zTrgt));

set(gca,'ZLim',[0,1],'ZTick',(0:0.2:1));
cb = colorbar;
set(cb,'YTick',(0:0.2:1));
set(cb,'YLim',[0,1]);

% set(cb,'YLim',[-60,0]);
% set(cb,'YTick',(-60:3:0));



set(gca,'FontSize',16);
grid on
view(2)

%% Save figure and .mat file
if(toSave)
    mtNm = [MtNm1,sprintf('_%d_%d_%d_%d',rhoTrgt,abs(zTrgt),LstDpth,SvInd)];
    rhozNm = [RhoZNm,sprintf('_%d_%d_%d_%d',rhoTrgt,abs(zTrgt),LstDpth,SvInd)];
    savefig(ChiFig1,fullfile(relPth,gglDrvPth,MatSvPth,mtNm));
%     save(fullfile(relPth,gglDrvPth,MatSvPth,[mtNm,'.mat']),'MF');
    save(fullfile(relPth,gglDrvPth,MatSvPth,[rhozNm,'.mat']),'RhoZ');
end
