% chi_comparisons
%
% compare methodologies for calculating chi and epsilon from the FLEAT data
%
% Jen, April 2024
%
%


%% load in data to use

% where data is located
dir0='/Users/jen/Library/CloudStorage/GoogleDrive-jmackinnon@ucsd.edu/My Drive/'; %jen's locatiom

load([dir0 'chi_testing.mat'])
% this has two stuctures
%   - M with mooring quantititues at a given depth
%   - S with the RBR solo raw temperature data (at a fast time) for this
%   depth. 



%% Calculate temperature frequency spectra over specified segement length

dt=10*60; %segment length, let's start with 10 minutes

dt0=nanmean(diff(S.dn))*24*3600; % sampling frequency of the Solo

lseg=round(dt/dt0);% number of points in each segment

nseg=floor(length(S.dn)/lseg); % number of segments in time series
warning off
clear Ptot time_seg
for iseg=1:nseg
    if rem(iseg,10000)==0; disp(['iseg = ' int2str(iseg) ' of ' int2str(nseg)]); end % to be able to see we are making progress!
    ii=[1:lseg]+(iseg-1)*lseg;
    temp=S.t(ii);
    [Pxx,om0]=fast_psd(despike(detrend(temp)),300,1/dt0);
    Ptot(iseg,:)=Pxx(:).';
    time_seg(iseg)=nanmean(S.dn(ii));
end
Om=ones(nseg,1)*om0(:).';

% Following Anna's paper, have a running mean of spectra over an hour
Ptot_sm=NaN*Ptot;
a=1; b=ones(1,round(3600/dt))/round(3600/dt);
for iom=1:length(om0);
    Ptot_sm(:,iom)=nanfilt(b,a,Ptot(:,iom));
end

% spectra slope criteria for rejecting data

%% now calculate Epsilon

% Here do the multiplication of the spectra by frequency **before**
% averaging, this is a type of pre-whitening and gets a good best estimate
% of the spectral level. 
iiom=find(om0>.02&om0<.05); % these are Anna's frequency limits
Specsum=nanmean(Ptot(:,iiom).*Om(:,iiom).^(5/3),2); Specsum=Specsum(:).';  

% interpolate mooring quantites onto the segment times. 
ui=abs(interp1(M.dn,M.u,time_seg)); dtdzi=interp1(M.dn,M.dtdz,time_seg); n2i=interp1(M.dn,M.n2,time_seg);

% some constants. 
gamma=0.2; Ct=0.4;

% calculate epsilon. See associated document for where this version of the formula
% comes from 
eps=(((2*pi)^(2/3)*Specsum.*(ui.^(-2/3)).*n2i)./(2*Ct*gamma*dtdzi.^2)).^(3/2);
epstest=

% let's have some standards for rejecting data
ib=find(eps<0); eps(ib)=NaN; eps=real(eps); % simple things

% more standards for rejecting data, from Anna's paper
ib=find((ui<.05)|(n2i<1e-6)|(dtdzi<1e-3)); eps(ib)=NaN;
%%
% spectra slope criteria for rejecting data
% first fit spectra in specified frequency range

% anna does the fit between 1/200 and 1/4 excluding 1/20-1/7
iiomb=find(om0>1/20&om0<1/7);
Ptot_sm(:,iiomb)=NaN;
iiom2=find(om0>1/200&om0<1/4);

for iseg=1:nseg
    ig=find(~isnan(Ptot_sm(iseg,iiom2)));
    p = polyfit(log10(om0(iiom2(ig))),log10(Ptot_sm(iseg,iiom2(ig))),1);
    Ptot_slope(iseg)=p(1);
end
% anna rejects slopes for temperature gradient spectra outside of [0 .66]
% here we have temperature spectra, not temperature gradient, so will
% reject slopes outside of [0 .66]-2
ib=find(Ptot_slope<-2 | Ptot_slope > -1.34);
eps0=eps;
eps(ib)=NaN;  % note this is rejecting something like 58% of the data.  Seems a bit harsh...?
% interpolate to fill in gaps
%ig=find(~isnan(eps)); 
%eps=interp1(time_seg(ig),eps(ig),time_seg);
% name this as our zeroth order epsilon estimate
% structure will be C for comparison
C.dn=time_seg;
C.eps_take1=eps; 



%% a couple possibilities for smoothing. Let's smooth over about an hour. 
a=1; b=ones(1,5)/5;
epsf=NaN*eps; ig=find(isfinite(log10(eps)));
epsf(ig)=10.^nanfilt(b,a,log10(eps(ig)));
epsf=10.^nanfilt(b,a,despike(log10(eps),1));

C.eps_take2=epsf; 



%% now let's try doing the same thing in wavenumber space


dt0=nanmean(diff(S.dn))*24*3600; % sampling frequency
dt=10*60; %segment length, let's start with 2 minutes
lseg=round(dt/dt0);% number of points in each segment
kk=1e-2:1e-2:2;
Ptot_k=NaN*ones(nseg,length(kk));
nseg=floor(length(S.dn)/lseg); % number of segments in time series
Ptot_k=NaN*ones(nseg,length(kk));
for iseg=1:nseg
    if rem(iseg,10000)==0; disp(['iseg = ' int2str(iseg) ' of ' int2str(nseg)]); end
    ii=[1:lseg]+(iseg-1)*lseg;
    temp=S.t(ii);
    time_seg(iseg)=nanmean(S.dn(ii));
    if ~isnan(ui(iseg))
        [Pxx,k0]=fast_psd(despike(detrend(temp)),300,1/(dt0*ui(iseg)));
        omtest=k0*ui(iseg); ig=find(omtest<.05); 
        Ptot_k(iseg,:)=interp1(k0(ig),Pxx(ig),kk);
    end
end

%% epsilon from wavenumber spectra
iik=find(kk>.1&kk<.6);
KK=ones(nseg,1)*kk(:).';

Specsum=nanmean(Ptot_k(:,iik).*KK(:,iik).^(5/3),2);
gamma=0.2; Ct=0.4;
Specsum=Specsum(:).';
%eps=(((2*pi)^(2/3)*Specsum.*(ui.^(4/3)).*n2i)./(2*Ct*gamma*dtdzi.^2)).^(3/2);
%eps=(((2*pi)^(2/3)*Specsum.*(ui.^(-2/3)).*n2i)./(2*Ct*gamma*dtdzi.^2)).^(3/2);
eps=2*pi/((2*Ct*gamma)^(3/2))*(n2i.^(3/2))./(dtdzi.^3).*Specsum.^(3/2);

ib=find(eps<0); eps(ib)=NaN; eps=real(eps);

% standards for rejecting data
ib=find((ui<.05)|(n2i<1e-6)|(dtdzi<1e-3)); eps(ib)=NaN;

C.eps_take3=eps;


%% make some plots
% time series 
figure(1); clf; fig_setup(1); set(gcf,'paperposition',[.5 .5 10 6]); wysiwyg
semilogy(M.dn,M.eps_anna1,C.dn,C.eps_take1,C.dn,C.eps_take2,'linewidth',1); shg 
hold on
eps_annaf=interp1(M.dn,M.eps_anna1,C.dn); %interpolated onto the same time grid
 ii=find(C.dn>datenum(2016,9,1)&C.dn<datenum(2017,3,1)); % after the weird biofouling period

plot([M.dn(1) M.dn(end)+10],nanmean(eps_annaf(ii))*[1 1],'k','linewidth',2)
plot([M.dn(1) M.dn(end)+10],nanmean(C.eps_take1(ii))*[1 1],'r','linewidth',2)
plot([M.dn(1) M.dn(end)+10],nanmean(C.eps_take2(ii))*[1 1],'b','linewidth',2)
datetick('x',12); shg
legend(gca,'Anna take 0','Jen take 1','Jen smoothed')
xlabel('date'); ylabel('Epsilon / W kg^{-1}')
set(gca,'ylim',[1e-12 1e-4])

%%
printdir=[dir0 'figs/'];
print(gcf,[printdir 'timeseries_comppare1'],'-dpng')

%% scatter plots
figure(1); clf; fig_setup(1); set(gcf,'paperposition',[.5 .5 6 5]); wysiwyg
loglog(C.eps_take1(ii),eps_annaf(ii),'.',[1e-15 1e-5],[1e-15 1e-5]);
axis equal
axis([1e-12 1e-4 1e-12 1e-4])
xlabel('Jen take 1, pre-whitened'); ylabel('Anna take 0')

%%
printdir=[dir0 'figs/'];
print(gcf,[printdir 'scatter1'],'-dpng')

%% frequency versus wavenumber calculations
figure(1); clf; fig_setup(1); set(gcf,'paperposition',[.5 .5 6 5]); wysiwyg
loglog(C.eps_take1(ii),C.eps_take3(ii),'.',[1e-15 1e-5],[1e-15 1e-5]);
axis equal
axis([1e-12 1e-4 1e-12 1e-4])
xlabel('Jen: frequency space'); ylabel('Jen: wavenumber space')

%%
printdir=[dir0 'figs/'];
print(gcf,[printdir 'scatter2'],'-dpng')




