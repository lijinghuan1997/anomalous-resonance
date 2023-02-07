%% load data
mms.db_init('local_file_db','C:\mms_db_EMIC')
ic = 1; 
timedur = '2020-01-14T19:23:20.00Z/2020-01-14T19:24:00.00Z';
Tint=irf.tint(timedur);

% here the electromagnetic fields are averaged over four spacecraft
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',Tint,?);',1:4);
c_eval('Bxyz?=Bxyz?.resample(Bxyz1);',2:4);
c_eval('Bgse?=mms.get_data(''B_gse_brst_l2'',Tint,?);',1:4);
c_eval('Bgse?=Bgse?.resample(Bgse1);',2:4);
c_eval('Exyz?=mms.get_data(''E_gse_edp_brst_l2'',Tint,?);',1:4);
c_eval('Exyz?=Exyz?.resample(Exyz1);',2:4);
c_eval('Egse?=mms.get_data(''E_gse_edp_brst_l2'',Tint,?);',1:4);
c_eval('Egse?=Egse?.resample(Egse1);',2:4);

Bxyz=Bxyz1;
Bxyz.data=(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4;
Exyz=Exyz1;
Exyz.data=(Exyz1.data+Exyz2.data+Exyz3.data+Exyz4.data)./4;

Bgse=Bgse1;
Bgse.data=(Bgse1.data+Bgse2.data+Bgse3.data+Bgse4.data)/4;
Egse=Egse1;
Egse.data=(Egse1.data+Egse2.data+Egse3.data+Egse4.data)/4;
% load the plasma moments from FPI-MOMS
c_eval('Vi = mms.get_data(''Vi_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('ni = mms.get_data(''Ni_fpi_brst_l2'',Tint,?);',ic);
c_eval('ne = mms.get_data(''Ne_fpi_brst_l2'',Tint,?);',ic);
% the magnetic field strength
magB = Bxyz.abs;
Bxyzmag = TSeries(Bxyz.time,[Bxyz.data magB.data]);
%%
% construct the FAC;  the difination is the same as the gyrophase code
dfB = 1/median(diff(Bgse.time.epochUnix));
Bgse3= Bgse.filt(0,0.05,dfB,3);
Bdir=irf.nanmean(Bgse3.data);
Bdir=Bdir/norm(Bdir);
Ydir=cross(Bdir,[0,0,1]);
Ydir=cross(Bdir,[1,0,0]);
Ydir=Ydir/norm(Ydir);
Xdir=cross(Ydir,Bdir);
A_=[Xdir;Ydir;Bdir];
b_=[1,0,0;0,1,0;0,0,1];
mat_coord =A_\b_;

%% waavelet
fmin = 0.09; fmax = 1;
Bgse1= Bgse.filt(fmin,fmax,dfB,3);
Bgse1.data=Bgse1.data*mat_coord;
nf = 100;  % number of frequency
Bwavelet = irf_wavelet(Bgse1,'nf',nf,'f',[fmin fmax]);
nc = 10;   % timepoint average
idx = nc/2:nc:length(Bwavelet.t)-nc/2;
Bwavelettimes = Bwavelet.t(idx);
Bwaveletx = zeros(length(idx),nf);
Bwavelety = zeros(length(idx),nf);
Bwaveletz = zeros(length(idx),nf);
for ii = 1:length(idx)
        Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
        Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
        Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specB=struct('t',Bwavelettimes);
specB.f=Bwavelet.f;
specB.p=Bwaveletx+Bwavelety; %+Bwaveletz; only care about the perpendicular components
specB.f_label='';
specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};

%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
% Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wci = e*B_SI/Mp;
% Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
% Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
% Fpp = Wpp/2/pi;
Fche = Fcp/4;
Fche2 = Fcp/2;
% Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
% Fce = irf.ts_scalar(magB.time,Fce);
% Flh = irf.ts_scalar(magB.time,Flh);
% Fpp = irf.ts_scalar(magB.time,Fpp);
Fcp = irf.ts_scalar(magB.time,Fcp);
Fche = irf.ts_scalar(magB.time,Fche);
Fche2 = irf.ts_scalar(magB.time,Fche2);

%%

% wave filter
fmin=0.25; fmax=0.65;
Bgse2= Bgse.filt(fmin,fmax,dfB,3);
Bgse2.data=Bgse2.data*mat_coord;
Bgse1=Bgse2;
dfE = 1/median(diff(Egse.time.epochUnix));
Egse2= Egse.filt(fmin,fmax,dfE,3);
Egse2.data=Egse2.data*mat_coord;
% Ew has the same time resolution as Bw.
Egse2=Egse2.resample(Bgse2);

EB=Bgse2;
EB.data=cross(Egse2.data,Bgse2.data); % poynting flux

Enorm=Egse2.abs;
Enorm.data=sqrt(Egse2.data(:,1).^2+Egse2.data(:,2).^2); % Ewave norm
Bnorm=Bgse2.abs;
Bnorm.data=sqrt(Bgse2.data(:,1).^2+Bgse2.data(:,2).^2); % Bwave norm
dotliang=dot(Egse2.data(:,1:2),Bgse2.data(:,1:2),2)./Bnorm.data; 
Enorm.data=sqrt(Enorm.data.^2-dotliang.^2); % only use the Ew component perpendicular to Bw
Enorm=TSeries(Enorm.time,smooth(Enorm.data,400));
Bnorm=TSeries(Bnorm.time,smooth(Bnorm.data,400));
speed=TSeries(Enorm.time,Enorm.data./Bnorm.data.*1e3); % estimate vphase  km/s


%%

h=irf_plot(6,'newfigure');
xSize=750; ySize=600;
set(gcf,'Position',[10 10 xSize ySize]);
xwidth = 0.86;
ywidth = 0.1;
set(h(1),'position',[0.10 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.10 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.10 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.10 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.10 0.97-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.10 0.97-6*ywidth xwidth ywidth]);


h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyzmag);
ylabel(h(1),{'B (nT)'},'Interpreter','tex','fontsize',12);
irf_zoom(h(1),'y',[-50 60]);
set(h(1),'ytick',[-30 0 30]);
set(h(1),'xlim',[0,60]);
set(h(1),'xtick',[10 15 20 25 30])
set(h(1),'xticklabel',[])

h(2)=irf_panel('Bspec');
irf_spectrogram(h(2),specB,'log');
hold(h(2),'on');
irf_plot(h(2),Fcp,'color','r','LineWidth',1.5)
irf_plot(h(2),Fche,'color','b','LineWidth',1.5)
hold(h(2),'off');
irf_legend(h(2),'f_{cp}',[0.3 0.60],'color','r','fontsize',12)
irf_legend(h(2),'f_{che}',[0.35 0.60],'color','b','fontsize',12)
caxis(h(2),[-2 2.3]);
line(h(2),[10,30],[0.25,0.25])
line(h(2),[10,30],[0.65,0.65])
set(h(2),'yscale','log');
set(h(2),'ytick',[1e-1 0.25 0.65]);
ylabel(h(2),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
set(h(2),'xtick',[10 15 20 25 30])
colormap(h(2),'jet');
set(h(2),'xticklabel',[])

h(3)=irf_panel('Bwave');
hold(h(3),'on');
irf_plot(h(3),TSeries(Bgse2.time,Bgse2.data(:,1)),'color','k','LineWidth',1.5)
irf_plot(h(3),TSeries(Bgse2.time,Bgse2.data(:,2)),'color','b','LineWidth',1.5)
irf_plot(h(3),TSeries(Bgse2.time,Bgse2.data(:,3)),'color','r','LineWidth',1.5)
hold(h(3),'off');
ylim(h(3),[-4,4])
set(h(3),'ytick',[-2,0,2])
set(h(3),'xtick',[10 15 20 25 30])
set(h(3),'xticklabel',[])

h(4)=irf_panel('Ewave');
hold(h(4),'on');
irf_plot(h(4),TSeries(Egse2.time,Egse2.data(:,1)),'color','k','LineWidth',1.5)
irf_plot(h(4),TSeries(Egse2.time,Egse2.data(:,2)),'color','b','LineWidth',1.5)
irf_plot(h(4),TSeries(Egse2.time,Egse2.data(:,3)),'color','r','LineWidth',1.5)
ylim(h(4),[-1.2,1.2])
set(h(4),'ytick',[-1,0,1])

hold(h(4),'off');
set(h(4),'xtick',[10 15 20 25 30])
set(h(4),'xticklabel',[])


h(5)=irf_panel('Poynting');
A=TSeries(EB.time,smooth(EB.data(:,1),400)*1e-9*1e-3/4/pi*1e7*1e6);
A1=TSeries(EB.time,smooth(EB.data(:,2),400)*1e-9*1e-3/4/pi*1e7*1e6);
A2=TSeries(EB.time,smooth(EB.data(:,3),400)*1e-9*1e-3/4/pi*1e7*1e6);
hold(h(5),'on')
irf_plot(h(5),A);
irf_plot(h(5),A1);
irf_plot(h(5),A2);
hold(h(5),'off')
set(h(5),'xtick',[10 15 20 25 30])
set(h(5),'ylim',[-1.3,0.5])
set(h(5),'ytick',[-1,0])
set(h(5),'xticklabel',[])

h(6)=irf_panel('Vphase');
irf_plot(h(6),speed)
% ylim(h(6),[70,350])
ylim(h(6),[45,250])
set(h(6),'ytick',[100,200])
set(h(6),'xtick',[10,15,20,25,30])
set(h(6),'xticklabel',[])


set(h(1:6),'xlim',[10,30])
irf_plot_axis_align(h(1:6));