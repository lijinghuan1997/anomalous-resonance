%% load data
mms.db_init('local_file_db','C:\mms_db_EMIC')
timedur = '2020-01-14T19:23:30.00Z/2020-01-14T19:23:50.00Z';
Tint = irf.tint(timedur);
ic=1;
c_eval('ePDist1=mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',1);
c_eval('ePDist2=mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',2);
ePDist2=ePDist2.resample(ePDist1.time);
c_eval('ePDist3=mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',3);
ePDist3=ePDist3.resample(ePDist1.time);
c_eval('ePDist4=mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',4);
ePDist4=ePDist4.resample(ePDist1.time);

c_eval('Bgse = mms.get_data(''B_gse_fgm_brst_l2'',Tint,?);',ic);

fmin=0.25;
fmax=0.65;
PA=0:5:175;
energy=ePDist1.depend{1}(1,:);
%% run
% define the FAC
Bgse1=Bgse;
magB = Bgse.abs;
dfB = 1/median(diff(Bgse.time.epochUnix));
Bgse3= Bgse.filt(0,0.05,dfB,3);
Bdir=irf.nanmean(Bgse3.data);
Bdir=Bdir/norm(Bdir);
Ydir=cross(Bdir,[1,0,0]);
Ydir=Ydir/norm(Ydir);
Xdir=cross(Ydir,Bdir);
A_=[Xdir;Ydir;Bdir]; % B0 as the z direction; y direction is the B0 cross [1,0,0]
b_=[1,0,0;0,1,0;0,0,1];
mat_coord =A_\b_; % transfor matrix

% filt Bw 
dfB = 1/median(diff(Bgse1.time.epochUnix));
Bgse2= Bgse1.filt(fmin,fmax,dfB,3);
Bgse2.data=Bgse2.data*mat_coord;

eID=16; % selected energy channel
pitch=90; % selected pitch angle


record=zeros(12,length(ePDist1.time));
for tID=1:length(ePDist1.time)-1 % the timepoint circulation
    result=plotgyro(ePDist1,ePDist2,ePDist3,ePDist4,Bgse,tID,eID,pitch);
    for gyro=1:12
         p1=result{gyro}>0;
         record(gyro,tID)=mean([result{gyro}(p1)]); % using average here
         if  isempty(find(p1>0))
             record(gyro,tID)=0;
         end          
    end
end
% "record" is the spectra with the same time resolution as the measurment

% average two timepoints; not necessary
record1=record(:,1:2:end);
record2=record(:,2:2:end);
width=size(record2,2);
recordave=(record1(:,1:width)+record2)/2;
% filter2 average 3*3
A=filter2(1/9*ones(3),recordave*1e25,'same');
A(1,:)=smooth(recordave(1,:)*1e25,3);
A(end,:)=smooth(recordave(end,:)*1e25,3);
datatotal{eID}=A;
A(A==0)=NaN;

data=A;
% spect.t=ePDist1.time(2:2:end).epochUnix;
% spect is the struct to use the irf_spectrogram function
c_eval('spect.t=ePDist?.time(1:2:end).epochUnix;',ic)
spect.f=15:30:345; % gyrophase angle for each grid
spect.p=data'; % the PSD data
Bphi1=TSeries(Bgse.time,transphase2(Bgse2,1)); % calculate the phase of the Bw; '1' represent the direction of Bw
Bphi1=Bphi1.tlim(irf.tint('2020-01-14T19:23:32.00Z/2020-01-14T19:23:46.50Z')); % limit the interval of Bphi
irf_spectrogram(spect,'lin'); % plot spectra
colormap('jet')
hold on
irf_plot(Bphi1,'k') % overplot the Bphi




function gyro=plotgyro(ePDist1,ePDist2,ePDist3,ePDist4,Bgse,tID,eID,pitch)

% define FAC
dfB = 1/median(diff(Bgse.time.epochUnix));
Bgse3= Bgse.filt(0,0.05,dfB,3);
Bdir=irf.nanmean(Bgse3.data);
Bdir=Bdir/norm(Bdir);
Ydir=cross(Bdir,[1,0,0]);
Ydir=Ydir/norm(Ydir);
Xdir=cross(Ydir,Bdir);
A_=[Xdir;Ydir;Bdir];
b_=[1,0,0;0,1,0;0,0,1];
mat_coord =A_\b_;
% the phase angle is an estimation, 
% you can also read the real phase angle from ePDist1.depend{2};
phi_edges = linspace(11.25/2,360-11.25/2,size(ePDist1.data,3));  % azimuthal angle bin edges, default
theta_edges = linspace(11.25/2,180-11.25/2,size(ePDist1.data,4)); % polar angle bin edges, default
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = -sind(THETA).*cosd(PHI); % '-' because the data shows which direction the particles were coming from
Y = -sind(THETA).*sind(PHI);
Z = -cosd(THETA);


gyro=cell(12,1);
for phi=15:30:345
for thetaid=1:length(theta_edges)
    for phiid=1:length(phi_edges)
        newvec=[X(thetaid,phiid),Y(thetaid,phiid),Z(thetaid,phiid)]*mat_coord; % particle velocity FAC
        targetvec=[sind(pitch)*cosd(phi),sind(pitch)*sind(phi),cosd(pitch)]; % target velocity
        angle=acosd(dot(newvec,targetvec)/norm(newvec)/norm(targetvec));
        if angle<=11.25
            gyro{floor((phi+14.99)/30)+1}= ...
            [gyro{floor((phi+14.99)/30)+1},ePDist1.data(tID,eID,phiid,thetaid),ePDist2.data(tID,eID,phiid,thetaid) ...
            ePDist3.data(tID,eID,phiid,thetaid),ePDist4.data(tID,eID,phiid,thetaid)];
        end
    end
end
end
end


function result=transphase2(Bgse,choice)
    [Bphi_,Bnorm]=cart2pol(Bgse.data(:,1),Bgse.data(:,2));
    if choice==1
        Bphi_=Bphi_*180/pi;
    else
        Bphi_=Bphi_*180/pi+180;
    end
    Pa=Bphi_;
    Pa(Pa<0)=Pa(Pa<0)+360;
    Pa(Pa<0)=Pa(Pa<0)+360;
    index=Pa<3;
    Pa(index)=NaN;
    index=Pa>357;
    Pa(index)=NaN;
    result=Pa;
end
