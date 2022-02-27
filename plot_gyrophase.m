%% load data from txt file
path='C:\Users\Àî¾©å¾\Documents\Visual Studio 2013\Projects\ceshiEMIC\ceshiEMIC\FSD.txt';
[t,vx,vy,vz,enem,Bx,By,Z]=textread(path,'%f %f %f %f %f %f %f %f','delimiter',' ');
E1=1/2*1.67*1e-27*(vx.^2+vy.^2+vz.^2)/1.6*1e19*1e10;
% change velocity into the spacecraft rest frame
vxpla=-20e3;
vypla=10e3;
vzpla=35e3;
vz=vz+vzpla/1e5;
vx=vx+vxpla/1e5;
vy=vy+vypla/1e5;
E=1/2*1.67*1e-27*(vx.^2+vy.^2+vz.^2)/1.6*1e19*1e10;
% calculate the phase angle of the particle and wave magnetic field
[Vphi,Vperp]=cart2pol(vx,vy);
[Bphi,Bperp]=cart2pol(Bx,By);
Vphi=Vphi*180/pi;
p=find(Vphi<0);
Vphi(p)=Vphi(p)+360;
Bphi=Bphi*180/pi;
p=find(Bphi<0);
Bphi(p)=Bphi(p)+360;
% obtain the phase difference
theta=Vphi-Bphi;
% make sure the theta is between (0,2pi)
poshigh=find(theta>=180);
theta(poshigh)=theta(poshigh)-360;
neglow=find(theta<=-180);
theta(neglow)=theta(neglow)+360;
%% record index of different particles
maxindextotal=[];
maxindextotal(1)=1;
for i=2:length(t)
    if t(i)>t(i-1)
        maxindextotal=[maxindextotal,i];
    else if t(i)==max(t)
            maxindextotal=[maxindextotal,i];
        else if (t(i)-t(i-1))<-0.1
                maxindextotal=[maxindextotal,i];
            end
        end
    end
end

%% 
% based on energy channels/ pitch angles in the fortran code
% adjust the ener and pt parameters
ener=1;
pt=1;
timepoint=65;  % spacecraft position number
gyro=12;       % particle gyrophase number
% the index adjustment (if ener or pt is not equal to 1)
xxx=pt*timepoint*gyro*(ener-1);
maxindex=maxindextotal(gyro*timepoint*(pt-1)+1+xxx:gyro*timepoint*pt+xxx); % the end channel +1

datasave574eV=cell(gyro,timepoint);
for i=1:timepoint*gyro 
    if i==timepoint*gyro && ener==1
        index=maxindex(i):length(t);
    else 
        index=maxindex(i):maxindex(i+1)-1;
    end
    % index is for one particle
    % to reduce the size of FSD.txt, we only output part of trajectory
    % 
    y=ceil(i/gyro);
    x=mod(i,gyro);
    if x==0
        x=gyro;
    end
    % save data for each particle
   datasave200eV{x,y}={-t(index),E(index),theta(index),Vperp(index)*1e5,vx(index)*1e5,vy(index)*1e5,vz(index)*1e5,Bx(index),By(index),E1(index)};
end
% save datasave200_4nT_xyshift_PA90_65point_pt1.mat datasave1000eV
save datasave200eV_4nT_xyshift_PA90_65point.mat datasave200eV
%% gyroplot
mp=1.67*1e-27;
qi=1.6*1e-19;
Tperp=130;
Tpara=50;
xwidth=0.75;
ywidth=0.1;
vxpla=-20e3;
vypla=10e3;
vzpla=35e3;
timepoint=65;
gyro=12;
const=8.5*1e6*(mp/2/pi/130/qi)^(2/2)*(mp/2/pi/50/qi)^(1/2)*1e13;
% the initial distribution function
f=@(vx,vy,vz)exp(-mp/2*(vx-vxpla).^2/Tperp/qi-mp/2*(vy-vypla).^2/Tperp/qi-mp/2*(vz-vzpla).^2/Tpara/qi);
psd=zeros(gyro,timepoint);
for ploti=1
%             c_eval('load(''datasave?_PA90_gyro15_345_amp=4nT_f0.5212_5phase.mat'');',ploti);
        dataplot=datasave200eV;  
for i=1:gyro
    for j=1:timepoint
        % obtain the particle initial velocity
        vx_=cell2mat(dataplot{i,j}(5)); 
        vx_=vx_(end);
        vy_=cell2mat(dataplot{i,j}(6));
        vy_=vy_(end);
        vz_=cell2mat(dataplot{i,j}(7));
        vz_=vz_(end);
        % calculate the corresponding PSD
       psd(i,j)=f(vx_,vy_,vz_);
    end
end
psd_=[psd;zeros(1,timepoint)];
psd_=[psd_,zeros(gyro+1,1)];
h=pcolor(psd_(:,1:63)*const);
set(h,'LineStyle','none');
% caxis([3000,5850])
set(gca,'ytick',[4 7 10]);
set(gca,'yticklabel',[90 180 270])
hold on
% calculate the phase angle of wave magnetic field
time=0:0.01:18.3;
Bphi_=360-2*pi*0.43*time/pi*180;
p=find(Bphi_<0);
while ~isempty(p)
    p=find(Bphi_<0);
    Bphi_(p)=Bphi_(p)+360;
    p=find(Bphi_<0);
end
q=find(Bphi_>357);
Bphi_(q)=NaN;
n=find(Bphi_<3);
Bphi_(n)=NaN;
% calculate the phase angle of wave electric field
Ephi_=360-2*pi*0.43*time/pi*180+90;
p=find(Ephi_<0);
while ~isempty(p)
    p=find(Ephi_<0);
    Ephi_(p)=Ephi_(p)+360;
    p=find(Ephi_<0);
end
q=find(Ephi_>357);
Ephi_(q)=NaN;
n=find(Ephi_<3);
Ephi_(n)=NaN;


yyaxis right
% plot the phase angle of wave magnetic field
plot(1.5:0.01*61/18.3:62.5,Bphi_,'linewidth',2)
ylim([0,360])

set(gca,'ytick',[])
colormap('jet')

end
