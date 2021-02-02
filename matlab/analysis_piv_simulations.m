clear all
close all
clc

%load PIV files
load('D:\simulations\example\PIV_data_corrected\large_stack_325.tif.mat')


%parameters
px2mic=3.25;
piv_window_size=16;%in pixels
overlap=0.5;%from the PIV settings
dt_PIV=5;%in minutes
dt_av=12;%dt for the slide averaging in time points: 12=1h
dt_plot=60;%dt between two measures or plots in time points: 24=2h 

dx_max_um=600;%dx max to compute the x correlations in um
dx_um=100;%width in um of the averaging window in x to get the x velocity profile along y

% analyze the PIV results
[u_mean_over_t, v_mean_over_t, u_middle, h, x_c_um, corr, corr_500 , mean_width_lanes]= analysis_simu(u,v,dt_av,dt_plot,px2mic,overlap,piv_window_size,dx_max_um,dx_um);% see function at the end


% plot the results
tend=size(u);%number of time frames
time=1:dt_plot:tend-dt_av;%frames at which we perform measurements 
time_hour=(time-1)*5/60;%corresponding time in hours
tf=size(time,2);%number of measurements time points 
px2mic_piv=overlap*piv_window_size*px2mic;
y_um=(1:size(u{1},1))*px2mic_piv;

for t=1:tf
  
    u_mean(t)=floor(nanmean(nanmean(abs(u_mean_over_t{t}))));
    v_mean(t)=floor(nanmean(nanmean(abs(v_mean_over_t{t}))));
    v_rms(t)=floor(nanmean(nanmean(sqrt(u_mean_over_t{t}.^2+v_mean_over_t{t}.^2))));
    
    u_mean_str=num2str(floor(nanmean(nanmean(abs(u_mean_over_t{t})))));
    v_mean_str=num2str(floor(nanmean(nanmean(abs(v_mean_over_t{t})))));
    v_rms_str=num2str(floor(nanmean(nanmean(sqrt(u_mean_over_t{t}.^2+v_mean_over_t{t}.^2)))));
    
    % Snapshot of x and y component of the velocity at different time
    % points
    clims=[-60 60];
    figure(1)
    sgtitle({"snapshot of x and y component of the velocity"},'FontSize',30)
    % x component u
    subplot(2,tf,t)
    imagesc(u_mean_over_t{t},clims)
    axis equal
    title(time_hour(t)+" h",'FontSize',20)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
     % y component v
    subplot(2,tf,t+tf)
    imagesc( v_mean_over_t{t},clims)
    axis equal
    title({"u mean="+u_mean_str+"um/h","v mean="+v_mean_str+"um/h","v rms="+v_rms_str+"um/h"},'FontSize',10)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    % histogram of x and y component of the velocity 
    figure(2)
    sgtitle({"histogram of x and y component of the velocity (um/h)"},'FontSize',30)
    % x component u
    subplot(2,tf,t)
    hist( u_mean_over_t{t}(:),100)
    xlim([-100 100])
    title(time_hour(t)+" h",'FontSize',20)
    set(gca,'ytick',[])
    % y component v
    subplot(2,tf,t+tf)
    title(time_hour(t)+' h','FontSize',20)
    hist(v_mean_over_t{t}(:),100)
    xlim([-100 100])
    title({"u mean="+u_mean_str+"um/h","v mean="+v_mean_str+"um/h","v rms="+v_rms_str+"um/h"},'FontSize',10)
    set(gca,'ytick',[])
    
    % vx correlation decrease along the x direction for different time
    % points
    figure(3)
    subplot(1,tf,t)
    plot(x_c_um,corr(t,:))
    hold on
    ylim([0 1])
    title(time_hour(t)+" h",'FontSize',20)
    
    % u profil in the middle of the field of view along y for different
    % time points
    figure(4)
    subplot(tf,1,t)
    plot(y_um,u_middle(:,t))
    hold on
    plot([0 y_um(end)], [0 0],'--r')
    ylim([-100 100])
    title(time_hour(t)+" h",'FontSize',20)

end


figure(3)
sgtitle({"Normalised v_x correlation C along the x direction"},'FontSize',30)
subplot(1,tf,1)
ylabel('normalized v_x correlation C','FontSize',15)
xlabel('dx in \mum','FontSize',15)


figure(4)
sgtitle({"v_x profile along the y direction"},'FontSize',30)
subplot(tf,1,1)
ylabel('v_x in \mum/h','FontSize',15)
subplot(tf,1,tf)
xlabel('y in \mum','FontSize',15)
ylim([-100 100])
title(time_hour(t)+" h",'FontSize',20)

% time evolution of the anisotropy ration h, the vx correlation at 500um x
% distance C(500 um) and the mean lanes width

figure(5)

% anisotropy ratio h 
subplot(1,3,1)
plot(time_hour,h)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('anisotropy ratio h','FontSize',15)
xlabel('time in hours','FontSize',15)

% C(500 um) 
subplot(1,3,2)
plot(time_hour,corr_500)
ylim([0 1])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('C(500 \mum)','FontSize',15)
xlabel('time in hours','FontSize',15)

%mean lane width
subplot(1,3,3)
plot(time_hour,mean_width_lanes)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('mean lanes width in \mum','FontSize',15)
xlabel('time in hours','FontSize',15)


% adjust the figures size
scrsz = get(0,'ScreenSize');
h=scrsz(4);
l=scrsz(3);
set(figure(1),'position',scrsz);
set(figure(2),'position',scrsz);
set(figure(3),'position',[1,h/4,l,h/4]);
set(figure(4),'position',[l/4,1,l/4,h]);
set(figure(5),'position',[l/4,h/4,l*2/3,h/3]);

%
%function analysing a PIV files

% parameters:
% u
% v
% dt_av: dt defining the time slide averaging in time frames, usually 12
% (=1h)
% dt_plot: dt between two time points to analyze, in time frames
% px2mic
% overlap (from the PIV settings)
% piv_window_size (from the PIV settings)
% dx_max_um: corresponding to the maximal x distance for computing x
% correlation in um
% dx_um: width of the averaging window to compute the v_x profile along y

% outputs:
% u_mean_over_t: u_mean_over_t{t} is the v_x component at time t after time
% averaging
% v_mean_over_t: same for v_y
% u_middle: u_middle(:,t) is the v_x profile averaged over dx_um
% h: h(t) is the anisotropy ratio at time t
% dx_c_um: is the dx at which the correlations are measured
% corr: corr(t,:) is the vx correlation along x
% corr_500: corr_500(t) is the value of the vx correlation for dx=500um
% mean_width_lanes: mean_width_lanes(t)

function [u_mean_over_t, v_mean_over_t, u_middle, h, x_c_um, corr, corr_500,mean_width_lanes]= analysis_simu(u,v,dt_av,dt_plot,px2mic,overlap,piv_window_size,dx_max_um,dx_um)

tend=size(u);
time=1:dt_plot:tend-dt_av;
tf=size(time,2);
px2mic_piv=overlap*piv_window_size*px2mic;
dx_max_px=floor(dx_max_um/px2mic_piv);

u_mean=zeros(1,tf);
v_mean=zeros(1,tf);
h=zeros(1,tf);
corr=zeros(tf,dx_max_px+1);
corr_500=zeros(1,tf);
mean_width_lanes=zeros(1,tf);
u_middle=zeros(size(u{1},1) ,tf);

for t=1:tf
    
    tplot=time(t);
    
    % time averaging of v_x and v_y
    [u_mean_over_t{t},v_mean_over_t{t}]=time_average(u,v,tplot,dt_av);% see function at the end
    
    % compute h
    u_mean(t)=floor(nanmean(nanmean(abs(u_mean_over_t{t}))));
    v_mean(t)=floor(nanmean(nanmean(abs(v_mean_over_t{t}))));
    h(t)=1-v_mean(t)/u_mean(t);
    
    % compute x correlations of v_x
    [x_c_um,corr(t,:)]=correlation_x(u_mean_over_t{t},px2mic_piv,dx_max_px);% see function at the end
    index_500um=floor(500/px2mic_piv)+1;
    corr_500(t)=corr(t,index_500um);
    
    
    % compute the v_x profile along y averaged over dx_um in the middle
    Lx=size(u_mean_over_t{t},2);
    Ly=size(u_mean_over_t{t},1);
    dx=floor((dx_um/px2mic_piv)/2)+1;
    u_middle(:,t)=mean(u_mean_over_t{t}(:,Lx/2-dx+Lx/2+dx),2);
    
    % compute mean lane width
    clearvars y_change_lane_um y_change_lane width_lanes
    y_change_lane=[];
    width_lanes=[];
    sign=u_middle(1,t)/abs(u_middle(1,t));%sign of first value
    for j=2:Ly
        if sign*u_middle(j,t)<0
            y_change_lane=[y_change_lane,j];
            sign=u_middle(j,t)/abs(u_middle(j,t));
        end
    end
    y_change_lane_um=y_change_lane*px2mic_piv;
    width_lanes=diff(y_change_lane_um);
    mean_width_lanes(t)=mean(width_lanes);
    
end
end

% time slide averaging over dt (in frames number)
function [u_mean_over_t, v_mean_over_t]=time_average(u,v,t,dt)
u_mean_over_t=zeros(size(u{1}));
v_mean_over_t=zeros(size(v{1}));

tmin=t;
tmax=t+dt;

for ti=tmin:1:tmax
    u_mean_over_t=u_mean_over_t+u{ti};
    v_mean_over_t=v_mean_over_t+v{ti};
end

u_mean_over_t=u_mean_over_t/(tmax-tmin+1);
v_mean_over_t=v_mean_over_t/(tmax-tmin+1);

end

% compute x correlations of an array 
% dj_max is the maximal dx at which we measure correlations (in pixels)
function [x_c_um,c]=correlation_x(array,px2mic_piv,dj_max)

Li=size(array,1);
Lj=size(array,2);
c=zeros(1,dj_max+1);

for dj=0:dj_max
    cdj=0;
    
    for j=1:Lj-dj
        prod=array(:,j).*array(:,j+dj);
        cdj=cdj+sum(prod./(abs(prod)));
    end
    
    c(dj+1)=cdj/(Li*(Lj-dj));
end

x_c_um=(0:dj_max)*px2mic_piv;

end