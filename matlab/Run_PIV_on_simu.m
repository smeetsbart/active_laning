clear all
close all
clc

% file information
pathname='D:\simulations\test\';
filename='stack-small.tif'; %the file is a tif sequence
nb_frames_per_hour=12;
px2mic=3.25;

%PIV settings
window_PIV_size=16;
overlap=0.5;
method='single';

%filters settings
threshold_snr=1.1;
threshold_global=5;
threshold_local=5;
window_local=3;

mkdir([pathname,'\','PIV_data']) %create a new folder to stock PIV results
mkdir([pathname,'\','PIV_data_corrected'])%create a new folder to stock PIV corrected results

write_PIV(pathname,filename,nb_frames_per_hour,px2mic,window_PIV_size,overlap,method)
correct_PIV(pathname,filename,threshold_snr,threshold_global,threshold_local,window_local)


function write_PIV(pathname,filename,nb_frames_per_hour,px2mic,window_PIV_size,overlap,method)

delta_t=1/nb_frames_per_hour;%in hours

% Write velocity field files
disp('writing...');
info = imfinfo([pathname,'\',filename]);
n_frames = numel(info);
t0=1;

%------- get piv data
for t=t0:n_frames-1
    image1=imread([pathname,'\',filename],t);
    image2=imread([pathname,'\',filename],t+1);
    [xt,yt,ut,vt,snrt,pkht] = matpiv(image1,image2,window_PIV_size,delta_t,overlap,method);
    
    x{t,1} = xt; y{t,1} = yt; u{t,1} = ut*px2mic; v{t,1} = vt*px2mic; snr{t,1}=snrt;pkh{t,1}=pkht;
end

% save piv data
save([pathname,'\','PIV_data\',filename,'.mat'],'x','y','u','v','snr','pkh');
end 


% filter PIV data

function correct_PIV(pathname,filename,threshold_snr,threshold_global,threshold_local,window_local)
disp('writing...');

load([pathname,'\PIV_data\',filename,'.mat'],'u','v','x','y','snr')
n_frames = size(u,1);
t0=1;
for t=t0:n_frames
    %------ filter results
    [su,sv]=snrfilt(x{t},y{t},u{t},v{t},snr{t},threshold_snr);
    [gu,gv]=globfilt(x{t},y{t},su,sv,threshold_global);
    [mu,mv]=localfilt(x{t},y{t},gu,gv,threshold_local,'median',window_local);
    [fu,fv]=naninterp(mu,mv,'linear');
    
    x{t,1} = x{t}; y{t,1} = y{t}; u{t,1} = fu; v{t,1} = fv;
end

%save piv data corrected 
save([pathname,'\','PIV_data_corrected\',filename,'.mat'],'x','y','u','v');

end



