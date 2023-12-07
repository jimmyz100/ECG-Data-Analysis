clc
clear all
close all

% Input signal file
MainProject = 'Surgery Test';
Folder = 'Cryosurgery Test Oct';
Day = '10.25.2021';
Group = '1_3';
input = load([ MainProject '\' Folder '\Day ' Day '\Grp' Group '.txt']);

% Normalizing data
time = (1:size(input, 1)) ./ 1000;
data = (input(:, 7)) .* 4.5 ./ (2^23 - 1);

% Plotting of signals
plot(time,data);
%xlim([56 59]);
%ylim([-1 2.5]);
set(gca,'TickLength',[0 0])
fs = 250;
 data = data(:); % vectorize
   ax = plot(data);axis tight;title('Raw Signal');
  xlim([56000 59000]);ylim([-1 2.5]);
%data=data(10000:80000); %%%%% do segmentation manually

% %%%%%Segmentation automatically
% [M,N]=size(data);
% k1 = 10;   %%% number of segments
% k2=N/k1;
% k3 = 1;
% for i=1:k1   
% segment1(i,:) = data(k3:i*k2);
% k3 = i*k2+1;
% end
% 
% for m=1:k1
% maxnum1(m)=max(abs(segment1(m,:)));
% end
%   
%  a1 = min(maxnum1);
%  b1 = max(maxnum1); 
%    for f=1:k1
%      for p=1:k2
%      if abs(segment1(f,p))>a1 
%          segment1(f,p)=a1;
%      end
%      end
%    end
%    
%    data2 = segment1(1,:);
%    for f=2:k1
%        data2 = cat(2, data2, segment1(f,:));
%    end
%    
%    figure
%    plot(data2);
%    title('cutting')
%    
%    
%    %%%%%%End of segmentation
% %   [qrs_i_raw,sign,en_thres] = qrs_detect(data(4000:6000),0.350,0.6,250,[],[],0);
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin_bounds(data,fs,1);
%%%R-R calculation
R = qrs_i_raw(2:2:end);
P = qrs_i_raw(1:2:end);
% 
% % D = [R;P];
% % PR=squareform(pdist(D'));
% % PR=R-P;
% [f,m]=size(R);
for i=1:m-1
   RR(i)= R(i+1)-R(i);
   PP(i)= P(i+1)-P(i);
end

% Determination of Means
PRavg = mean(PR)
RRavg = mean(RR)
PPavg = mean(PP)
PPStd = std(PP)
RRStd = std(RR)

%   figure
%   plot(data); 
%   hold on
%   for i=1:length(qrs_i_raw)
%     x=qrs_i_raw(i);
% 
%     plot(x,1,'r--o','LineWidth',2)
%     hold on
%  end