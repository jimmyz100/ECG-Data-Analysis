function [finalData] = extractNoise(input,fs, coef )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
plot_flag = 1;
%===================Remove baseline first==================================
fDomain=fft(input); 
fDomain(1 : round(length(fDomain)*1.5/fs))=0; 
fDomain(end - round(length(fDomain)*1.5/fs) : end)=0; 
final_input = real(ifft(fDomain));
if plot_flag ==1
    figure, subplot(211); plot(input);
    subplot(212); plot(final_input,'r');
end
%=========================================================================%
maxWindow = coef * fs;
max_index = zeros(round(length(final_input)/maxWindow),1);
min_index = zeros(round(length(final_input)/maxWindow),1);
j = 1;
for i = 1 : maxWindow: length(final_input) - maxWindow
    max_index(j) = find(final_input(i: i + maxWindow) == max(final_input(i: i + maxWindow)));
    min_index(j) = find(final_input(i: i + maxWindow) == min(final_input(i: i + maxWindow)));
    maxPoint(j) = max_index(j) + i - 1;
    minPoint(j) = min_index(j) + i - 1;
    j = j + 1;
    
end
aveMax = mean(final_input(maxPoint));
aveMin = mean(final_input(minPoint));
j = 1;
flag1 = 2; flag2 = 2; flag3 = 0;
% find points which is beyond 2*Average Min
for i = 1: length(minPoint)
    if abs(final_input(minPoint(i))) >= 2* abs(aveMin)
        bottomPoint(j) = minPoint(i);
        j = j+1;
        flag1 = 0;
    end
end
if (j == 1)
    flag1 = 10;
    bottomPoint(j) = 0;
    bottomPoint(j) = [];
end
% find points which is beyond 2*Average Max
j = 1;
for i = 1: length(maxPoint)
    if abs(final_input(maxPoint(i))) >=  2 * aveMax
        topPoint(j) = maxPoint(i);
        j = j+1;
        flag2 = 1;
    end
end
 if (j == 1)
     flag2 = 10;
     topPoint(j) = 0;   
     topPoint(j) = [];
 end
j = 0;
%if the signal is good, we don't need to trim off everything
if flag1 ~= 10 || flag2 ~= 10
    for i = 1: length(bottomPoint) + length(topPoint)
        if (isempty(bottomPoint))
            combinedPoint(i) = topPoint(i);     
        elseif(isempty(topPoint))     
            combinedPoint(i) = bottomPoint(i);
        elseif (isempty(bottomPoint) && isempty(topPoint))
           break; 
        end
        if (i <= length(bottomPoint))
            combinedPoint(i) = bottomPoint(i);
        else
           combinedPoint(i) = topPoint(i - length(bottomPoint)); 
        end
    end
    combinedPoint = unique(combinedPoint);
    %=========================================================================%
    subData = cell(length(combinedPoint) + 1,1);
    subData{1} = final_input(1: combinedPoint(1)-fs);
    if combinedPoint(end) < length(final_input) - fs -1
       subData{end} = final_input(combinedPoint(end) + fs : length(final_input)); 
    else
        subData{end} = [];
    end
    for i = 1: length(combinedPoint)-1
        subData{i+1} = final_input(combinedPoint(i) + fs : combinedPoint(i+1) - fs);
    end
    for i = 1:  length(subData)
       if length(subData{i}) < round(4.2*fs)
         subData{i} = [];
         j = j + 1;
       end
    end
    finalData = cell(length(combinedPoint) - j+1,1);
    k = 0;
    for i = 1 : length(subData)
       if (~isempty(subData{i}))
           k = k+1;
           finalData{k} = subData{i};
       end
    end
elseif flag1 == 10 && flag2 ==10
    finalData = cell(1,1);
    finalData{1} = final_input; 
end
if plot_flag == 1
    figure;
    plot(final_input)
    hold on
    plot(maxPoint, final_input(maxPoint),'r');
    plot(minPoint, final_input(minPoint),'r');
    plot(1: length(final_input), aveMax*ones(length(final_input),1),'g');
    plot(1: length(final_input), aveMin*ones(length(final_input),1),'g');
    if (flag1 == 0)
       plot(bottomPoint, final_input(bottomPoint),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
    end
    if(flag2 == 1)
       plot(topPoint, final_input(topPoint),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
    end 
    if(flag3 == 1)
        plot(cutPoint(:,1), final_input(cutPoint(:,1)),'r*','MarkerEdgeColor','g','MarkerFaceColor','r','MarkerSize',6); 
        plot(cutPoint(:,2), final_input(cutPoint(:,2)),'r*','MarkerEdgeColor','g','MarkerFaceColor','r','MarkerSize',6); 
    end
end

 

end


