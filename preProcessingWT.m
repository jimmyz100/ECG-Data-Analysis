function [finalOutput] = preProcessingWT(noOfData, input, u)
%===========================Read me========================================
%Firstly, the signal is devided into equal segment using function
%partitionData
%Sencondly, these segments are applied WT to filter, then they are
%recombined again under newFish name
%===========================Display signal time domain=====================
coef = 4096;
fishData = cell(noOfData,1); %noOfData is the number of segment partitioned
%by trim off process for each signal
cout = 0;
for i = 1 : noOfData
    fishData{i}(:,1) = 1: coef*(round(length(input{i}))/coef - 1);
    fishData{i}(:,2) = input{i}(1: length(fishData{i}(:,1)));
%     if length(fishData{i}(:,1)) < coef
%         cout = cout +1;
%     end
end
dataForWT = partitionData(fishData, coef, noOfData);% from each segment,it 
%is devided in 1024, 4096 and etc for WT process
filtered_fishData = cell(noOfData,1);
finalOutput = cell(noOfData,1);
for i = 1: length(dataForWT)
   
    [row col] = size(dataForWT{i});
    swa = cell(row,1);
    swd = cell(row,1);
    nswd = cell(row,1);
    nswa = cell(row,1);
    for j = 1: row
        [swa{j}, swd{j}]= swt(dataForWT{i}(j,:), 10,'coif5'); 
        ThreshML(j,:) = wthrmngr('sw1ddenoLVL','sqrtbal_sn',swd{j},'mln');
         %bal_sn  sqrtbal_sn
        % sqtwolog: universal thereshold
        % thr = wthrmngr('sw1ddenoLVL','sqtwolog',swtdec,scale)
        % thr = wthrmngr('sw1ddenoLVL','rigrsure',swtdec,scale)
        % thr = wthrmngr('sw1ddenoLVL','heursure',swtdec,scale)
        % thr = wthrmngr('sw1ddenoLVL','minimaxi',swtdec,scale)
        % mln is scale: one, sln
        selectThrFunc = '1';
        for k = 1: 9
           nswd{j}(k,:) =  thrFunction(swd{j}(k,:), 0.4*ThreshML(j,k), selectThrFunc);
        end
        nswd{j}(10,:) = zeros(1, length(dataForWT{i}));
        nswa{j}(10,:) = zeros(1, length(dataForWT{i}));
        filtered_fishData{i}(j,:) = iswt(nswa{j}, nswd{j},'coif5');
    end
end
finalOutput = combineData(filtered_fishData);
% dlmwrite('testECG_data.txt',filter_ECG_test,'delimiter','\t');
%%
% save processed data to a new folder and rename them.
for i = 1 : length(dataForWT)
    name = sprintf('procssedSA_%d_%d.txt',u,i);
    newdataPath = fullfile('processedData\',name);
    dlmwrite(newdataPath,finalOutput{i},'delimiter','\t');
    figure, plot(finalOutput{i});
    title(['fish',num2str(u),'-', num2str(i)]);
    hold on
    plot(input{i});
end
end