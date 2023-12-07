clear all;  % Clear workspace
clc;        % Clear the command window.
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.

no = 396;
start = 395;
Trimtxtfast(start,no)
    A_HR = zeros(1,no);
    A_QTc = zeros(1,no);
    A_Ratio = zeros(1,no);
    A_RR = zeros(1,no);
    A_RTc = zeros(1,no);
    A_QT = zeros(1,no);
    A_RT = zeros(1,no);
    A_ST = zeros(1,no);

for noSignal = start:no
% noSignal = 13;
 Average_HR = [];
 QTc = [];
 Ratio = [];
 Average_RR = [];
 RTc = [];
 Average_QT = [];
 Average_RT = [];
 Average_ST = [];

for u = noSignal : noSignal
    name = sprintf('%d.txt', u);
    dataPath = fullfile('C:\Users\jimmy\OneDrive\School Work\PhD\PhD Research\Methamphetamine Study\Large Exp\Data for Processing',name);
   
    input = load(dataPath);
    finalData = extractNoise( input(:,2) ,1000,5 );
    noOfData = length(finalData);
    finalOutput = preProcessingWT(noOfData,finalData,u);
    finalData = finalOutput';
    noOfData = length(finalData);
    %===========================Loading all data sets==========================
  
    fishPackage = cell(noOfData,1);

    %========================Declare constant==================================
    Fs = 1000;
    RR_interval = 60;% 60 for fs = 1000; 8 for fs = 100
    slopeLevel = 5;% 5 for fs = 1000; 3 for fs =100
    windowLevel = 100;% 100 for fs = 1000; 8 for fs = 100
    QS_window = 6;% 6 for fs = 1000; 2 for fs = 100
    threshold = zeros(1,1);
    Q_wave = cell(noOfData,1);
    S_wave = cell(noOfData,1);
    P_wave = cell(noOfData,1);
    T_wave = cell(noOfData,1);
    peak = cell(noOfData,1);
    x = cell(noOfData,1);
    R_peak = cell(noOfData,1);
    R_peak_final = cell(noOfData,1);
    %=======================Debug==============================================
    no = 1;
    flag1 = 1;
    flag2 = 0;
    flag3 = 0;
    %==========================================================================
    for i = 1 : noOfData

        fishPackage{i}(:,2) = finalData{i};
        fishPackage{i}(:,1) = 1: length(finalData{i});
        %===========================Spectrum Power Density ========================
        % NFFT=4096;                 
        % X1=fftshift(fft(ecg_zebrafish(:,2),NFFT));         
        % Px1=X1.*conj(X1)/(NFFT*length(ecg_zebrafish)); %Power of each freq components         
        % fVals=Fs*(-NFFT/2:NFFT/2-1)/NFFT;  
        % figure;
        % plot(fVals,10*log10(Px1),'b');         
        % title('Power Spectral Density of Geless signal');         
        % xlabel('Frequency (Hz)')         
        % ylabel('Power(dB)');
        %===========================Slope calculation==============================
        x{i} = differential(fishPackage{i}(:,2),slopeLevel).^2;

        %=========================Peak Detection===================================
        %%
        peak{i} = max_detection(x{i},windowLevel);

        %%
        %============================= Peak remove=================================
        % Nhung dinh nao ma co Amplitude > 4* the average of amplitude
        % The peaks which are less than 0.8 * the average of amplitude  is removed.
        threshold(i) = mean(peak{i}(:,2));
        %=================================Debug====================================
        if (flag3 ==3)
            figure
            plot(x{no});
            title('Peak of Slope signal');
            xlabel('Sample')
            ylabel('Amplitude')
            hold on
            plot(peak{no}(:,1),x{no}(peak{no}(:,1)),'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);
            threshold(i) = mean(peak{no}(:,2));
            peak{no}(peak{no}(:,2) > 20*threshold(i),:)=[];
            length(peak{no}(:,1));
            plot(peak{no}(:,1),x{no}(peak{no}(:,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
            threshold(i) = mean(peak{no}(:,2)); % Tinh lai threshold sau khi loai nhung dinh lon hon 40* old Threshold cu
            figure
            title('Tin Hieu Do Doc Final');
            xlabel('Sample')
            ylabel('Amplitude')
            hold on
            plot(peak{no}(:,1),x{no}(peak{no}(:,1)),'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);
            peak{no}(peak{no}(:,2)< threshold(i),:)=[];
            length(peak{no}(:,1))
            plot(peak{no}(:,1),x{no}(peak{no}(:,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
            legend('DinhLoaiBo','DinhGiuLai');
            plot(x{no});
            flag3 = 1;
        end
        %==========================================================================
        if (flag1 ==1)
            peak{i}(peak{i}(:,2) > 20*threshold(i),:)=[];
            threshold(i) = mean(peak{i}(:,2)); % Tinh lai threshold sau khi loai nhung dinh lon hon 40* old Threshold cu
            peak{i}(peak{i}(:,2)< threshold(i),:)=[];
        end
        if (flag1 == 10)
            figure
            plot(fishPackage{i}(:,2));
            title(['ECG signal part ', num2str(i)]);
            xlabel('Sample')
            ylabel('Amplitude')
            hold on
            xlim([0 length(fishPackage{i}(:,2))])
            plot(peak{i}(:,1),fishPackage{i}(peak{i}(:,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
        end
       
        %%=================================Finalize R =============================

        for j = 1:length(peak{i}(:,1));

        temp(j,:)=((peak{i}(j,1) - RR_interval : peak{i}(j,1) + RR_interval)); % tim trong khoang 37 truoc va 37 sau

        R_peak{i}(j,1) = max_find(fishPackage{i}(:,2),temp(j,:));
        end
        %%
        
        %==========================Peak Final=====================================
        meanDistance = mean(diff(R_peak{i}(:,1)));
        Average_HR(i) = 60/meanDistance*Fs;
        R_peak_final{i}(:,1)= R_peak{i}(:,1);
        R_peak_final{i}(:,2) = fishPackage{i}(R_peak{i}(:,1),2);
        %%
        if (flag1 == 10)
            figure
            plot(fishPackage{i}(:,2));
            title(['ECG signal part PEAK FINAL ', num2str(i)]);
            xlabel('Sample')
            ylabel('Amplitude')
            hold on
            xlim([0 length(fishPackage{i}(:,2))])
            plot(R_peak_final{i}(:,1),fishPackage{i}(R_peak_final{i}(:,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
        end
        %%
        %%====================Trim off meaningless data==============================
        if(flag1 == 1)
          distance = diff(R_peak_final{i}(:,1));
            for k = 1:length(distance)-3
            flagx = 1;      
            for j = 1:3
                if((distance(k+j)<0.3*distance(k))||(distance(k+j)>2.3*distance(k)))
                    flagx = 0;
                    break;
                end
            end
            %find a good average
            if(flagx == 1)
                averageRPeak = (distance(k)+distance(k+1)+distance(k+2)+distance(k+3))/4;
            end
            end  
        end
        kj = 1;     
        for kj=1:length(R_peak_final{i}(:,1))-2
            if kj==length(R_peak_final{i}(:,1))-2
                break;
            end
            len = length(R_peak_final{i}(:,1));
            if (kj < len - 1)
                if((R_peak_final{i}(kj+1,1) - R_peak_final{i}(kj,1))<0.6*averageRPeak)
                    R_peak_final{i}(kj+1,:) = [];
                    %kj = kj-1;
                end
             end
        end
        %%
%         if (flag1 ==0)
%         figure
%         plot(fishPackage{no}(:,2));
%         title('ECG signal');
%         xlabel('Sample')
%         ylabel('Amplitude')
%         hold on
%         xlim([0 length(fishPackage{no}(:,2))])
%         plot(R_peak_final{no}(:,1),fishPackage{no}(R_peak_final{no}(:,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
%         end

        %%
        %============================SA detection=================================%
        if (flag1 == 2)
            j = 0;
            SA_point = zeros(length(R_peak_final(:,1)),2);
            for i = 1: length(R_peak_final(:,1)) - 3
                post = R_peak_final(i +2,1) - R_peak_final(i + 1,1);
                prior = R_peak_final(i +1,1) - R_peak_final(i,1);
               if (post >= 2*prior)
                   j = j+1;
                   SA_point(j,1) = R_peak_final(i+1,1);
                   SA_point(j,2) = R_peak_final(i+2,1);

               end
            end
            if (j == 0)
                disp('No SA symtom');
                k = 0;
            else
                k = 0;
                disp ('Yes SA symtom ')
                % Delete 0 value in SA_point
                for i = 1: length( SA_point(:,1))
                   if (SA_point(i,1) ==0)
                       if (k == 0)
                       j = i;
                       k = 1;
                       end
                   end
                   j2 = i;  
                end
                SA_point(j:j2,:) =[];
            end
        end

        %%
        if (flag1 == 2)
            figure
            plot(ecg_zebrafish(:,2));
            title('Tin Hieu ECG');
            xlabel('Sample')
            ylabel('Amplitude')
            hold on
            xlim([0 length(ecg_zebrafish)])
            plot(R_peak_final(:,1),ecg_zebrafish(R_peak_final(:,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
            if (k == 1)
                plot(SA_point(:,1),ecg_zebrafish(SA_point(:,1),2),'ro','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);
                plot(SA_point(:,2),ecg_zebrafish(SA_point(:,2),2),'ro','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);
            end
            xlim([0 length(ecg_zebrafish)]);
        end

        %%
        %============================Q & S detection===============================   
        [Q_wave{i}, S_wave{i}] = findQS(fishPackage{i} ,R_peak_final{i}, QS_window);
        if (flag1 == 2)
            figure
            plot(fishPackage{i}(:,2));
            title(['Plot Q,R,S wave on ECG signal part ',num2str(i)]);
            xlabel('Sample')
            ylabel('Amplitude')
            hold on
            plot(R_peak_final{i}(:,1),fishPackage{i}(R_peak_final{i}(:,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
            plot(Q_wave{i}(1:length(Q_wave{i}) - 1,1),fishPackage{i}(Q_wave{i}(1:length(Q_wave{i}) - 1,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',6);
            plot(S_wave{i}(1:length(S_wave{i}) - 1,1),fishPackage{i}(S_wave{i}(1:length(S_wave{i}) - 1,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerSize',6);
        end


        %%
        %====================== Find P wave & T wave===============================
         [P_wave{i}, T_wave{i}, T_wave_offset{i}] = findPT(fishPackage{i}(:,2),R_peak_final{i}(:,1), Q_wave{i},S_wave{i});
        if (flag1 == 1)
            figure
            plot(fishPackage{i}(:,2));
            title(['QRS and P, T Waves part ', num2str(i)]);
            hold on;

            plot(R_peak_final{i}(:,1),fishPackage{i}(R_peak_final{i}(:,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
            plot(Q_wave{i}(1:length(Q_wave{i}) - 1,1),fishPackage{i}(Q_wave{i}(1:length(Q_wave{i}) - 1,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',6);
            plot(S_wave{i}(1:length(S_wave{i}) - 1,1),fishPackage{i}(S_wave{i}(1:length(S_wave{i}) - 1,1),2),'ro','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerSize',6);
            plot(P_wave{i}(2:end-1,1),fishPackage{i}(P_wave{i}(2:end-1,1),2),'g*','MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',6); 
            plot(T_wave{i}(1:length(T_wave{i}) - 1,1),fishPackage{i}(T_wave{i}(1:length(T_wave{i}) - 1,1),2),'g*','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6); 
            plot(T_wave_offset{i}(:,1),fishPackage{i}(T_wave_offset{i}(:,1),2),'g+','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6); 
        end
        Average_Pwave(i) = sum(abs(P_wave{i}(:,2)))/(length(P_wave{i}) - 2);
        
        Average_Rwave(i) = sum(abs(R_peak_final{i}(2:end-1,2)))/(length(R_peak_final{i}) - 2);
        
        Ratio(i) =  Average_Pwave(i)/Average_Rwave(i);
        Average_RR(i) = sum(abs(diff(R_peak_final{i}(2:end-1,1)))/Fs)/(length(R_peak_final{i}) - 2);
        %Average_RT(i) = sum(abs(R_peak_final{i}(1:length(R_peak_final{i}) - 1,1) - T_wave{i}(1:length(T_wave{i}) - 1,1))/Fs)/(length(T_wave{i}) - 1);
        Average_RT(i) = sum(abs(R_peak_final{i}(1:length(R_peak_final{i}) - 1,1) - T_wave{i}(1:length(T_wave{i}) - 1,1))/Fs)/(length(T_wave{i}) - 1);
        Average_QT(i) = sum(abs(Q_wave{i}(1:length(R_peak_final{i}) - 1,1) - T_wave_offset{i}(:,1)))/(length(T_wave{i}) - 1);
        QTc(i) =(Average_QT(i))./sqrt(Average_RR(i));
        RTc(i) =(Average_RT(i))./sqrt(Average_RR(i));
        Average_ST(i) = sum(abs(S_wave{i}(1:length(R_peak_final{i}),1) - T_wave{i}(:,1)))/(length(T_wave{i}) - 1);
        %NoRPeak(i) = length(P_wave{i});
        %%
    end
    mean_Average_Pwave(u) = mean(Average_Pwave);
    mean_Average_Rwave(u) = mean(Average_Rwave);
    mean_Ratio(u) = mean(Ratio);
    mean_Average_HR(u) = mean(Average_HR);
    %mean_Average_RT(u) = mean(Average_RT);
    mean_Average_RR(u) = mean (Average_RR);
    mean_QTc(u,1) = mean(QTc);
    mean_QTc(u,2) = std(QTc);
   % mean_QTc(u,2) = mean(QTc);
    mean_RTc(u) = mean(RTc);
    mean_ST(u) = mean(Average_ST);
end
Ave_Ratio = mean_Ratio';
% Ave_RT = Average_RT';
% Ave_QT = Average_QT';
% NoR = NoRPeak';
% A_QTc = QTc';
% A_RTc = RTc';
Ave_RR = mean_Average_RR';
Ave_HR = mean_Average_HR';
Ave_R = mean_Average_Rwave';
Ave_RTc = mean_RTc';
Ave_QTc = mean_QTc';
Ave_ST = mean_ST';

A_HR(1,noSignal) = sum(Average_HR)/length(Average_HR);
A_QTc(1,noSignal) = sum(QTc)/length(QTc);
A_Ratio(1,noSignal) = sum(Ratio)/length(Ratio);
A_RR(1,noSignal) = sum(Average_RR)/length(Average_RR);
A_RTc(1,noSignal) = sum(RTc)/length(RTc);
A_QT(1,noSignal) = sum(Average_QT)/length(Average_QT);
A_RT(1,noSignal) = sum(Average_RT)/length(Average_RT);
A_ST(1,noSignal) = sum(Average_ST)/length(Average_ST);
%A_QRS(1,noSignal) = sum(Average_QRS)/length(Average_QRS);

DataA = readtable('C:\Users\jimmy\OneDrive\School Work\PhD\PhD Research\Methamphetamine Study\Large Exp\Data Analysis.xlsx');
DataA(1,noSignal) = array2table(A_HR(1,noSignal));
DataA(2,noSignal) = array2table(A_QT(1,noSignal));
DataA(3,noSignal) = array2table(A_QTc(1,noSignal));
DataA(4,noSignal) = array2table(A_Ratio(1,noSignal));
DataA(5,noSignal) = array2table(A_RR(1,noSignal));
DataA(6,noSignal) = array2table(A_RT(1,noSignal));
DataA(7,noSignal) = array2table(A_RTc(1,noSignal));
DataA(8,noSignal) = array2table(A_ST(1,noSignal));
%DataA(9,noSignal) = array2table(A_QRS(1,noSignal));
writetable(DataA,'C:\Users\jimmy\OneDrive\School Work\PhD\PhD Research\Methamphetamine Study\Large Exp\Data Analysis.xlsx')
end



