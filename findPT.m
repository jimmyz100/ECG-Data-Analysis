function [ P_wave ,T_wave, T_wave_offset] = findPT( input, RR_index, Q_index, S_index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here   , S_index, S_window, n
% Find P wave
P_index_temp = zeros(length(Q_index),1);
P_wave = zeros(length(Q_index),2);
T_wave = zeros(length(S_index),2);
index = 0;
RR_interval_T = diff(RR_index);
RR_interval_Q = diff(RR_index(1: end));
inputT = cell(length(S_index)-1,1);
inputQ = cell(length(Q_index)-1,1);
Tcoef1 = 0.5;%0.5
Tcoef12 = 0.2;
Qcoef = 0.2;%0.25
for i = 2: length(Q_index)-1;
    inputQ{i} =  input(Q_index(i) - round(Qcoef*RR_interval_Q(i-1)) : Q_index(i));
    max_temp = -20;
    for j = 1: round(Qcoef * RR_interval_T(i-1))
       if  (inputQ{i}(j) > max_temp)
           max_temp = inputQ{i}(j);
           index = j;
       end  
    end
    P_wave(i,1) = Q_index(i) + index - round(Qcoef*RR_interval_T(i-1)) - 1;
    P_wave(i,2) = input( P_wave(i,1));
%    P_index_temp(i) = high_detection(input(Q_index(i) - Q_window : Q_index(i)), 25);
%    P_index(i) = Q_index(i) + P_index_temp(i) - Q_window - 1;
end
%============================Find T wave===================================


for i = 1: length(S_index)-1
    
    m = 1;
    inputT{i} =  input(S_index(i)+round(Tcoef12*RR_interval_T(i)): S_index(i)+ round(Tcoef1*RR_interval_T(i)));
    max_temp = -20;
    min_temp = 10;
    for j = 1: round(Tcoef1*RR_interval_T(i))-round(Tcoef12*RR_interval_T(i))+1
       if  (inputT{i}(j) > max_temp)
           max_temp = inputT{i}(j);
           index_max = j;
       end  
       if (inputT{i}(j) < min_temp)
           min_temp = inputT{i}(j);
           index_min = j;
       end
    end
    Tmax_wave(i,1) = S_index(i) + index_max + round(Tcoef12*RR_interval_T(i));
    Tmax_wave(i,2) = input(Tmax_wave(i,1));
    Tmin_wave(i,1) = S_index(i) + index_min + round(Tcoef12*RR_interval_T(i));
    Tmin_wave(i,2) = input( Tmin_wave(i,1));
    if (abs(Tmax_wave(i,2)- input(S_index(i))) > abs(Tmin_wave(i,2)-input(S_index(i))))
        T_wave(i,1) = S_index(i) + index_max + round(Tcoef12*RR_interval_T(i));
        T_wave(i,2) = input(T_wave(i,1));
    else
        T_wave(i,1) = S_index(i) + index_min + round(Tcoef12*RR_interval_T(i));
        T_wave(i,2) = input(T_wave(i,1));
    end
    % T_offset detection
    if (input(T_wave(i,1)) > 0 )
       while (input(T_wave(i,1) + m) > 0)
        m = m +1; 
       end
    else
        while (input(T_wave(i,1) + m) < 0)
        m = m +1; 
       end
    end
    T_wave_offset(i,1) = T_wave(i,1) + m;
    T_wave_offset(i,2) = input(T_wave(i,1) + m);
end
end

