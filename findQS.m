function [Q , S] = findQS(input,R,window)
Q = zeros(length(R),1);
S = zeros(length(R),1);
for i= 1: length(R)-1
    point = find(input(:,1)==R(i,1));
    % find Q point
    for j = 1: 1: window
        while( input(point, 2) >= input(point - j,2))
            point =  point -1;
        end
        Q(i) = point;
        
    end
    %%
%     while 1
%         flag=1;
%        for j=1:1:window
%         if((input(point,2)>input(point+j,2)) || (input(point,2)>(input(point-j,2))))
%             point = point-1;
%             flag=0;                    
%             break;
%         end                    
%        end
%        if(flag==1)
%              Q(i)=point;
%              break;
%        end
% end 
%     point = find(input(:,1)==R(i,1));
%   while 1
%         flag=1;
%        for j=1:1:window
%         if((input(point,2)>input(point+j,2)) || (input(point,2)>(input(point-j,2))))
%             point=point+1;
%             flag=0;
%             break;
%         end 
%                     
%        end
%        if(flag==1)
%              S(i)=point;
%              break;
%        end
%   end 
%%
point = find(input(:,1)==R(i,1));
point = point + 5;
 % find Q point
   % for j = 1: 1: window
        while( input(point, 2) >= input(point +1,2))
            point =  point +1;
        end
        S(i) = point;
       
    %end
%%
end
end