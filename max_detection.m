function output = max_detection(input,parameter)

[x,y] = size(input);
if x==1
    if y<2*parameter+1
        
        warning('input is short');
    else
        
        input_final = input';
    end
else
    if x<2*parameter+1
        
        warning('input is short');
    else
        
        if y==1
            
            input_final = input;
        else
            
            warning('input must be vector');
        end
    end
        
end
Window = cell(length(input_final)-2*parameter,1);
peak = zeros(length(input_final)-2*parameter,2);
position_max_Window = cell(length(input_final)-2*parameter,1);

for i = 1:length(input_final)-2*parameter
    
    Window{i,1}(:,1) = input_final(i:2*parameter+i,1);
    max_Window = max(Window{i,1}(:,1));
    position_max_Window{i,1}(:,1) = find(Window{i,1}(:,1)==max_Window);
    if length(position_max_Window{i,1}(:,1))==1
        
        if Window{i,1}(parameter+1,1) == max_Window
            
            position_max = parameter+i;
            max_value = max_Window;
        else
            position_max = 0;
            max_value = 0;
        end
    else
        
        if Window{i,1}(parameter+1,1) == max_Window
            for j = 1:parameter+i
                if input_final(parameter+i-j)==Window{i,1}(parameter+1,1);
                    continue;
                else
                    if input_final(parameter+i-j)> Window{i,1}(parameter+1,1)
                        position_max = 0;
                        max_value = 0;
                        break;
                    else
                        position_max = parameter+i-j+1;
                        max_value = max_Window;
                        break;
                    end
                end
            end
        else
            position_max = 0;
            max_value = 0;
        end
    end
peak(i,1) = position_max;
peak(i,2) = max_value;
end
find_zero = find(peak(:,1)==0);
output = removerows(peak,'ind',find_zero);
output = unique(output,'rows');
if isempty(output)
    warning('signal has not peak');
end

end

%==========================================================================
% step 1: tao cua so chay
% step 2: so sanh voi max cua cua so. neu diem chinh giua la max thi nhan
% step 3: loai cac truong hop nhieu peak lien nhau vi cung bang gia tri
% step 4: chon gia tri bang nhau cuoi cung ben tay trai