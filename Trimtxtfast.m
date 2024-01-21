function output = Trimtxtfast(FileStart,FileEnd,MainProject,Folder)
% Trim multiple ECG txt files, move to Processing folder for Data Analysis
close all


DataStore = ['C:\Users\JimmyZhang\OneDrive - personalmicrosoftsoftware.uci.edu\School Work\PhD\PhD Research\' MainProject '\' Folder '\Data for Processing'];
TrimTable = readtable(['C:\Users\JimmyZhang\OneDrive - personalmicrosoftsoftware.uci.edu\School Work\PhD\PhD Research\' MainProject '\' Folder '\TrimTable.xlsx']);
Day = table2array(TrimTable(:,1));
Group = table2array(TrimTable(:,2));
Start = table2array(TrimTable(:,3));
End = table2array(TrimTable(:,4));
for a = FileStart:FileEnd
    data = [];
    new = [];
    FileNum = num2str(a);
    if isnumeric(Day(a))
        DayS = num2str(Day(a));
    else
        DayS = Day(a);
    end
    GroupS = Group{a};
    StartS = Start(a);
    EndS = End(a);
input = ['C:\Users\JimmyZhang\OneDrive - personalmicrosoftsoftware.uci.edu\School Work\PhD\PhD Research\' MainProject '\' Folder '\Day ' DayS '\Grp' GroupS '.txt'];
input1 = ['C:\Users\JimmyZhang\OneDrive - personalmicrosoftsoftware.uci.edu\School Work\PhD\PhD Research\' MainProject '\' Folder '\Day ' DayS '\Grp' GroupS 'U.txt'];
NewFileName = [DataStore '\' FileNum '.txt'];
data = load(input);
for i = StartS:EndS
    new(i-StartS+1,:) = data(i,:);
    new(i-StartS+1,2) = new(i-StartS+1,2);
end
if isfile(input1)
    input2 = ['C:\Users\JimmyZhang\OneDrive - personalmicrosoftsoftware.uci.edu\School Work\PhD\PhD Research\' MainProject '\' Folder '\Day ' DayS '\Grp' GroupS 'UU.txt'];
    writematrix(new,input2,'Delimiter','tab');
    movefile(input2,DataStore);
    NewLocationName = [DataStore '\Grp' GroupS 'UU.txt'];
    movefile(NewLocationName,NewFileName);
else
    writematrix(new,input1,'Delimiter','tab');
    movefile(input1,DataStore);
    NewLocationName = [DataStore '\Grp' GroupS 'U.txt'];
    movefile(NewLocationName,NewFileName);
end
end
