function [outputfile] = delete_overlap_atoms(file,bond_density,len,wid,stacking)
%   deletes overlapping atoms from the structure

file_to_edit = file;  %% STC

%% %%% Section 3.1:get all files in current directory and select our desired file %%


%% %%% Section 3.2:location of file %%
% file_loc = [directory '\' file_picked];

temp_file1 = 'tmp1_raw.data';
temp_file2 = 'tmp2_initials.lmp';
output_temp = 'overlap_del_tmp3.data';
output_temp2 = 'overlap_del_tmp4.data';
outputfile = sprintf('gr%dx%d%s_%.2fBD.data',len,wid,stacking,bond_density);

A = regexp(fileread(file_to_edit),'\n','split');
whichline = find(contains(A,'Atoms'));
whichline = whichline+1;

str = sprintf('awk "NR > %d {print}" %s > %s',whichline,file,temp_file1);
command = (str);
system(command);


fid = fopen(temp_file1);
tline = fgetl(fid);
tlines = cell(0,1);

while ischar(tline)
    
    disp(tline);
    tlines{end+1,1} = tline;   %listing each line as a separate array
    tline = fgetl(fid);
    
end

fclose(fid);

%% %%% Section 3.3: deleting overlapping atoms (atoms with same x, y, z values) %%% %%

rows = length(tlines);
str = split(tlines);
lines = str;
x1 = zeros(rows,1);
y1 = zeros(rows,1);
z1 = zeros(rows,1);

for t = 1: rows
    
    x1(t) = str2double(str{t,3});
    y1(t) = str2double(str{t,4});
    z1(t) = str2double(str{t,5});
    
end

list = zeros(rows,1);

for y = 1 : length(tlines)
    
    str = tlines(y);
    k_values = ['NaN'];
    indices = zeros(1,1);
    isK = cellfun(@(x)isequal(x,k_values),str);
    list(y) = isK;
    
end

list1=[];

for i = 1: rows
    
    for j = i+1: rows
        
        if x1(i) == x1(j) && y1(i) == y1(j) && z1(i) == z1(j)
            
            tlines{i,1} = [];
            lines{i,1}(1:end) = [];
            lines{i,2}(1:end) = [];
            lines{i,3}(1:end) = [];
            lines{i,4}(1:end) = [];
            lines{i,5}(1:end) = [];
            list1(i) = i;
            
        end
        
    end
    
end

list1(list1 == 0) = [];
list1=list1';


fid = fopen(temp_file1);
s = textscan(fid,'%f %f %f %f %f');
fclose(fid);

ind = s{1};
type=s{2};
x= s{3};
y=s{4};
z=s{5};

for h = [list1]
    
    ind(h) = [];
    type(h)=[];
    x(h)= [];
    y(h)=[];
    z(h)=[];
    
end

%% %%% Section 3.4: writing output file after deletion %%% %%

fid2 = fopen(output_temp, 'w');

for o = 1: length(x)
    
    fprintf(fid2, '%d\t %f\t %f\t %f\n',type(o),x(o),y(o),z(o));
    
end

fclose(fid2);

onetimerun= true;
if onetimerun == true
    
    % to not re-execute command 3 again and again (for test purpose)
    str0 = sprintf("rm %s'",outputfile);
    command0 = (str0);
    system(command0);
    
    str1 = sprintf('awk "{print NR, $1, $2, $3, $4, $5}" %s > %s',output_temp,output_temp2);
    command1 = (str1);
    system(command1);
    
    str2 = sprintf('awk "NR < %d {print}" %s > %s',whichline+1, file, temp_file2);
    command2 = (str2);
    system(command2);
    
    str3 = sprintf('cat %s %s > %s',temp_file2,output_temp2,outputfile);
    command3 = (str3);
    system(command3);
    
    
    TOTAL_ATOMS = length(type);
    str4 = sprintf('sed -i "2 c %d atoms" %s',TOTAL_ATOMS,outputfile);
    command4 = str4;
    system(command4);
    
    str5 = sprintf("rm %s %s %s %s'",temp_file1, temp_file2, output_temp, output_temp2);
    command5 = (str5);
    system(command5);
end
end

