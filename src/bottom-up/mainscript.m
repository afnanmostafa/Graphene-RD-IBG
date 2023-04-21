%%% Afnan Mostafa
%%% 06/25/2022 

%% %%% Section 1: clear variables %%% %%
clear
clc
close
rng('shuffle');

%% %%% Section 2: Reading the data file and storing spatial coordinate values %%% %%
bond_density = 5;
% bond_density = [1.25,2.50,3.75,5.00,6.25,7.50,8.75,10.00,11.25,12.50,13.75,15.00,16.25,17.50,18.75,20.00,25.00,30.00,35.00,40.00,45.00]; % (in %)

for gg = 1:length(bond_density)
    [file,outputfile,rand_bot_atoms,rand_top_atoms,len,wid,stacking] = ibbgwH_generator(bond_density(gg));
end

for wp = 1:length(bond_density)
    [outputfile2] = delete_overlap_atoms(outputfile,bond_density,len,wid,stacking);
end

%% function to get ibbg at discrete bond densities
function [file,outputfile,rand_bot_atoms,rand_top_atoms,len,wid,stacking] = ibbgwH_generator(bond_density)

[file,len,wid,stacking,~] = generate_graphene(10,10,'ab',2);

%% check layers == 2
if layers ~= 2
    error('ERROR: layers: use 2 layers only; this does not work for multi-layer structures');
end

%%
lmp_input = sprintf('in%.2fbd.lmp',bond_density);
inputfile = file;
hydrogenfile = 'gh.data';
outputfile = sprintf('ibbgHab10x10_%.2f.data',bond_density);

[cell_column,~] = readtextfile(file,5,16,'','#');

index=cell_column{1};
atom_type=cell_column{2};
x = cell_column{3};
y = cell_column{4};
z = cell_column{5};

x_bot = x(1:length(x)/2);
x_top = x(length(x)/2+1:end);
y_bot = y(1:length(y)/2);
y_top = y(length(y)/2+1:end);
z_bot = z(1:length(z)/2);
z_top = z(length(z)/2+1:end);

edge_cases = 0;

[bot_align_atoms,top_align_atoms] = get_aligned_atoms(index,x,y,stacking);

%% %%% Section 2.1: random bonds' selection depending on bonding density %%% %%
%% bond density calculation
bonding_density_in_fraction = bond_density/100;
random_C_atoms = round(length(atom_type)*bonding_density_in_fraction);

%% edge case 1
if edge_cases == 1
    if length(bot_align_atoms)*2 < random_C_atoms
        error('fix your int_bond_density; lower it slightly unless random_C_atoms => length(bot_align_atoms');
    else
    end
end

%% check: max aligned atoms not less than randonly selected C atoms
max_align_atoms = length(bot_align_atoms)*2;
if max_align_atoms < random_C_atoms
    old_rand = random_C_atoms;
    difference = abs(random_C_atoms - max_align_atoms);
    random_C_atoms = random_C_atoms - difference;
    maxBD = (random_C_atoms/length(atom_type))*100;
    userBD = bond_density;
    BD_change = userBD - maxBD;
    fprintf(2,'MAX C atoms for interlayer bonding is: %d; \nreducing random_C_atoms (current: %d) to match this number for this specific structure\n',max_align_atoms,old_rand);
    fprintf(2,'Max BD: %0.4f%%, User input = %.4f%%, diff = %.4f%%\n',maxBD,userBD,BD_change);
else
    % do nothing or remove the "else"
end

[indx_rand] = gen_rand_bot(bot_align_atoms,bond_density);
rand_bot_atoms = bot_align_atoms(indx_rand(1:end));
rand_top_atoms = top_align_atoms(indx_rand(1:end));

%% %%% Section 2.2: select atoms within range %%% %%
[hydro_C,hydro_C_top,cout,cout2] = get_hydro(rand_bot_atoms,rand_top_atoms,x_bot,x_top,y_top,x,y);

%% %%% Section 2.3: writing LAMMPS input script for C-C atoms to be moved towards each other in order to be bonded %%% %%
%%% generating LAMMPS input file which will, ultimately, give us the
%%% interlayer-bonded bilayer graphene structure (post_equil.data)

write_equil_lmps(lmp_input,rand_bot_atoms,rand_top_atoms,len,wid,stacking,bond_density);

%% %%% Section 2.4: data file including the H atoms (No C atoms in this data file) %%% %%

write_hydro_lmps(hydrogenfile,atom_type,x_bot,y_bot,z_bot,z_top,x,y,hydro_C,hydro_C_top,cout,cout2)

onetimerun= true;

if onetimerun == true
    
    TOTAL_ATOMS = length(atom_type)+cout+cout2;
    
    % to not re-execute command 3 again and again (for test purpose)
    str0 = sprintf("rm %s ",outputfile);
    command0 = (str0);
    system(command0);
    
    str1 = sprintf("cat %s %s > %s",inputfile, hydrogenfile, outputfile);
    command1 = (str1);
    system(command1);
    
    str2 = sprintf('sed -i "4 c 3 atom types" %s',outputfile);
    command2 = str2;
    system(command2);
    
    str3 = sprintf('sed -i "13 a 3 1.00784" %s',outputfile);
    command3 = str3;
    system(command3);
    
    str4 = sprintf('sed -i "2 c %d atoms" %s',TOTAL_ATOMS,outputfile);
    command4 = str4;
    system(command4);
    
end
end

%% %%%                                            End of script                                                             %%% %%