function [bot_align_atoms,top_align_atoms_added] = get_aligned_atoms(index,x,y,stacking)
%   picks atoms having same x and y coordinates but z
%   this function takes indices, x-coordinate, and y-coordinates of atoms
%   as input and picks the aligned atoms, returns their indices

x_bot = x(1:length(x)/2);
x_top = x(length(x)/2+1:end);
y_bot = y(1:length(y)/2);
y_top = y(length(y)/2+1:end);




bot_align_atoms = zeros(length(x_top),1);
top_align_atoms = zeros(length(x_top),1);

counter1 = 0;
for i = 1:length(x_bot)
    for j = 1:length(x_top)
        if (y_bot(i) == y_top(j)) && (x_bot(i) == x_top(j))
            counter1 = counter1 + 1;
            bot_align_atoms(counter1,1) = index(i);
            top_align_atoms(counter1,1) = index(j);
        end
    end
end

off_limit = 1.5;

x_boundary_max = max(x)-off_limit;
x_boundary_min = min(x)+off_limit;

y_boundary_max = max(y)-off_limit;
y_boundary_min = min(y)+off_limit;

bot_align_atoms = bot_align_atoms(1:counter1);
top_align_atoms = top_align_atoms(1:counter1);

top_align_atoms_added = zeros(1,1);

atoms_not_at_boundary = (x_bot(bot_align_atoms) < x_boundary_max)...
    & (x_bot(bot_align_atoms) > x_boundary_min) &...
    (y_bot(bot_align_atoms) < y_boundary_max)...
    & (y_bot(bot_align_atoms) > y_boundary_min);

bot_align_atoms = bot_align_atoms(atoms_not_at_boundary);

top_align_atoms = top_align_atoms(atoms_not_at_boundary);

%% getting index differences

diff = unique(abs(bot_align_atoms - top_align_atoms));

for gg = 1:length(top_align_atoms)
    if (abs((top_align_atoms(gg)) - (bot_align_atoms(gg))) == diff(1))
        top_align_atoms_added(gg) = bot_align_atoms(gg) + length(x_bot)-diff(1);
    else
        top_align_atoms_added(gg) = bot_align_atoms(gg) + length(x_bot)-diff(2);
    end
    
end

top_align_atoms_added = top_align_atoms_added';


if strcmp(stacking,'aa')
    top_align_atoms_added = top_align_atoms_added(1:2:end);
    bot_align_atoms = bot_align_atoms(1:2:end);
end


end

