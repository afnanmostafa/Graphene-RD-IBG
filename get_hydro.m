function [hydro_C,hydro_C_top,cout,cout2] = get_hydro(rand_bot_atoms,rand_top_atoms,x_bot,y_bot,x_top,x,y)
%   returns indices of passivating H atoms
%   gets 3 nearest neighbor C atoms of interlayer bond forming C atom in
%   each layer and passivates them with H

mid_atom_list = rand_bot_atoms;
hydro_C = zeros(1,1);
mid_x = x(mid_atom_list);
mid_y = y(mid_atom_list);

rad = 1.5; %% radius of circle centering C atoms that are forming int. bond
cout = 0;
for ww = 1:length(mid_atom_list)
    for ia = 1:length(x_bot)
        if mid_atom_list(ww) ~= ia
            if (((x(ia)-mid_x(ww))^2)+(y(ia)-mid_y(ww))^2) < rad^2
                cout = cout+1;
                hydro_C(cout,1) = ia;
            end
        end
        
    end
end



mid_atom_list_top = rand_top_atoms;
hydro_C_top = zeros(1,1);

mid_x_top = x(mid_atom_list_top);
mid_y_top = y(mid_atom_list_top);
cout2 = 0;

for vv = 1: length(mid_atom_list_top)
    for ib=length(x_top)+1:length(x)
        if mid_atom_list_top(vv) ~= ib
            if (((x(ib)-mid_x_top(vv))^2)+(y(ib)-mid_y_top(vv))^2) < rad^2
                cout2 = cout2+1;
                hydro_C_top(cout2,1) = ib;
                
            end
        end
        
    end
end
end

