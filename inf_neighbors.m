function [no_inf] = inf_neighbors(node, statev, g)
% A function that counts the number of infectious neighbours 'node' has, as
% a member of graphlet 'g', in state 'statev'.

no_inf = 0;
for i = 1:length(g)
    if g(node,i)== 1 && statev(i) == 0
        no_inf = no_inf  + 1;
    end
end
end

