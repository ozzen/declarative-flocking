function C = connectivity(p_net)
    N = size(p_net,1);
    C = rank(p_net)/(N-1); 
end