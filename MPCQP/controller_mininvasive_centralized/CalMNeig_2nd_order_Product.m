function [MNeig,H,f,A,b,lb,ub,alp] = CalMNeig_2nd_order_Product(pos,vel,uhat,params)
% 2nd order CalMNeig calculates CBF constraints based on second order dynamics
% v2 changes CBF to 1/h
alp = 0;
N = length(pos(:,1));
MNeig = zeros(N, N); % construct Neigbor matrix
damax = 2*params.amax; %delta amax
H = 2*eye(2*N);
f = -2*reshape(uhat',2*N,1);

A = []; b = []; %compute A, b, H, f matries for quadprog
% h_hist = zeros(1, nchoosek(params.n,2)); 
% Lfh_hist = h_hist;
% Lgh_hist = zeros( nchoosek(params.n,2), 2*params.n);
idx = 1;
h_hist = []; Lfh_hist =[]; Lgh_hist = [];
for i = 1:(N-1)
    for j = (i+1):N
        if norm(pos(i,:)-pos(j,:)) < 3.07
        MNeig(i,j) = 1;
        D = norm(pos(i,:) - pos(j,:));
        dx12 = (pos(i,:) - pos(j,:))*(vel(i,:) - vel(j,:))';
        dx22 = (vel(i,:) - vel(j,:))*(vel(i,:) - vel(j,:))';
        h = sqrt(2*damax*(D-params.Ds)) + dx12/D;
        
        Lfh = sqrt(damax/2/(D-params.Ds))*dx12/D - dx12^2/D^3 + dx22/D;
        Lgh = zeros(1,2*N);
        Lgh((2*i-1):(2*i)) = (pos(i,:) - pos(j,:))/D;
        Lgh((2*j-1):(2*j)) = -(pos(i,:) - pos(j,:))/D;
            h_hist = [h_hist h];
            Lfh_hist = [Lfh_hist Lfh];
            Lgh_hist = [Lgh_hist; Lgh];
%         h_hist(idx) = h;
%         Lfh_hist(idx) = Lfh;
%         Lgh_hist(idx,:) = Lgh;
%         idx = idx + 1;
        end
    end
end

Nb = length(h_hist); % #of barriers
h_bar = prod(h_hist);
if Nb ~= 0
    A = zeros(1,2*N);
    b = alphaFunction(h_bar, params.gamma); %divide h_bar on both sides
    alp = alphaFunction(h_bar, params.gamma);
end
for i = 1:Nb
    A = A - Lgh_hist(i,:)/h_hist(i);
    b = b + Lfh_hist(i)/h_hist(i);
end

lb = -1.2 * params.amax*ones(2*params.n,1);
ub = 1.2 * params.amax*ones(2*params.n,1);

end
