function [ resA, resB, alp ] = GenAandb( pos, vel, u, params)

damax = 2 * params.amax;
vmax = params.vmax;
Dmax = params.Dmax;

dt = params.dt;

Ds = params.Ds;
gamma = params.gamma;
hor = params.h;

pos = pos';
vel = vel';

N = params.n;
alp = 0;

resA = zeros(hor, 2*N*hor+1);
resB =zeros(hor, 1);

for k = 0:hor-1
    h_hist = zeros(1, nchoosek(params.n,2)); 
    Lfh_hist = h_hist;
    Lgh_hist = zeros( nchoosek(params.n,2), 2*params.n);
%     h_hist = []; Lfh_hist =[]; Lgh_hist = [];
    idx = 1;
    for i = 1:(N-1)
        for j = (i+1):N
            D = norm(pos(i,:) - pos(j,:));
            dx12 = (pos(i,:) - pos(j,:))*(vel(i,:) - vel(j,:))';
            dx22 = (vel(i,:) - vel(j,:))*(vel(i,:) - vel(j,:))';

            h = sqrt(2*damax*(D-Ds)) + dx12/D;

            Lfh = sqrt(damax/2/(D-Ds))*dx12/D - dx12^2/D^3 + dx22/D;

            Lgh = zeros(1,2*N);
            Lgh((2*i-1):(2*i)) = (pos(i,:) - pos(j,:))/D;
            Lgh((2*j-1):(2*j)) = -(pos(i,:) - pos(j,:))/D;

%                 h_hist = [h_hist h];
%                 Lfh_hist = [Lfh_hist Lfh];
%                 Lgh_hist = [Lgh_hist; Lgh];

            h_hist(idx) = h;
            Lfh_hist(idx) = Lfh;
            Lgh_hist(idx,:) = Lgh;
            idx = idx + 1;
        end
    end
    
    h_bar = prod(h_hist);

    A = sum(-Lgh_hist./repmat(h_hist', 1, 2*N), 1);
    b = sum(Lfh_hist./h_hist) + alphaFunction( h_bar, gamma );
    if k == 0
        alp = alphaFunction( h_bar, gamma );
    end
    
    resA(k + 1, k*N + 1: (k+1)*N) = A(1:2:end);
    resA(k + 1, (hor+k)*N + 1: (hor+k+1)*N) = A(2:2:end);
    resB(k+1) = b;
    
    %% Update state
    acc = [u(k*N+1:(k+1)*N) u((hor+k)*N + 1: (hor+k+1)*N)];
    vel = vel + params.ct*acc;
        for jj = 1:N
            if norm(vel(jj,:)) > vmax
                vel(jj,:) = (vmax/norm( vel(jj,:) )) * vel(jj,:);
            end
        end
    pos = pos + params.ct*vel;
    
end
[A_constrained, b_constrained] = GenConstAcc(params, vel');

resA = [resA; A_constrained];
resB = [resB; b_constrained];

end


