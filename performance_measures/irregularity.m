% from Olfati-saber, but it's the same as J_p of Zhang2015
function e = irregularity(q,~,d,pnet)
N = size(q,1);
e = 0;
for i=1:N
%     distances = sqrt(sum((q - repmat(q(i,:), N, 1)).^2,2));
%     % neighbourhood
% 	N_i = distances < r; N_i(i)=0;
%     e = e + sum((distances(N_i) - d).^2);

    q_Ni = q(logical(pnet(i,:)),:);
    if ~isempty(q_Ni)
        distances = sqrt(sum((q_Ni - repmat(q(i,:),size(q_Ni,1), 1)).^2,2));
        e = e + sum((distances - d).^2);
    end
     
end
% e = e/d/sum(sum(pnet));
e = e/(sum(sum(pnet))+1);
end