figure;
%%
A1 = 1;
A2 = 8;

v1 = [vx(:,A1) vy(:,A1)];
v2 = [vx(:,A2) vy(:,A2)];

p1 = [x(:,A1) y(:,A1)];
p2 = [x(:,A2) y(:,A2)];

v12 = v1 - v2;
p12 = p1 - p2;

bd = traj.bd;

%%
% subplot(2,1,1);

distance = sqrt(sum(p12.^2, 2));
vRel = -sum(v12 .* p12, 2) ./ distance;

steps = floor((vRel) / (params.amax * params.dt));
d_br_max = params.dt * (steps .* vRel - 0.5 * (steps - 1) .* steps * params.amax * params.dt);

v_rel = v12;
p_rel = p12;
p_rel_n = p_rel ./ repmat(sqrt(sum(p_rel.^2,2)), 1, 2);
delta_v = - sum(p_rel_n .* v_rel, 2);
h_pairwise = 2 * sqrt(params.amax * (sqrt(sum(p_rel.^2, 2)) - params.Ds)) -...
    delta_v;


% h_condition = plot(d_br_max + dmin, 'r', 'LineWidth', 1.2);
% hold on
% h_dist = plot(distance, 'b', 'LineWidth', 1.2);
% H1 = [h_dist, h_condition];
% legend(H1, {'d_{ij}', 'd_{br}+d_{min}'},'FontSize',14);
hold on
h_barrierfunction_pw = plot(real(h_pairwise), 'r', 'LineWidth', 1.2);
lie_derivative_pw = plot(bd(A2, :, A1), 'b', 'LineWidth', 1.2);

%%
% subplot(2,1,2);

h_policy1 = plot(4*policy(:,A1)-4,'--', 'Color', [0.5 0.5 0.5],'LineWidth', 1.2);
h_policy2 = plot(4*policy(:,A2)-4,':', 'Color', [0.5, 0.5, 0.5],'LineWidth', 1.8);

% H2 = [h_policy1, h_policy2];
% legend(H2, {['policy agent ' num2str(A1)], ['policy agent ' num2str(A2)]},'FontSize',14);

H = [h_barrierfunction_pw, lie_derivative_pw, h_policy1, h_policy2];
% H = [h_barrierfunction_pw, h_policy1, h_policy2];
legend(H, {'h_{ij}', 'Lie-derivative',['policy agent ' num2str(A1)], ['policy agent ' num2str(A2)]},'FontSize',14);

set(gcf, 'Position', [0 300 1440 200]);
saveas(gcf, [destPath dirName '/' num2str(A1), '_' num2str(A2) '.jpg']);

%%
plot(sqrt(sum(p12.^2, 2)), 'g', 'LineWidth', 1.5)





