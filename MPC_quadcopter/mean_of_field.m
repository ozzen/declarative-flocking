function [mean_data] = mean_of_field(traj, i)
    trajs_array = squeeze(struct2cell(traj));
    fnams = fieldnames(traj);

    fprintf('%s\n', fnams{i});
    cur_data = cell2mat(trajs_array(i,:));
    mean_data = mean(cur_data,2);
end

