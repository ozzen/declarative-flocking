function [x, y, vx, vy, ax, ay] = read_shouvik(filename, rept)
    load(filename);

    x = pred_traj(rept).q(:,:,1);
    y = pred_traj(rept).q(:,:,2);
    
    vx = pred_traj(rept).p(:,:,1);
    vy = pred_traj(rept).p(:,:,2);
    
    ax = pred_traj(rept).a(:,:,1);
    ay = pred_traj(rept).a(:,:,2);   
end

