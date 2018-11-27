function [norm_Pts, T] = normalize2DPoints(pts)

[rows,colmn] = size(pts);
rows = rows - 1;

trans_Pts = pts(1:rows,:);  

% disp(pts);
% disp(trans_Pts);

cntr_Pts = mean(trans_Pts,2);
% disp(cntr_Pts);
cntr_X = cntr_Pts(1,1);
% disp(cntr_Pts(1,1));
cntr_Y = cntr_Pts(2,1);
% disp(cntr_Pts(2,1));

for i = 1:colmn
    
    trans_Pts(:,i) = trans_Pts(:,i) - cntr_Pts;  % subtracting mean

end

d2 = sum(trans_Pts.^2); % adding the new points
d1 = sqrt(d2);
d = mean(d1);

scl = 1.414/d;
% disp(scl);

Sx = -scl*cntr_X;
Sy = -scl*cntr_Y;
% disp(Sx);

T = [scl, 0, Sx; 0, scl, Sy; 0, 0, 1]; % T matrix

norm_Pts = T*pts; % Normalized points

end