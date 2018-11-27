clc;
close all;

img1 = imread('../images/img1.png');
img2 = imread('../images/img2.png');

img1_gray = rgb2gray(img1);
img2_gray = rgb2gray(img2);

point1 = detectSURFFeatures(img1_gray,'MetricThreshold',10);
point2 = detectSURFFeatures(img2_gray,'MetricThreshold',10);

[f1,vpts1] = extractFeatures(img1_gray,point1);
[f2,vpts2] = extractFeatures(img2_gray,point2);

indexPairs = matchFeatures(f1,f2);
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));

% disp(matchedPoints1);

mch_trans1 = transpose(matchedPoints1.Location);
mch_trans2 = transpose(matchedPoints2.Location);

% disp(mch_trans1);

ones_trans1 = ones(1,size(mch_trans1,2));
ones_trans2 = ones(1,size(mch_trans2,2));

hom_pts1 = [mch_trans1;ones_trans1];
hom_pts2 = [mch_trans2;ones_trans2];

% disp(hom_pts1);


[F, idxs] = estimateFundamentalMatrixRANSAC(hom_pts1, hom_pts2);
fprintf("Fundamental Matrix \n");
disp(F);
% disp(idxs);

K = [558.7087,  0.0,  310.3210; 0.0,  558.2827,  240.2395; 0.0,  0.0,  1.0];
 
E = K' * F * K;
% disp(E);
[u, d, v] = svd(E);

new_d = diag([(d(1,1) + d(2,2))/2, (d(1,1) + d(2,2))/2, 0]);

E = u * new_d * v';
fprintf("Essential Matrix \n");
disp(E);

hom_pts1 = hom_pts1(:,idxs);
hom_pts2 = hom_pts2(:,idxs);

[R, t] = decomposeEssentialMatrix(E, hom_pts1, hom_pts2, K);
fprintf("Rotational Matrix \n");
disp(R);
fprintf("Translation Matrix \n");
disp(t);

ProjMat_1 = K*[eye(3,3) [0 0 0]'];
ProjMat_2 = K*[R,t];
pts_3D = algebraicTriangulation(hom_pts1, hom_pts2, ProjMat_1, ProjMat_2);
pts_3D_hom = pts_3D./repmat(pts_3D(4,:), 4, 1);
ones_arr = ones(1,4);
Rot_trans = [R, t; ones_arr];
% disp(Rot_trans);

figure;
plot3(pts_3D_hom(1,:),pts_3D_hom(2,:),pts_3D_hom(3,:),'*');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D structure');
hold on;
T = [R, t; 0 0 0 1];
plotCameraFrustum(T, 'r', 0.3);
hold off;

figure; 
showMatchedFeatures(img1_gray, img2_gray, matchedPoints1, matchedPoints2,'montage','Parent',axes );
legend('matched points 1','matched points 2');

% imshow(img1_gray);
% imshow(img2_gray);