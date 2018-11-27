function [F, idxs] = estimateFundamentalMatrixRANSAC(pts1 ,pts2)

function f = estimateFundamentalMatrix(pts1, pts2)

    sz = size(pts1,2);
%     disp(sz);
    arr = zeros(sz,9);
    
    [pts1, T1] = normalize2DPoints(pts1);
    [pts2, T2] = normalize2DPoints(pts2);
    
    for ii=1:sz
        
        trans_pts1 = transpose(pts1(:,ii));
        trans_pts2 = transpose(pts2(:,ii));
        arr(ii,:) = kron(trans_pts1, trans_pts2);
    
    end

    [~,~,V] = svd(arr); %Singular value decomposition
    f1 = transpose(reshape(V(:, end), 3, 3)); % writing in the form of 3*3 matrix
    [U, S, V] = svd(f1); % Again SVD
    new_S = diag([S(1, 1), S(2, 2), 0]);
    
    f2 = U * new_S * transpose(V);
    f = transpose(T2) * f2 * T1; % calculating F
    f = f/f(end);
    
end

s = 8;
[~,N] = size(pts1);
e =0.65;
p = 0.99;
Max_Inlier = 0;
idxs =[];
thrs = 0.025;

Max_Num_Iter = log(1-p)/log(1-(1-e).^s);

for i=1:Max_Num_Iter
    
    tmpNumCrctMatch = 0;    % Temporary number of correct matches
    tmpMtchdIdx = []; % Temporary matched index
    idx = randperm(N,8); % randomly selecting the permutation
    tmp_F = estimateFundamentalMatrix(pts1(:,idx),pts2(:,idx)); % initially estimating F

    for j=1:N
        tmp_val = transpose(pts2(:,j)) * tmp_F * pts1(:,j);
        % disp(tmp_val);
        if(abs(tmp_val) < thrs) % setting threshold
            
            tmpNumCrctMatch = tmpNumCrctMatch + 1;
            % disp(tmpNumCrctMatch);
            tmpMtchdIdx = [tmpMtchdIdx;j];
            % disp(tmpMtchIdx);
            
        end
    end

    if (tmpNumCrctMatch == N) % number of correct matches = max inliers
        
        Max_Inlier = N;     % max number of inliers
        idxs = tmpMtchdIdx;     % indices
        break;
        
    end

    if (tmpNumCrctMatch > Max_Inlier) % To get maximum number of inliers 
        
        idxs = tmpMtchdIdx;
        Max_Inlier = tmpNumCrctMatch;
    
    end
end

inlr1 = pts1(:,idxs);   % Inliers in image1
inlr2 = pts2(:,idxs);   % Inliers in image2
F = estimateFundamentalMatrix(inlr1,inlr2); % Final F matrix

end