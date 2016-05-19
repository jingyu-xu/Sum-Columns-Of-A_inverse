fid = fopen('breastinfo_simple.txt','r');
breastID = str2double(fgetl(fid));
s1 = str2double(fgetl(fid));
s2 = str2double(fgetl(fid));
s3 = str2double(fgetl(fid));
class = str2double(fgetl(fid));
fclose(fid);

load mtype.mat;
load pval.mat;

muscle_wall = 153;
skin_start = 138;

% Convert vector into cube
mtype_cube = zeros(s1,s2,s3); % each voxel is .5mmx.5mmx.5mm
pval_cube = zeros(s1,s2,s3);
cur_pos = 1;
for k=1:s3
    for j=1:s2
        for i= 1:s1
            mtype_cube(i,j,k) = mtype(cur_pos);
            pval_cube(i,j,k) = pval(cur_pos);
            cur_pos = cur_pos + 1;
        end 
    end
end

% subsample cubes in order to solve sparse matrix
s1_ss = floor(s1/2); % voxels are now 1mmx1mmx1mm
s2_ss = floor(s2/2);
s3_ss = floor(s3/2);
xi = 1; yi = 1; zi = 1;
mtype_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
pval_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
for z=1:2:s3-1
    for y = 1:2:s2
        for x = 1:2:s1
            mid = mtype_cube(x,y,z);
            pid = pval_cube(x,y,z);
            mtype_cube_subsamp(xi,yi,zi) = mid;
            pval_cube_subsamp(xi,yi,zi) = pid;
            xi = xi+1;
        end
        xi = 1;
        yi = yi + 1;
    end
    yi = 1;
    zi = zi + 1;
end
% some voxels not converted to muscle during subsampling
% so do that now
% still need to figure out how to get pval converted for fdtd
for x=1:s1_ss
    for y=1:s2_ss
        for z=1:s3_ss
            if x > 153
               mtype_cube_subsamp(x,y,z) = -4;
                pval_cube_subsamp(x,y,z) = 1;
            end
        end
    end
end

% save mtype_cube_subsamp; save pval_cube_subsamp.mat;
% load mtype_cube_subsamp.mat;
[s1_ss, s2_ss, s3_ss] = size(mtype_cube_subsamp);
model = mtype_cube_subsamp; air_id = -1;
% sliceomatic(model)
% model = model(18:end,:,:);
% [s1_ss, s2_ss, s3_ss] = size(model);

tumor_on = 0; tumor_depth = 10;  tumor_radius = 10; Tambient = 27; Tart = 37; 
[A_inverse_nom_sum,tum_pos_cen_nom] = gen_columns_of_A_inverse(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start);

% Generate temperature anomalies with radius = 10
tumor_on = 1; tum_y_cen = 90; tum_z_cen = floor(s3_ss/2);

tumor_radius = 10; tumor_depth = 10; tum_x_cen = 20 + tumor_depth;% tumor dept is 1cm
[A_inverse_abn1_sum,tum_pos_cen_abn] = gen_columns_of_A_inverse(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);

% convert the vector of sum of columns of A inverse to 3d
A_inverse_abn1_sum_3d = convert_1d_to_3d(A_inverse_abn1_sum,s1_ss,s2_ss,s3_ss);
% [m,n] = size(A_inverse_abn1_sum);
% Plot the sum of columns of A inverse

% plot 2d surface
A_inverse_2d = reshape(A_inverse_abn1_sum_3d(:,:,floor(s3_ss/2)),s1_ss,s2_ss);
figure
surf(A_inverse_2d);
figure
plot(A_inverse_abn1_sum_3d(20:100,90,floor(s3_ss/2)));
j = 2;
for i = floor(s3_ss/2)-25:floor(s3_ss/2)+25
    figure(j)
    plot(A_inverse_abn1_sum_3d(20:100,90,i));
    j = j + 1;
end


% Plot cube
% cm = colormap; % return the current color map
% colorID = zeros(m,1);
% myColor = zeros(m,3);
% for i = 1:m
%     colorID(i,1) = max(1,sum( A_inverse_abn1_sum(i,1)> (0:1/length(cm(:,1)):1)));
%     myColor(i,:) = cm(colorID(i,1),:);
% end
% myColor_3d = convert_1d_to_3d(myColor,s1_ss,s2_ss,s3_ss);
% Plot cube use plotcube,but can not assign color to the voxel
figure
for x = 1:10
    for y = 1:10
        for z = 1:10
            CubeLength = 1;
            plotcube([CubeLength CubeLength CubeLength],[x,y,z],0.9,[0 0 1]);
        end
    end
end













