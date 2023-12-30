%% Poincare Map 

% Project: A novel geometrical approach toward obtaining Poincare sections
% The Ohio State University

% All rights reserved. To cite this code please refer to the paper:

% Shahhosseini, A., Tien, MH. & D’Souza, K. 
% Poincare maps: a modern systematic approach toward obtaining effective sections. 
% Nonlinear Dyn 111, 529–548 (2023). 
% https://doi.org/10.1007/s11071-022-07864-y

% In case of any issues with running the code, please feel free to 
% contact the authors

%% Clearing MATLAB

clc;clear;close all

%% How to use?

% To use this code, provide the code with an n x m matrix where 
% n is the dimension of the system and represents the number of 
% variables and m is the data point of each variable. For example,
% for the case of Lorenz system, n = 3 and m can be any number based
% on the simulation time and time step of the simulation. For a 
% smooth performance, at least 5000 points for each dimension is 
% needed. An increase in the number of dimensions necessitates 
% a larger sample for adequate cluster analysis.

%% Important Note

% Please make sure to fully read the "Limitation" section of the paper to
% avoid numerical issues in running this code. It is critically important
% to provide "data of adequate quality" to this algorithm since its whole foundation is
% built upon geometry extraction from raw data and if the raw data do
% not meet its lenient criteria, the corresponding results, specially for the case
% of high-dimensional systems will be poor.


%% Trajectory Data Aacquisition

% Provide a n x m matrix here that represents the trajectory of your
% dynamics. For the case of the Lorenz attractor, the matrix consists of
% the variables, x-y-z and is provided below.

% 
load('Lorenz.mat')
X = [x;y;z];


global Dimension;

[n,m] = size(X);

Trajectory_Data = X;
if n > m % makes sure that your data are in the correct format!
    Trajectory_Data = transpose(Trajectory_Data);
end

[n,m] = size(Trajectory_Data);

Dimension = n;
Data_length = m;

Trajectory_Numerical = double(Trajectory_Data);


%% Algorithms Intrinsic Parameters 

% DO NOT CHANGE UNLESS YOU ABSOLUTELY! KNOW WHAT YOU ARE DOING!

Offset_Inclusion_Parameter = 0.3;              % The limits of the consideration of COC data relative to the trajectory data.
Maximum_Number_of_CoC_Clusters = 2;            % The maximum number of rotary flows allowed. Keep small or otherwise, computational costs will become high.
Sphere_To_Partition_Radius_Ratio = 0.3;        % The radius of the data in the vicinity of secondary centriods to determine the directionality vector
Maximum_Cluster_In_Partition = 3;              % The maximum number of clusters allowed in each partition. Keep small or otherwise, computational costs will become high.
Meshing_Quality_for_GR = 0.4;                  % The increment for plotting the Poincare sections
Coefficient_of_Contraction = 1;                % Coefficient that contracts COC data to let them be within the bounds of trajectory

%% Finding the Centers of Curvatures (CoC) 

global Data_Centroids

Data_Centroids = zeros(Dimension,length(Trajectory_Numerical)-Dimension);

for i = 1:length(Trajectory_Numerical)-Dimension-1
    i
    if Dimension == 2
        Temp_Trajectory = Trajectory_Numerical(:,i:i+Dimension);
    else
        Temp_Trajectory = Trajectory_Numerical(:,i:i+Dimension-1);     
    end
    [Data_Centroids_Local] = Curvature_Identifier2(Temp_Trajectory,i);
    Data_Centroids(:,i) = Data_Centroids_Local;
end

Data_Centroids = Data_Centroids*Coefficient_of_Contraction;
Data_Centroids = double(Data_Centroids);


%% Applying Bounds to the Obtained Data

counter_COCs = 1;

for j = 1:Data_length-Dimension
    counter_variable = 0;
    for i = 1:Dimension
        if Data_Centroids(i,j) >= (min(Trajectory_Numerical(i,:)) - Offset_Inclusion_Parameter*(max(Trajectory_Numerical(i,:)) - min(Trajectory_Numerical(i,:)))) && Data_Centroids(i,j) <= (max(Trajectory_Numerical(i,:)) + Offset_Inclusion_Parameter*(max(Trajectory_Numerical(i,:)) - min(Trajectory_Numerical(i,:))))
            counter_variable = counter_variable + 1;
        end
    end
    if counter_variable == Dimension
        % It is within the bounds
        Refined_COCs(:,counter_COCs) = Data_Centroids(:,j);
        counter_COCs = counter_COCs + 1;
    end
end


%% CoC Clustering Process

X = Refined_COCs';

% idx is the index for the cluster number associated to each data
% C is the location of the centroids.

%%%%%% Obtaining the optimal number of clusters

eva = evalclusters(X,'kmeans','silhouette','KList',[1:Maximum_Number_of_CoC_Clusters]);


%%%%%%%%%%%%%%%%%%%%    Clustering      %%%%%%%%%%%%%%%%%%%%%%%

opts = statset('Display','final');
number_of_clusters = eva.OptimalK;
[idx,C] = kmeans(X,number_of_clusters,'Distance','cityblock','Replicates',5,'Options',opts);

%% Partitioning Process

Centroid_Data = zeros(number_of_clusters,Dimension);
Centroid_Delta = zeros(number_of_clusters-1,Dimension);

for i = 1:number_of_clusters
   Centroid_Data(i,:) = C(i,:);
end

for i = 1:number_of_clusters-1
    for j = 1:Dimension
        Centroid_Delta(i,j) = abs(Centroid_Data(i+1,j) - Centroid_Data(i,j));
    end
end

if number_of_clusters > 2
    Partition_Decision_Variable = max(Centroid_Delta)./min(Centroid_Delta);
else
    Partition_Decision_Variable = 1./Centroid_Delta;
end


% The orientation of the partitions can then be easily identified!

Partition = Partition_Decision_Variable;

flag = zeros(1,Dimension);

index_orientation = find(min(Partition) == Partition);

flag(index_orientation) = 1;

%%%%%%%%%% Calculating the partitioning hyperplanes

Partitioning_Hyperplanes = C(:,index_orientation);

% Please have in mind that the "index_orientation" determines the
% orientation of the hyperplanes where "Partitioning_Hyperplanes"
% determines the value. For example if "index_orientation"= 1 then the
% direction of the hyperplanes are parallel to the x-axis and their
% corresponding values are obtained as "Partitioning_Hyperplanes"=c

% Now that the partitioning hyperplanes and their orientation are obtained,
% it is possible to partition the space into several subspaces and begin
% the internal clusterings

%% Finding Data Centroids of Each Partitioned Data

% Sorting the trajectory data based on the orientation of the hyperplane

Position = Trajectory_Numerical';
sorted_Position = sortrows(Position,index_orientation);
sorted_Coc_Position = sortrows(C,index_orientation);

sorted_Position_crossing = sorted_Position(:,index_orientation);
sorted_Coc_Position_crossing = sorted_Coc_Position(:,index_orientation);

% Now a search is required to find the crossings of the trajectory data
% with the partitioning hyperplanes
Index_crossing = zeros(1,number_of_clusters);

for i = 1:number_of_clusters
    Index_crossing(i) = find(sorted_Position_crossing>sorted_Coc_Position_crossing(i),1);
end

% Now, the inner data must be partitioned

for i = 1:number_of_clusters+1
    if i == 1
        
        % Obtaining the Data in the First Zone
        Zone_Data = sorted_Position(1:Index_crossing(i),:);
        
        %Clustering the Data in the First Zone
        temp_eva = evalclusters(Zone_Data,'kmeans','silhouette','KList',[1:Maximum_Cluster_In_Partition]);
        opts = statset('Display','final');
        number_of_local_clusters = temp_eva.OptimalK;
        [idx_temp,Centroid_Location] = kmeans(Zone_Data,number_of_local_clusters,'Distance','cityblock','Replicates',5,'Options',opts);
        Mesh_Centroids = Centroid_Location;
        
    elseif i > 1 && i~=(number_of_clusters+1)
        
        % Obtaining the Data in the ith Zone
        Zone_Data = sorted_Position(Index_crossing(i-1):Index_crossing(i),:);
        
        %Clustering the Data in the ith Zone
        temp_eva = evalclusters(Zone_Data,'kmeans','silhouette','KList',[1:Maximum_Cluster_In_Partition]);
        opts = statset('Display','final');
        number_of_local_clusters = temp_eva.OptimalK;
        [idx_temp,Centroid_Location] = kmeans(Zone_Data,number_of_local_clusters,'Distance','cityblock','Replicates',5,'Options',opts);
        Mesh_Centroids = [Mesh_Centroids;Centroid_Location];
        
    else
        Zone_Data = sorted_Position(Index_crossing(i-1):end,:);
        
        temp_eva = evalclusters(Zone_Data,'kmeans','silhouette','KList',[1:Maximum_Cluster_In_Partition]);
        opts = statset('Display','final');
        number_of_local_clusters = temp_eva.OptimalK;
        [idx_temp,Centroid_Location] = kmeans(Zone_Data,number_of_local_clusters,'Distance','cityblock','Replicates',5,'Options',opts);
        Mesh_Centroids = [Mesh_Centroids;Centroid_Location];
    end
end

%% Spherical Data Sampling and Trajectory Orientation Identification
global m1
[m,n] = size(Mesh_Centroids);
[m1,n1] = size(C);

Sorted_C = sortrows(C,index_orientation);
Sorted_Mesh_Centroids = sortrows(Mesh_Centroids,index_orientation);
Index_Bound = [1, Index_crossing, length(Position)];
Spherical_Data_Storage = zeros(length(Position),Dimension*m);
global Vector_Data
Vector_Data = zeros(length(Position),Dimension*m);
Threshold = zeros(1,m1+1);


for i = 1:m1+1
    Threshold(i) = Sphere_To_Partition_Radius_Ratio*(sorted_Position(Index_Bound(i+1),index_orientation) - sorted_Position(Index_Bound(i),index_orientation));
end


for i = 1:m
    
    temp = find(Sorted_Mesh_Centroids(i,index_orientation) > Sorted_C(:,index_orientation));
    if isempty(temp)
        Partition_Number = 1;
    else
        Partition_Number = max(temp)+1;
    end

    for j = 1:length(Position)
        if sqrt(sum((Sorted_Mesh_Centroids(i,:) - Position(j,:)).^2)) < Threshold(Partition_Number)
            Spherical_Data_Storage(j,Dimension*i-(Dimension-1):Dimension*i) = Position(j,:);
        end
    end
end

for i = 1:m
    for j = 1:length(Position)-1 
        temp1 = Spherical_Data_Storage(j+1,Dimension*i-(Dimension-1):Dimension*i);
        temp2 = Spherical_Data_Storage(j,Dimension*i-(Dimension-1):Dimension*i);
        if (sum(abs(temp1)) ~= 0 || sum(abs(temp2)) ~= 0) && (sum(temp1 ~= temp2) ~= 0)
            Vector_Data(j,Dimension*i-(Dimension-1):Dimension*i) = Spherical_Data_Storage(j+1,Dimension*i-(Dimension-1):Dimension*i) - Spherical_Data_Storage(j,Dimension*i-(Dimension-1):Dimension*i);
            Vector_Data(j,Dimension*i-(Dimension-1):Dimension*i) = Vector_Data(j,Dimension*i-(Dimension-1):Dimension*i)/(sqrt(sum(Vector_Data(j,Dimension*i-(Dimension-1):Dimension*i).*Vector_Data(j,Dimension*i-(Dimension-1):Dimension*i))));
        end
    end
end
global oss Average_Directionality_Vector SphereCenters Temp_Partition_Number AttractorCenteres
global Real_Directionality_Vector

SphereCenters = Sorted_Mesh_Centroids;
AttractorCenteres = Sorted_C;
Average_Directionality_Vector = zeros(m,Dimension);

Real_Directionality_Vector = zeros(m,Dimension+1);
x0 = 0.1 + zeros(1,Dimension);

for oss = 1:m
    [U,resnorm] = lsqnonlin(@Maximization_Function,x0);
    Average_Directionality_Vector(oss,:) = U;
end

x0 = [Average_Directionality_Vector(1,:) 0.01];

for oss = 1:m
    i = oss;
    temp = find(Sorted_Mesh_Centroids(i,index_orientation) > Sorted_C(:,index_orientation));
    if isempty(temp)
        Temp_Partition_Number = 1;
    else
        Temp_Partition_Number = max(temp)+1;
    end
    [U,resnorm] = lsqnonlin(@Optimization_Function,x0);
    Real_Directionality_Vector(oss,:) = U;
end

%% Solver to Draw Poincare Map

% The idea here is to find the crossings between the poincare sections and
% the trajectory in their partition of interest so the
% idea is to obtain m poincare maps for the trajectory of interest and
% attempt to identify the regime of motion!


global Poincare_Map_Data jo io
Poincare_Map_Data = zeros(length(Position),Dimension*m);
Poincare_Map_Data_Planar = zeros(length(Position),Dimension*m);

for io = 1:m % To identify the crossings of all m sections!
    sign_plane_now = 0;
    sign_plane_past = 0;
    % The first part is to identify the partition of interest to be able to
    % determine the subspace of interest!
    temp = find(Sorted_Mesh_Centroids(io,index_orientation) > Sorted_C(:,index_orientation));
    if isempty(temp)
        Partition_Number = 1;
    else
        Partition_Number = max(temp)+1;
    end
    
    % The second step is to identify the crossing points (to the highest
    % accuracy). After this detection, the corresponding data will be
    % stored in the Poincare_map data!

    if Partition_Number == 1
        % Do nothing 
    elseif Partition_Number == (m1+1)
        Partition_Number = Partition_Number-1;
    else
        distance1 = sqrt( sum (  (SphereCenters(io,1:Dimension) - AttractorCenteres(Partition_Number,1:Dimension)).^2));
        distance2 = sqrt( sum (  (SphereCenters(io,1:Dimension) - AttractorCenteres(Partition_Number-1,1:Dimension)).^2));
    
        if distance1 > distance2
            Partition_Number = Partition_Number-1;
        end
    end

    for jo = 1:length(Position)
        sign_plane_past = sign_plane_now;
        if Position(jo,index_orientation) < Sorted_C(Partition_Number,index_orientation)
            % The above condition means that the point is in the desired
            % partition!
            sign_plane_now = sign(sum(Real_Directionality_Vector(io,1:Dimension).*Position(jo,:)) - Real_Directionality_Vector(io,Dimension+1));
            if sign_plane_now*sign_plane_past == -1
                
                 x0 = Position(jo,1:Dimension);
                [U,resnorm] = lsqnonlin(@Interpolation_Function,x0);
                Poincare_Map_Data(jo,Dimension*io-(Dimension-1):Dimension*io) = U(1:Dimension);
            end
        end
    end
    
    % The last step is to rotate the data to be in their corresponding
    % section coordinate to have a n-1 dimensional mapping for the nD system
    % of interest!
    
    wh = Real_Directionality_Vector(io,1:Dimension);
    w2 = null(wh);
    for lll = 1:Dimension-1
        rotation_column(:,lll) = w2(:,lll);
    end

    R = [rotation_column';wh];

    Poincare_Map_Data_Planar(:,Dimension*io-(Dimension-1):Dimension*io) = (R*Poincare_Map_Data(:,Dimension*io-(Dimension-1):Dimension*io)')';
    
    % The last part is to graphically illustrate the Poincare map data and
    % observe the existence of chaos or periodic solutions!
        
end



%% Graphical Demonstration


if Dimension == 2
    

figure();
plot(Trajectory_Data(1,:),Trajectory_Data(2,:));
hold on
plot(Sorted_C(:,1),Sorted_C(:,2),'o','linewidth',3)
plot(Sorted_Mesh_Centroids(:,1),Sorted_Mesh_Centroids(:,2),'x','linewidth',7)
xlim([min(Trajectory_Data(1,:)) max(Trajectory_Data(1,:))])
ylim([min(Trajectory_Data(2,:)) max(Trajectory_Data(2,:))])

for i = 1:m
    temp = find(Sorted_Mesh_Centroids(i,index_orientation) > Sorted_C(:,index_orientation));
    if isempty(temp)
        Partition_Number = 1;
    else
        Partition_Number = max(temp)+1;
    end
    
    if Partition_Number == 1
        X =  min(Trajectory_Data(1,:)):Meshing_Quality_for_GR:Sorted_C(1,1);
    elseif Partition_Number == (m1+1)
        X = Sorted_C(end,1):Meshing_Quality_for_GR:max(Trajectory_Data(1,:));
    else
        X = Sorted_C(Partition_Number-1,1):Meshing_Quality_for_GR:Sorted_C(Partition_Number,1);
    end
        Y = (1/(Real_Directionality_Vector(i,2)))*(Real_Directionality_Vector(i,3)-Real_Directionality_Vector(i,1)*X);
    plot(X,Y)
end



elseif Dimension == 3
    % 3D Trajectory and 2D Poincare Map
    

figure();
plot3(Trajectory_Data(1,:),Trajectory_Data(2,:),Trajectory_Data(3,:));
hold on
plot3(Sorted_C(:,1),Sorted_C(:,2),Sorted_C(:,3),'o','linewidth',3)
plot3(Sorted_Mesh_Centroids(:,1),Sorted_Mesh_Centroids(:,2),Sorted_Mesh_Centroids(:,3),'x','linewidth',7)
xlim([min(Trajectory_Data(1,:)) max(Trajectory_Data(1,:))])
ylim([min(Trajectory_Data(2,:)) max(Trajectory_Data(2,:))])
zlim([min(Trajectory_Data(3,:)) max(Trajectory_Data(3,:))])

for i = 1:m
    temp = find(Sorted_Mesh_Centroids(i,index_orientation) > Sorted_C(:,index_orientation));
    if isempty(temp)
        Partition_Number = 1;
    else
        Partition_Number = max(temp)+1;
    end
    
    if index_orientation == 2
        if Partition_Number == 1
            [X,Y]=meshgrid(min(Trajectory_Data(1,:)):Meshing_Quality_for_GR:max(Trajectory_Data(1,:)),min(Trajectory_Data(2,:)):Meshing_Quality_for_GR:Sorted_C(1,2));
        elseif Partition_Number == (m1+1)
            [X,Y]=meshgrid(min(Trajectory_Data(1,:)):Meshing_Quality_for_GR:max(Trajectory_Data(1,:)),Sorted_C(end,2):Meshing_Quality_for_GR:max(Trajectory_Data(2,:)));
        else
            [X,Y]=meshgrid(min(Trajectory_Data(1,:)):Meshing_Quality_for_GR:max(Trajectory_Data(1,:)),Sorted_C(Partition_Number-1,2):Meshing_Quality_for_GR:Sorted_C(Partition_Number,2));
        end
        Z = (1/(Real_Directionality_Vector(i,3)))*(Real_Directionality_Vector(i,4)-Real_Directionality_Vector(i,1)*X-Real_Directionality_Vector(i,2)*Y);
    elseif index_orientation == 1
        if Partition_Number == 1
            [X,Y]=meshgrid(min(Trajectory_Data(1,:)):Meshing_Quality_for_GR:Sorted_C(1,1),min(Trajectory_Data(2,:)):Meshing_Quality_for_GR:max(Trajectory_Data(2,:)));
        elseif Partition_Number == (m1+1)
            [X,Y]=meshgrid(Sorted_C(end,1):Meshing_Quality_for_GR:max(Trajectory_Data(1,:)),min(Trajectory_Data(2,:)):Meshing_Quality_for_GR:max(Trajectory_Data(2,:)));
        else
            [X,Y]=meshgrid(Sorted_C(Partition_Number-1,1):Meshing_Quality_for_GR:Sorted_C(Partition_Number,1),min(Trajectory_Data(2,:)):Meshing_Quality_for_GR:max(Trajectory_Data(2,:)));
        end
        Z = (1/(Real_Directionality_Vector(i,3)))*(Real_Directionality_Vector(i,4)-Real_Directionality_Vector(i,1)*X-Real_Directionality_Vector(i,2)*Y);
    else
        if Partition_Number == 1
            [Y,Z]=meshgrid(min(Trajectory_Data(2,:)):Meshing_Quality_for_GR:max(Trajectory_Data(2,:)),min(Trajectory_Data(3,:)):Meshing_Quality_for_GR:Sorted_C(1,3));
        elseif Partition_Number == (m1+1)
            [Y,Z]=meshgrid(min(y):Meshing_Quality_for_GR:max(y),Sorted_C(end,3):Meshing_Quality_for_GR:max(Trajectory_Data(3,:)));
        else
            [Y,Z]=meshgrid(min(Trajectory_Data(2,:)):Meshing_Quality_for_GR:max(Trajectory_Data(2,:)),Sorted_C(Partition_Number-1,3):Meshing_Quality_for_GR:Sorted_C(Partition_Number,3));
        end
        X = (1/(Real_Directionality_Vector(i,1)))*(Real_Directionality_Vector(i,4)-Real_Directionality_Vector(i,3)*Z-Real_Directionality_Vector(i,2)*Y);   
    end
    
    s = surf(X,Y,Z,'FaceColor',[0.1*rand() 0.1*rand() 0.1*rand()]);
    s.EdgeColor = 'none';
end
    
    
A = zeros(3,3*m);

temp_Poincare_Map_Data_Planar = Poincare_Map_Data_Planar;

for i = 1:m
    counter = 1;
    for j = length(temp_Poincare_Map_Data_Planar):-1:1
        if (sum(temp_Poincare_Map_Data_Planar(j,3*i-2:3*i) == [0 0 0])) ~= 3
            A(counter,3*i-2:3*i) = temp_Poincare_Map_Data_Planar(j,3*i-2:3*i);
            counter = counter +1;
            temp_Poincare_Map_Data_Planar(j,3*i-2:3*i) = [1e5 1e5 1e5];
            if counter == 4
                break
            end
        end
    end
end

figure();
for io = 1:m
    subplot(ceil(m/2),2,io)
    plot(Poincare_Map_Data_Planar(:,3*io-2),Poincare_Map_Data_Planar(:,3*io-1),'o');
    hold on
    plot(A(:,3*io-2),A(:,3*io-1),'o')
end
    
    
    
    
    
    
    
    
    
elseif Dimension == 4
    % 4D Trajectory and 3D Poincare Map
    
A = zeros(3,4*m);

temp_Poincare_Map_Data_Planar = Poincare_Map_Data_Planar;

for i = 1:m
    counter = 1;
    for j = length(temp_Poincare_Map_Data_Planar):-1:1
        if (sum(temp_Poincare_Map_Data_Planar(j,4*i-3:4*i) == [0 0 0 0])) ~= 4
            A(counter,4*i-3:4*i) = temp_Poincare_Map_Data_Planar(j,4*i-3:4*i);
            counter = counter +1;
            temp_Poincare_Map_Data_Planar(j,4*i-3:4*i) = [1e5 1e5 1e5 1e5];
            if counter == 4
                break
            end
        end
    end
end

figure();
for io = 1:m
    subplot(ceil(m/2),2,io)
    plot3(Poincare_Map_Data_Planar(:,4*io-3),Poincare_Map_Data_Planar(:,4*io-2),Poincare_Map_Data_Planar(:,4*io-1),'o');
    hold on
    plot3(A(:,4*io-3),A(:,4*io-2),A(:,4*io-1),'o')
end
    
    
elseif Dimension > 4
    fprintf('\n\n\n Systems with a dimension higher than 4 cannot be graphically illustrated!')
end