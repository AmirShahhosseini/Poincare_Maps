function [Data_Centroids_Local] = Curvature_Identifier2(Temp_Trajectory,i)

% The x vector represnts the x components of the trajectory positions while
% the y vector represnts the y components of the trajectory positions. E.G.
% p1 = (x(1),y(1)), p2 = (x(2),y(2)) and p3 = (x(3),y(3));

% The center of the curvature is obtained for each three points as to facilitate 
% the progress of the algorithm.
global Data_Centroids Position_Data Dimension
 
if i == 1
    if Dimension == 2
        Initialization_point = [2*sum(Temp_Trajectory,2)/(Dimension+1); 3];  
    else
        Initialization_point = 2*sum(Temp_Trajectory,2)/Dimension;
    end
else 
    if Dimension == 2
        Initialization_point = [Data_Centroids(:,i-1); 3]; 
    else
        Initialization_point = Data_Centroids(:,i-1);
    end
end
 

x0 = Initialization_point;


Position_Data = Temp_Trajectory;


[U,resnorm] = lsqnonlin(@Fitting_Function2,x0);

if Dimension == 2
    Data_Centroids_Local = [U(1) U(2)];
else
    Data_Centroids_Local = U;
end
end

