function F = Optimization_Function(U)

global Average_Directionality_Vector oss SphereCenters AttractorCenteres Temp_Partition_Number m1 Dimension

F(1) =  abs(sum(U(1:Dimension).^2) - 1);

F(2) = 5*sum((U(1:Dimension) - Average_Directionality_Vector(oss,1:Dimension)).^2);

F(3) = 15*((sum(U(1:Dimension).*SphereCenters(oss,1:Dimension)) - U(Dimension+1))^2);

if Temp_Partition_Number == 1
    % Do nothing
    
elseif Temp_Partition_Number == (m1+1)
    Temp_Partition_Number = Temp_Partition_Number-1;
else
    distance1 = sqrt( sum (  (SphereCenters(oss,1:Dimension) - AttractorCenteres(Temp_Partition_Number,1:Dimension)).^2));
    distance2 = sqrt( sum (  (SphereCenters(oss,1:Dimension) - AttractorCenteres(Temp_Partition_Number-1,1:Dimension)).^2));
    if distance1 > distance2
        Temp_Partition_Number = Temp_Partition_Number -1;
    end
end

F(4) = (sum(U(1:Dimension).*AttractorCenteres(Temp_Partition_Number,:)) - U(Dimension+1))^2;
end
