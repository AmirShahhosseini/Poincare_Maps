function F = Maximization_Function(U)

global Vector_Data oss Dimension

F(1) =  abs(sum(U.^2) - 1);
Sum1 = 0;

for i = 1:length(Vector_Data)
    Sum1 = Sum1 + abs(sum(U.*Vector_Data(i,Dimension*oss-(Dimension-1):Dimension*oss)));
end

F(2) = 1/(1+(Sum1)^2);

end
