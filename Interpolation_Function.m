function F = Interpolation_Function(U)

global Poincare_Map_Data jo Real_Directionality_Vector io Dimension

A_components = U - Poincare_Map_Data(jo,Dimension*io-(Dimension-1):Dimension*io);
A_magnitude = sqrt( sum( A_components.^2));
A = A_components/A_magnitude;

B_components = -U + Poincare_Map_Data(jo-1,Dimension*io-(Dimension-1):Dimension*io);
B_magnitude = sqrt( sum( B_components.^2));
B = B_components/B_magnitude;

C = sum(A.*B);
F(1) = abs(C-1);

F(2) = abs(sum(Real_Directionality_Vector(io,1:Dimension).*U(1:Dimension)) - Real_Directionality_Vector(io,Dimension+1));

end
