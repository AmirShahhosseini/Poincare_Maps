function F = Fitting_Function2(U)

global Position_Data Dimension


if Dimension == 2

    F(1) = abs( (Position_Data(1,1)-U(1))^2 + (Position_Data(2,1)-U(2))^2 - (U(3))^2);
    F(2) = abs( (Position_Data(1,2)-U(1))^2 + (Position_Data(2,2)-U(2))^2 - (U(3))^2);
    F(3) = abs( (Position_Data(1,3)-U(1))^2 + (Position_Data(2,3)-U(2))^2 - (U(3))^2);
    
else


    vector = zeros(Dimension-1,Dimension);

    for ii = 1:Dimension-1
        vector(ii,:) = Position_Data(:,ii+1)-Position_Data(:,1);
    end

    CROSS_MATRIX = vector;

    n = zeros(1,length(CROSS_MATRIX));

    for ok =1:length(CROSS_MATRIX)
        temp = CROSS_MATRIX;
        temp(:,ok) = [];
        n(ok) = ((-1)^(ok+1))*det(temp);
    end

    mag_n = sqrt(sum(n.^2));
    if mag_n ==0
        % pass
    else
        n = n/mag_n;
    end

    d = sum(n.*Position_Data(:,1)');

    DU = zeros(1,Dimension);

    for kk = 1:Dimension
        S = 0;
        for jj = 1:length(Position_Data)
            S = S + ((Position_Data(jj,kk) - U(jj)).^2);
        end
        DU(kk) = sqrt(S);
    end

    temp = Dimension-1; 
    SSS = 1;
    counter = 1;
    for opo = 1:nchoosek(Dimension,2)
        F(opo) = abs(DU(SSS) - DU(SSS+counter));
        temp = temp - 1;
        counter = counter +1;
        if temp == 0
            SSS = SSS +1;
            counter = 1;
            temp = Dimension - SSS;
        end
    end
    F(opo+1) = abs(sum(n.*U') - d);
end
end

