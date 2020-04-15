function Area = integ(Y,X)
% Y è vettore con il valore della funzione in funzione del vettore X -->
% Y=f(X)
    for ii=1:199
        An(ii) = ( (Y(ii) + Y(ii+1)) * (X(ii+1)-X(ii)) )/2;
    end
    Area = sum(An);
end