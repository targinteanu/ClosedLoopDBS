function yf = myFastForecastAR(arMdl, yPast, k)
% a faster version of forecast when arMdl is a no-input, discrete AR model
% with the same sampling rate as yPast, and yPast is a column vector. 

y = [yPast; nan(k,1)]; % [yPast; yf]

A = arMdl.A;
A1 = A(1); A = A(2:end);
N = length(A);

for t = (length(yPast)+1):length(y)
    yp = y((t-1):-1:(t-N));
    y(t) = (-A*yp)/A1;
end

yf = y((length(yPast)+1):end);

end