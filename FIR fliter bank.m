function [H] = FIR(M)

b = fir1(62,1/M,'low');
a = b;
for i = 1:M-2
    c = fir1(62,[i/M,(i+1)/M],'bandpass');
    a = [a; c];
end
d = fir1(62,(M-1)/M,'high');

a = [a; d];

H = a;