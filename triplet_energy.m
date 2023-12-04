clear all; close all; clc;

% Constants
Ec = 1;
Ej = 15 * Ec;
lambda = 2;

% Vars
n = 10;
N = -n:n;
L = length(N);

ng = 2;
Ng = -ng:0.01:ng;
l = length(Ng);

E = zeros(L, l);

for j = 1:l
        H = zeros(L);
    
        for i = 1:L
            H(i, i) = 4 * Ec * (N(i) - Ng(j))^2;
        
            if (i ~= L)
                H(i + 1, i) = Ej;
                H(i, i + 1) = Ej;
            end
    
            if (i + 2 < L)
                C = Ej / 2 / (1 + lambda^2)^2;
                H(i + 2, i) = C;
                H(i, i + 2) = C;
            end
        end
        
        e = eig(H);
        E(:, j) = e;
end

f = 5;
filter = E(1:f, :);

figure(1);
hold on;
for i = 1:f
    plot(Ng, filter(i, :) /4 / Ec, 'DisplayName', string(N(i) + n));
end

ti = ['Ej = ' num2str(Ej)  'Ec,  lambda = ' num2str(lambda)];
lgd = legend;
title(lgd, 'N');
xlabel('Ng (charge offset)');
ylabel('E');
title(ti);


