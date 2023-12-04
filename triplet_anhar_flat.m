clear all; close all; clc;

% loops
ratios = 1:0.5:50;
ls =0:0.1:5;

R = length(ratios);
P = length(ls);

% Constants
Ec = 1;
Ejs = Ec * ratios;
eta = 0.5;

% Vars
n = 10;
N = -n:n;

ng = 5;
dg = 0.01;
Ng = -ng:dg:ng;

% Grids
flats = zeros(R, P);
anhars = flats;

for i =1:R
    Ej = Ejs(i);
    for j = 1:P
        lambda = ls(j);
        energies = triplet_hamiltonian(eta, Ec, Ej, lambda, N, Ng);

        E0 = energies(1, :);
        E1 = energies(2, :);
        E2 = energies(3, :);

        w1 = max(E1) - min(E1);

        E0 = mean(E0);
        E1 = mean(E1);
        E2 = mean(E2);

        flats(i, j) = flatness(E0, E1, w1);
        anhars(i, j) = anharmonicity(E0, E1, E2);
    end
end

figure(1);
imagesc(ratios, ls, anhars );
colorbar;
colormap("gray");
set(gca,'YDir','normal');
xlabel('$E_J/E_C$');
ylabel('$\lambda$');
set(1,'DefaultTextInterpreter', 'latex') 

figure(2);
imagesc(ratios, ls, flats);
colorbar;
colormap("bone");
set(gca,'YDir','normal');
xlabel('$E_J/E_C$');
ylabel('$\lambda$');
set(2,'DefaultTextInterpreter', 'latex') 



function energies = triplet_hamiltonian(eta, Ec, Ej, lambda, N, Ng)
    L = length(N);
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

    f = 3;
    energies = E(1:f, :);
end

function f = flatness(E0, E1, w1)
    w10 = E1 - E0;
    f = w1 / w10;
    f = -log(f/2);
end

function a = anharmonicity(E0, E1, E2)
    w21 = E2 - E1;
    w10 = E1 - E0;
    w20 = E2 - E0;
    a = (w21 - w10) / w20;
end

