clear all
clc

% Función de Wigner en (0,0)

load MUBRAYOS.mat
load MUBSK.mat

S0 = [1;0]; % basis state |0>
S1 = [0;1]; % basis state |1>

%E = kron(S0, kron(S0, kron(S0,kron(S0,S0))));
%E1 = kron(S1, kron(S1, kron(S1,kron(S1,S1))));

E1 = MUBK(32*2+1:(32*3),1);
%kron(S0,kron(S1,kron(S1,kron(S1,S1))));
RoE = E1*ctranspose(E1);%+E*ctranspose(E))/(2);

W1 = zeros(32);
W2 = W1;

for i1=1:32
    for i2=1:32
        [ NR, NK ] = Nucleo(i1-1,i2-1, MUB, MUBK);
        WR = trace(NR*RoE)/32;
        WK = trace(NK*RoE)/32;
        W1(i1,i2) = WR;
        W2(i1,i2) = WK;
    end
end

figure(1)
Z = W1;
b = bar3(Z);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%title('Wigner \mid\rightarrow\rightarrow\rightarrow\rightarrow\rightarrow\rangle')
title('Wigner Est�ndar para un eigen-estado de MUB- Est�ndar')
xlabel('\alpha');
ylabel('\beta');
colormap(jet)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
Z = W2;
b = bar3(Z);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%title('Wigner \mid\rightarrow\rightarrow\rightarrow\rightarrow\rightarrow\rangle')
title('Wigner Kantor para un eigen-estado de MUB-Est�ndar')
xlabel('\alpha');
ylabel('\beta');
colormap(jet)   
