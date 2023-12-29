function [ NR,NK ] = Nucleo(x,y,MUB,MUBK)
   %Función que calcula el núcleo de Wigner en el punto x,y

   NR = zeros(32);
   NK = NR;
   
   m = 0:31;
   n = 5;
   m = gf(m,n);
   x1 = gf(x,5,37);
   T = x1+x1.^2+x1.^4+x1.^8+x1.^16;
   y1 = gf(y,5,37);
   
   Gam = y1+m*x1;
   GamK = y1+(m.^2).*x1+m.*T+m.*x1+(m.*x1).^2+(m.*x1).^4+(m.*x1).^8+(m.*x1).^16;
   G = gf2dec(Gam,5,37);
   GK = gf2dec(GamK,5,37);
   C = 1;

   for k1=1:32:1025
      if k1 ~= 33
           NR = NR+MUB(k1:k1+31,G(C)+1)*ctranspose(MUB(k1:k1+31,G(C)+1));
           NK = NK+MUBK(k1:k1+31,GK(C)+1)*ctranspose(MUBK(k1:k1+31,GK(C)+1));
           C = C+1;
      else
           NR = NR+MUB(k1:k1+31,x+1)*transpose(MUB(k1:k1+31,x+1));
           NK = NK+MUBK(k1:k1+31,x+1)*transpose(MUBK(k1:k1+31,x+1));
       end
    end
   NR = NR-eye(32);
   NK = NK-eye(32);
end
