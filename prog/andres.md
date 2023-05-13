### Original Andrés code

R<w> := GaloisRing(2^2, 5);
R;

x := [* Zero(R), One(R) *];

F := [* *];
ya := F;

for i1 := 1 to 30 by 1 do
    x := Append(x, (R.1^2 + 2*R.1 + 1)^i1);
end for;

Tra := function(a,b);
    T := 0;
    for k in [0,1,2,3,4] do
       T := T + (a)^(2^k) + 2*(b^(2^k));
    end for;
    return T;
end function;

for i1 := 1 to 32 by 1 do
    for i2 := 1 to 32 by 1 do
        for i3 := 1 to 32 by 1 do  
            y1:=x[i3]*((x[i1]^2)*x[i3]+x[i1]*Tra(x[i3],x[1])+Tra(x[i1]*x[i3],x[1])+2*x[i2]);
            F := Append(F, y1);
        end for;
    end for;
end for;

for i1 := 1 to 32*32*32 by 1 do
    for i2 := 1 to 32 by 1 do
        for i3 := 1 to 32 by 1 do  
            if F[i1] eq (x[i2]+2*x[i3]) then
                ya:=Append(ya,Tra(x[i2],x[i3]));
                break i2;
            end if;
        end for;
    end for;
end for;

kk := 0;
for i2 := 1 to 32 by 1 do 
   for i3 := 1 to 32 by 1 do  
        print "\n", ya[i3+kk], ya[i3+kk+32],ya[i3+kk+32*2],ya[i3+kk+32*3],ya[i3+kk+32*4],ya[i3+kk+32*5],ya[i3+kk+32*6],ya[i3+kk+32*7],ya[i3+kk+32*8],ya[i3+kk+32*9],ya[i3+kk+32*10],ya[i3+kk+32*11],ya[i3+kk+32*12],ya[i3+kk+32*13],ya[i3+kk+32*14],ya[i3+kk+32*15],ya[i3+kk+32*16],ya[i3+kk+32*17],ya[i3+kk+32*18],ya[i3+kk+32*19],ya[i3+kk+32*20],ya[i3+kk+32*21],ya[i3+kk+32*22],ya[i3+kk+32*23],ya[i3+kk+32*24],ya[i3+kk+32*25],ya[i3+kk+32*26],ya[i3+kk+32*27],ya[i3+kk+32*28],ya[i3+kk+32*29],ya[i3+kk+32*30],ya[i3+kk+32*31];
   end for;
   kk := 32*32*i2;
end for;

---
### Cleaned up 2 qubit version

R<w> := GaloisRing(2^2, 2);
R;

Teic := [* Zero(R), One(R) *];

for i1 := 1 to 2 by 1 do
    Teic := Append(Teic, (R.1^2 + 2*R.1 + 1)^i1);
end for;
Teic;

F := [* *];
ya := F;

Tra := function(a,b);
    T := 0;
    for k in [0,1] do
       T := T + (a)^(2^k) + 2*(b^(2^k));
    end for;
    return T;
end function;

for i1 := 1 to 4 by 1 do
    for i2 := 1 to 4 by 1 do
        for i3 := 1 to 4 by 1 do  
            y1 := Teic[i3]^2*Teic[i1] + 2*Teic[i3]*Teic[i2];
            F := Append(F, y1);
        end for;
    end for;
end for;

for i1 := 1 to 4*4*4 by 1 do
    for i2 := 1 to 4 by 1 do
        for i3 := 1 to 4 by 1 do  
            if F[i1] eq (Teic[i2] + 2*Teic[i3]) then
                ya := Append(ya, Tra(Teic[i2], Teic[i3]));
                break i2;
            end if;
        end for;
    end for;
end for;

kk := 0;
for i2 := 1 to 4 by 1 do 
   for i3 := 1 to 4 by 1 do  
        print "\n", ya[i3+kk], ya[i3+4+32], ya[i3+kk+4*2], ya[i3+kk+4*3];
   end for;
   kk := 4 * 4 * i2;
end for;

---
### Cleaned up 5 qubit version



