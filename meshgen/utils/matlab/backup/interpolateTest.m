airx = pts(1:50,1);
airy = pts(1:50,2);

len = length(airx);

count = 0;
for i = 1:len-1
   delta = airx(i+1) - airx(i);
  
   count = count + 1;   
   xi(count) = airx(i);
   
   count = count + 1;   
   xi(count) = airx(i) + 0.25*delta;
   
   count = count + 1;   
   xi(count) = airx(i) + 0.50*delta;
   
   count = count + 1;   
   xi(count) = airx(i) + 0.75*delta;
   
   
end

xi = xi';

yi = pchip(airx,airy,xi);

xnaca = xi;
ynaca = -naca(xnaca);

close all;clc
figure(1)
plot(xi,yi,'r',airx,airy,'b.',xnaca,ynaca,'ks')
legend('PCHIP','DATA','NACA')