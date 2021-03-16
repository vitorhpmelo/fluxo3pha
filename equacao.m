clear all;
tetap=[0 -120 120];
tetas=[-0.339 -120.344 119.629];
Vp=[1 1 1]*12470/sqrt(3);
Vs=[0.98708 0.99169 0.98905]*12470/sqrt(3);
tetap=tetap*pi/180;
tetas=tetas*pi/180;

Vpcom=transpose(Vp.*exp(i*tetap));
Vscom=transpose(Vs.*exp(i*tetas));

Za=[0.4576+i*1.078,0.1559+i*0.5017,0.1535+i*0.3849];
Zb=[0,0.4666+i*1.0482,0.1580+i*0.4236];
Zc=[0,0,0.4615+i*1.0651];

Z=0.3788*[Za;Zb;Zc];
D=diag(Z);
Z=Z+transpose(Z)-diag(D);

Ypp=inv(Z);
Yss=inv(Z);
Yps=-inv(Z);
Ysp=-inv(Z);

Gpp=real(Ypp);
Gps=real(Yps);
Gss=real(Yss);
Gsp=real(Ysp);
Bpp=imag(Ypp);
Bps=imag(Yps);
Bss=imag(Yss);
Bsp=imag(Ysp);



%matriz impedancia 

Ips=Ypp*Vpcom+Yps*Vscom;

%% fluxo direto

Pps_a=(Vp(1)^2)*Gpp(1,1);

for i=2:3
  Pps_a=Pps_a+Vp(1)*Vp(i)*(Gpp(1,i)*cos(tetap(1)-tetap(i))+Bpp(1,i)*sin(tetap(1)-tetap(i)));
end

for i=1:3
  Pps_a=Pps_a+Vp(1)*Vs(i)*(Gps(1,i)*cos(tetap(1)-tetas(i))+Bps(1,i)*sin(tetap(1)-tetas(i)));
end

Pps_b=(Vp(2)^2)*Gpp(2,2);

for i=[1,3]
  Pps_b=Pps_b+Vp(2)*Vp(i)*(Gpp(2,i)*cos(tetap(2)-tetap(i))+Bpp(2,i)*sin(tetap(2)-tetap(i)));
end

for i=1:3
  Pps_b=Pps_b+Vp(2)*Vs(i)*(Gps(2,i)*cos(tetap(2)-tetas(i))+Bps(2,i)*sin(tetap(2)-tetas(i)));
end


Pps_c=(Vp(3)^2)*Gpp(3,3);

for i=[1,2]
  Pps_c=Pps_c+Vp(3)*Vp(i)*(Gpp(3,i)*cos(tetap(3)-tetap(i))+Bpp(3,i)*sin(tetap(3)-tetap(i)));
end

for i=1:3
  Pps_c=Pps_c+Vp(3)*Vs(i)*(Gps(3,i)*cos(tetap(3)-tetas(i))+Bps(3,i)*sin(tetap(3)-tetas(i)));
end


%% fluxo inverso

Psp_a=(Vs(1)^2)*Gss(1,1);

for i=2:3
  Psp_a=Psp_a+Vs(1)*Vs(i)*(Gss(1,i)*cos(tetas(1)-tetas(i))+Bss(1,i)*sin(tetas(1)-tetas(i)));
end

for i=1:3
  Psp_a=Psp_a+Vs(1)*Vp(i)*(Gsp(1,i)*cos(tetas(1)-tetap(i))+Bsp(1,i)*sin(tetas(1)-tetap(i)));
end

Psp_b=(Vs(2)^2)*Gss(2,2);

for i=[1,3]
  Psp_b=Psp_b+Vs(2)*Vs(i)*(Gss(2,i)*cos(tetas(2)-tetas(i))+Bss(2,i)*sin(tetas(2)-tetas(i)));
end

for i=1:3
  Psp_b=Psp_b+Vs(2)*Vp(i)*(Gsp(2,i)*cos(tetas(2)-tetap(i))+Bsp(2,i)*sin(tetas(2)-tetap(i)));
end


Psp_c=(Vs(3)^2)*Gss(3,3);

for i=[1,2]
  Psp_c=Psp_c+Vs(3)*Vs(i)*(Gss(3,i)*cos(tetas(3)-tetas(i))+Bss(3,i)*sin(tetas(3)-tetas(i)));
end

for i=1:3
  Psp_c=Psp_c+Vs(3)*Vp(i)*(Gsp(3,i)*cos(tetas(3)-tetap(i))+Bsp(3,i)*sin(tetas(3)-tetap(i)));
end


%% reativo

% fluxo direto

Qps_a=-(Vp(1)^2)*Bpp(1,1);

for i=2:3
  Qps_a=Qps_a-Vp(1)*Vp(i)*(Bpp(1,i)*cos(tetap(1)-tetap(i))-Gpp(1,i)*sin(tetap(1)-tetap(i)));
end

for i=1:3
  Qps_a=Qps_a-Vp(1)*Vs(i)*(Bps(1,i)*cos(tetap(1)-tetas(i))-Gps(1,i)*sin(tetap(1)-tetas(i)));
end

Qps_b=-(Vp(2)^2)*Bpp(2,2);

for i=[1,3]
  Qps_b=Qps_b-Vp(2)*Vp(i)*(Bpp(2,i)*cos(tetap(2)-tetap(i))-Gpp(2,i)*sin(tetap(2)-tetap(i)));
end

for i=1:3
  Qps_b=Qps_b-Vp(2)*Vs(i)*(Bps(2,i)*cos(tetap(2)-tetas(i))-Gps(2,i)*sin(tetap(2)-tetas(i)));
end


Qps_c=-(Vp(3)^2)*Bpp(3,3);

for i=[1,2]
  Qps_c=Qps_c-Vp(3)*Vp(i)*(Bpp(3,i)*cos(tetap(3)-tetap(i))-Gpp(3,i)*sin(tetap(3)-tetap(i)));
end

for i=1:3
  Qps_c=Qps_c-Vp(3)*Vs(i)*(Bps(3,i)*cos(tetap(3)-tetas(i))-Gps(3,i)*sin(tetap(3)-tetas(i)));
end

% fluxo inverso

Qsp_a=-(Vs(1)^2)*Bss(1,1);

for i=2:3
  Qsp_a=Qsp_a-Vs(1)*Vs(i)*(Bss(1,i)*cos(tetas(1)-tetas(i))-Gss(1,i)*sin(tetas(1)-tetas(i)));
end

for i=1:3
  Qsp_a=Qsp_a-Vs(1)*Vp(i)*(Bsp(1,i)*cos(tetas(1)-tetap(i))-Gsp(1,i)*sin(tetas(1)-tetap(i)));
end

Qsp_b=-(Vs(2)^2)*Bss(2,2);

for i=[1,3]
  Qsp_b=Qsp_b-Vs(2)*Vs(i)*(Bss(2,i)*cos(tetas(2)-tetas(i))-Gss(2,i)*sin(tetas(2)-tetas(i)));
end

for i=1:3
  Qsp_b=Qsp_b-Vs(2)*Vp(i)*(Bsp(2,i)*cos(tetas(2)-tetap(i))-Gsp(2,i)*sin(tetas(2)-tetap(i)));
end


Qsp_c=-(Vs(3)^2)*Gss(3,3);

for i=[1,2]
  Qsp_c=Qsp_c-Vs(3)*Vs(i)*(Bss(3,i)*cos(tetas(3)-tetas(i))-Gss(3,i)*sin(tetas(3)-tetas(i)));
end

for i=1:3
  Qsp_c=Qsp_c-Vs(3)*Vp(i)*(Bsp(3,i)*cos(tetas(3)-tetap(i))-Gsp(3,i)*sin(tetas(3)-tetap(i)));
end

