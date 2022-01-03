clear all;
clc;
syms tetap_a tetap_b tetap_c;
syms tetas_a tetas_b tetas_c;
syms Vp_a Vp_b Vp_c;
syms Vs_a Vs_b Vs_c;
Gpp=sym("Gpp",[3 3]);
Gps=sym("Gps",[3 3]);
Gsp=sym("Gsp",[3 3]);
Gss=sym("Gss",[3 3]);

Bpp=sym("Bpp",[3 3]);
Bps=sym("Bps",[3 3]);
Bsp=sym("Bsp",[3 3]);
Bss=sym("Bss",[3 3]);

syms Gt Bt;
tetap=[tetap_a,tetap_b,tetap_c];
tetas=[tetas_a,tetas_b,tetas_c];
Vp=[Vp_a, Vp_b, Vp_c];
Vs=[Vs_a, Vs_b, Vs_c];

teta3ph_p=[0 -120 120]*pi/180;
teta3ph_s=[-30 -150 90]*pi/180;

%% fluxo de potência ativa de p para s 

%fase a 
Pps_a=(Vp(1)^2)*Gpp(1,1);

for i=2:3
  Pps_a=Pps_a+Vp(1)*Vp(i)*(Gpp(1,i)*cos(tetap(1)-tetap(i))+Bpp(1,i)*sin(tetap(1)-tetap(i)));
end


for i=1:3
  Pps_a=Pps_a+Vp(1)*Vs(i)*(Gps(1,i)*cos(tetap(1)-tetas(i))+Bps(1,i)*sin(tetap(1)-tetas(i)));
end

%fase b
Pps_b=(Vp(2)^2)*Gpp(2,2);

for i=[1,3]
  Pps_b=Pps_b+Vp(2)*Vp(i)*(Gpp(2,i)*cos(tetap(2)-tetap(i))+Bpp(2,i)*sin(tetap(2)-tetap(i)));
end

for i=1:3
  Pps_b=Pps_b+Vp(2)*Vs(i)*(Gps(2,i)*cos(tetap(2)-tetas(i))+Bps(2,i)*sin(tetap(2)-tetas(i)));
end

%fase c
Pps_c=(Vp(3)^2)*Gpp(3,3);

for i=[1,2]
  Pps_c=Pps_c+Vp(3)*Vp(i)*(Gpp(3,i)*cos(tetap(3)-tetap(i))+Bpp(3,i)*sin(tetap(3)-tetap(i)));
end

for i=1:3
  Pps_c=Pps_c+Vp(3)*Vs(i)*(Gps(3,i)*cos(tetap(3)-tetas(i))+Bps(3,i)*sin(tetap(3)-tetas(i)));
end


%% fluxo de potência ativa de s para p
%fase a 
Psp_a=(Vs(1)^2)*Gss(1,1);

for i=2:3
  Psp_a=Psp_a+Vs(1)*Vs(i)*(Gss(1,i)*cos(tetas(1)-tetas(i))+Bss(1,i)*sin(tetas(1)-tetas(i)));
end

for i=1:3
  Psp_a=Psp_a+Vs(1)*Vp(i)*(Gsp(1,i)*cos(tetas(1)-tetap(i))+Bsp(1,i)*sin(tetas(1)-tetap(i)));
end
%fase b
Psp_b=(Vs(2)^2)*Gss(2,2);

for i=[1,3]
  Psp_b=Psp_b+Vs(2)*Vs(i)*(Gss(2,i)*cos(tetas(2)-tetas(i))+Bss(2,i)*sin(tetas(2)-tetas(i)));
end

for i=1:3
  Psp_b=Psp_b+Vs(2)*Vp(i)*(Gsp(2,i)*cos(tetas(2)-tetap(i))+Bsp(2,i)*sin(tetas(2)-tetap(i)));
end

%fase c
Psp_c=(Vs(3)^2)*Gss(3,3);

for i=[1,2]
  Psp_c=Psp_c+Vs(3)*Vs(i)*(Gss(3,i)*cos(tetas(3)-tetas(i))+Bss(3,i)*sin(tetas(3)-tetas(i)));
end

for i=1:3
  Psp_c=Psp_c+Vs(3)*Vp(i)*(Gsp(3,i)*cos(tetas(3)-tetap(i))+Bsp(3,i)*sin(tetas(3)-tetap(i)));
end


%% fluxo de potência reativa de p para s

%fase a 

Qps_a=-(Vp(1)^2)*Bpp(1,1);

for i=2:3
  Qps_a=Qps_a-Vp(1)*Vp(i)*(Bpp(1,i)*cos(tetap(1)-tetap(i))-Gpp(1,i)*sin(tetap(1)-tetap(i)));
end

for i=1:3
  Qps_a=Qps_a-Vp(1)*Vs(i)*(Bps(1,i)*cos(tetap(1)-tetas(i))-Gps(1,i)*sin(tetap(1)-tetas(i)));
end
%fase b
Qps_b=-(Vp(2)^2)*Bpp(2,2);

for i=[1,3]
  Qps_b=Qps_b-Vp(2)*Vp(i)*(Bpp(2,i)*cos(tetap(2)-tetap(i))-Gpp(2,i)*sin(tetap(2)-tetap(i)));
end

for i=1:3
  Qps_b=Qps_b-Vp(2)*Vs(i)*(Bps(2,i)*cos(tetap(2)-tetas(i))-Gps(2,i)*sin(tetap(2)-tetas(i)));
end

%fase c
Qps_c=-(Vp(3)^2)*Bpp(3,3);

for i=[1,2]
  Qps_c=Qps_c-Vp(3)*Vp(i)*(Bpp(3,i)*cos(tetap(3)-tetap(i))-Gpp(3,i)*sin(tetap(3)-tetap(i)));
end

for i=1:3
  Qps_c=Qps_c-Vp(3)*Vs(i)*(Bps(3,i)*cos(tetap(3)-tetas(i))-Gps(3,i)*sin(tetap(3)-tetas(i)));
end
%% fluxo de potência reativa de s para p
%fase a

Qsp_a=-(Vs(1)^2)*Bss(1,1);

for i=2:3
  Qsp_a=Qsp_a-Vs(1)*Vs(i)*(Bss(1,i)*cos(tetas(1)-tetas(i))-Gss(1,i)*sin(tetas(1)-tetas(i)));
end

for i=1:3
  Qsp_a=Qsp_a-Vs(1)*Vp(i)*(Bsp(1,i)*cos(tetas(1)-tetap(i))-Gsp(1,i)*sin(tetas(1)-tetap(i)));
end

%fase b
Qsp_b=-(Vs(2)^2)*Bss(2,2);

for i=[1,3]
  Qsp_b=Qsp_b-Vs(2)*Vs(i)*(Bss(2,i)*cos(tetas(2)-tetas(i))-Gss(2,i)*sin(tetas(2)-tetas(i)));
end

for i=1:3
  Qsp_b=Qsp_b-Vs(2)*Vp(i)*(Bsp(2,i)*cos(tetas(2)-tetap(i))-Gsp(2,i)*sin(tetas(2)-tetap(i)));
end

%fase c
Qsp_c=-(Vs(3)^2)*Gss(3,3);

for i=[1,2]
  Qsp_c=Qsp_c-Vs(3)*Vs(i)*(Bss(3,i)*cos(tetas(3)-tetas(i))-Gss(3,i)*sin(tetas(3)-tetas(i)));
end

for i=1:3
  Qsp_c=Qsp_c-Vs(3)*Vp(i)*(Bsp(3,i)*cos(tetas(3)-tetap(i))-Gsp(3,i)*sin(tetas(3)-tetap(i)));
end

% vetor h simbolico
hp=[Pps_a,Pps_b,Pps_c,Psp_a,Psp_b,Psp_c];
hq=[Qps_a,Qps_b,Qps_c,Qsp_a,Qsp_b,Qsp_c];

H=jacobian(hp,[tetap tetas Vp Vs]);

H=simplify(H);


%% colunas isoladas da jacobiana transposta de p para s
GradPps_a=gradient(Pps_a,[tetap tetas Vp Vs]);
GradPps_a=subs(GradPps_a,[tetap tetas Vp Vs],[teta3ph_p teta3ph_s 1 1 1 1 1 1]);



%elementos da diagonal dos quadripolos
auxg=Gt;
auxb=Bt;
%ss e pp
GradPps_a=subs(GradPps_a,[Bss(1,1) Gss(1,1) Bpp(1,1) Gpp(1,1)],[auxb auxg auxb*2/3 auxg*2/3]);
%sp e ps
auxg=Gt*sqrt(3)/3;
auxb=Bt*sqrt(3)/3;
GradPps_a=subs(GradPps_a,[Bsp(1,1) Gsp(1,1) Bps(1,1) Gps(1,1)],[-auxb -auxg -auxb -auxg]);
%elementos fora da diagonal dos quadripolos 
auxg=Gt*sqrt(3)/3;
auxb=Bt*sqrt(3)/3;
%a b
GradPps_a=subs(GradPps_a,[Bsp(1,2) Gsp(1,2) Bps(1,2) Gps(1,2)],[0 0 auxb auxg]);
%a c
GradPps_a=subs(GradPps_a,[Bsp(1,3) Gsp(1,3) Bps(1,3) Gps(1,3)],[auxb auxg 0 0]);

auxg=Gt/3;
auxb=Bt/3;

GradPps_a=subs(GradPps_a,[Bss(1,2) Gss(1,2) Bpp(1,2) Gpp(1,2)],[0  0 -auxb -auxg]);
GradPps_a=subs(GradPps_a,[Bss(1,3) Gss(1,3) Bpp(1,3) Gpp(1,3)],[0  0 -auxb -auxg]);


%% colunas isoladas da jacobiana transposta de s para p


GradPsp_a=gradient(Psp_a,[tetap tetas Vp Vs]);
GradPsp_a=subs(GradPsp_a,[tetap tetas Vp Vs],[teta3ph_p teta3ph_s 1 1 1 1 1 1]);



%elementos da diagonal dos quadripolos
auxg=Gt;
auxb=Bt;
%ss e pp
GradPsp_a=subs(GradPsp_a,[Bss(1,1) Gss(1,1) Bpp(1,1) Gpp(1,1)],[auxb auxg auxb*2/3 auxg*2/3]);
%sp e ps
auxg=Gt*sqrt(3)/3;
auxb=Bt*sqrt(3)/3;
GradPsp_a=subs(GradPsp_a,[Bsp(1,1) Gsp(1,1) Bps(1,1) Gps(1,1)],[-auxb -auxg -auxb -auxg]);
%elementos fora da diagonal dos quadripolos 
auxg=Gt*sqrt(3)/3;
auxb=Bt*sqrt(3)/3;
%a b
GradPsp_a=subs(GradPsp_a,[Bsp(1,2) Gsp(1,2) Bps(1,2) Gps(1,2)],[0 0 auxb auxg]);
%a c
GradPsp_a=subs(GradPsp_a,[Bsp(1,3) Gsp(1,3) Bps(1,3) Gps(1,3)],[auxb auxg 0 0]);

auxg=Gt/3;
auxb=Bt/3;

GradPsp_a=subs(GradPsp_a,[Bss(1,2) Gss(1,2) Bpp(1,2) Gpp(1,2)],[0  0 -auxb -auxg]);
GradPsp_a=subs(GradPsp_a,[Bss(1,3) Gss(1,3) Bpp(1,3) Gpp(1,3)],[0  0 -auxb -auxg]);