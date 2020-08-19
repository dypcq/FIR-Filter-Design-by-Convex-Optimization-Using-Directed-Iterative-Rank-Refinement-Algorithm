function [h,ph_dev,s_en,Rat,status,high_sol,r_save] = IRR_v13(L,f_s,f_p,s_v,p_l_v,...
    p_u_v,irr_in,irr_c,phase1,phase2,p_add,type,Fs,stop_lim,err)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% L          : Order of the filter,
% f_s        : Stopband frequency beginning sample,
% f_p        : Passband frequency finishing sample,
% s_v        : Stopband value (on linear real line),
% p_l_v      : Passband lower mask value (on linear real line),
% p_u_v      : Passband upper mask value (on linear real line),
% phase1     : Slope of the group delay mask,
% phase2     : DC of the group delay mask,
% type       : Type of the group delay mask (constant, linear, quadratic),
% Fs         : Sampling frequency,
% stop_lim   : Stopping ratio limit for the algorithm,
% err        : Error margin on frequency domain constraints,
% p_add      : Iteration number without group delay constraints are added,
% irr_in     : Initial upper limit on rank constraint,
% irr_c      : Multiplicative constant on irr_in,

% Outputs:
% h          : Filter coefficients,
% ph_dev     : Phase deviation from best line, (UNUSED)
% s_en       : Stopband energy,
% Rat        : Ratio of first left singular vector to total energy of
%              singular vectors.

% Description:
% Given filter specifications, this function finds the filter coefficients.
% Both magnitude and group delay constraints can be applied in a flexible manner.
% Algorithm uses iterative SDP with rank refinement algorithm. Filter 
% coefficients are real.

% Modification:
% 1. Group delay constraints can be added later.
% 2. Positive side of the group delay constraints is used.
% 3. If optimization fails to be feasible in iteration i, for i>1, then
% rank constrained is smoothed DETERMINISTICALLY and problem is resolved.
% 4. cvx_optval can be NaN in any iteration.
% 5. Not all group delay constraints are added to optimization problem. One
% group delay constraint is added in 10 samples. (Can be changed.)
% 6. Objective function now maximizes the energy of first left-singular
% vector.
% 7. Comments are removed.(Look at IRR_v7 for comments!)
% 8. Tries all solutions in all directions.
% 9. Modified over IRR_v10.

% Authors    : Mehmet Dedeoðlu
% Date       : 05/16/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
close all;
high_sol = 0;

% Initializations
sample = Fs*1+1;
fgrid = linspace(-Fs/2,Fs/2,sample);

% Produce frequency domain pulses
freq = fgrid(Fs/2+1:end);
Hf = zeros(L,length(fgrid));
Hf_r = zeros(L,length(fgrid));
Hf_i = zeros(L,length(fgrid));
gs = zeros(L,length(fgrid));
gc = zeros(L,length(fgrid));
dgs = zeros(L,length(fgrid));
dgc = zeros(L,length(fgrid));
for a = 0:L-1
    Hf(a+1,:) = exp(-1i*2*pi*(a-L/2+0.5)*fgrid/length(fgrid));
    gs(a+1,:) = -sin(2*pi*(a-L/2+0.5)*fgrid/length(fgrid));
    gc(a+1,:) = cos(2*pi*(a-L/2+0.5)*fgrid/length(fgrid));
    dgs(a+1,:) = -(a-L/2+0.5)/Fs*2*pi*cos(2*pi*(a-L/2+0.5)*fgrid/length(fgrid));
    dgc(a+1,:) = -(a-L/2+0.5)/Fs*2*pi*sin(2*pi*(a-L/2+0.5)*fgrid/length(fgrid));
    Hf_s = Hf(:,Fs/2+1:end);
    Hf_r = real(Hf);
    Hf_rs = Hf_r(:,Fs/2+1:end);
    Hf_i = imag(Hf);
    Hf_is = Hf_i(:,Fs/2+1:end);
end;

% Produce magnitude mask
publl = -f_s;
pubul = f_s;
pubv = p_u_v;
plbll = -f_p;
plbul = f_p;
plbv = p_l_v;
sll = -f_s;
sul = f_s;
sv = s_v;
[M1,M2] = filtermask_v1(fgrid,publl,pubul,pubv,plbll,plbul,plbv,sll,sul,sv);
M1_s = M1(Fs/2+1:end);
M2_s = M2(length(M2)/2+0.5:end);
M3_s = M1_s(M1_s>=0.99);
figure;plot(freq/Fs,10*log10(M1_s));hold on;plot(freq(freq<plbul)/Fs,...
    10*log10(M2_s));ylim([10*log10(s_v/sqrt(10)) 1]);
ylabel('Magnitude (dB)');xlabel('Frequency');
handaxes2 = axes('position', [0.18 0.2 0.3 0.5]);
plot(freq/Fs,10*log10(M1_s));hold on;plot(freq(freq<plbul)/Fs,...
    10*log10(M2_s));ylim([-1/2 0.05]);xlim([0 pubul/Fs]);
set(handaxes2, 'box', 'on');

% Produce phase mask
if type == 0
    val = phase1;
elseif type == 1
    val = [phase1,phase2];
else
    val = phase1;
end;
[PM_u,PM_l] = phasemask_v2(type,fgrid,val);
figure;plot(fgrid/Fs,((PM_u)));hold on;plot(fgrid/Fs,(zeros(length(PM_u),1)));
ylabel('Phase');xlabel('Frequency (Hz)');
figure;plot(fgrid(fgrid>sll & fgrid<sul)/Fs,atan(sqrt(PM_u(fgrid>sll & fgrid<sul))));
hold on;plot(fgrid(fgrid>sll & fgrid<sul)/Fs,-atan(sqrt(PM_u(fgrid>sll & fgrid<sul))));
ylabel('Phase');xlabel('Frequency (Hz)');grid on;

% Produce matrices for phase constraint
I = zeros(L,L,length(fgrid));
R = zeros(L,L,length(fgrid));
for c = 1:length(fgrid)
    hf_i = Hf_i(:,c);
    I(:,:,c) = hf_i*hf_i';
    
    hf_r = Hf_r(:,c);
    R(:,:,c) = PM_u(c)*(hf_r*hf_r');
end;
P = I-R;
P = P(:,:,(fgrid<plbul & fgrid>=0));
I = I(:,:,(fgrid<plbul & fgrid>=0));
R = R(:,:,(fgrid<plbul & fgrid>=0));

% Produce group delay matrices
Dk = zeros(L,L,length(fgrid));
Cg = zeros(L,L,length(fgrid));
for c = 1:length(fgrid)
    gs_i = gs(:,c);
    gc_i = gc(:,c);
    dgs_i = dgs(:,c);
    dgc_i = dgc(:,c);
    Dk(:,:,c) = gs_i*dgc_i'-gc_i*dgs_i';
    chf = Hf(:,c);
    Cg(:,:,c) = chf*chf';
end;

% Produce C matrices (positive semidefinite matrices)
f_sample = length(freq);
C = zeros(L,L,f_sample);
for c = 1:f_sample
    hf = Hf_s(:,c);
    C(:,:,c) = hf*hf';
end;

Clowconstraint = C(:,:,(freq<plbul));
% Cpass = -(sum(Clowconstraint,3)-sum(C(:,:,(freq>plbul)),3));

Cpass = sum(C(:,:,(freq>sul)),3);
% C_equ = sum(C,3);

% Produce group delay constraints
glk = -phase1;
guk = -glk;

glkvec = glk*ones(1,length(fgrid(fgrid<plbul & fgrid>-plbul)));
gukvec = guk*ones(1,length(fgrid(fgrid<plbul & fgrid>-plbul)));

Lk = glk*Cg(:,:,(fgrid<plbul & fgrid>=0))+Dk(:,:,(fgrid<plbul & fgrid>=0));
Uk = guk*Cg(:,:,(fgrid<plbul & fgrid>=0))+Dk(:,:,(fgrid<plbul & fgrid>=0));

% Optimization
Rat = 0;
figure;
order_r = L;
delSing=irr_in;
delta = 1;
X_old = ones(order_r,order_r);
ratio_old = 0;
temp_mul = irr_c;
Uvec=zeros(order_r,order_r-1);
UU1=zeros(order_r,order_r);
rfig=[];
for n=1:1
    disp(n);
%--------------------------------------------------------------------------
%solve semidefinite program
cvx_begin
cvx_precision best
cvx_solver SDPT3
variable X(order_r,order_r) symmetric
variable t(1)
maximize delta*trace(UU1*X)-(1-delta)*t
subject to
    for i=1:length(M3_s)
        trace(C(:,:,i)*X) <= M3_s(i)-err;
    end
    for i=length(M3_s)+1:f_sample
        trace(C(:,:,i)*X) <= t;
    end
    for i=1:length(M2_s)
        trace(Clowconstraint(:,:,i)*X) >= M2_s(i)+err;
    end
%     for i=1:order_r-1
%         Uvec(:,i)'*X*Uvec(:,i) <= delSing;
%     end
%     if n>p_add
%         for i=[1:length(fgrid(fgrid<plbul & fgrid>=0)),length(fgrid(fgrid<plbul & fgrid>=0))]
% %             trace(P(:,:,i)*X) <= 0;
%             trace(Lk(:,:,i)*X) <= 0;
%             trace(Uk(:,:,i)*X) >= 0;
%         end
%     end
X==semidefinite(order_r);
cvx_end

X = full(X);
nan_check = isnan(X);
if sum(nan_check(:)) == 0 && isnan(cvx_optval) == 0
    temp_mul = irr_c;
    high_sol = 1;
    [U1,S1,V1]=svd(X);
    UU1 = U1(:,1)*(U1(:,1))';%%%%%%%
    Uvec=U1(:,2:end);
    Rat = S1(1,1)/(sum(diag(S1)));
    rfig = [rfig,Rat];
    if Rat > ratio_old
        X_old = X;
        ratio_old = Rat;
    end;
    plot(diag(S1),'.-');axis tight;ylim([-1 6]);xlim([-1 order_r]);grid;
    title(['Iteration number: ',num2str(n),' ratio: ',...
        num2str(Rat),' rank limit: ',num2str(delSing)]);drawnow;
    ylabel('Magnitude of singular values');xlabel('Singular value number');
    delSing=delSing*irr_c;
    status = 0;
    if Rat>stop_lim
        status = 1000000;
        break;
    end;
    
else
    temp_mul = [temp_mul,((1+1/prod(temp_mul))/2)];
    delSing = delSing*temp_mul(end);
    if n == 1
        X = ones(order_r,order_r);
    else
        X = X_old;
    end;
    disp('Convex optimization failed!');
    status = n;
end;

end

% Make distributed optimization
% Decide on which singular vectors will be used
singular = diag(S1)/(sum(diag(S1)));
aim_num = sum(singular>0);
% Make opimization for all aim_num
XX = X;
r_save = ratio_old;
for ii=1:aim_num
[X,high_sol,status,Rat,rfig,figs1,tt] = optimizer_v6(ii,30,U1,S1,C,Cpass,Clowconstraint,M1_s,...
    M2_s,M3_s,err,delSing,p_add,fgrid,plbul,Lk,Uk,order_r,irr_c,ratio_old,...
    stop_lim,temp_mul,X_old,f_sample,rfig,delta);
if Rat > r_save
    r_save = Rat;
    XX = X;
end;
if Rat > stop_lim
    break;
end;
end;
X = XX;

% Construction
[UU,SS,VV]=svd(X);
g_vec_r=UU(:,1)*sqrt(SS(1,1));

% % Resulting coefficients
h = g_vec_r;

% Plot time domain and frequency domain response from time domain samples.
Ts = 1/Fs;
tgrid = [0:Ts:(length(h)-1)*Ts]-Ts*(L/2-0.5);
fir = h;
figure;stem(tgrid,fir);xlabel('Time (seconds)');ylabel('Magnitude');
title(['Designed filter in time domain with # of coefficients = ',num2str(L)]);
[H,I] = meshgrid(tgrid,fgrid);
v = exp(1i.*2.*pi.*H.*I);
v = v';
Ffir = fir'*v;

Mfig = [10*log10(tt(end)*ones(length(fgrid(-f_s>=fgrid)),1));...
    10*log10(M1(-f_s<fgrid & fgrid<f_s)).';...
    10*log10(tt(end)*ones(length(fgrid(f_s<=fgrid)),1))];
figure;plot(fgrid/Fs,20*log10(abs(Ffir)));hold on;
plot(fgrid/Fs,Mfig,'r');hold on;
plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,10*log10(M2),'r');ylim([10*log10(s_v/100) 1]);
xlabel('Frequency');ylabel('Magnitude (dB)');
title(['Ratio = ',num2str(r_save)]);

figure;plot(fgrid/Fs,20*log10(abs(Ffir)));hold on;
plot(fgrid(-f_s<fgrid & fgrid<f_s)/Fs,10*log10(M1(-f_s<fgrid & fgrid<f_s)),'r');hold on;
plot(fgrid(fgrid<=-f_s)/Fs,10*log10(tt(end)*ones(length(fgrid(-f_s>=fgrid)),1)),'r');hold on;
plot(fgrid(fgrid>=f_s)/Fs,10*log10(tt(end)*ones(length(fgrid(f_s<=fgrid)),1)),'r');hold on;
plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,10*log10(M2),'r');ylim([10*log10(s_v/100) 1]);
xlabel('Frequency');ylabel('Magnitude (dB)');
title(['Ratio = ',num2str(r_save)]);

handaxes2 = axes('position', [0.38 0.2 0.3 0.5]);
plot(fgrid/Fs,20*log10(abs(Ffir)));hold on;plot(fgrid/Fs,10*log10(M1),'r');hold on;
plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,10*log10(M2),'r');
ylim([-1/2 0.05]);xlim([-pubul/Fs pubul/Fs]);
set(handaxes2, 'box', 'on');

toc;

Ffir_angle = unwrap(angle(Ffir));
Ffir_angle = Ffir_angle-Ffir_angle(fgrid==0);
Ffir_angle_pass = Ffir_angle(fgrid>plbll & fgrid<plbul);
freq_pass = fgrid(fgrid>plbll & fgrid<plbul);
est = leastSqrEst(Ffir_angle_pass,freq_pass);
est_err = (Ffir_angle_pass-est)*(Ffir_angle_pass-est)';
lin_phase_step = (Ffir_angle_pass(end)-Ffir_angle_pass(1))/(freq_pass(end)-freq_pass(1));
lin_phase = Ffir_angle_pass(1):lin_phase_step:Ffir_angle_pass(end);
lin_phase_err = (Ffir_angle_pass-lin_phase)*(Ffir_angle_pass-lin_phase)';
figure;plot(fgrid/Fs,Ffir_angle);hold on;plot(freq_pass/Fs,est,'g');hold on;
plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,atan(sqrt(PM_u(fgrid>plbll & fgrid<plbul))),'r');
hold on;plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,-atan(sqrt(PM_u(fgrid>plbll & fgrid<plbul))),'r');%hold on;plot(freq_pass/Fs,lin_phase,'g');
title(['estimation error: ',num2str(est_err),' linear phase error: ',num2str(lin_phase_err)]);grid on;

figure;plot(figs1,'Marker','x');

groupdelay = zeros(1,501);
for i=1:length(fgrid)
    groupdelay(i) = -(trace(h'*Dk(:,:,i)*h))/(trace(h'*Cg(:,:,i)*h));
end;

figure;plot(fgrid/Fs,groupdelay);hold on;
plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,glkvec,'r');
hold on;plot(fgrid(fgrid>plbll & fgrid<plbul)/Fs,gukvec,'r');

figure;plot(tt);

s_en = Ffir(fgrid<sll | fgrid>sul)*(Ffir(fgrid<sll | fgrid>sul))'/length(Ffir(fgrid<sll | fgrid>sul));
ph_dev = est_err;

% Save the results
save(['IRR_v13-L-',num2str(L),'-f_s-',num2str(f_s),'-f_p-',num2str(f_p),...
    '-s_v-',num2str(s_v),'-p_l_v-',num2str(p_l_v),'-p_u_v-',num2str(p_u_v),...
    '-irr_in-',num2str(irr_in),'-irr_c-',num2str(irr_c),'-phase1-',...
    num2str(phase1),'-phase2-',num2str(phase2),'-p_add-',num2str(p_add),...
    '-type-',num2str(type),'-Fs-',num2str(Fs),'-stop_lim-',num2str(stop_lim),...
    '-err-',num2str(err),'.mat'],'h','ph_dev','s_en','Rat','status','high_sol');

end