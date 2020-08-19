function [X,high_sol,status,ratio_old,rfig,figs1] = optimizer_v5(aim,n,U1,S1,C,Cpass,...
    Clowconstraint,M1_s,M2_s,err,delSing,p_add,fgrid,plbul,Lk,Uk,order_r,...
    irr_c,ratio_old,stop_lim,temp_mul,X_old,f_sample,rfig,delta)

Rat = ratio_old;
high_sol = 0;
ratio_temp=0;
UU1 = U1(:,aim)*(U1(:,aim))';
Uvec = [U1(:,1:aim-1),U1(:,aim+1:end)];
S1(1,1)=0;
figs1=[];
for ii=1:n
    disp(ii);
    disp(S1(1,1));
    delta = 1-Rat^(1);
%--------------------------------------------------------------------------
%solve semidefinite program
cvx_begin
cvx_precision best
cvx_solver SDPT3
variable X(order_r,order_r) symmetric
maximize delta*trace(UU1*X)-(1-delta)*(real(trace(Cpass*X)))
subject to
    for i=1:f_sample
        trace(real(C(:,:,i)*X)) <= M1_s(i)-err;
    end
    for i=1:length(M2_s)
        trace(real(Clowconstraint(:,:,i)*X)) >= M2_s(i)+err;
    end
%     for i=1:order_r-1
%         Uvec(:,i)'*X*Uvec(:,i) <= delSing;
%     end
    if n>p_add
        for i=[1:length(fgrid(fgrid<plbul & fgrid>=0)),length(fgrid(fgrid<plbul & fgrid>=0))]
%             trace(P(:,:,i)*X) <= 0;
            real(trace(Lk(:,:,i)*X)) <= 0;
            real(trace(Uk(:,:,i)*X)) >= 0;
        end
    end
X==semidefinite(order_r);
cvx_end

X = full(X);
nan_check = isnan(X);
if sum(nan_check(:)) == 0 && isnan(cvx_optval) == 0
    temp_mul = irr_c;
    high_sol = 1;
    [U1,S1,V1]=svd(X);
    figs1=[figs1,S1(1,1)];
    UU1 = U1(:,1)*(U1(:,1))';%%%%%%%
    Uvec=U1(:,2:end);
    Rat = S1(1,1)/(sum(diag(S1)));
    rfig = [rfig,Rat];
%     if abs(Rat-ratio_temp)<0.00001%%%%%%%%%%%
%         break;
%     end;%%%%%%%%%%%
    if Rat > ratio_old
        X_old = X;
        ratio_old = Rat;
    end;
    ratio_temp=Rat;%%%%%%%%%%
    semilogy(diag(S1)/sum(diag(S1)),'.-');axis tight;ylim([-1 6]);xlim([-1 order_r]);grid;
    title([num2str(aim),'th. Iteration number: ',num2str(ii),' ratio: ',...
        num2str(Rat),' rank limit: ',num2str(delSing)]);drawnow;
    ylabel('Magnitude of singular values');xlabel('Singular value number');
    delSing=delSing*irr_c;
    status = 0;
    if Rat>stop_lim
        status = 1000000;
        break;
    end;
    
else
    rfig=[rfig,Rat];
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

end;
X = X_old;
end