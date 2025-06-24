function [delta_J,d,K] = backward_pass_iLQR(Nx,Nt,Nu, xtraj, xgoal,utraj,QN,Q,R,dfdx,dfdu,dAdx,dAdu,dBdx,dBdu)
    import casadi.*

    delta_J = 0.0;
    p(:,Nt) = QN*(xtraj(:,Nt)-xgoal);
    P(:,:,Nt) = QN;
    
    A = cell(1, Nt-1);
    B = cell(1, Nt-1);
    %Ax = cell(1, Nt-1);
    %Bx = cell(1, Nt-1);
    %Au = cell(1, Nt-1);
    %Bu = cell(1, Nt-1);

    parfor k = 1:(Nt-1)
        A{k} = full(dfdx(xtraj(:,k), utraj(:,k)));
        B{k} = full(dfdu(xtraj(:,k), utraj(:,k)));
        
        %Ax{k} = full(dAdx(xtraj(:,k), utraj(:,k)));
        %Bx{k} = full(dBdx(xtraj(:,k), utraj(:,k)));
        %Au{k} = full(dAdu(xtraj(:,k), utraj(:,k)));
        %Bu{k} = full(dBdu(xtraj(:,k), utraj(:,k)));

    end

    for k = (Nt-1):-1:1
        %% Calculate derivatives
        gx = Q*(xtraj(:,k)-xgoal) + A{k}'*p(:,k+1);
        gu = R*utraj(:,k) + B{k}'*p(:,k+1);
    
        %Gxx = Q + A{k}'*P(:,:,k+1)*A{k} + kron(p(:,k+1)',eye(Nx))*comm(Nx,Nx)*Ax{k};
        %Guu = R + B{k}'*P(:,:,k+1)*B{k} + kron(p(:,k+1)',eye(Nu))*comm(Nx,Nu)*Bu{k};
        %Gxu = A{k}'*P(:,:,k+1)*B{k} + kron(p(:,k+1)',eye(Nx))*comm(Nx,Nx)*Au{k};
        %Gux = B{k}'*P(:,:,k+1)*A{k} + kron(p(:,k+1)',eye(Nu))*comm(Nx,Nu)*Bx{k};

        Gxx = Q + A{k}'*P(:,:,k+1)*A{k};
        Guu = R + B{k}'*P(:,:,k+1)*B{k};
        Gxu = A{k}'*P(:,:,k+1)*B{k};
        Gux = B{k}'*P(:,:,k+1)*A{k};
        
        %% regularzation
        beta = 0.1;
        while ~isPositiveDe([Gxx, Gxu; Gux, Guu])
            Gxx = Gxx+ A{k}'*beta*eye(size(A{k},1))*A{k};
            Guu = Guu+ B{k}'*beta*eye(size(B{k},1))*B{k};
            Gxu = Gxu+ A{k}'*beta*eye(size(B{k},1))*B{k};
            Gux = Gux+ B{k}'*beta*eye(size(A{k},1))*A{k};
            beta = 2*beta;
            disp("regularizing G");
        end
        

    %for k = (Nt-1):-1:1
        %A = full(dfdx(xtraj(:,k), utraj(:,k)));
        %B = full(dfdu(xtraj(:,k), utraj(:,k)));
        
        %Ax = full(dAdx(xtraj(:,k), utraj(:,k)));
        %Bx = full(dBdx(xtraj(:,k), utraj(:,k)));
        %Au = full(dAdu(xtraj(:,k), utraj(:,k)));
        %Bu = full(dBdu(xtraj(:,k), utraj(:,k)));

        %% Calculate derivatives
        %gx = Q*(xtraj(:,k)-xgoal) + A'*p(:,k+1);
        %gu = R*utraj(:,k) + B'*p(:,k+1);
    
        %Gxx = Q + A'*P(:,:,k+1)*A + kron(p(:,k+1)',eye(Nx))*comm(Nx,Nx)*Ax;
        %Guu = R + B'*P(:,:,k+1)*B + kron(p(:,k+1)',eye(Nu))*comm(Nx,Nu)*Bu;
        %Gxu = A'*P(:,:,k+1)*B + kron(p(:,k+1)',eye(Nx))*comm(Nx,Nx)*Au;
        %Gux = B'*P(:,:,k+1)*A + kron(p(:,k+1)',eye(Nu))*comm(Nx,Nu)*Bx;

        %Gxx = Q + A'*P(:,:,k+1)*A;
        %Guu = R + B'*P(:,:,k+1)*B;
        %Gxu = A'*P(:,:,k+1)*B;
        %Gux = B'*P(:,:,k+1)*A;
        
        %% regularzation
        %beta = 0.1;
        %while ~isPositiveDe([Gxx, Gxu; Gux, Guu])
            %Gxx = Gxx+ A'*beta*eye(size(A,1))*A;
            %Guu = Guu+ B'*beta*eye(size(B,1))*B;
            %Gxu = Gxu+ A'*beta*eye(size(B,1))*B;
            %Gux = Gux+ B'*beta*eye(size(A,1))*A;
            %beta = 2*beta;
            %disp("regularizing G");
        %end
        %% calculate controll
        d(:,k) = Guu\gu;
        K(:,:,k)= Guu\Gux;
    
        p(:,k)= gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(:,k) - Gxu*d(:,k);
        P(:,:,k)= Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux;
    
        delta_J = delta_J+ gu'*d(:,k);
    end
end