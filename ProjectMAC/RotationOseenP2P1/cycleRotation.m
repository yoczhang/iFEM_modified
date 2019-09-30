function bigu = cycleRotation(uh, vh, ph,  A, B, f1, f2, g, Sp, ProU, ResU, ProP, ResP, level, node, elem)


%% Cycle
error = 1; step = 0; s1 = 5; s2 = 5; residual = []; 
Nu = zeros(level,1); Np = zeros(level,1); 
Nu(level) = size(uh,1); Np(level) = size(ph,1);

bigA = [A{level} B{level}'; B{level} sparse(Np(level),Np(level))];
bigF = [f1; f2; g]; bigu = [uh; vh; ph];
r = bigF - bigA*bigu;  
figure, trisurf(elem{level}, node{level}(:,1), node{level}(:,2), r(1:Np(level))); shading interp; colorbar; view(3);


while (error > 1e-8)
    eh = zeros(2*Nu(level)+Np(level),1); 
    eh = Vcycle(eh, r, level);
    bigu = bigu + eh;
    
    r = bigF - bigA*bigu; error = norm(r) / norm(bigF)
    step = step + 1; residual(step) = error;
end

%% Vcycle 
    function bige = Vcycle(bige, rhs, J)
        if J == 2
            Nu(J) = size(A{J},1)/2; Np(J) = size(B{J},1);
            temp = 1e-10*spdiags(ones(Np(J),1), 0, Np(J), Np(J));
            bigA3 = [A{J} B{J}'; B{J} temp];
            bige = bigA3 \ rhs; 
            return;
        end
        
        % Pre-smooth
        for i = 1:s1
            bige = LSC_DGS_rotation_smoother(bige, A{J}, B{J}, rhs, Sp{J}, node{J}, elem{J});
        end
        
        % Get residual
        Nu(J) = size(A{J},1)/2; Np(J) = size(B{J},1);
        bigA2 = [A{J} B{J}'; B{J} sparse(Np(J),Np(J))];
        resh = rhs - bigA2*bige;
%         figure, trisurf(elem{J}, node{J}(:,1), node{J}(:,2), resh(1:Np(J))); shading interp; colorbar; view(2);
        
        % Restrict residual to the coarse grid
        resH = [ResU{J-1}*resh(1:Nu(J)); 
                ResU{J-1}*resh(Nu(J)+1:2*Nu(J)); 
                ResP{J-1}*resh(2*Nu(J)+1:end)];
        
        % Cycle
        Nu(J-1) = size(A{J-1},1)/2; Np(J-1) = size(B{J-1},1);
        bigeH = zeros(2*Nu(J-1)+Np(J-1),1);
        bigeH = Vcycle(bigeH, resH, J-1);
%         figure, trisurf(elem{J-1}, node{J-1}(:,1), node{J-1}(:,2), bigeH(1:Np(J-1))); shading interp; colorbar; view(2);

        % Prolongate the correction to the fine grid
        bigeh = [ProU{J-1}*bigeH(1:Nu(J-1)); 
                 ProU{J-1}*bigeH(Nu(J-1)+1:2*Nu(J-1)); 
                 ProP{J-1}*bigeH(2*Nu(J-1)+1:end)];
        bige = bige + bigeh;
        
        % Post-smooth
        for i = 1:s2
            bige = LSC_DGS_rotation_smoother(bige, A{J}, B{J}, rhs, Sp{J}, node{J}, elem{J});
        end
    end
end