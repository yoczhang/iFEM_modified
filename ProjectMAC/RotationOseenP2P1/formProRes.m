function [ProP, ResP] = formProRes(Nc, Nf, HB, N2c, N2f, edgef)

% [ProU, ResU, ProP, ResP] = formProRes(Nc, Nf, HB, N2c, N2f, edgef)
% nodec, elemc, elem2dofc, edgec, bdDofc
% nodef, elemf, elem2doff, edgef, bdDoff

% nodec = [0,0; 1,0; 1,1; 0,1]; elemc = [2,3,1; 4,1,3];
% [elem2dofc, edgec, ~] = dofP2(elemc); 
% [nodef,elemf,~,HB] = uniformrefine(nodec,elemc);
% [elem2doff, edgef, ~] = dofP2(elemf); 
% 
% Nc = size(nodec,1); Nf = size(nodef,1);
% N2c = Nc + size(edgec,1); N2f = Nf + size(edgef,1);
%% P1
Pro1 = speye(Nc);

i = double([HB(:,1); HB(:,1)]);
j = double([HB(:,2); HB(:,3)]);
s = 0.5*ones(2*size(HB,1),1);
Pro2 = sparse(i, j, s, Nf, Nc);

ProP = [Pro1; Pro2(Nc+1:end,:)]; 
% ResP = sparse((1:Nc)', (1:Nc)', ones(Nc,1), Nc, Nf);
ResP = ProP';
temp = sum(ResP,2); temp = repmat(temp, 1, Nf);
ResP = ResP./temp; 

% ProP = ProP + zeros(Nf, Nc);
% ResP = ResP + zeros(Nc, Nf);
clear Pro1 i j s Pro2 temp

%% P2
% Pro1 = speye(N2c);
% 
% ii = double([(N2c+1:N2f)'; (N2c+1:N2f)']);
% jj = double([edgef(:,1); edgef(:,2)]);
% ss = 0.5*ones(2*(N2f-N2c),1);
% Pro2 = sparse(ii, jj, ss, N2f, N2c);
% 
% ProU = [Pro1; Pro2(N2c+1:N2f,:)]; 
% % ResU = sparse((1:N2c)', (1:N2c)', ones(N2c,1), N2c, N2f); 
% ResU = ProU';
% temp = sum(ResU,2); temp = repmat(temp, 1, N2f);
% ResU = ResU./temp;
% 
% % ProU = ProU + zeros(N2f, N2c);
% % ResU = ResU + zeros(N2c, N2f);
% 
% clear ii jj ss Pro1 Pro2 temp

end