function [L2errPGamma, H1errPGamma] = getErrorFracturePart_2(pGammah,pGamma,pGammax,gcrl, gprl,Tb,basis_type)

L2errPGamma = 0; H1errPGamma = 0; h = 2/size(Tb,2);

for k = 1:size(Tb,2)
    beginPoint = [-1+(k-1)*h, 0]; endPoint = [-1+k*h, 0];
    pGammahLocal = pGammah(Tb(:,k));
    [gclx,gplx] = generate_Gauss_local_1D(gcrl, gprl, beginPoint, endPoint, 0);
    for i = 1:length(gclx)
        L2errPGamma = L2errPGamma + gclx(i)*(pGamma(gplx(i)) - ...
                              feSolution1D(gplx(i), 0, basis_type, beginPoint, endPoint, pGammahLocal))^2;
        H1errPGamma = H1errPGamma + gclx(i)*(pGammax(gplx(i)) - ...
                              feSolution1D(gplx(i), 1, basis_type, beginPoint, endPoint, pGammahLocal))^2;
    end
end
L2errPGamma = sqrt(L2errPGamma); H1errPGamma = sqrt(H1errPGamma);

end