function L2errPGamma = getErrorFracturePart_P0DG_2(pGammah, pGamma, gcrl, gprl, Tb)

L2errPGamma = 0; h = 1/size(Tb,2);

for k = 1:size(Tb,2)
    beginPoint = [(k-1)*h, 0]; endPoint = [k*h, 0];
    pGammahLocal = pGammah(k);
    [gclx,gplx] = generate_Gauss_local_1D(gcrl, gprl, beginPoint, endPoint, 0);
    for i = 1:length(gclx)
        L2errPGamma = L2errPGamma + gclx(i)*(pGamma(gplx(i)) - pGammahLocal)^2;
    end
end
L2errPGamma = sqrt(L2errPGamma); 

end