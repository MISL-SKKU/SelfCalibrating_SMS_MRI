function res=gnrHSL_fast(kdataMB,Vs,Nsc,kdataSMS3D,mask,param);

[sx,sy,nc] = size(kdataMB);
% mb = size(kdataSMS3D,4);
kSize=param.kSize; iter=param.iter;
mb=param.mb; ppa = param.ppa;
vs = mb*ppa;
tmp = im2row(kdataMB,kSize); [tsx,tsy,tsz] = size(tmp);
tmpA = reshape(tmp,tsx,tsy*tsz); A = tmpA;
mu=1/1000000.0; mu2=0.0000000000001;
% mu=1/1000000000.0; mu2=0.0000000000001;

ind = (squeeze(kdataSMS3D(:,10,1,:))~=0);

%% calculating pre-defined function
disp(sprintf('Calculating pre-defined functions'));
for m=1:mb
    NNH(:,:,m)=Nsc{1,m}*Nsc{1,m}';
    invNNH(:,:,m)=inv(mu2.*eye(size(NNH(:,:,m),2))+NNH(:,:,m));
    invVHV(:,:,m)=inv(mu.*eye(size(Vs{1,m},2))+Vs{1,m}'*Vs{1,m});
    ANNH{1,m}=A*NNH(:,:,m); VVHV{1,m}=Vs{1,m}*invVHV(:,:,m);
    VNNH{1,m}=Vs{1,m}'*invNNH(:,:,m); resA_Base{1,m}=ANNH{1,m}*invNNH(:,:,m);
end
clear Nsc NNH invNNH invVHV ANNH
%% Aliasing Sepration: SMS-HSL
for l=1:iter
    if (l==iter) disp(sprintf('Aliasing Separation: # of iterations=%d',l)); end
    for m=1:mb
        if l==1
            resA=resA_Base{1,m};
        else
            resA=resA_Base{1,m}+mu2.*resU(:,:,m)*VNNH{1,m};
        end
        res(:,:,:,m) = row2im(reshape(resA,tsx,tsy,tsz),[sx,sy,nc],kSize);
        
    end
    % hybrid data fidelity
    res(:,:,:,vs)=0;
    tmp = shift_by_mask(res,mask,0);
    
    mod_res=fft(repmat(tmp(:,:,:,1:mb),1,1,1,ppa),[],4);
    
% mod_res=fftshift(fft(shift_cycles(res,caipi,-1*cycles,ctr),[],4),4);
    for m=1:vs
        mod_res(ind(:,m)~=0,:,:,m)=kdataSMS3D(ind(:,m)~=0,:,:,m);
        mod_res(:,:,:,2:ppa:end)=0;

    end
    res = shift_by_mask(ifft(mod_res,[],4),mask,1);
    res = res(:,:,:,1:mb);
%         res = shift_cycles(ifft(ifftshift(mod_res,4),[],4),caipi,cycles,ctr);

    for m=1:mb
        tmp = im2row(res(:,:,:,m),kSize); [tsx,tsy,tsz] = size(tmp);
        resA = reshape(tmp,tsx,tsy*tsz);
        resU(:,:,m)=resA*VVHV{1,m};
    end
end

end
