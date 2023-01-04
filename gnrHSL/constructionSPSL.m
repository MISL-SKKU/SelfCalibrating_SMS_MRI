function [Vs Nsc] = constructionSPSL(kdata,param,dip);
if nargin<3
    dip=0;
end


[sx,sy,nc,vs]=size(kdata);
mb=param.mb;
Vs=cell(1,mb); Nsc=cell(1,mb);
kSize=param.kSize; r=param.r; nul=param.nul;

if length(nul)<mb
    nul=repmat(nul(1),1,mb);
end

for m=1:mb
    
    tmp = im2row(kdata(:,:,:,m),kSize); [tsx,tsy,tsz] = size(tmp);
    A = reshape(tmp,tsx,tsy*tsz); [U,S,V] = svd(A,'econ');
    if (dip==1)
%     figure(101),subplot(2,ceil(mb/2),m), plot(diag(S))
    figure(700+m),plot(diag(S))
    end
    Vs{1,m}=V(:,1:r);

    tmpkdata=zeros(sx,sy,nc);
    tmpsdata=[];

    for n=1:vs
        if n ~= m
            tmpkdata = tmpkdata+kdata(:,:,:,n); 
%             
            tmp = im2row(kdata(:,:,:,n),kSize); [tsx,tsy,tsz] = size(tmp);
            B = reshape(tmp,tsx,tsy*tsz);
            tmpsdata = [tmpsdata; B]; 
        end
    end    
    tmp = im2row(tmpkdata(:,:,:)/vs,kSize); [tsx,tsy,tsz] = size(tmp);
    C= reshape(tmp,tsx,tsy*tsz); 
    
    X = [tmpsdata;C];
    
    
    
    
    [U,S,V] = svd(X,'econ');
    if (dip==1)
%         figure(111),subplot(2,ceil(mb/2),m), plot(diag(S))
    figure(710+m),plot(diag(S))

    end
    Nsc{1,m}=V(:,nul(m)+1:end);
end

end