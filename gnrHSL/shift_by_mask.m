function [ kdata ] = shift_by_mask( kdata,mask,inv )
%UNTITLED11 이 함수의 요약 설명 위치
%   자세한 설명 위치
[sx,sy,nc,ns] = size(kdata);

% fft!!!!!
LPhase = fft(mask,[],2);
LAngle = angle(LPhase);

if inv==1
%     mask = flip(mask,2);
% mask = circshift(mask,2,2);
% LPhase = flipud(LPhase);
LAngle = angle(LPhase)*-1;

end

% figure,imshow(angle(LPhase),[-pi pi],'initialMagnification',300,'Colormap',jet)
%

LinPhase=[];
for m=1:ns
    LinPha=repmat(exp(-1i*LAngle(:,m)),[1 sy]);
    LinPhase(:,:,:,m)=repmat(LinPha,[1 1 nc]);
end

for m=1:ns
kdata(:,:,:,m) = squeeze(LinPhase(:,:,:,m)).*kdata(:,:,:,m);
end

end

