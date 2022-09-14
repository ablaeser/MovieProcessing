function [R,C,reg] = RectifyStack(vol, start) 
[~,~,Nz] = size(vol); % Nx,Ny
R = zeros(Nz,1);
C = zeros(Nz,1);    
reg = zeros(size(vol));

metric = registration.metric.MeanSquares;
optimizer = registration.optimizer.RegularStepGradientDescent;
%figure('WindowState','maximized');
target = vol(:,:,start); %imgaussfilt(vol(:,:,start), 3); % 
%forward
for i = start:Nz
    source = vol(:,:,i); %imgaussfilt(vol(:,:,i), 3); % 
    tform = imregtform(source, target, 'translation', optimizer, metric);
    R(i) = tform.T(3,2); C(i) = tform.T(3,1);
    %{
    subplot(1,2,1); imshow(target, []); title(sprintf('z = %i: target',i));
    subplot(1,2,2); imshow(source, []); title(sprintf('source. [Row, Col] = [%2.1f, %2.1f]',R(i), C(i))); pause;
    %}
    target = imtranslate(source,[C(i),R(i)]);
    reg(:,:,i) = target;
end

%backwards
target = vol(:,:,start);
for i = flip(1:start)
    source = vol(:,:,i); %imgaussfilt(vol(:,:,i), 3); %
    tform = imregtform(source, target, 'translation', optimizer, metric);
    R(i) = tform.T(3,2); C(i) = tform.T(3,1);
    %{
    subplot(1,2,1); imshow(target, []); title(sprintf('z = %i: target',i));
    subplot(1,2,2); imshow(source, []); title(sprintf('source. Translation [Row, Col] = [%2.1f, %2.1f]',R(i), C(i))); pause;
    %}
    target = imtranslate(source,[C(i),R(i)]);
    reg(:,:,i) = target;
end
end