 function [ n , Pn ] = ProbDensityHistV3(image)
<<<<<<< HEAD
crop=[75,235,150,150];
bgcrop = [100 100 50 50];
folder='\\Elder-pc\j\Elder Backup Raw Images\2016\2016-02\2016-02-20\';
Nimg=imcrop(image,crop);
Nimg = imresize(Nimg,1/2);
=======
crop=[70,225,120,120];
bgcrop = [100 100 50 50];
folder='\\Elder-pc\j\Elder Backup Raw Images\2016\2016-02\2016-02-20\';
Nimg=imcrop(image,crop);
Nimg = imresize(Nimg,1/2,'Method','Box');
>>>>>>> roop
figure(1)
subplot(1,2,1)
imagesc(Nimg);
axis image
colormap gray
Nimg(isnan(Nimg))=0;
Nimg(Nimg==inf)=0;
Nimg(Nimg==-inf)=0;
BGimg = imcrop(image,bgcrop);
BGimg = imresize(BGimg,1/2);
BGcount=sum(BGimg(:))/length(BGimg(:));
Nimg=real(Nimg(:))-BGcount;
N=size(Nimg,1);
Nmax=max(Nimg)*1.1; % the maximum of the grid
Nmin=0; % the minimum of the grid
Gsize=100; % size of the grid
n=linspace(Nmin,Nmax,Gsize);
dn=(Nmax-Nmin)/(Gsize-1); % interval of the grid
Pn=n*0; %P(n)



for i=1:N
    K=round(((Nimg(i))-Nmin)/dn)+1;
    if K>=0 && K<Gsize
        Pn(K+1)=Pn(K+1)+Nimg(i)/dn;
    end
end
[~,ix] = max(Pn);
n= n/n(ix);
subplot(1,2,2)
<<<<<<< HEAD
plot(n,Pn,'.','MarkerSize',20)
=======
plot(n,Pn,'k-o',...
    'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
>>>>>>> roop
xlabel('n/n_0')
ylabel('Counts')
set(gca,'FontSize',14)
xlim([0,max(n)])
<<<<<<< HEAD
=======
ylim([-max(Pn)/10 max(Pn)])
>>>>>>> roop
