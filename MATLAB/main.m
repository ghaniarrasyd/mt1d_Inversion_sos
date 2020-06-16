%Magnetotelluric 1-D inversion using original SOS algorithm(Cheng & Prayogo)
clear all; clc;
%%
% Create synthetic data
% In this example I used 3 layers model with parameter stated below
rsin = [100 10 1000]; %resistivities of each layer in synthetic model
tsin = [500 1500]; %thickness of each layer in synthetic model
freq = logspace(-3, 3, 50); %frequency from 10^-3 - 10^3
nsin = length(rsin); %number of layers
[rho_sin,dsin] = plot_rho_h(rsin,tsin);

% Generate synthetic data
[app_data, phase_data] = FWDMT1D(rsin, tsin, freq);

%%
% Model space
npop=50; %number of populations/models
nlayer=3; %number of predicted layers
ngen=50; %number of iteration

rmin=1; %minimum value of resistivity
rmax=2000; %maximum value of resistivity
tmin=200; %minimum value of thickness
tmax=2000; %maximum value of thickness

% Generate initial model
for ipop=1:npop
    for imod=1:nlayer
        rho(ipop,imod)=rmin+rand*(rmax-rmin);
    end
    for imod=1:nlayer-1
        thick(ipop,imod)=tmin+rand*(tmax-tmin);
    end
end

% Calculate apparent resistivity and phase for each model, and generate
% misfit of each model
for ipop=1:npop
    [apparentResistivity, phase]=FWDMT1D(rho(ipop,:),thick(ipop,:),freq);
     app_mod(ipop,:)=apparentResistivity;
     phase_mod(ipop,:)=phase;
     
    [misfit]=misfitMT(app_data,phase_data,app_mod(ipop,:),phase_mod(ipop,:));
    E(ipop)=misfit;
end

%%
% Inversion process
% We generate the inversion for 50 iterations
for igen=1:ngen % loop over iteration
for i=1:npop % loop over model
% Identify best model from the whole model populations
idx=find(E==min(E));
rbest=rho(idx(1),:);
tbest=thick(idx(1),:);

% First phase, Mutualisms
j=randi(npop,1);
k=randi(npop,1);
while j==i || k==i
    j=randi(npop,1);
    k=randi(npop,1);
end
    % Model tested
    rmut=[rho(i,:);rho(j,:)];
    tmut=[thick(i,:);thick(j,:)];
    % Mutual vector
    for n=1:nlayer
        mv_r(n)=(rho(i,n)+rho(j,n))/2;
    end
    for n=1:nlayer-1
        mv_t(n)=(thick(i,n)+thick(j,n))/2;
    end
    % Benefit factor
    bf=1;
    
    % Calculation of new model with mutualisms proscedure
    for l=1:2
        for n=1:nlayer
            rum(l,n)=rmut(l,n)+rand*(rho(k,n)-mv_r(n)*bf);
            if rum(l,n)<rmin
                rum(l,n)=rmin;
            end
            if rum(l,n)>rmax
                rum(l,n)=rmax;
            end
        end
        for n=1:nlayer-1
            tum(l,n)=tmut(l,n)+rand*(thick(k,n)-mv_t(n)*bf);
            if tum(l,n)<tmin
                tum(l,n)=tmin;
            end
            if tum(l,n)>tmax
                tum(l,n)=tmax;
            end
        end
    end

    % Data calculation and misfit
    for l=1:2
        [app_um, phase_um]=FWDMT1D(rum(l,:),tum(l,:),freq);     
        [error]=misfitMT(app_data,phase_data,app_um,phase_um);
        Em(l)=error;
        
        % Update model if misfit of new model better than previous one
        if l==1
            if Em(l)<E(i)
                rho(i,:)=rum(l,:);
                thick(i,:)=tum(l,:);
                app_mod(i,:)=app_um;
                phase_mod(i,:)=phase_um;
                E(i)=Em(l);
            end
        else
            if Em(l)<E(j)
                rho(j,:)=rum(l,:);
                thick(j,:)=tum(l,:);
                app_mod(j,:)=app_um;
                phase_mod(j,:)=phase_um;
                E(j)=Em(l);
            end
        end
    end

% Second phase: Commensalisms
j=randi(npop,1);
while j==i
    j=randi(npop,1);
end
    % Calculation of new model with commensalisms proscedure
    for n=1:nlayer
        rcom(1,n)=rho(i,n)+(0.4+0.5*rand)*(rbest(n)-rho(j,n));
        if rcom(1,n)<rmin
            rcom(1,n)=rmin;
        end
        if rcom(1,n)>rmax
            rcom(1,n)=rmax;
        end
    end
    for n=1:nlayer-1
        tcom(1,n)=thick(i,n)+(0.4+0.5*rand)*(tbest(n)-rho(j,n));
        if tcom(1,n)<tmin
            tcom(1,n)=tmin;
        end
        if tcom(1,n)>tmax
            tcom(1,n)=tmax;
        end
    end
    
    % Data calculation and misfit
    [app_com, phase_com]=FWDMT1D(rcom(1,:),tcom(1,:),freq);
    [Ec]=misfitMT(app_data,phase_data,app_com,phase_com);
    
    % Update model if misfit of new model better than previous one
    if Ec<E(i)
        rho(i,:)=rcom(1,:);
        thick(i,:)=tcom(1,:);
        app_mod(i,:)=app_com;
        phase_mod(i,:)=phase_com;
        E(i)=Ec;
    end
    
%Third phase, parasitism
j=randi(npop,1);
while j==i
    j=randi(npop,1);
end
    % Calculation of parasite model with parasitism proscedure
    rpar=rho(i,:);
    tpar=thick(i,:);
    p1=randi(2,1);
    if p1==1
        p2=randi(nlayer,1);
        rpar(p2)=rmin+rand*(rmax-rmin);
    else
        p3=randi(nlayer-1,1);
        tpar(p3)=tmin+rand*(tmax-tmin);
    end
    
    % Data calculation and misfit
    [app_par, phase_par]=FWDMT1D(rpar,tpar,freq);
    [Ep]=misfitMT(app_data,phase_data,app_par,phase_par);
    
    % Update model if misfit of new model better than previous one
    if Ep<E(j)
        rho(j,:)=rpar(1,:);
        thick(j,:)=tpar(1,:);
        app_mod(j,:)=app_par;
        phase_mod(j,:)=phase_par;
        E(j)=Ep;
    end    
end

%%
% Update best model for each iterations
Emin=100;
for ipop=1:npop
    if E(ipop)<Emin
        Emin=E(ipop);
        rmodel=rho(ipop,:);
        tmodel=thick(ipop,:);
        app_best=app_mod(ipop,:);
        phase_best=phase_mod(ipop,:);
    end
end
% Best misfit for each iterations
Egen(igen)=Emin;

end

%%
% Visualizations
[rho_sin,dsin]=plot_rho_h(rsin,tsin);
[rho_mod,dmod]=plot_rho_h(rmodel,tmodel);
dsin(6)=dmod(2*nlayer);

% Plot model
figure (1)
subplot(2, 2, 1)
loglog(1./freq,app_data,'.b',1./freq,app_best,'r','MarkerSize',12,'LineWidth',2.5);
axis([10^-3 10^3 1 10^3]);
legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold');
xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
ylabel('App. Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');

subplot(2, 2, 3)
loglog(1./freq,phase_data,'.b',1./freq,phase_best,'r','MarkerSize',12,'LineWidth',2.5);
axis([10^-3 10^3 0 90]);
set(gca,'YTick',[0:15:90]);
set(gca, 'YScale', 'linear');
legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold');
xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
ylabel('Phase (deg)','FontSize',12,'FontWeight','Bold');

subplot(2, 2, [2 4])
stairs(rho_sin,dsin,'--b','Linewidth',2);
hold on
stairs(rho_mod,dmod,'-r','Linewidth',2.5);
hold off
legend({'Synthetic Model','Calculated Model'},'EdgeColor','none','Color','none','FontWeight','Bold');
axis([1 10^4 0 5000]);
xlabel('Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',12,'FontWeight','Bold');
set(gca,'YDir','Reverse');
set(gca, 'XScale', 'log');
set(gcf, 'Position', get(0, 'Screensize'));

% Plot misfit
figure (2)
plot(1:ngen, Egen);
xlabel('Number of iterations');
ylabel('Misfit');

%% Save the best model into external files
cal=[freq; app_best; phase_best; Egen];
best=[rmodel tmodel];
dlmwrite('SOScal.txt',cal,'delimiter','\t');
dlmwrite('SOSmod.txt',best,'delimiter','\t');