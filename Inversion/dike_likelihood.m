function nll = dike_likelihood(m,~,params)%[nll, dhat]
%LK modified July 2022

nChains = params.Pm.nChains;
%keyboard
%% update parameters to include the ones that vary between iterations

params.Pm.Tbackground = 10.^m(1,:);
params.Pm.twidth = 10.^m(2,:);
params.Pm.tcenter = 10.^m(3,:);
params.Pm.k=m(4,:);

if params.Pm.DataFlag == 1 || params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    params.Pm.EaZr = m(6,:);
end
if params.Pm.DataFlag == 1 || params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 6
    params.Pm.EaAp = m(5,:);
end
if params.Pm.DataFlag == 3|| params.Pm.DataFlag == 4 %|| params.Pm.DataFlag == 5
    params.Pm.EaBt=m(7,:);
end
if params.Pm.DataFlag == 5
    params.Pm.EaBt=m(6,:);
end
if params.Pm.DataFlag == 3|| params.Pm.DataFlag == 4 %|| params.Pm.DataFlag == 6
    params.Pm.AFTrmr0=m(8,:);
    %params.Pm.ZFTrmr0=m(9,:);
end
if params.Pm.DataFlag == 6
    params.Pm.AFTrmr0=m(6,:);
end
    

params.Pm.kappa = params.Pm.k/(params.Pm.cp*params.Pm.rho); %thermal diffusivity

%construct finite difference stencil, might vary between iterations
% N=params.Pm.N;

%     params.Pm.beta = params.Pm.dxdZinv.*params.Pm.k.*ones(length(params.Pm.x),1); %coefficients in transformed coordinates
% 
%     e2=(params.Pm.beta(2:N)+params.Pm.beta(3:N+1));
%     e1=(params.Pm.beta(2:N)+params.Pm.beta(1:N-1));
%     e0=-(2*params.Pm.beta(2:N)+params.Pm.beta(3:N+1)+params.Pm.beta(1:N-1));
%     Dint=1./(2.*params.Pm.dZ^2)*(diag(e1(2:N-1),-1)+diag(e2(1:N-2),1)+diag(e0,0));
%     
%     D=zeros(N+1);
% %sparse representation - doesnt speed things up unless N is large    
% %    e11=[e1(2:Pm.N-1);0];e22=[0;e1(2:Pm.N-1)];
% %    Dint2 = 1./(2.*Pm.dZ^2) .*spdiags([e11 e0 e22],-1:1,Pm.N-1,Pm.N-1);
% %     
%     D(2:N,2:N)=Dint;
%     D(2,1)=1./(2.*params.Pm.dZ^2)*(params.Pm.beta(2)+params.Pm.beta(1));
%     D(N,N+1)=1./(2.*params.Pm.dZ^2)*(params.Pm.beta(N+1)+params.Pm.beta(N));
%     
%     params.Pm.D=D;
    
%    clear D Dint e0 e1 e2 
%% execute models
%run thermal diffusion model
[~,time,theta] = diffusion_1d_dike(params.Pm);

if params.Pm.DataFlag == 1 || params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    %run chemical model for each chronometer 
    %zircon
    Zr_agemodel = thermochron_calc(params.data.ZPts,theta,params.Pm.x,time,params.Tm.D0a2_Zr,params.Tm.EaZr,params.Tm.trZ,params.Tm.tiZ,params.Tm,params.Pm,1);
end
if params.Pm.DataFlag == 1 || params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 6
    %apatite
    Ap_agemodel = thermochron_calc(params.data.APts,theta,params.Pm.x,time,params.Tm.D0a2_Ap,params.Tm.EaAp,params.Tm.trA,params.Tm.tiA,params.Tm,params.Pm,1);
    
end
if params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    %biotite
    Bt_agemodel = thermochron_calc(params.data.BPts,theta,params.Pm.x,time,params.Tm.D0a2_Bt,params.Tm.EaBt,params.Tm.trB,params.Tm.tiB,params.Tm,params.Pm,2);
end
if params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 6
    %AFT fission track age model
    AFT_agemodel = thermochron_calc(params.data.AFTPts,theta,params.Pm.x,time,params.Pm.AFTrmr0,params.Tm.EaAp,params.Tm.trAFT,params.Tm.tiAFT,params.Tm,params.Pm,3);
end
%paleomagnetic resetting: compute the distance from dike contact at which
%temperature is predicted to reach Curie. 

MaxT=max(theta,[],2);

if floor(max(MaxT))>params.Pm.Tliqd+1
    IT=find(floor(MaxT)>params.Pm.Tliqd);
    disp(['max temp is ' num2str(max(MaxT))]);
    disp(['at these spatial indices: ' num2str(IT)']);
    disp('MCMC params: ')
    disp([params.Pm.k, params.Pm.tcenter, params.Pm.twidth, params.Pm.Tbackground])
    error('Thermal diffusion model is likely unstable, ending simulation')
end
%keyboard

if params.Pm.DataFlag == 2 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    X_MGT = interp1(MaxT(params.Pm.x>params.Pm.DikeThick),params.Pm.x(params.Pm.x>params.Pm.DikeThick),params.Pm.CurieT);
end

%% calculate residuals
if  params.Pm.DataFlag == 1 || params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    nll_Z = nlpCalc(Zr_agemodel,params.data.ZrAge,params.data.Zr_covinv,1,nChains);
end
if  params.Pm.DataFlag == 1 || params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 6
    nll_A = nlpCalc(Ap_agemodel,params.data.ApAge,params.data.Ap_covinv,1,nChains);
end
if  params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    nll_B = nlpCalc(Bt_agemodel,params.data.BtAge,params.data.Bt_covinv,1,nChains);
end
if  params.Pm.DataFlag == 3 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 6
    nll_AFT = nlpCalc(AFT_agemodel,params.data.AFTAge,params.data.AFT_covinv,1,nChains);
end
if  params.Pm.DataFlag == 2 || params.Pm.DataFlag == 4 || params.Pm.DataFlag == 5
    nll_MGT = nlpCalc(X_MGT,params.data.MGTdist,params.data.MGT_covinv,1,nChains);
end

%construct log(likelihood) for different data use cases
if  params.Pm.DataFlag == 1
    nll = nll_Z + nll_A ;
elseif params.Pm.DataFlag == 2
    nll = nll_MGT;
elseif params.Pm.DataFlag == 3
    nll = nll_Z + nll_A + nll_B + nll_AFT;
elseif params.Pm.DataFlag == 4
    nll = nll_Z + nll_A + nll_B + nll_AFT + nll_MGT;
elseif params.Pm.DataFlag == 5
    nll = nll_Z + nll_B + nll_MGT; 
elseif params.Pm.DataFlag == 6
    nll = nll_A + nll_AFT;     
end
%disp(nll_G)
%disp(nll_A)

% if nargout>1
%     dhat.Ve = Ve;
%     dhat.Vg = Vg;
% end

%% Helper Functions
%==========================================================================

    function nlp = nlpCalc(dhat,d,cov_inv,gamma,nChains)
        
        d = repmat(d,1,nChains);
        residual = d - dhat;
        chisq = zeros(1,nChains);
        for cc = 1:nChains
            chisq(cc) = residual(:,cc)' * cov_inv * residual(:,cc);
        end
        nlp = length(d)*log(gamma) + 0.5./gamma.^2.*chisq;
    end



end
