
function[effects]=MResmodel_for_GSA(Option,EnPriceScenario,N)

%MResmodel: calculates first and total order effects for 9 variables used
%to calculate returns for Clients and ESCOs.
%inputs  - Option: sets case study project (ProjectA, ProjectB, ProjectC or ProjectD)
%        - EnPriceScenario: sets energy price scenario 1= High, 2 = Low, 3
%           = Reference
%        - N: number of runs
%outputs - effects: table setting out first order and total effects
%contributions to variance in outcomes for ESCOs and Clients

load 'inputsforlhsusamplegen'


%% Set Fixed Inputs
ConSig = datenum(2014,08,01); % date of contract signature
Inflation = 0.02;    % set inflation rate (fixed over model period)
BaseDate = datenum (2014,04,01); % model base date for indexation
PC = datenum (2014,11,01); % Practical Completion Date (works complete)
ProcStart = datenum (2014,01,01); % procurement start date
RoR = 0.05; % rate of return used to calculate client and ESCO NPVs
Modelperiod = 25*365; % model period set to 25 years, NPV is calculated over this period
Modelend = PC+Modelperiod; % model end date
k = 9; %number of variables

% set up vectors for dates
indata = datenum(ProcStart):datenum(Modelend);
dvec   = datevec(indata);
duniq  = unique(dvec(:, 1:2), 'rows');
cashflowdates = datenum(duniq(:,1), duniq(:,2), 1); %vector of 1st of month dates for period to be considered
modelmonths = months(ProcStart, Modelend, 1)+1 ;% adjusted to match no. months in cashflow dates



    % Get Energy Prices
    %set dates for which prices are needed
    aa = cashflowdates(1);
    ab = cashflowdates(end);
    ac = EnergyPrices(:,1) >= aa & EnergyPrices(:,1) <=ab;

    %create vectors of prices
    ElecPriceHigh = EnergyPrices(:,3) .* ac;
    ElecPriceLow = EnergyPrices (:,5) .* ac;
    ElecPriceRef = EnergyPrices(:,7).* ac;
    GasPriceHigh = EnergyPrices(:,2) .* ac;
    GasPriceLow = EnergyPrices(:,4) .* ac;
    GasPriceRef = EnergyPrices (:,6) .* ac;

    %remove non-zero elements from price vectors
    ElecPriceHigh(ElecPriceHigh==0) = [];
    ElecPriceLow(ElecPriceLow==0) = [];
    ElecPriceRef(ElecPriceRef==0) = [];
    GasPriceHigh(GasPriceHigh==0) = [];
    GasPriceLow(GasPriceLow==0) = [];
    GasPriceRef(GasPriceRef==0) = [];


    a = ProcStart.*ones(modelmonths,1);
    b = ConSig.*ones(modelmonths,1);
%% generate samples     
%   generation of uniform samples in (0;1)
  for i=1:N
   
    T(i,:)= LPTAU51(i,2*k);
   
  end
   
       
        
%   transformation from uniform in (0,1) to actual distribution (all
%    distributions assumed to be uniform)
       transvec1 = Option(2,:)- Option(1,:);
       col1 = transvec1(1)* T(:,1)+ Option(1,1);
       col2 = transvec1(2)* T(:,2)+ Option(1,2);
       col3 = transvec1(3)* T(:,3)+ Option(1,3);
       col4 = transvec1(4)* T(:,4)+ Option(1,4);
       col5 = transvec1(5)* T(:,5)+ Option(1,5);
       col6 = transvec1(6)* T(:,6)+ Option(1,6);
       col7 = transvec1(7)* T(:,7)+ Option(1,7);
       col8 = transvec1(8)* T(:,8)+ Option(1,8);
       col9 = transvec1(9)* T(:,9)+ Option(1,9);
       col10 = transvec1(1)* T(:,10)+ Option(1,1);
       col11 = transvec1(2)* T(:,11)+ Option(1,2);
       col12 = transvec1(3)* T(:,12)+ Option(1,3);
       col13 = transvec1(4)* T(:,13)+ Option(1,4);
       col14 = transvec1(5)* T(:,14)+ Option(1,5);
       col15 = transvec1(6)* T(:,15)+ Option(1,6);
       col16 = transvec1(7)* T(:,16)+ Option(1,7);
       col17 = transvec1(8)* T(:,17)+ Option(1,8);
       col18 = transvec1(9)* T(:,18)+ Option(1,9);
 
       A=[col1,col2,col3,col4,col5,col6,col7,col8,col9];
       B=[col10,col11,col12,col13,col14,col15,col16,col17,col18];
       
%   PREPARATION OF THE RADIAL SAMPLE MATRIX X(k+2,k) 
      
C = [];
for i = (1:k);
    ab = A;
    ab(:,i)=B(:,i);
    C(:,end+1:end+k)=ab;
end
inputs = [A,B,C];
%%  Run Model
for i = 1:k+2; % this is intended to loop through model for k+2 sets of input values

ERC = inputs(:,9*(i-1)+2);
EIGPC = inputs(:,9*(i-1)+3);
EIC=inputs(:,9*(i-1)+4);
EM=inputs(:,9*(i-1)+5);
CTC = inputs(:,9*(i-1)+1);
GS = inputs(:,9*(i-1)+6);
ES = inputs(:,9*(i-1)+7);
GGS = inputs(:,9*(i-1)+8);
GES = inputs(:,9*(i-1)+9);

% initial calculations

%calculate guarantee period based on mean savings and prices current at
% start of procurement.  Ignore transaction costs in calculation

%checkcalc1 = mean(ES)
%checkcalc2 = ElecPriceRef(1,1)
%checkcalc3 = mean(GS)
%checkcalc4 = GasPriceRef(1,1)
Expectedannualsaving = mean(ES)*ElecPriceRef(1,1) +mean(GS)*GasPriceRef(1,1);
guaranteeperiod = mean(EIC)/Expectedannualsaving;
guaranteeend = PC + (365*guaranteeperiod);
%datestr(guaranteeend)%check end of guarantee

% calculate ESCO Cash Flows

% procurement cost condition
c = cashflowdates>=a;% date is after procurement start
d = cashflowdates<b;% date is before contract signature
e = c.*d;
% procurement cost
ESCOProcurementCost = e*((-ERC'-EIGPC')/months(ProcStart,ConSig,1));

%ESCOProccostcheck = sum(ESCOProcurementCost(:,10000))

% financial close bullet payment
f = (cashflowdates)==b; 
FCbullet = f*(ERC'+EIGPC');

%FCbulletcheck = sum(FCbullet(:,10000))

%installation costs
g = PC.*ones(modelmonths,1);
h = cashflowdates>=b; % date is after contract signature
j = cashflowdates<g; % date is before practical completion
kk =h.*j;
ESCOinstallationcosts = kk*(-EIC')/(months(ConSig,PC,1));
%ESCOinstallcheck = sum(ESCOinstallationcosts(:,10000))
%ESCOinstallcheck2 = EIC(1,10000)

% Practical Completion Bullet Payment
l = (cashflowdates)==g; 
PCbullet = l*(EIC'.*(1+EM'));
%PCbulletcheck = sum(PCbullet(:,10000))

% guarantee period calcs
q = cashflowdates>=PC;
ah = cashflowdates<guaranteeend;% date is before guaranteeend
aj = ah.*q;
diag_aj = diag(aj);

% Electricity Shortfall
m = ES' - GES'; % calculate difference
n = ES'<GES'; % only use if negative difference
p = 1/12*(n.*m); % need monthly amount
y = q; % applies from PC
Elecshortfall = y*p;
ESCOelecshortfall =diag_aj* Elecshortfall;% only applies before guarantee end

% Gas Shortfall
s = GS' - GGS'; %calculate difference
t = GS' < GGS';% only applies if negative
u = 1/12*(t.*s);% need monthly amount

Gasshortfall = y*u;% only applies from PC
ESCOgasshortfall = diag_aj*Gasshortfall; % only applies before guarantee end

% calculate Client Cash Flows 

% procurement cost
ClientProcurementCost = e*((-CTC')*30/(ConSig-ProcStart));
%ClientProcCostcheck = sum(ClientProcurementCost(:,1))

% Electricity Shortfall
ClientElecshortfall = -ESCOelecshortfall;


% Gas Shortfall
ClientGasshortfall = -ESCOgasshortfall;

% Calculate Energy costs
% Convert Energy Price vectors to diagonal matrices

EPH = diag (ElecPriceHigh);
EPL = diag (ElecPriceLow);
EPR = diag (ElecPriceRef);
GPH = diag (GasPriceHigh);
GPL = diag (GasPriceLow);
GPR = diag (GasPriceRef);

if EnPriceScenario == 1;
    EP = EPH;
    GP = GPH;
else if EnPriceScenario == 2;
        EP = EPL;
        GP = GPL;
    else if EnPriceScenario == 3;
            EP=EPR;
            GP=GPR;
    end
    end
end
        

ESCOelecshortfallcost = EP * ESCOelecshortfall;
ESCOGasshortfallcost = GP * ESCOgasshortfall;

Clientelecshortfallcost = EP *ClientElecshortfall;
ClientGasshortfallcost = GP* ClientGasshortfall;

ClientElecSaving = EP* 1/12*(y*ES' ) ; %client receives all of actual saving
ClientGasSaving = GP * 1/12*(y*GS' );

% Cashflow

ESCOCashflowreal = ESCOProcurementCost + FCbullet + ESCOinstallationcosts +...
    PCbullet + ESCOelecshortfallcost + ESCOGasshortfallcost;

%checkescoinstal = sum(EIC(:,10000))
%checkPCbullet = sum(PCbullet(:,10000))
%Checkmargin = EM(1,10000)

ClientCashflowreal = ClientProcurementCost - FCbullet - PCbullet + ...
    Clientelecshortfallcost + ClientGasshortfallcost + ClientElecSaving + ...
    ClientGasSaving;

% apply inflation
% define indexation period
indexationperiod = floor((cashflowdates - BaseDate)/365);

% calculate inflation factor
ad = (1 + Inflation).^indexationperiod;
add = diag(ad); %% transform vector to diagonal matrix

% nominal cashflows
ESCOCashflownom = add * ESCOCashflowreal;
ClientCashflownom = add* ClientCashflowreal;

% Calculate NPVs

ESCONPV(:,i)= pvvar(ESCOCashflownom,(1/12*RoR));
ClientNPV(:,i) = pvvar(ClientCashflownom, (1/12*RoR));

end

%% Variance calculations
        yAESCO=ESCONPV(:,1);
        yBESCO=ESCONPV(:,2);
        yAbESCO=ESCONPV(:,3:end);
        yAClient=ClientNPV(:,1);
        yBClient=ClientNPV(:,2);
        yAbClient=ClientNPV(:,3:end);
        nn=1/N;
        f02ESCO = (nn*sum(yBESCO))^2;
        f02Client = (nn*sum(yBClient))^2;
        yBESCO2 = nn*sum(yBESCO.*yBESCO);
        yBClient2 = nn*sum(yBClient.*yBClient);
      
        
 for j= 1:k
     
     yaycESCO(:,j) = sum(yAESCO.*yAbESCO(:,j));
     ybycESCO(:,j) = sum(yBESCO.*yAbESCO(:,j));
     yaycClient(:,j) = sum(yAClient.*yAbClient(:,j));
     ybycClient (:,j)= sum(yBClient.*yAbClient(:,j));
     
 end 
     siESCO = ((nn*ybycESCO)-f02ESCO)/(yBESCO2-f02ESCO);
     sTESCO = 1-((nn*yaycESCO)-f02ESCO)/(yBESCO2-f02ESCO);
     siClient = ((nn* ybycClient)-f02Client)/(yBClient2-f02Client);
     sTClient = 1-((nn*yaycClient)-f02Client)/(yBClient2-f02Client);
 
 
 
 
%% create results table

variables = {'CTC','ERC','EIGPC','EIC','EM','GS','ES','GGS','GES'};
FirstorderESCO = siESCO';
TotaleffectsESCO = sTESCO';
FirstorderClient = siClient';
TotaleffectsClient= sTClient';
effects = table(FirstorderESCO, TotaleffectsESCO, FirstorderClient,TotaleffectsClient,...
    'RowNames',variables);

 
end
 
        
