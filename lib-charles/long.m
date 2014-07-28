% Long Correction
% Charles Le Losq
% CIW Washington 2013

% Long's correction of Raman spectra and normalisation
% last rev. Oct 2010, converted to Matlab
% ensures strictly increasing values of wavenumber
% calc. e.s.e. as Long cor. norm. sqrt(n_raw) 3d output col.
% exp. format to avoid null ese.
% program long3;

%Original program in pascal J. Roux, modified for matlab C. Le Losq.


function spectre = long(data,temp,wave)


h = 6.62606896e-34;   % J.s    Plank constant
k = 1.38066e-23;      % J/K    Boltzman
c = 2.9979e8;         % m/s    Speed of light
v = wave;             % cm-1   Excitating laser line
nu0 = 1.0./v*1e9;
T = temp + 273.15;    % K temperature


x = data(:,1);
y = data(:,2);

test = find(y<=0);
test2 = find(y>0);

% Calculate the error on data as sqrt(y). If y <= 0, then error = 1.

ese = ones(1,length(y))';
ese(test2) = sqrt(y(test2));
ese(test) = 1;

% For retrieving the good errors after long correction, one simple way is
% to work with %...

error = ese./y;


% then we proceed to the correction (Neuville and Mysen, 1996; Le Losq et
% al., 2012)

nu = 100.0.*x; % cm-1 -> m-1 Raman shift
rnu = nu0-nu; % nu0 is in m-1
t0 = nu0.*nu0.*nu0.*nu./rnu./rnu./rnu./rnu;
t1 = -h.*c.*nu./k./T; % c in m/s  : t1 dimensionless
t2 = 1 - exp(t1);
long = y.*t0.*t2;% pour les y
long2 = ese.*t0.*t2;% pour les erreurs

norm = max(long);
long = long./max(long);

eselong = error.*long;
spectre = [x,long,eselong];


end
