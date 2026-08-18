// https://space.stackexchange.com/questions/18289/how-to-get-semi-major-axis-from-tle
n = tle.n; // rad/s
mu = %CL_mu;// 3.986004418e+14 m^3.s^-2
sma = (mu/n^2.)^(1./3.); // semi-major axis, m
ecc = tle.ecc; // eccentricity
inc = tle.inc; // inclination, rad
argp = tle.argp; // argument of perigee, rad
raan = tle.raan; // RAAN, rad
ma = tle.ma; // mean anomaly, rad

kep0 = [sma; ecc; inc; argp; raan; ma]; // Keplerian MEAN orbital elements

// Initial date (cjd, time scale: TREF)
cjd0 = tle.epoch_cjd;
// Final dates
cjd = cjd0 + (0:1:3653); // 10 years

// STELA model parameters
// http://www.astronautix.com/m/meteor-3m.html
params = CL_stela_params();
// https://www.telematik-zentrum.de/projects/qube2/
params.mass = 16; // kg, platform + payload, max
params.drag_area = 0.04; // m^2 (min 0.02 - max 0.06)

/*
Gemini: Average Density by Altitude
200 km (Very Low Earth Orbit / VLEO): ~ 2e-10 to 5e-10 kg/m^3 (High drag; rapid orbital decay)
400 km (ISS Altitude): 2e-12 to 15e-12 kg/m^3 (Moderate drag; requires periodic reboosts)
600 km: 1e10-14 to 30e10-14 kg/m^3 (Low drag; varying significantly by solar cycle phases)
800+ km: 1e-18 kg/m^3 or less (Very minimal drag; multi-century orbital lifespans)
*/

bal = CL_tle_getBalCoef(tle);
params.drag_coef = params.mass/params.drag_area*bal;

// Propagation
[kep, info] = CL_stela_extrap('kep', cjd0, kep0, cjd, params, ['m', 'i']); // [sma,ecc,inc,aper,raan,ma] [m,rad]

// Plot inclination, deg
scf();
xgrid;
plot(cjd, kep(3,:)*180/%pi);
xtitle('','Time, mjd','Inclination, deg');
[y] = settings(gca());

// Plot semimajor axis, km
scf();
xgrid;
plot(cjd, kep(1,:)/1000, 'linewidth', 2);
xtitle('','Time, mjd','Semimajor axis, km');
[y] = settings(gca());

apogee  = sma*(1 + kep(2,:)) - %CL_eqRad;
perigee = sma*(1 - kep(2,:)) - %CL_eqRad;

// Plot apogee, perigee, km
scf();
xgrid;
plot(cjd, apogee/1000, 'linewidth', 2);
set(gca(), 'auto_clear', 'off');
plot(cjd, perigee/1000, 'linewidth', 2);
xtitle('','Time, mjd','Apogee and Perigee, km');
[y] = settings(gca());

// Plot RAAN, rad
scf();
xgrid;
plot(cjd, kep(5,:)*%CL_rad2deg, 'linewidth', 2);
xtitle('','Time, mjd','RAAN, deg');
[y] = settings(gca());
