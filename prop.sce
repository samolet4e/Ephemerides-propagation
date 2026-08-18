CL_init;
// =====================================================
// Satellite ephemeris
// =====================================================

// Date/time of orbital elements (TREF)

UTC_offset = 3; // Sofia, BG = UTC + 3; You must enter the accurate UTC offset of your location
t = datevec(now());
t(4) = t(4) - UTC_offset;
cjd0 =  CL_dat_convert('cal', 'cjd', t');// - 1/12/2;
cjd = cjd0 + (0 : 1/1440 : 1/6); // 4 hours

// Balkan-2
// https://www.n2yo.com/satellite/?s=69024
str = [..
'1 69024U 26100AW  26229.50886231  .00003134  00000-0  14604-3 0  9993'; ..
'2 69024  97.3977 126.2108 0004478 125.0504 235.1152 15.20618553 16120'];

tle = CL_tle_parse(str);
[pos_ecf, vel_ecf] = CL_tle_genEphem(tle, cjd, 'ECF', %CL_UT1_TREF, %CL_TT_TREF);

// =====================================================
// Ground stations geodetic coordinates
// =====================================================
//       East lon, rad    North lat, rad.. Altitude,m
Sofia = [23.3219*%pi/180; 42.6977*%pi/180; 500];
Nuremberg = [11.0767*%pi/180; 49.4521*%pi/180; 302];

// Earth->Sun and Earth->Moon in ECI 
Sun_eci = CL_eph_sun(cjd); 
Moon_eci = CL_eph_moon(cjd);

// =====================================================
// Check Sun-synchronicity
// =====================================================
// Mean local time of ascending node
[pos_eci, vel_eci] = CL_fr_convert('ECF','ECI', cjd, pos_ecf, vel_ecf);
kep = CL_oe_car2kep(pos_eci, vel_eci); // [sma,ecc,inc,aper,raan,ma] [m,rad]
mlh = CL_op_locTime(cjd, 'ra', kep(5,:), 'mlh');

// Plot
scf();
plot(cjd - cjd0, mlh, 'linewidth', 2); 
xtitle('','Days','Mean local time of ascending node, hours');
CL_g_stdaxes();
[y] = settings(gca());

// =====================================================
// Ground stations visibility
// =====================================================
// Min elevation for visibility
min_elev = 10*%pi/180;

// Computation of visibility intervals
tvisi1 = CL_ev_stationVisibility(cjd, pos_ecf, Sofia, min_elev);

// Plot 
scf();
plot2d3(tvisi1(1,: )- cjd0, (tvisi1(2,:) - tvisi1(1,:))*1440, style=2);
xtitle('','Days','Visibility duration, min');
h = CL_g_select(gca(), 'Polyline');
h.thickness = 2;
CL_g_stdaxes();
[y] = settings(gca());
//CL_g_legend(gca(), 'Sofia');//["sta1", "sta2"]);

// Plot East azimuth <-> elevation for Sofia
scf();
for k = 1:size(tvisi1,2)
  tk = linspace(tvisi1(1,k), tvisi1(2,k), 100);
  [az, el] = CL_gm_stationPointing(Sofia, CL_interpLagrange(cjd, pos_ecf, tk));
  plot(-az*180/%pi, el*180/%pi);
end
xtitle('','East azimuth, deg','Satellite elevation from Sofia, deg');
CL_g_stdaxes();
[y] = settings(gca());

// =====================================================
// (geocentric) Longitude/latitude plot
// =====================================================
scf(); 
// Plot Earth map
CL_plot_earthMap(color_id = color('seagreen'), res = 'high', thickness = 1);
// Plot ground tracks
CL_plot_ephem(pos_ecf,color_id = 2, thickness = 2);

// Plot visibility circles
// min_elev: min elevation for visibility
// rmin/rmax: min/max orbital radius (from Earth center)
min_elev = 10*%CL_deg2rad;
rmin = min(CL_norm(pos_ecf));
rmax = max(CL_norm(pos_ecf));
az = linspace(0,2*%pi,100);

CL_plot_ephem(CL_gm_stationVisiLocus(Sofia, az, min_elev, rmin), color_id = 2);
CL_plot_ephem(CL_gm_stationVisiLocus(Sofia, az, min_elev, rmax), color_id = 2);
xtitle('','Longitude, deg','Latitude, deg');
[y] = settings(gca());

// Plot sub-satellite point
ecef.x = interp1(cjd, pos_ecf(1,:), cjd0, 'linear');
ecef.y = interp1(cjd, pos_ecf(2,:), cjd0, 'linear');
ecef.z = interp1(cjd, pos_ecf(3,:), cjd0, 'linear');

//[geoSat] = wgs84EcefToGeo(ecef);
//[geoSat] = ecef2geodetic(ecef);
geo = CL_co_car2ell([ecef.x,ecef.y,ecef.z]', %CL_eqRad, %CL_obla);
geo(1) = modulo(geo(1) + %pi, 2*%pi) - %pi; // wrap longitude into [-pi, pi)
geoSat.lon = geo(1)*%CL_rad2deg; geoSat.lat = geo(2)*%CL_rad2deg; geoSat.alt = geo(3);

set(gca(), 'auto_clear', 'off'); // Matlab/Octave 'hold on'
scatter(geoSat.lon, geoSat.lat, 64, 'marker', 11, 'markerEdgeColor', 'red', 'markerBackground', 'yellow');

// =====================================================
// Sun and Moon directions in satellite frame
// Satellite frame supposed to be "qsw"
// =====================================================
M_eci2sat = CL_fr_qswMat(pos_eci, vel_eci);

// Satellite->Sun and satellite->Moon directions
Sun_dir  = M_eci2sat*CL_unitVector(Sun_eci  - pos_eci);
Moon_dir = M_eci2sat*CL_unitVector(Moon_eci - pos_eci);

// Plot angles
f = scf();
plot(cjd - cjd0, CL_vectAngle(Sun_dir, [0;0;-1])*180/%pi, 'r', 'thickness', 2);
plot(cjd - cjd0, CL_vectAngle(Moon_dir, [0;0;1])*180/%pi, 'b', 'thickness', 2);
xtitle('','Days','Angle with satellite frame axis, deg'); 
CL_g_legend(gca(), ['Sun <-> -Z', 'Moon <-> Z']);
CL_g_stdaxes();
[y] = settings(gca());

// =====================================================
// Eclipse periods of Sun by Earth
// =====================================================
// Eclipse intervals (umbra)
interv = CL_ev_eclipse(cjd, pos_eci, Sun_eci, typ = "umb");

dur = (interv(2,:) - interv(1,:))*1440; // min
x = [interv(1,:); interv; interv(2,:); %nan*interv(1,:)] - cjd0; 
y = [zeros(dur); dur; dur; zeros(dur); %nan*ones(dur)]; 

// Plot
scf();
plot(x(:)', y(:)', 'thickness', 2);
xtitle('','Days','Time in Earth shadow, min');
CL_g_stdaxes();
[y] = settings(gca());
