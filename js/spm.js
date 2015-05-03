var spm = function() {

    var out = {};

    out.calcFusion = function( Te, ne, T_frac, He3_frac, R, voltage, thermal_eff ) {
    // Calculates the fusion parameters for the SPMX reactor for given
    // reactor parameters
    // Te = electron temperature
    // ne = electron density
    // T_frac: n_T/n_D
    // He3_frac: n_He3/n_D
    // R = radius in meters
    // effective plugging voltage in volts
    // thermal/electric efficiency


        // constants
        var e = 1.602E-19;
        var me = 9.109E-31;
        var mi = 1.67E-27;
        var mu0 = 4.0*Math.PI*1E-7;
        var eps0 = 8.85E-12;

        // plasma energy/temperature (J)
        var kT = Te*1000*1.602E-19;
        var betac = 0.5;

        // density of ion species (fuels)
        var n_D = ne/(1.0 + T_frac + 2*He3_frac);
        var n_T = T_frac*n_D;
        var n_He3 = He3_frac*n_D;
        var ni = n_D + n_T + n_He3;

        // numerical model for volume and surface area of plasma at betac
        var Vp = (0.253*betac + 0.439*betac*betac)*R*R*R;
        var Ap = (3.15 + 0.466*betac)*R*R;



        // fusion cross sections from NRL Plasma Formulary, Revised 2009 p44-45 (m^3/s)
        var sigma_DD = [1.5E-28, 5.4E-27, 1.8E-25, 1.2E-24, 5.2E-24, 2.1E-23, 4.5E-23, 8.8E-23, 1.8E-22];
        var sigma_DT = [5.5E-27, 2.6E-25, 1.3E-23, 1.1E-22, 4.2E-22, 8.7E-22, 8.5E-22, 6.3E-22, 3.7E-22];
        var sigma_DHe3 = [1.0E-32, 1.4E-29, 6.7E-27, 2.3E-25, 3.8E-24, 5.4E-23, 1.6E-22, 2.4E-22, 2.3E-22];
        var sigma_TT = [3.3E-28, 7.1E-27, 1.4E-25, 7.2E-25, 2.5E-24, 8.7E-24, 1.9E-23, 4.2E-23, 8.4E-23];
        var sigma_THe3 = [1.0E-34, 1.0E-31, 2.1E-28, 1.2E-26, 2.6E-25, 5.3E-24, 2.7E-23, 9.2E-23, 2.9E-22];

        // temperatures at which sigma is tabulated (keV)
        var Te_sigma = [ 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0];

        // avg fusion products per fusion event: [neutron, proton, D, T, He3, He4]
        var DD_product = [0.5, 0.5, 0.0, 0.5, 0.5, 0.0];
        var DT_product = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        var DHe3_product = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        var TT_product = [2.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        var THe3_product = [0.57, 0.57, 0.43, 0.0, 0.0, 1.0];

        // fusion energy per product (J): [neutron, proton, D, T, He3, He4]
        var DD_E = [2.45*1.602E-13, 3.02*1.602E-13, 0.0, 1.01*1.602E-13, 0.82*1.602E-13, 0.0];
        var DT_E = [14.1*1.602E-13, 0.0, 0.0, 0.0, 0.0, 3.5*1.602E-13];
        var DHe3_E = [0.0, 14.7*1.602E-13, 0.0, 0.0, 0.0, 3.6*1.602E-13];
        var TT_E = [11.3*1.602E-13/3.0, 0.0, 0.0, 0.0, 0.0, 11.3*1.602E-13/3.0];
        var THe3_E = [12.1*1.602E-13/3.0, 12.1*1.602E-13/3.0, 9.5*1.602E-13, 0.0, 0.0, (0.57*(12.1/3.0) + 0.43*4.8)*1.602E-13];


        var Te_1=0;
        var Te_2=0;
        var frac_1=0;
        var frac_2=0;

        // find Te range for sigma data points
        for (var i = 1; i < Te_sigma.length; i++) {
            if (Te_sigma[i] >= Te) {
                Te_1 = i-1;
                Te_2 = i;

                frac_1 = (Te_sigma[i] - Te)/(Te_sigma[i] - Te_sigma[i-1]);
                frac_2 = 1 - frac_1;

                break;
            }
        }

        if (Te_1 == 0) {
            console.log('Te not found');
            return;
        }

        // interpolate sigma data points for fusion cross section at Te
        var sigma_DD_Te = frac_1*sigma_DD[Te_1] + frac_2*sigma_DD[Te_2];
        var sigma_DT_Te = frac_1*sigma_DT[Te_1] + frac_2*sigma_DT[Te_2];
        var sigma_DHe3_Te = frac_1*sigma_DHe3[Te_1] + frac_2*sigma_DHe3[Te_2];
        var sigma_TT_Te = frac_1*sigma_TT[Te_1] + frac_2*sigma_TT[Te_2];
        var sigma_THe3_Te = frac_1*sigma_THe3[Te_1] + frac_2*sigma_THe3[Te_2];



        // plasma pressure, Z, required magnetic field
        var pressure = (ne + ni)*kT;
        var Z = ne/ni;
        var Bc = Math.sqrt(2*mu0*pressure/betac);

        // fusion rate density
        var DD_rate = 0.5*sigma_DD_Te*n_D*n_D;
        var DT_rate = sigma_DT_Te*n_D*n_T;
        var DHe3_rate = sigma_DHe3_Te*n_D*n_He3;
        var TT_rate = 0.5*sigma_TT_Te*n_T*n_T;
        var THe3_rate = sigma_THe3_Te*n_T*n_He3;

        // fusion power as neutrons
        var neutron_power = Vp*(
            DD_rate*DD_E[0]*DD_product[0] +
            DT_rate*DT_E[0]*DT_product[0] +
            DHe3_rate*DHe3_E[0]*DHe3_product[0] +
            TT_rate*TT_E[0]*TT_product[0] +
            THe3_rate*THe3_E[0]*THe3_product[0]);

        var charged_power = 0;

        for(var j = 1; j < 6; j++) {
            charged_power += DD_rate*DD_E[j]*DD_product[j] +
                DT_rate*DT_E[j]*DT_product[j] +
                DHe3_rate*DHe3_E[j]*DHe3_product[j] +
                TT_rate*TT_E[j]*TT_product[j] +
                THe3_rate*THe3_E[j]*THe3_product[j];
        }

        charged_power *= Vp;

        // total fusion power
        var total_power = neutron_power + charged_power;


        // energy loss by Bremsstrahlung radiation
        var brem_loss = Vp*4.22E-29*Math.sqrt(1.6E-19*Te*1000)*ne*(n_D + n_T + 4*n_He3);

        console.log(brem_loss);

        // plasma parameter
        var log_Lambda = 23.5 - 0.5*Math.log(1E-6*ne) + (5/4)*Math.log(Te*1000) - Math.sqrt(1E-5 + (1/16)*Math.pow(Math.log(Te*1000)-2, 2));

        // diffusion loss
        var diff_loss = Ap*611*(Math.pow(e,3))*Math.sqrt(me/mi)*(Math.pow(Z,2))*(Math.pow(ne,(3/2)))*log_Lambda/(2520*Math.sqrt(2*mu0)*(Math.pow(Math.PI,(3/2)))*(Math.pow(eps0,2))*Math.sqrt(kT*(1 + 1/Z)));

        // energy loss through cusps

        // cusp half-width
        var a = Math.sqrt(2.0*kT*mi)/(Z*e*Bc);

        // effective cusp area
        var Ac = 2*Math.PI*a*(2*R + a);

        // plugging voltages

        var phi_e = -voltage;
        var phi_i = voltage;

        // plasma potential determined by ambipolar cusp losses
        var phi_p = Te*1000*Math.log(Math.pow(mi/me, (1/4)));

        var cuspIonFlux = 0.1*(ne/Z)*Math.sqrt(2*kT/(Math.PI*mi))*Math.exp(-Z*e*(phi_i - phi_p)/kT);

        // assume ion loss = electron loss through cusp due to ambipolar
        // condition.
        var cusp_loss = cuspIonFlux*Ac*(3*kT/2 - e*(phi_e - phi_p) + 3*kT/2 + Z*e*(phi_i - phi_p));

        //total energy losses and confinement time
        var total_loss = cusp_loss + brem_loss + diff_loss;

        //T_tot = (3/2)*(ne + ni)*kT/P_loss;
        var Q = total_power/total_loss;

        // assuming heating efficiency is around 50%, and some of the lost power will be recovered as thermal energy
        var net_power = (total_power + total_loss)*thermal_eff - 2.0*total_loss;

        return {
            'B': Bc,
            'total_loss': total_loss,
            'total_power': total_power,
            'net_power': net_power,
            'Q': Q,
            'neutron_frac': neutron_power/total_power,
            'charged_frac': charged_power/total_power
        };
    };

    out.make = function (){

        var Te = 20;
        var ne = 8E20;
        var T_frac = 1.0;
        var He3_frac = 0;
        var radius = 2.0;
        var voltage = 300000;

        var recalc = function(){
            var reactor_out = out.calcFusion(Te, ne, T_frac, He3_frac, radius, voltage, 0.4);

            document.getElementById("total_power").innerText = reactor_out.total_power.toExponential(2) + " W";
            document.getElementById("total_loss").innerText = reactor_out.total_loss.toExponential(2) + " W";
            document.getElementById("Q").innerText = reactor_out.Q.toFixed(1);
            document.getElementById("net_power").innerText = reactor_out.net_power.toExponential(2) + " W";
            document.getElementById("neutron_frac").innerText = reactor_out.neutron_frac.toFixed(2);
            document.getElementById("charged_frac").innerText = reactor_out.charged_frac.toFixed(2);
            document.getElementById("B").innerText = reactor_out.B.toFixed(1) + " T";
        };


        $("#Te").slider({
            value: Te,
            orientation: "horizontal",
            range: "max",
            min: 1,
            max: 100,
            animate: true,
            change: function(event, ui) {
                Te = ui.value;
                recalc();
            }
        });

        $("#ne").slider({
            value: 8,
            orientation: "horizontal",
            range: "max",
            min: 1,
            max: 10,
            animate: true,
            change: function(event, ui) {
                ne = ui.value*1E20;
                recalc();
            }
        });

        $("#T_frac").slider({
            value: T_frac,
            orientation: "horizontal",
            range: "max",
            min: 0,
            max: 2,
            animate: true,
            change: function(event, ui) {
                T_frac = ui.value;
                recalc();
            }
        });


        $("#He3_frac").slider({
            value: He3_frac,
            orientation: "horizontal",
            range: "max",
            min: 0,
            max: 2,
            animate: true,
            change: function(event, ui) {
                He3_frac = ui.value;
                recalc();
            }
        });

        $("#radius").slider({
            value: radius,
            orientation: "horizontal",
            range: "max",
            min: 1,
            max: 10,
            animate: true,
            change: function(event, ui) {
                radius = ui.value;
                recalc();
            }
        });

        $("#voltage").slider({
            value: voltage/1000,
            orientation: "horizontal",
            range: "max",
            min: 0,
            max: 1000,
            animate: true,
            change: function(event, ui) {
                voltage = ui.value*1000;
                recalc();
            }
        });
    };

    out.make();

    return out;

};
