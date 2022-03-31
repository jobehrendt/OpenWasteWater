within OpenWasteWater;

package ADM_P

  record Solubles_ADM
    Real su "gCOD/m3 monosaccharides/sugar";
    Real aa "gCOD/m3 amino acids";
    Real fa "gCOD/m3 long chain fatty acids";
    Real va "gCOD/m3 valerate";
    Real bu "gCOD/m3 butyrate";
    Real pro "gCOD/m3 propionate";
    Real ac "gCOD/m3 acetate";
    Real H2 "gCOD/m3 hydrogen";
    Real ch4 "gCOD/m3 methane";
    Real IC "mole C/m3 inorganic carbon";
    Real IN "gN/m3 inorganic nitrogen";
    Real I "gCOD/m3 soluble inerts";
    Real P "gP/m3";
    Real Fe2 "gFe2/m3";
    Real Fe3 "gFe3/m3";
    Real Al3 "gAl3/m3";
  end Solubles_ADM;

  record Particulates_ADM
    Real c "gCOD/m3 composites";
    Real ch "gCOD/m3 carbohydrates";
    Real pr "gCOD/m3 proteins";
    Real li "gCOD/m3 lipids";
    Real su "gCOD/m3 sugar degraders";
    Real aa "gCOD/m3 amino acids degraders";
    Real fa "gCOD/m3 fatty acids degraders";
    Real c4 "gCOD/m3 butyrate & valerate degraders";
    Real pro "gCOD/m3 propionate degraders";
    Real ac "gCOD/m3 acetate degraders";
    Real H2 "gCOD/m3 hydrogen degraders";
    Real I "gCOD/m3 inert";
    Real FeP "gFePO4/m3";
    Real FeOH "gFe(OH)3/m3";
    Real AlP "gAlPO4/m3";
    Real AlOH "gAl(OH)3/m3";
  end Particulates_ADM;

  record Gaseous_ADM
    Real CO2 "mole/m3 ";
    Real H2 "gCOD/m3 ";
    Real NH3 "gN/m3 ";
    Real CH4 "gCOD/m3";
  end Gaseous_ADM;

  record BiogasValues
    Real Qnorm "m3/d";
    Real xCO2 "% ";
    Real xH2 "% ";
    Real xNH3 "% ";
    Real xCH4 "%";
  end BiogasValues;


  connector InFlow
    Real T "C";
    flow Real Q "m3/d";
    input Solubles_ADM S;
    input Particulates_ADM X;
   annotation(
      defaultComponentName = "InFlow",
       Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {143,89,2}, fillColor = {143,89,2}, fillPattern = FillPattern.Solid), Text(extent = {{-88, 92}, {88, -94}}, lineColor = {255, 255, 255}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "%name")}),
      Documentation(info = "<html>
    <H1 style=\"font-size:20px\">Connector (for ADM1 componets and Flow. Potential variable is the Temerature </H1>
    <p style=\"font-size:20px\">
    This connector defines a connector for the IWA Anaerobic Digestion Model No. 1
    for all 24 compounds as particulete and dissolved variables and the flowrate. To fullfill the
    modelica requirement the temperature is added as potential variable.
    </p>
  </html>"));
  end InFlow;




  connector OutFlow
    Real T "C";
    flow Real Q "m3/d";
    output Solubles_ADM S;
    output Particulates_ADM X;
  annotation(
      defaultComponentName = "OutFlow",
       Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {233,185,110}, fillColor = {233,185,110}, fillPattern = FillPattern.Solid), Text(extent = {{-88, 92}, {88, -94}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "%name")}),
      Documentation(info = "<html>
    <H1 style=\"font-size:20px\">Connector (for ADM1 componets and Flow. Potential variable is the Temerature </H1>
    <p style=\"font-size:20px\">
    This connector defines a connector for the IWA Anaerobic Digestion Model No. 1
    for all 24 compounds as particulete and dissolved variables and the flowrate. To fullfill the
    modelica requirement the temperature is added as potential variable.
    </p>
  </html>"));
  end OutFlow;


  class Csum
  protected       //hides the values from plot
    extends biochemical;
    //Carbon Content
    parameter Real C_su = 6 / 192;
    parameter Real C_aa = 0.03;
    parameter Real C_fa = 0.0217;
    parameter Real C_va = 5 / 208;
    parameter Real C_bu = 4 / 160;
    parameter Real C_pro = 3 / 112;
    parameter Real C_ac = 2 / 64;
    parameter Real C_ch4 = 1 / 64;
    parameter Real C_SI = 0.03;
    parameter Real C_xc = 0.0279;
    parameter Real C_ch = 0.0313;
    parameter Real C_pr = 0.03;
    parameter Real C_li = 0.0217;
    parameter Real C_biom = 5 / 160;
    parameter Real C_XI = 0.03;

//Csum

    Real Cr1 = f_SI_xc * C_SI - C_xc + f_ch_xc * C_ch + f_pr_xc * C_pr + f_li_xc * C_li + f_XI_xc * C_XI;
    Real Cr2 = C_su - C_ch;
    Real Cr3 = C_aa - C_pr;
    Real Cr4 = (1 - f_fa_li) * C_su + f_fa_li * C_fa - C_li;
    Real Cr5 = (-C_su) + (1 - Y_su) * (f_bu_su * C_bu + f_pro_su * C_pro + f_ac_su * C_ac) + Y_su * C_biom;
    Real Cr6 = (-C_aa) + (1 - Y_aa) * (f_va_aa * C_va + f_bu_aa * C_bu + f_pro_aa * C_pro + f_ac_aa * C_ac) + Y_aa * C_biom;
    Real Cr7 = (-C_fa) + (1 - Y_fa) * 0.7 * C_ac + Y_fa * C_biom;
    Real Cr8 = (-C_va) + (1 - Y_c4) * (0.54 * C_pro + 0.31 * C_ac) + Y_c4 * C_biom;
    Real Cr9 = (-C_bu) + (1 - Y_c4) * 0.8 * C_ac + Y_c4 * C_biom;
    Real Cr10 = (-C_pro) + (1 - Y_pro) * 0.57 * C_ac + Y_pro * C_biom;
    Real Cr11 = (-C_ac) + (1 - Y_ac) * C_ch4 + Y_ac * C_biom;
    Real Cr12 = (1 - Y_H2) * C_ch4 + Y_H2 * C_biom;
    annotation(
      defaultComponentName = "Csum",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Carbon Content factors</H1>
      <p style=\"font-size:20px\">
      Calculation factors v_i,j of the Peterson Matrix for ADM1, for the calculation of S.IC.
      Numeration according to the corresponding process.
      Values for carbon content are taken from
        'Y. Feng. Calibration and verication of a mathematical model for the simulation of
        blackwater/biowaste digestion. Master's thesis, Hamburg University of Technology, 2004.'
      </p>
    </html>"));
  end Csum;


  class biochemical
  protected       //hides the values from plot
    //values from Masterthesis Appendix E ADM1 and ADM1
    //Nitrogen Content
    parameter Real N_aa = 0.098 "g N/gCOD";
    parameter Real N_SI = 0.028 "g N/gCOD";
    parameter Real N_xc = 0.028 "g N/gCOD";
    parameter Real N_pr = 0.098 "g N/gCOD";
    parameter Real N_biom = 0.0875 "g N/gCOD";
    parameter Real N_XI = 0.028 "g N/gCOD";
    //stoichiometric parameters for mass flux
    parameter Real f_SI_xc = 0.1 "soluble inert from composites";
    parameter Real f_XI_xc = 0.25 "particulate inert from composites";
    parameter Real f_ch_xc = 0.2 "carbohydrates from composites";
    parameter Real f_pr_xc = 0.2 "proteins from composites";
    parameter Real f_li_xc = 0.25 "lipids from composites";
    parameter Real f_fa_li = 0.95 "LCFA from lipids";
    parameter Real f_H2_su = 0.1906 "hydrogen from sugars";
    parameter Real f_bu_su = 0.1328 "butyrate from sugars";
    parameter Real f_pro_su = 0.2690 "propionate from sugars";
    parameter Real f_ac_su = 0.4076 "acetate from sugars";
    parameter Real f_H2_aa = 0.06 "hydrogen from amino acids";
    parameter Real f_va_aa = 0.23 "valerate from amino acids";
    parameter Real f_bu_aa = 0.26 "butyrate from amino acids";
    parameter Real f_pro_aa = 0.05 "propionate from amino acids";
    parameter Real f_ac_aa = 0.4 "acetate from amino acids";
    //disintegration & hydrolysis
    parameter Real k_dis = 0.5 "1/d disintegration rate of composites";
    parameter Real k_hyd_ch = 10 "1/d hydrolysis rate of carbohydrates";
    parameter Real k_hyd_pr = 10 "1/d hydrolysis rate of proteins";
    parameter Real k_hyd_li = 10 "1/d hydrolysis rate of lipids";
    //biomass decay
    parameter Real k_dec = 0.02 "1/d";
    //max. uptake rate
    parameter Real k_m_su = 30 "1/d, sugar degraders";
    parameter Real k_m_aa = 50 "1/d, amino acid degraders";
    parameter Real k_m_fa = 6 "1/d, LCFA degraders";
    parameter Real k_m_c4 = 20 "1/d, valerate and butyrate degraders";
    parameter Real k_m_pro = 13 "1/d, proprionate degraders";
    parameter Real k_m_ac = 8 "1/d, acetate degraders";
    parameter Real k_m_H2 = 35 "1/d, hydrogen degraders";
    //half-saturation concentration
    parameter Real K_S_su = 500 "g COD/m³";
    parameter Real K_S_aa = 300 "g COD/m³";
    parameter Real K_S_fa = 400 "g COD/m³";
    parameter Real K_S_c4 = 200 "g COD/m³";
    parameter Real K_S_pro = 100 "g COD/m³";
    parameter Real K_S_ac = 150 "g COD/m³";
    parameter Real K_S_H2 = 0.007 "g COD/m³";
    //yield from substrates to degraders
    parameter Real Y_su = 0.1 "gCOD/gCOD, sugar degraders";
    parameter Real Y_aa = 0.08 "gCOD/gCOD, amino acid degraders";
    parameter Real Y_fa = 0.06 "gCOD/gCOD, LCFA degraders";
    parameter Real Y_c4 = 0.06 "gCOD/gCOD, valerate and butyrate degraders";
    parameter Real Y_pro = 0.04 "gCOD/gCOD, proprionate degraders";
    parameter Real Y_ac = 0.05 "gCOD/gCOD, acetate degraders";
    parameter Real Y_H2 = 0.06 "gCOD/gCOD, hydrogen degraders";
    //Inhibition and limitation
    parameter Real K_I_H2_fa = 0.005 "g COD/m³, hydrogen to LCFA uptake";
    parameter Real K_I_H2_c4 = 0.01 "g COD/m³, hydrogen to bu and va uptake";
    parameter Real K_I_H2_pro = 0.0035 "g COD/m³, hydrogen to pro uptake";
    parameter Real K_I_NH3_ac = 25.2 "g NH3-N/m³, free NH3 to acetate uptake";
    parameter Real K_S_NH3_all = 1.4 "g N/m³ Inorganic nitrogen limitation";
    parameter Real pH_ul_su_pro = 5.5 "pH upper limit, su-pro degraders";
    parameter Real pH_ll_su_pro = 4.0 "pH lower limit, su-pro degraders";
    parameter Real pH_ul_ac = 8.0 "pH upper limit, ac degraders";
    parameter Real pH_ll_ac = 6.0 "pH lower limit, ac degraders";
    parameter Real pH_ul_H2 = 6.0 "pH upper limit, H2 degraders";
    parameter Real pH_ll_H2 = 5.0 "pH lower limit, H2 degraders";
   annotation(
      defaultComponentName = "biochemical",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Biochemical Parameters for the ADM</H1>
      <p style=\"font-size:20px\">
      This class 'biochemical' contains all biochemical coefficients regarding the IWA
    Anaerobic Digestion Model No.1 with the suggested default values.

     Values for nitrogen content are taken from
        'Y. Feng. Calibration and verication of a mathematical model for the simulation of
        blackwater/biowaste digestion. Master's thesis, Hamburg University of Technology, 2004.'
      </p>
    </html>"));

  end biochemical;




  class physicochemical
  protected       //hides the values from plot
    //acid-base equilibrium coefficients
    parameter Real K_W = 10 ^ ((-14) + 6) "mM²";
    parameter Real K_a_va = 10 ^ ((-4.80) + 3) "mole/m³, mM";
    parameter Real K_a_bu = 10 ^ ((-4.84) + 3) "mole/m³, mM";
    parameter Real K_a_pro = 10 ^ ((-4.88) + 3) "mole/m³, mM";
    parameter Real K_a_ac = 10 ^ ((-4.76) + 3) "mole/m³, mM";
    parameter Real K_a_CO2 = 10 ^ ((-6.35) + 3) "mole/m³, mM";
    parameter Real K_a_HCO3_neg = 10 ^ ((-10.33) + 3) "mole/m³, mM";
    parameter Real K_a_NH4_pos = 10 ^ ((-9.25) + 3) "mole/m³, mM";
      //free energy
    parameter Real H0_a_CO2 = 7646 "J/mole";
    parameter Real H0_a_HCO3_neg = 14850 "J/mole";
    parameter Real H0_a_W = 55900 "J/mole";
    parameter Real H0_a_NH4_pos = 51965 "J/mole";
    //Liquid-gas transfer parameter values
    //Henry's law constant
    parameter Real K_H_H2 = 0.78 "mMliq/bargas at 298K";
    parameter Real K_H_ch4 = 1.4 "mMliq/bargas at 298K";
    parameter Real K_H_CO2 = 35 "mMliq/bargas at 298K";
    parameter Real K_H_NH3 = 57540 "mMliq/bargas at 298K";
      //free energy
    parameter Real H0_H_H2 = -4180 "J/mole";
    parameter Real H0_H_ch4 = -14240 "J/mole";
    parameter Real H0_H_CO2 = -19410 "J/mole";
      //
    parameter Real R = 8.314472 "universal gas law constant, J/mole K";
    parameter Real P_atm = 1 "atmospheric pressure, bar";
    parameter Real k_p = 10000 "pipe resistance coefficent, m³/d bar";
  annotation(
    defaultComponentName = "physicochemical",
    Documentation(info = "<html>
    <H1 style=\"font-size:20px\">Physicochemical Parameters for the ADM</H1>
    <p style=\"font-size:20px\">
    This class 'physicochemical' contains all physicochemical coefficients regarding the IWA
  Anaerobic Digestion Model No.1 with the suggested default values.

   Additional values were are taken from
      'Y. Feng. Calibration and verication of a mathematical model for the simulation of
      blackwater/biowaste digestion. Master's thesis, Hamburg University of Technology, 2004.'
    </p>
  </html>"));
  end physicochemical;


  model inflow_ASM
    parameter Real T = 38 "operational temperature, C";
    Real TSS_In "gTSS/m3";
    OpenWasteWater.ADM_P.OutFlow Out1 annotation(
      Placement(visible = true, transformation(origin = {61, -1}, extent = {{-25, -25}, {25, 25}}, rotation = 0), iconTransformation(origin = {61, -1}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
    OpenWasteWater.ASM1P.InPipe In1 annotation(
      Placement(visible = true, transformation(origin = {-56, -2}, extent = {{-26, -26}, {26, 26}}, rotation = 0), iconTransformation(origin = {-68, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    Out1.Q + In1.Q = 0.0;
    Out1.T = T;
//SOLUBLES
    Out1.S.su = 0;
    Out1.S.aa = 0;
    Out1.S.fa = 0;
    Out1.S.va = 0;
    Out1.S.bu = 0;
    Out1.S.pro = 0;
    Out1.S.ac = 0;
    Out1.S.H2 = 0;
    Out1.S.ch4 = 0;
    Out1.S.IC = In1.S.ALK "mole/m³, S_ALK";
    Out1.S.IN = In1.S.NH "gN/m³, S_NH";
    Out1.S.I = In1.S.I "gCOD/m³, S_I";
    Out1.S.P = In1.S.P;
    Out1.S.Fe2 = In1.S.Fe2;
    Out1.S.Fe3 = In1.S.Fe3;
    Out1.S.Al3 = In1.S.Al3;
  
//PARTICULATES
    Out1.X.c = In1.X.S + In1.X.H + In1.X.A + In1.X.P "gCOD/m³, All X, exluding X_I";
    Out1.X.ch = 0;
    Out1.X.pr = 0;
    Out1.X.li = 0;
    Out1.X.su = 0;
    Out1.X.aa = 0;
    Out1.X.fa = 0;
    Out1.X.c4 = 0;
    Out1.X.pro = 0;
    Out1.X.ac = 0;
    Out1.X.H2 = 0;
    Out1.X.I = In1.X.I "gCOD/m³, X_I";
    Out1.X.FeP = In1.X.FeP;
    Out1.X.FeOH = In1.X.FeOH;
    Out1.X.AlP = In1.X.AlP;
    Out1.X.AlOH = In1.X.AlOH;
  
    TSS_In = 0.75 * (In1.X.I + In1.X.S + In1.X.H + In1.X.A + In1.X.P) + In1.X.FeP + In1.X.FeOH + In1.X.AlP + In1.X.AlOH;
    
    annotation(
      defaultComponentName = "inflow_ADM1",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Inteface ASM1 to ADM1 </H1>
      <p style=\"font-size:20px\">
      Inflow for combination with the ASM1.
        Conversion according to the documentation of the ADM1 from the IWA, Appendix C.3, Table C.2

      Parameter T is important as it defines the operational temperature of the digester.
      </p>
    </html>"),
  Icon(graphics = {Rectangle(origin = {2, 0}, lineThickness = 0.75, extent = {{-22, 62}, {22, -62}})}));
  end inflow_ASM;



  model inflow "without WWTP"
    parameter Real Q = 0.4 "m³/d";
    parameter Real T = 38 "operational temperature, C";
    OutFlow Out1;
  equation
    Out1.Q = -Q;
    Out1.T = T;
//SOLUBLES
    Out1.S.su = 0;
    Out1.S.aa = 0;
    Out1.S.fa = 0;
    Out1.S.va = 0;
    Out1.S.bu = 0;
    Out1.S.pro = 0;
    Out1.S.ac = 0;
    Out1.S.H2 = 0;
    Out1.S.ch4 = 0;
    Out1.S.IC = 4 "mole/m³, S_ALK";
    Out1.S.IN = 0.4 "gN/m³, S_NH";
    Out1.S.I = 37 "gCOD/m³, S_I";
    Out1.S.P = 1.0 "gP/m3";
    Out1.S.Fe2 = 1.0 "gFe2/m3";
    Out1.S.Fe3 = 1.0 "gFe3/m3";
    Out1.S.Al3 = 1.0 "gAl3/m3";
  
  //PARTICULATES
    Out1.X.c = 107000 "gCOD/m³, All X, exluding X_I";
    Out1.X.ch = 0;
    Out1.X.pr = 0;
    Out1.X.li = 0;
    Out1.X.su = 0;
    Out1.X.aa = 0;
    Out1.X.fa = 0;
    Out1.X.c4 = 0;
    Out1.X.pro = 0;
    Out1.X.ac = 0;
    Out1.X.H2 = 0;
    Out1.X.I = 2800 "gCOD/m³, X_I";
    Out1.FeP = 1.0 "gFePO4/m3";
    Out1.FeOH = 1.0 "gFe(OH)3/m3";
    Out1.AlP = 1.0 "gAlPO4/m3";
    Out1.AlOH = 1.0 "gAl(OH)3/m3";
  
    annotation(
      defaultComponentName = "inflow",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">single inflow for ADM1 </H1>
      <p style=\"font-size:20px\">
  Inflow for the ADM1 as single model without combination with the ASM

      Parameter T is important as it defines the operational temperature of the digester.

  The chosen values come from the sludge flow of the simple WWTP model
      </p>
    </html>"));
  end inflow;




  function pH_lower_Inhibition
    input Real pH_ll, pH_ul, pH;
    output Real I_pH;
  algorithm
    if pH < pH_ul then
      I_pH := exp(-3 * ((pH - pH_ul) / (pH_ul - pH_ll)) ^ 2);
    else
      I_pH := 1;
    end if;
    annotation(
      defaultComponentName = "pH_Inhibition",
      Documentation(info = "<html>
    <H1 style=\"font-size:20px\">pH Inhibition Factors</H1>
    <p style=\"font-size:20px\">
    Inhibtion of the substrate uptake by low pH.
    Inhibition equation according to the IWA Anaerobic Digestion Model No.1
    </p>
  </html>"));
  end pH_lower_Inhibition;



  function pH_Inhibition
    input Real pH_ll, pH_ul, pH;
    output Real I_pH;
  algorithm
    I_pH := (1 + 2 * 10 ^ (0.5 * (pH_ll - pH_ul))) / (1 + 10 ^ (pH - pH_ul) + 10 ^ (pH_ll - pH));
    annotation(
      defaultComponentName = "pH_Inhibition",
      Documentation(info = "<html>
    <H1 style=\"font-size:20px\">pH Inhibition Factors</H1>
    <p style=\"font-size:20px\">
    Inhibtion of the substrate uptake by lower and upper pH limits.
    Inhibition equation according to the IWA Anaerobic Digestion Model No.1.
    </p>
  </html>"));
  end pH_Inhibition;


  function K
    input Real K1, H, T2;
    output Real K2;
  protected
    Real T1 = 298 "K";
    extends physicochemical;
  algorithm
    K2 := K1 * exp(H / R * (1 / T1 - 1 / (T2 + 273.15)));
  annotation(
      defaultComponentName = "K",
      Documentation(info = "<html>
    <H1 style=\"font-size:20px\">Temperature adaption for the K-values</H1>
    <p style=\"font-size:20px\">
    Temperature adaption for K_H (Henry's law constants) and K_a (acid dissociation coefficients)
    Conversion equation according to the IWA Anaerobic Digestion Model No.1.
    </p>
  </html>"));
  end K;



  model CSTR_ADM
    //General
    parameter Real V_R = 1200 "reactor volume, m³";
    parameter Real kLa = 20 "1/d";
    parameter OpenWasteWater.ADM_P.Solubles_ADM Sini(
      su = 12,
      aa = 5.36,
      fa = 102.7,
      va = 10.75, 
      bu = 14.3,
      pro = 16.85,
      ac = 39,
      H2 = 0.00024,
      ch4 = 152,
      IC = 52.4,
      IN = 527,
      I = 6400,
      P = 1.0,
      Fe2 = 0,
      Fe3 = 0,
      Al3 = 0);
    parameter OpenWasteWater.ADM_P.Particulates_ADM Xini(
      c = 5500,
      ch = 53,
      pr = 54,
      li = 67,
      su = 835,
      aa = 630,
      fa = 560,
      c4 = 273,
      pro = 129,
      ac = 820,
      H2 = 390,
      I = 34000,
      FeP = 0,
      FeOH = 0,
      AlP = 0,
      AlOH = 0);
  
  
    InFlow In1;
    OutFlow Out1;
    InPgas InP;
    GasFlowOut OutG;
    Solubles_ADM S;
    Particulates_ADM X;
    //Real SRT "d, optional";
    //parameters
    extends biochemical;
    extends physicochemical;
    extends Csum;
    Real COD "gCOD/m³";
    //Inhibitors
    Real I_pH, I_IN, I_pH_ac, I_pH_H2, I_H2_fa, I_H2_c4, I_H2_pro, I_NH3_ac;
    //equilibrium
    Real pH "= 7.0";
    //Real pH1 = 7.0;
    Real S_NH3, S_NH4_pos "g N/m³";
    Real S_CO2, S_CO3_2neg, S_HCO3_neg, S_OH_neg "mole/m³";
    Real S_ac_neg, S_pro_neg, S_bu_neg, S_va_neg ", S_H_pos g COD/m³";
    output Real S_H_pos(start = 0.0001) "mole/m³";
    Real S_cat = 10.0 "mole/m³";
    //For finding of cat/an: pH and pH1 = desired pH change equation with log(S_H_pos) to pH1
    Real S_an = 0.0 "mole/m³";
    //Gas transfer rates
    Real r_lg_H2 "g COD/m³*d";
    Real r_lg_ch4 "g COD/m³*d";
    Real r_lg_CO2 "mole/m³*d";
    Real r_lg_NH3 "g N/m³*d";
  protected
    //biochemical rates "g/COD m³"
    Real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19;
  initial equation
//Solubles
  S = Sini;
  
  //PARTICULATES
    X = Xini;
  
  equation
//Derivatives
////Solubles
    der(S.su) = In1.Q / V_R * (In1.S.su - S.su) + r2 + (1 - f_fa_li) * r4 - r5;
    der(S.aa) = In1.Q / V_R * (In1.S.aa - S.aa) + r3 - r6;
    der(S.fa) = In1.Q / V_R * (In1.S.fa - S.fa) + f_fa_li * r4 - r7;
    der(S.va) = In1.Q / V_R * (In1.S.va - S.va) + (1 - Y_aa) * f_va_aa * r6 - r8;
    der(S.bu) = In1.Q / V_R * (In1.S.bu - S.bu) + (1 - Y_su) * f_bu_su * r5 + (1 - Y_aa) * f_bu_aa * r6 - r9;
    der(S.pro) = In1.Q / V_R * (In1.S.pro - S.pro) + (1 - Y_su) * f_pro_su * r5 + (1 - Y_aa) * f_pro_aa * r6 + (1 - Y_c4) * 0.54 * r8 - r10;
    der(S.ac) = In1.Q / V_R * (In1.S.ac - S.ac) + (1 - Y_su) * f_ac_su * r5 + (1 - Y_aa) * f_ac_aa * r6 + (1 - Y_fa) * 0.7 * r7 + (1 - Y_c4) * 0.31 * r8 + (1 - Y_c4) * 0.8 * r9 + (1 - Y_pro) * 0.57 * r10 - r11;
    der(S.H2) = In1.Q / V_R * (In1.S.H2 - S.H2) + (1 - Y_su) * f_H2_su * r5 + (1 - Y_aa) * f_H2_aa * r6 + (1 - Y_fa) * 0.3 * r7 + (1 - Y_c4) * 0.15 * r8 + (1 - Y_c4) * 0.2 * r9 + (1 - Y_pro) * 0.43 * r10 - r12 - r_lg_H2;
    der(S.ch4) = In1.Q / V_R * (In1.S.ch4 - S.ch4) + (1 - Y_ac) * r11 + (1 - Y_H2) * r12 - r_lg_ch4;
    der(S.IC) = In1.Q / V_R * (In1.S.IC - S.IC) - (Cr1 * r1 + Cr2 * r2 + Cr3 * r3 + Cr4 * r4 + Cr5 * r5 + Cr6 * r6 + Cr7 * r7 + Cr8 * r8 + Cr9 * r9 + Cr10 * r10 + Cr11 * r11 + Cr12 * r12) - r_lg_CO2;
    der(S.IN) = In1.Q / V_R * (In1.S.IN - S.IN) - Y_su * N_biom * r5 + (N_aa - Y_aa * N_biom) * r6 - Y_fa * N_biom * r7 - Y_c4 * N_biom * r8 - Y_c4 * N_biom * r9 - Y_pro * N_biom * r10 - Y_ac * N_biom * r11 - Y_H2 * N_biom * r12 - r_lg_NH3;
    der(S.I) = In1.Q / V_R * (In1.S.I - S.I) + f_SI_xc * r1;

    der(S.P) = In1.Q / V_R * (In1.S.P - S.P) + 0.015 * r1 - 0.01 *(r5+r6+r7+r8+r9+r10+r11+r12-r13-r14-r15-r16-r18-r19);
    der(S.Fe2) = In1.Q / V_R * (In1.S.Fe2 - S.Fe2);
    der(S.Fe3) = In1.Q / V_R * (In1.S.Fe3 - S.Fe3);
    der(S.Al3) = In1.Q / V_R * (In1.S.Al3 - S.Al3);
  
  
  
  
  ////Particulates
    der(X.c) = In1.Q / V_R * (In1.X.c - X.c) - r1 + r13 + r14 + r15 + r16 + r17 + r18 + r19;
    der(X.ch) = In1.Q / V_R * (In1.X.ch - X.ch) + f_ch_xc * r1 - r2;
    der(X.pr) = In1.Q / V_R * (In1.X.pr - X.pr) + f_pr_xc * r1 - r3;
    der(X.li) = In1.Q / V_R * (In1.X.li - X.li) + f_li_xc * r1 - r4;
    der(X.su) = In1.Q / V_R * (In1.X.su - X.su) + Y_su * r5 - r13;
    der(X.aa) = In1.Q / V_R * (In1.X.aa - X.aa) + Y_aa * r6 - r14;
    der(X.fa) = In1.Q / V_R * (In1.X.fa - X.fa) + Y_fa * r7 - r15;
    der(X.c4) = In1.Q / V_R * (In1.X.c4 - X.c4) + Y_c4 * r8 + Y_c4 * r9 - r16;
    der(X.pro) = In1.Q / V_R * (In1.X.pro - X.pro) + Y_pro * r10 - r17;
    der(X.ac) = In1.Q / V_R * (In1.X.ac - X.ac) + Y_ac * r11 - r18;
    der(X.H2) = In1.Q / V_R * (In1.X.H2 - X.H2) + Y_H2 * r12 - r19;
    der(X.I) = In1.Q / V_R * (In1.X.I - X.I) + f_XI_xc * r1;

    der(X.FeP) = In1.Q / V_R * (In1.X.FeP - X.FeP);
    der(X.FeOH) = In1.Q / V_R * (In1.X.FeOH - X.FeOH);
    der(X.AlP) = In1.Q / V_R * (In1.X.AlP - X.AlP);
    der(X.AlOH) = In1.Q / V_R * (In1.X.AlOH - X.AlOH);
  
  
  //Inhibition and Limitation
////pH-Inhibition
    I_pH = pH_lower_Inhibition(pH_ll_su_pro, pH_ul_su_pro, pH);
    I_pH_ac = pH_Inhibition(pH_ll_ac, pH_ul_ac, pH);
    I_pH_H2 = pH_lower_Inhibition(pH_ll_H2, pH_ul_H2, pH);
////non-competetive Inhibition
    I_H2_fa = K_I_H2_fa / (K_I_H2_fa + S.H2);
    I_H2_c4 = K_I_H2_c4 / (K_I_H2_c4 + S.H2);
    I_H2_pro = K_I_H2_pro / (K_I_H2_pro + S.H2);
    I_NH3_ac = K_I_NH3_ac / (K_I_NH3_ac + S_NH3);
////IN Limitation
    I_IN = S.IN / (K_S_NH3_all + S.IN);
//biochemical process rates
////disintegration and hydrolysis
    r1 = k_dis * X.c;
    r2 = k_hyd_ch * X.ch;
    r3 = k_hyd_pr * X.pr;
    r4 = k_hyd_li * X.li;
////uptake
    r5 = k_m_su * (S.su / (K_S_su + S.su)) * X.su * I_pH * I_IN;
    r6 = k_m_aa * (S.aa / (K_S_aa + S.aa)) * X.aa * I_pH * I_IN;
    r7 = k_m_fa * (S.fa / (K_S_fa + S.fa)) * X.fa * I_pH * I_IN * I_H2_fa;
    r8 = k_m_c4 * (S.va / (K_S_c4 + S.va)) * X.c4 * (S.va / (S.va + S.bu)) * I_pH * I_IN * I_H2_c4;
    r9 = k_m_c4 * (S.bu / (K_S_c4 + S.bu)) * X.c4 * (S.bu / (S.bu + S.va)) * I_pH * I_IN * I_H2_c4;
    r10 = k_m_pro * (S.pro / (K_S_pro + S.pro)) * X.pro * I_pH * I_IN * I_H2_pro;
    r11 = k_m_ac * (S.ac / (K_S_ac + S.ac)) * X.ac * I_pH_ac * I_IN * I_NH3_ac;
    r12 = k_m_H2 * (S.H2 / (K_S_H2 + S.H2)) * X.H2 * I_pH_H2 * I_IN;
////decay
    r13 = k_dec * X.su;
    r14 = k_dec * X.aa;
    r15 = k_dec * X.fa;
    r16 = k_dec * X.c4;
    r17 = k_dec * X.pro;
    r18 = k_dec * X.ac;
    r19 = k_dec * X.H2;
//physico-chemical processes
////liquid-liquid equilibrium processes
    S_cat + S_NH4_pos / 14 + S_H_pos - S_HCO3_neg - 2 * S_CO3_2neg - S_ac_neg / 64 - S_pro_neg / 112 - S_bu_neg / 160 - S_va_neg / 208 - S_OH_neg - S_an = 0;
    S_OH_neg - K(K_W, H0_a_W, In1.T) / S_H_pos = 0;
    S_va_neg - K_a_va * S.va / (K_a_va + S_H_pos) = 0;
    S_bu_neg - K_a_bu * S.bu / (K_a_bu + S_H_pos) = 0;
    S_pro_neg - K_a_pro * S.pro / (K_a_pro + S_H_pos) = 0;
    S_ac_neg - K_a_ac * S.ac / (K_a_ac + S_H_pos) = 0;
    S_HCO3_neg - K(K_a_CO2, H0_a_CO2, In1.T) * S.IC / (K(K_a_CO2, H0_a_CO2, In1.T) + S_H_pos) = 0;
    S_CO3_2neg - K(K_a_HCO3_neg, H0_a_HCO3_neg, In1.T) * S_HCO3_neg / (K(K_a_HCO3_neg, H0_a_HCO3_neg, In1.T) + S_H_pos) = 0;
    S.IC - S_CO2 - S_HCO3_neg - S_CO3_2neg = 0;
    S_NH4_pos - S_H_pos * S.IN / (K(K_a_NH4_pos, H0_a_NH4_pos, In1.T) + S_H_pos) = 0;
    S.IN - S_NH3 - S_NH4_pos = 0;
    pH = (-log10(S_H_pos)) + 3;
////liquid-gas
    r_lg_H2 = kLa * (S.H2 - K(K_H_H2, H0_H_H2, In1.T) * InP.P_H2 * 16);
    r_lg_ch4 = kLa * (S.ch4 - K(K_H_ch4, H0_H_H2, In1.T) * InP.P_ch4 * 64);
    r_lg_CO2 = kLa * (S_CO2 - K(K_H_CO2, H0_H_CO2, In1.T) * InP.P_CO2);
    r_lg_NH3 = kLa * (S_NH3 - K_H_NH3 * InP.P_NH3 * 14);
    COD = S.su + S.aa + S.fa + S.va + S.bu + S.ac + S.H2 + S.ch4 + S.I + X.c + X.ch + X.pr + X.li + X.su + X.aa + X.fa + X.c4 + X.pro + X.ac + X.H2 + X.I "gCOD/m³";
//OUT
    Out1.Q + In1.Q = 0;
    Out1.T = In1.T;
    Out1.S = S;
    Out1.X = X;
    OutG.r_lg_H2 = r_lg_H2;
    OutG.r_lg_ch4 = r_lg_ch4;
    OutG.r_lg_CO2 = r_lg_CO2;
    OutG.r_lg_NH3 = r_lg_NH3;
    OutG.T = In1.T;
    OutG.V_R = V_R;

  annotation(
      defaultComponentName = "ADM1_CSTR",
      Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Digestion Tank of the ADM1</H1>
  <p style=\"font-size:20px\">
  This Model is a model for an CSTR hat work like described for the ADM no. 1. The CSTR works in combination with the ASM1 if the 'inflow_ASM' is connected in between.

  Input parameters are the oxygen mass transfer coefficient kLa and the
  reactor volume V_R as well as the chosen scenario. The reactor volume is by default chosen that the SRT is 20 days. The SRT calculation is therefore commented out and can be activated optional.
  The operational temperature can be defined and adjusted by the parameter T in the inflow model.

  The CSTR does not include the headspace for the gas compounds. The headspace is modelled separately and connected via the gas connectors. The gas transfer rates are transmitted to the headspace in exchange with the present gas pressure in the headspace.

  The equilibrium equations determine the pH in the reactor. As the system is very sensitive errors can occur if the initial values in the reactor do not fit.

  The initial values are given in the class 'initial_values_ADM' and can be determined by the chosen scenario.

  </p>
  </html>"));
  end CSTR_ADM;







  model Headspace_ADM
    parameter Real V_G = 2 "Headspace volume, m3";
    parameter Gaseous_ADM Gini(CO2=14.6, H2=0.0037, NH3=0.005, CH4=1646.0) "for initial values";
    Real S_H2g "g COD/m³";
    Real S_ch4g "g COD/m³";
    Real S_CO2g "mole/m³";
    Real S_NH3g "g N/m³";
    Real P_H2, P_ch4, P_CO2, P_NH3 "bar";
    Real P_tot, P_H2O "bar";
    Real q_gas "m³/d";
    Real q_gas_norm "m³/d";
    extends physicochemical;
   // extends initial_values_ADM;
    GasFlowIn InG;
    OutPgas OutP;
  initial equation
    S_H2g = Gini.H2;
    S_ch4g = Gini.CH4;
    S_CO2g = Gini.CO2;
    S_NH3g = Gini.NH3;
  equation
//pressure
    P_H2 = S_H2g / 16 * R * 10 ^ (-5) * (InG.T + 273.15);
    P_ch4 = S_ch4g / 64 * R * 10 ^ (-5) * (InG.T + 273.15);
    P_CO2 = S_CO2g * R * 10 ^ (-5) * (InG.T + 273.15);
    P_NH3 = S_NH3g / 14 * R * 10 ^ (-5) * (InG.T + 273.15);
    P_tot = P_H2 + P_ch4 + P_CO2 + P_H2O + P_NH3;
    P_H2O = 0.0313 * exp(5290 * (1 / 298 - 1 / (InG.T + 273.15)));
//gas flow
    q_gas = k_p * (P_tot - P_atm);
    q_gas_norm = q_gas * (P_tot - P_H2O) * (273.15 / (InG.T + 273.15));
//derivatives
    der(S_H2g) =  (-q_gas * S_H2g / V_G)  + InG.r_lg_H2  * (InG.V_R / V_G);
    der(S_ch4g) = (-q_gas * S_ch4g / V_G) + InG.r_lg_ch4 * (InG.V_R / V_G);
    der(S_CO2g) = (-q_gas * S_CO2g / V_G) + InG.r_lg_CO2 * (InG.V_R / V_G);
    der(S_NH3g) = (-q_gas * S_NH3g / V_G) + InG.r_lg_NH3 * (InG.V_R / V_G);
//Out
    OutP.P_H2 = P_H2;
    OutP.P_ch4 = P_ch4;
    OutP.P_CO2 = P_CO2;
    OutP.P_NH3 = P_NH3;
  annotation(
      defaultComponentName = "ADM1_Headspace",
      Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Headspace of the anaerobic digestion Tank of the ADM1</H1>
  <p style=\"font-size:20px\">
  This Model is a model for the headspace of the CSTR.

  Input parameters are the oxygen mass transfer coefficient kLa and the
  reactor volume V_R as well as the chosen scenario. The headspace volume is by default chosen to be 0.25 * V_Digester.

  The headspace is connected to the CSTR via the gas connectors. The gas transfer rates are transmitted from the CSTR to the headspace in exchange with the present gas pressure in the headspace.

  The initial values are given in the class 'initial_values_ADM' and can be determined by the chosen scenario.

  </p>
  </html>"));
  end Headspace_ADM;




  connector GasFlowIn
    input Real T "C";
    input Real r_lg_H2 "g COD/m³*d";
    input Real r_lg_ch4 "g COD/m³*d";
    input Real r_lg_CO2 "mole/m³*d";
    input Real r_lg_NH3 "g N/m³*d";
    input Real V_R "m3";
  annotation(
    defaultComponentName = "GasFlowIn",
    Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Gas flow from CSTR to Headspace</H1>
  <p style=\"font-size:20px\">
  This connector defines a connector between the anaerobic digestion tank of the ADM1 and the belonging headspace. The connector includes the gas transfer rates of the 4 gas compounds, the operational temperature of the CSTR and Volume of the digester as required for the derivatives of the gas compounds.
  </p>
  </html>"));
  end GasFlowIn;



  connector GasFlowOut
    output Real T;
    output Real r_lg_H2  "g COD/m³*d";
    output Real r_lg_ch4 "g COD/m³*d";
    output Real r_lg_CO2 "mole/m³*d";
    output Real r_lg_NH3 "g N/m³*d";
    output Real V_R;
  annotation(
  defaultComponentName = "GasFlowOut",
  Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Gas flow from CSTR to Headspace</H1>
  <p style=\"font-size:20px\">
  This connector defines a connector between the anaerobic digestion tank of the ADM1 and the belonging headspace. The connector includes the gas transfer rates of the 4 gas compounds, the operational temperature of the CSTR and Volume of the digester as required for the derivatives of the gas compounds.
  </p>
  </html>"));
  end GasFlowOut;


  connector InPgas
    input Real P_H2;
    input Real P_ch4;
    input Real P_CO2;
    input Real P_NH3;
  annotation(
  defaultComponentName = "InPgas",
  Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Present pressure in the Headspace</H1>
  <p style=\"font-size:20px\">
  This connector defines a connector between the anaerobic digestion tank of the ADM1 and the belonging headspace. The connector transmits the information about the present gas pressure in the headspace to the CSTR for the calculation of the gas transfer rates of the 4 gas compounds.
  </p>
  </html>"));
  end InPgas;



  connector OutPgas
    output Real P_H2;
    output Real P_ch4;
    output Real P_CO2;
    output Real P_NH3;
  annotation(
  defaultComponentName = "OutPgas",
  Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Present pressure in the Headspace</H1>
  <p style=\"font-size:20px\">
  This connector defines a connector between the anaerobic digestion tank of the ADM1 and the belonging headspace. The connector transmits the information about the present gas pressure in the headspace to the CSTR for the calculation of the gas transfer rates of the 4 gas compounds.
  </p>
  </html>"));
  end OutPgas;

    connector inBiogasSensor
      input BiogasValues Value;
      annotation(
        defaultComponentName = "InBiogasSensor",
        Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {148, 0, 204}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {0, 0, 15}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
        Documentation(info = "<html>
        <H1 style=\"font-size:20px\">Connector (Biogas parameter)  </H1>
        <p style=\"font-size:20px\">
        Used for control of flowrates to enhance the performance of the WWTP and/or operate economic.
        </p>
      </html>"));
    end inBiogasSensor;

    connector outBiogasSensor
      output BiogasValues Value;
      annotation(
        defaultComponentName = "OutBiogasSensor",
        Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {255, 182, 203}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {155, 155, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
        Documentation(info = "<html>
        <H1 style=\"font-size:20px\">Connector (Biogas parameter)  </H1>
        <p style=\"font-size:20px\">
        Used for control of flowrates to enhance the performance of the WWTP and/or operate economic.
        </p>
      </html>"));
    end outBiogasSensor;






  model outflow "without WWTP"
    Real Q, T;
    OpenWasteWater.ADM_P.Solubles_ADM S;
    OpenWasteWater.ADM_P.Particulates_ADM X;
    InFlow In1;
  equation
    Q = In1.Q;
    T = In1.T;
    S = In1.S;
    X = In1.X;
  annotation(
    defaultComponentName = "outflow",
    Documentation(info = "<html>
    <H1 style=\"font-size:20px\">single outflow for ADM1 </H1>
    <p style=\"font-size:20px\">
  Outflow for the ADM1 as single model without combination with the ASM
    </p>
  </html>"));
  end outflow;


  model outflow_ASM
    parameter Real T = 15;
  OpenWasteWater.ADM_P.InFlow In1 annotation(
      Placement(visible = true, transformation(origin = {-76, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-56, 0}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
  OpenWasteWater.ASM1P.OutPipe Out1 annotation(
      Placement(visible = true, transformation(origin = {56, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {59, -8.88178e-16}, extent = {{-27, -28}, {27, 28}}, rotation = 0)));
  equation
    Out1.S.I = In1.S.I;
    Out1.S.S = In1.S.su + In1.S.aa + In1.S.fa + In1.S.va + In1.S.bu + In1.S.pro + In1.S.ac + In1.S.H2 + In1.S.ch4;
    Out1.X.I = In1.X.I;
    Out1.X.S = In1.X.c + In1.X.ch + In1.X.pr + In1.X.li + In1.X.su + In1.X.aa + In1.X.fa + In1.X.c4 + In1.X.pro + In1.X.ac + In1.X.H2;
    Out1.X.H = 0;
    Out1.X.A = 0;
    Out1.X.P = 0;
    Out1.S.O2 = 0;
    Out1.S.NO = 0;
    Out1.S.NH = In1.S.IN;
    Out1.S.ND = 0;
    Out1.X.ND = 0;
    Out1.S.ALK = In1.S.IC;
    Out1.S.P = In1.S.P;
    Out1.S.Fe2 = In1.S.Fe2;
    Out1.S.Fe3 = In1.S.Fe3;
    Out1.S.Al3 = In1.S.Al3;
    Out1.X.FeP = In1.X.FeP;
    Out1.X.FeOH = In1.X.FeOH;
    Out1.X.AlP = In1.X.AlP;
    Out1.X.AlOH = In1.X.AlOH;
  
    Out1.Q + In1.Q = 0.0;
    Out1.T = T;
  annotation(
    defaultComponentName = "outflow_ASM",
    Documentation(info = "<html>
    <H1 style=\"font-size:20px\">Interface ADM1 to ASM1</H1>
    <p style=\"font-size:20px\">
    Outflow for combination with the ASM1.

    Conversion according to the documentation of the ADM1 from the IWA, Appendix C.3, Table C.2

    Temperature T is set back to the wastewater temperature
    </p>
  </html>"),
  Icon(graphics = {Rectangle(origin = {1, -1}, lineThickness = 0.75, extent = {{-23, 59}, {23, -59}})}));
  end outflow_ASM;




class initial_values_ADM
protected
  Real CSTR[24, 2];
  Real GAS[4, 2];
equation
//Scenario1
//WWTP_ADM CSTR and Headspace
    CSTR[1, 1] = 12 "S.su";
  CSTR[2, 1] = 5.36 "S.aa";
  CSTR[3, 1] = 102.7 "S.fa";
  CSTR[4, 1] = 10.75 "S.va";
  CSTR[5, 1] = 14.3 "S.bu";
  CSTR[6, 1] = 16.85 "S.pro";
  CSTR[7, 1] = 44 "S.ac";
  CSTR[8, 1] = 0.00024 "S.H2";
  CSTR[9, 1] = 164 "S.ch4";
  CSTR[10, 1] = 63.85 "S.IC";
  CSTR[11, 1] = 721 "S.IN";
  CSTR[12, 1] = 6320 "S.I";
  CSTR[13, 1] = 6250 "X.c";
  CSTR[14, 1] = 62.5 "X.ch";
  CSTR[15, 1] = 62.5 "X.pr";
  CSTR[16, 1] = 78 "X.li";
  CSTR[17, 1] = 948 "X.su";
  CSTR[18, 1] = 714 "X.aa";
  CSTR[19, 1] = 632 "X.fa";
  CSTR[20, 1] = 308 "X.c4";
  CSTR[21, 1] = 146 "X.pro";
  CSTR[22, 1] = 933 "X.ac";
  CSTR[23, 1] = 448 "X.H2";
  CSTR[24, 1] = 40250 "X.I";
//GAS
    GAS[1, 1] = 0.0043 "S_H2g";
  GAS[2, 1] = 2310 "S_ch4g";
  GAS[3, 1] = 20.1 "S_CO2g";
  GAS[4, 1] = 0.0047 "S_NH3g";
//Scenario2
//PC_WWTP_ADM CSTR and Headspace
    CSTR[1, 2] = 12 "S.su";
  CSTR[2, 2] = 5.36 "S.aa";
  CSTR[3, 2] = 102.7 "S.fa";
  CSTR[4, 2] = 10.75 "S.va";
  CSTR[5, 2] = 14.3 "S.bu";
  CSTR[6, 2] = 16.85 "S.pro";
  CSTR[7, 2] = 39 "S.ac";
  CSTR[8, 2] = 0.00024 "S.H2";
  CSTR[9, 2] = 152 "S.ch4";
  CSTR[10, 2] = 52.4 "S.IC";
  CSTR[11, 2] = 527 "S.IN";
  CSTR[12, 2] = 4470 "S.I";
  CSTR[13, 2] = 4400 "X.c";
  CSTR[14, 2] = 44 "X.ch";
  CSTR[15, 2] = 44 "X.pr";
  CSTR[16, 2] = 55 "X.li";
  CSTR[17, 2] = 669 "X.su";
  CSTR[18, 2] = 504 "X.aa";
  CSTR[19, 2] = 445 "X.fa";
  CSTR[20, 2] = 217 "X.c4";
  CSTR[21, 2] = 103 "X.pro";
  CSTR[22, 2] = 657 "X.ac";
  CSTR[23, 2] = 316 "X.H2";
  CSTR[24, 2] = 25400 "X.I";
//GAS
    GAS[1, 2] = 0.0059 "S_H2g";
  GAS[2, 2] = 2750 "S_ch4g";
  GAS[3, 2] = 23.5 "S_CO2g";
  GAS[4, 2] = 0.0022 "S_NH3g";

  annotation(
    defaultComponentName = "Scenarios_ADM",
    Documentation(info = "<html>
<H1 style=\"font-size:20px\">Scenarios_ADM</H1>
<p style=\"font-size:20px\">
This class contains the initial values for all the two models with ADM in order to avoid changing the values frequently in the 'CSTR_ADM' and 'Headspace_ADM' models.
 New scenarios can be added by expanding the number of columns n in the matrix CSTR[m,n].

Scenario 1 - WWTP_ADM model with ES = 385 m3/d
Scenario 2 - PC_WWTP_ADM model with ES = 385 m3/d
</p>
</html>"));
end initial_values_ADM;

package adm_models

model ADM1
      extends OpenWasteWater.Icons.Digester;
      parameter Real V_D = 1200.0 "m3 volume digester";
      parameter Real V_H = 120.0 "m3 volume headspace";
      parameter Real kLa = 20 "1/d";
      parameter OpenWasteWater.ADM_P.Solubles_ADM Sini(
       su = 12,
       aa = 5.36,
       fa = 102.7,
       va = 10.75, 
       bu = 14.3,
       pro = 16.85,
       ac = 39,
       H2 = 0.00024,
       ch4 = 152,
       IC = 52.4,
       IN = 527,
       I = 6400,
       P = 1.0,
       Fe2 = 0,
       Fe3 = 0,
       Al3 = 0);
      parameter OpenWasteWater.ADM_P.Particulates_ADM Xini(
       c = 4400,
       ch = 44,
       pr = 44,
       li = 55,
      su = 669,
      aa = 504,
      fa = 445,
      c4 = 217,
      pro = 103,
      ac = 657,
      H2 = 316,
      I = 34000,
      FeP = 0,
      FeOH = 0,
      AlP = 0,
      AlOH = 0);

      parameter OpenWasteWater.ADM_P.Gaseous_ADM Gini(CO2=50.5, H2=0.0078, NH3=0.00156, CH4=6260.0);      

      OpenWasteWater.ADM_P.CSTR_ADM R1(V_R = V_D, kLa = kLa, Sini = Sini, Xini = Xini);
      OpenWasteWater.ADM_P.Headspace_ADM Gas(V_G = V_H, Gini=Gini);
      OpenWasteWater.ADM_P.InFlow In1 annotation(
        Placement(visible = true, transformation(origin = {-92, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, -7}, extent = {{-14, -15}, {14, 15}}, rotation = 0)));
      OpenWasteWater.ADM_P.OutFlow Out1 annotation(
        Placement(visible = true, transformation(origin = {86, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, -3}, extent = {{-14, -15}, {14, 15}}, rotation = 0)));
      OpenWasteWater.ADM_P.outBiogasSensor OutBiogasSensor annotation(
        Placement(visible = true, transformation(origin = {22, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {22, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      OutBiogasSensor.Value.Qnorm = Gas.q_gas_norm "dry biogas flow norma conditions";
      OutBiogasSensor.Value.xCO2  = Gas.P_CO2/(Gas.P_CO2+Gas.P_ch4+Gas.P_H2+Gas.P_NH3)*100.0;
      OutBiogasSensor.Value.xH2   = Gas.P_H2 /(Gas.P_CO2+Gas.P_ch4+Gas.P_H2+Gas.P_NH3)*100.0;
      OutBiogasSensor.Value.xCH4  = Gas.P_ch4/(Gas.P_CO2+Gas.P_ch4+Gas.P_H2+Gas.P_NH3)*100.0;
      OutBiogasSensor.Value.xNH3  = Gas.P_NH3/(Gas.P_CO2+Gas.P_ch4+Gas.P_H2+Gas.P_NH3)*100.0;
      connect(In1, R1.In1);
      connect(R1.OutG, Gas.InG);
      connect(Gas.OutP, R1.InP);
      connect(R1.Out1, Out1);
      annotation(
        defaultComponentName = "ADM1",
        Documentation(info = "<html>
<H1 style=\"font-size:20px\">ADM1</H1>
<p style=\"font-size:20px\">
This Model is the single ADM1 model based on the Anearobic Sludge Model No.1 by the IWA
</html>"),
        experiment(StartTime = 0, StopTime = 14, Tolerance = 1e-08, Interval = 0.01));

end ADM1;

end adm_models;

end ADM_P;