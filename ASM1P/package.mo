within OpenWasteWater;

package ASM1P "Component models for the Activated Sludge Model No.1 with P-precipitation"
  partial model kinetic
    parameter Real muemH15 = 4.0 "d-1";
    parameter Real Y_H = 0.67 "gXH/gSS";
    parameter Real K_S = 10.0 "gCOD/m3";
    parameter Real K_O2H = 0.2 "gO2/m3";
    parameter Real K_NO = 0.5 "gO2/m3";
    parameter Real K_P = 0.01 "gP/m3";
    parameter Real eta_D = 0.8 "-";
    parameter Real b_H15 = 0.3 "d-1";
    parameter Real muemA15 = 0.5 "d-1";
    parameter Real Y_A = 0.24 "gXA/gSNH";
    parameter Real fp = 0.08 "gXP/gXB";
    parameter Real iXB = 0.08 "gSNH/gXB";
    parameter Real iXP = 0.06 "-";
    parameter Real b_A15 = 0.05 "d-1";
    parameter Real eta_H = 0.8 "-";
    parameter Real k_h15 = 3.0 "d-1";
    parameter Real K_X15 = 0.1 "-";
    parameter Real K_NH = 1.0 "gN/m3";
    parameter Real K_O2A = 0.4 "gO2/m3";
    parameter Real K_ALK = 0.5 "mol/m3";
    parameter Real k_a15 = 0.05 "m3/(g d)";
    parameter Real iPBM = 0.01 "gP/gX.BM";
    parameter Real iPb = 0.01 "gP/gX.BM";
    parameter Real k_PRE = 1.2 "1/(gMe(OH)3 d)";
    parameter Real k_RED = 0.6 "1/d";
    parameter Real k_Feox = 400.0 "1/d";
    parameter Real k_POH = 1000.0 "1/d";
    Real S_O2sat, TT;
    Real muemH, b_H, muemA, b_A, k_a, k_h, K_X;
  algorithm
    S_O2sat := fS_O2sat(TT);
    muemH := muemH15 * exp(0.069 * (TT - 15));
    b_H := b_H15 * exp(0.069 * (TT - 15));
    muemA := muemA15 * exp(0.098 * (TT - 15));
    b_A := b_A15 * exp(0.08 * (TT - 15));
    k_a := k_a15 * exp(0.069 * (TT - 15));
    k_h := k_h15 * exp(0.11 * (TT - 15));
    K_X := K_X15 * exp(0.11 * (TT - 15));
    annotation(
      defaultComponentName = "kinetic",
      Documentation(info = "<html>
            <H1 style=\"font-size:20px\">Kinetic coefficient for ASM1 </H1>
            <p style=\"font-size:20px\">
            This class ''kinetic'' contains all coefficients regarding the IWA
            Activated Sludge Model No.1 with the suggested default values. For practical
            application these values have do be adoted with model kalibration.
            </p>
        </html>"));
  end kinetic;

  function fS_O2sat
    input Real T;
    output Real S_O2sat;
  algorithm
    S_O2sat := exp(7.7117 - 1.31403 * log(T + 45.93));
    annotation(
      defaultComponentName = "fS_O2sat",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Saturation concnetration of oygen </H1>
      <p style=\"font-size:20px\">
      This function determines the saturation concentration of
      oxygen in water aerated with air as function
      of the temperature<br />
      (T in °C S_O2sat in g/m3)
      </p>
    </html>"));
  end fS_O2sat;

  function fJS "Sedimentation velocity function"
    input Real X "g/m3";
    //input Real ISV "ml/g";
    output Real JS "g/(m2 d)";
    //    parameter Real ISV = 120.0 "ml/g";
    //    parameter Real h_in = 2.0 "m";
  protected
    Real v0str, v0 "maximum settling velocity";
    Real nv "exponent as part of the Vesilind equation";
    Real XTSS;
    Real rh, rp;
  algorithm
    v0 := 474.0 "m/d";
    v0str := 250.0 "m/d";
    rh := 0.000576 "m3/(g SS)";
    rp := 0.00286 "m3/(g SS)";
    XTSS := X * 0.75;
    JS := min(v0str, v0 * exp(-rh * XTSS) - v0 * exp(-rp * XTSS)) * XTSS / 0.75 "g/(m2 d)";
    annotation(
      defaultComponentName = "fJS",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Sludge Sedimentation Flux </H1>
      <p style=\"font-size:20px\">
      This function determines the sludge sedimentation flux according to
      Takacs. Input ist the concetration of the particulat compounents in
      g COD / m3. Output in the sludge flux in g COD (m2 d).
      </p>
    </html>"));
  end fJS;

  function fkLa
    output Real kLa "specific mass transfer coefficient, 1/d";
    input Real Q "air flow rate, m³/d";
    input Real V_R = 1333 "reactor volume, m³";
    input Real H = 4.5 "fluid height, m";
  protected
    Real kLa_stern "-";
    Real v = 10 ^ (-6) "kinematic viscosity of water, m²/s";
    Real g = 9.81 "gravity, m/s²";
    Real w "m/d";
    Real w_stern "-";
  algorithm
    w := Q * H / V_R;
    w_stern := w / (84600 * (g * v) ^ (1 / 3));
    kLa_stern := 1.17 * 10 ^ (-4) * w_stern ^ (-0.1);
    kLa := kLa_stern * w / (v ^ 2 / g) ^ (1 / 3);
    annotation(
      defaultComponentName = "klacalc",
      Documentation(info = "<html>
  <H1 style=\"font-size:20px\">kLa calculation</H1>
  <p style=\"font-size:20px\">
 This function calculates the specific mass transfer coefficient kLa in '1/d' for a defined air flow rate Q_air.
 The division by 84600 in the term 'w_stern' is the conversion from m/d to m/s in order to fit with the units of gravity [m/s²] and viscosity [m²/s].
 The equations are valid for a gas distribution systems with perforated bottom, sintered plate or frit.
 Based on
  'M. Zlokarnik. Verfahrenstechnische Grundlagen der reaktorgestaltung. Acta Biotechnologica 1, 1981.'
  </p> </html>"));
  end fkLa;

  record Soluble
    Real O2 "gO2/m3 dissolved oxygen";
    Real I "gCOD/m3 inert soluble organic material";
    Real S "gCOD/m3 readily biodegradable organic substances";
    Real NH "gN/m3 ammonium + ammonia N";
    Real NO "gN/m3 nitrite + nitrate N";
    Real ND "gN/m3 dissolved organic N";
    Real ALK "mol/m3 alkalinity";
    Real P "gP/m3";
    Real Fe2 "gFe2/m3";
    Real Fe3 "gFe3/m3";
    Real Al3 "gAl3/m3";
  end Soluble;

  record Particulate
    Real H "gCOD/m3 heterotrophic bacteria";
    Real A "gCOD/m3 autotrophic bacteria";
    Real S "gCOD/m3 particulate slowly degradable substrates";
    Real I "gCOD/m3 inert particulate organic material";
    Real P "gCOD/m3 inert particulate organic matter resulting from decay";
    Real ND "gN/m3 particulate organic nitrogen";
    Real FeP "gFePO4/m3";
    Real FeOH "gFe(OH)3/m3";
    Real AlP "gAlPO4/m3";
    Real AlOH "gAl(OH)3/m3";
  end Particulate;

  connector InPipe
    Real T;
    flow Real Q;
    input OpenWasteWater.ASM1P.Soluble S;
    input OpenWasteWater.ASM1P.Particulate X;
    annotation(
      defaultComponentName = "InPipe",
      Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Text(extent = {{-88, 92}, {88, -94}}, lineColor = {255, 255, 255}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "%name")}),
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Connector (for ASM1 componets and Flow. Potential variable is the Temerature </H1>
      <p style=\"font-size:20px\">
      This connector defines a connector for the IWA Activated Sludge Model No. 1
      for all 13 components as particulate and dissolved variables and the flowrate. To fullfill the
      modelica requirement the temperature is added as potential variable.
      </p>
    </html>"));
  end InPipe;

  connector OutPipe
    Real T;
    flow Real Q;
    output OpenWasteWater.ASM1P.Soluble S;
    output OpenWasteWater.ASM1P.Particulate X;
    annotation(
      defaultComponentName = "OutPipe",
      Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Text(extent = {{-88, 92}, {94, -92}}, textString = "%name")}),
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Connector (for ASM1 componets and Flow. Potential variable is the Temerature </H1>
      <p style=\"font-size:20px\">
      This connector defines a connector for the IWA Activated Sludge Model No. 1
      for all 13 compounents as particulete and dissolved variables and the flowrate. To fullfill the
      modelica requirement the temperature is added as potential variable.
      </p>
    </html>"));
  end OutPipe;

  connector inQ
    input Real Q;
    annotation(
      defaultComponentName = "InQ",
      Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
      Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {173, 255, 47}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {0, 0, 15}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Connector (for a Flow only (Air and Water)  </H1>
      <p style=\"font-size:20px\">
      Used for Airflow to reactors and pump flows (return sludge, waste sludge and recirculation.
      </p>
    </html>"));
  end inQ;

  connector outQ
    output Real Q;
    annotation(
      defaultComponentName = "OutQ",
      Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
      Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {248, 248, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {155, 155, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Connector (for a Flow only (Air and Water)  </H1>
      <p style=\"font-size:20px\">
      Used for Airflow to reactors and pump flows (return sludge, waste sludge and recirculation.
      </p>
    </html>"));
  end outQ;

  partial model WWParameters
    Real TSS, COD, BOD, NH4_N, NO3_NO2_N, TKN, NT, PO4_P, TP;
    parameter Real wfp = 0.08 "gXP/gXB";
    parameter Real wiXB = 0.08 "gSNH/gXB";
    parameter Real wiXP = 0.06 "-";
    parameter Real wiPBM = 0.01 "gP/gX.BM";
    Soluble S;
    Particulate X;
  equation
    TSS = 0.75 * (X.I + X.S + X.H + X.A + X.P) + X.FeOH + X.FeP + X.AlOH + X.AlP;
    COD = X.I + X.S + X.H + X.A + X.P + S.I + S.S;
    BOD = 0.25 * (S.S + X.S + (1 - wfp) * (X.H + X.A));
    NH4_N = S.NH;
    NO3_NO2_N = S.NO;
    TKN = X.ND + S.ND + S.NH + wiXB * (X.H + X.A) + wiXP * (X.I + X.P);
    NT = TKN + S.NO;
    PO4_P = S.P;
    TP = S.P + wiPBM * (X.I + X.S + X.H + X.A + X.P) + 0.2054 * X.FeP + 0.254 * X.AlP;
  end WWParameters;

  model Mixer2
    extends OpenWasteWater.Icons.mixer2;
    InPipe In1 annotation(
      Placement(transformation(extent = {{-110, 15}, {-90, 35}}, rotation = 0))), In2 annotation(
      Placement(transformation(extent = {{-110, -25}, {-90, -5}}, rotation = 0)));
    OutPipe Out1 annotation(
      Placement(transformation(extent = {{90, -5}, {110, 15}}, rotation = 0)));
  equation
    0 = Out1.Q + In1.Q + In2.Q;
  //  Out1.T = In1.T;
    0 = Out1.T   * Out1.Q + In1.T   * In1.Q + In2.T   * In2.Q;
    0 = Out1.S.I * Out1.Q + In1.S.I * In1.Q + In2.S.I * In2.Q;
    0 = Out1.S.S * Out1.Q + In1.S.S * In1.Q + In2.S.S * In2.Q;
    0 = Out1.X.I * Out1.Q + In1.X.I * In1.Q + In2.X.I * In2.Q;
    0 = Out1.X.S * Out1.Q + In1.X.S * In1.Q + In2.X.S * In2.Q;
    0 = Out1.X.H * Out1.Q + In1.X.H * In1.Q + In2.X.H * In2.Q;
    0 = Out1.X.A * Out1.Q + In1.X.A * In1.Q + In2.X.A * In2.Q;
    0 = Out1.X.P * Out1.Q + In1.X.P * In1.Q + In2.X.P * In2.Q;
    0 = Out1.X.FeP * Out1.Q + In1.X.FeP * In1.Q + In2.X.FeP * In2.Q;
    0 = Out1.X.AlP * Out1.Q + In1.X.AlP * In1.Q + In2.X.AlP * In2.Q;
    0 = Out1.X.FeOH * Out1.Q + In1.X.FeOH * In1.Q + In2.X.FeOH * In2.Q;
    0 = Out1.X.AlOH * Out1.Q + In1.X.AlOH * In1.Q + In2.X.AlOH * In2.Q;
    0 = Out1.S.O2 * Out1.Q + In1.S.O2 * In1.Q + In2.S.O2 * In2.Q;
    0 = Out1.S.NO * Out1.Q + In1.S.NO * In1.Q + In2.S.NO * In2.Q;
    0 = Out1.S.NH * Out1.Q + In1.S.NH * In1.Q + In2.S.NH * In2.Q;
    0 = Out1.S.ND * Out1.Q + In1.S.ND * In1.Q + In2.S.ND * In2.Q;
    0 = Out1.X.ND * Out1.Q + In1.X.ND * In1.Q + In2.X.ND * In2.Q;
    0 = Out1.S.ALK * Out1.Q + In1.S.ALK * In1.Q + In2.S.ALK * In2.Q;
    0 = Out1.S.Fe2 * Out1.Q + In1.S.Fe2 * In1.Q + In2.S.Fe2 * In2.Q;
    0 = Out1.S.Fe3 * Out1.Q + In1.S.Fe3 * In1.Q + In2.S.Fe3 * In2.Q;
    0 = Out1.S.Al3 * Out1.Q + In1.S.Al3 * In1.Q + In2.S.Al3 * In2.Q;
    0 = Out1.S.P * Out1.Q + In1.S.P * In1.Q + In2.S.P * In2.Q;
    annotation(
      defaultComponentName = "Mixer2",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Mixer2 </H1>
      <p style=\"font-size:20px\">
      This model mixes 2 flows of connector ( in1 and in2) for the IWA Activated Sludge Model No. 1
      for all 13 compounents and connect the result in the connector Pipe Out1.
      </p>
    </html>"));
  end Mixer2;

  model Mixer3
    extends OpenWasteWater.Icons.mixer3;
    Mixer2 M1 annotation(
      Placement(visible = true, transformation(origin = {-26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.Mixer2 M2 annotation(
      Placement(visible = true, transformation(origin = {8, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.InPipe In1 annotation(
      Placement(visible = true, transformation(origin = {-100, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.InPipe In2 annotation(
      Placement(visible = true, transformation(origin = {-98, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.InPipe In3 annotation(
      Placement(visible = true, transformation(origin = {-100, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OutPipe Out1 annotation(
      Placement(visible = true, transformation(origin = {96, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {96, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(In2, M1.In2) annotation(
      Line(points = {{-98, -4}, {-36, -4}, {-36, 38}}, color = {0, 0, 255}));
    connect(In3, M2.In2) annotation(
      Line(points = {{-100, -46}, {-2, -46}, {-2, 20}}, color = {0, 0, 255}));
    connect(In1, M1.In1) annotation(
      Line(points = {{-100, 36}, {-68, 36}, {-68, 42}, {-36, 42}}, color = {0, 0, 255}));
    connect(M2.Out1, Out1) annotation(
      Line(points = {{18, 22}, {90, 22}, {90, 0}, {96, 0}}));
    connect(M1.Out1, M2.In1) annotation(
      Line(points = {{-16, 40}, {-2, 40}, {-2, 24.5}}));
  end Mixer3;

  model Divider2
    extends OpenWasteWater.Icons.divider2;
    OpenWasteWater.ASM1P.InPipe In1 annotation(
      Placement(visible = true, transformation(extent = {{-110, -7}, {-90, 13}}, rotation = 0), iconTransformation(origin = {-96, 3}, extent = {{-14, -16}, {14, 16}}, rotation = 0)));
    OpenWasteWater.ASM1P.OutPipe Out2 annotation(
      Placement(visible = true, transformation(extent = {{90, -25}, {110, -5}}, rotation = 0), iconTransformation(origin = {90, -15}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    OpenWasteWater.ASM1P.OutPipe Out1 annotation(
      Placement(visible = true, transformation(extent = {{90, 16}, {110, 36}}, rotation = 0), iconTransformation(origin = {90, 26}, extent = {{-12, -14}, {12, 14}}, rotation = 0)));
  equation
    0 = Out1.Q + Out2.Q + In1.Q;
    Out1.T = In1.T;
    Out2.T = In1.T;
    Out1.S = In1.S;
    Out1.X = In1.X;
    Out2.S = In1.S;
    Out2.X = In1.X;
    annotation(
      defaultComponentName = "Split2",
      Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Split2 (for Water and Sludge) </H1>
      <p style=\"font-size:20px\">
      This model splits a flows of connector in1 for the IWA Activated Sludge Model No. 1
      for all 13 compounents and connect the result in the connector 2 connectors (out1 and ou2).
      The model needs the parameter Q2 for the flowrate of the connector out2.
      </p>
    </html>"));
  end Divider2;

  model Inflow_simple "inflow ASM1"
    extends OpenWasteWater.Icons.WWSource;
    parameter Real InO2 = 0;
    parameter Real InQ = 18000.0;
    parameter Real InCOD = 570.0;
    parameter Real InN = 49.1;
    parameter Real InALK = 7.0;
    parameter Real InP = 9.0;
    parameter Real InFe3 = 0.0;
    parameter Real T = 20.0;
    Real Q;
    OutPipe Out1 annotation(
      Placement(visible = true, transformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    Q = InQ;
    Out1.S.I = 0.089 * InCOD;
    Out1.S.S = 0.300 * InCOD;
    Out1.X.I = 0.000 * InCOD;
    Out1.X.S = 0.600 * InCOD;
    Out1.X.H = 0.015 * InCOD;
    Out1.X.A = 0.005 * InCOD;
    Out1.X.P = 0.000 * InCOD;
    Out1.X.FeP = 0.000 * InCOD;
    Out1.X.AlP = 0.000 * InCOD;
    Out1.X.FeOH = 0.000 * InCOD;
    Out1.X.AlOH = 0.000 * InCOD;
    Out1.S.O2 = InO2;
    Out1.S.NH = 0.623 * InN;
    Out1.S.NO = 0.000;
    Out1.S.ND = 0.142 * InN;
    Out1.X.ND = 0.235 * InN;
    Out1.S.ALK = InALK;
    Out1.S.P = InP;
    Out1.S.Fe2 = 0.0;
    Out1.S.Fe3 = InFe3;
    Out1.S.Al3 = 0.0;
    Out1.Q = -abs(Q);
    Out1.T = T;
  end Inflow_simple;

  model InflowSludgeTest "inflow ASM1"
    extends OpenWasteWater.Icons.WWSource;
    parameter Real T = 38.0;
    Real Q;
    OutPipe Out1 annotation(
      Placement(visible = true, transformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    Q = 60.0 + 10.0 * sin(2 * Modelica.Constants.pi * time) "m3/d";
    Out1.S.I = 30.0;
    Out1.S.S = 26.0;
    Out1.X.I = 20000.0;
    Out1.X.S = 21900.0;
    Out1.X.H = 27980.0;
    Out1.X.A = 1394.0;
    Out1.X.P = 8300.0;
    Out1.S.O2 = 0.0;
    Out1.S.NH = 11.6;
    Out1.S.NO = 9.1;
    Out1.S.ND = 2.94;
    Out1.X.ND = 1100.0;
    Out1.S.ALK = 4.92;
    Out1.S.P = 0.5;
    Out1.S.Fe2 = 0.0;
    Out1.S.Fe3 = 0.0;
    Out1.S.Al3 = 0.0;
    Out1.X.FeP = 20.0;
    Out1.X.FeOH = 20.0;
    Out1.X.AlP = 20.0;
    Out1.X.AlOH = 20.0;
    Out1.Q = -abs(Q);
    Out1.T = T;
  end InflowSludgeTest;

  model Dosing "inflow ASM1"
    extends OpenWasteWater.Icons.WWSource;
    parameter Real InFe2 = 0.0;
    parameter Real InFe3 = 3000.0;
    parameter Real InAl3 = 0.0;
    parameter Real betaP = 1.5;
    Real QFe;
    OpenWasteWater.ASM1P.OutPipe Out1 annotation(
      Placement(visible = true, transformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.TechUnits.inWWSensor SensorInfDos annotation(
      Placement(visible = true, transformation(origin = {-90, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
// control Fe Pump
    QFe = SensorInfDos.Value.P * SensorInfDos.Value.Q / ((InFe2 + InFe3) / 1.8 + InAl3 / 0.871) * betaP;
    Out1.S.I = 0.0;
    Out1.S.S = 0.0;
    Out1.X.I = 0.0;
    Out1.X.S = 0.0;
    Out1.X.H = 0.0;
    Out1.X.A = 0.0;
    Out1.X.P = 0.0;
    Out1.X.FeP = 0.0;
    Out1.X.AlP = 0.0;
    Out1.X.FeOH = 0.0;
    Out1.X.AlOH = 0.0;
    Out1.S.O2 = 0.0;
    Out1.S.NH = 0.0;
    Out1.S.NO = 0.0;
    Out1.S.ND = 0.0;
    Out1.X.ND = 0.0;
    Out1.S.ALK = 0.0;
    Out1.S.P = 0.0;
    Out1.S.Fe2 = InFe2;
    Out1.S.Fe3 = InFe3;
    Out1.S.Al3 = InAl3;
    Out1.Q = -abs(QFe);
    Out1.T = 15.0;
  end Dosing;

  model Inflow "inflow ASM1"
    extends OpenWasteWater.Icons.WWSource;
    parameter Real T = 15 "°C";
    parameter String Inf_File = "Resources/ASM1P/Inf_strm_P.txt";
    Modelica.Blocks.Sources.CombiTimeTable T_data(tableOnFile = true, fileName = Inf_File, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, tableName = "t_data", columns = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, startTime = 0.0, timeScale = 1, offset = {0}, verboseRead = true);
    OutPipe Out1 annotation(
      Placement(visible = true, transformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.TechUnits.outWWSensor OutSensor1 annotation(
      Placement(visible = true, transformation(origin = {56, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {81, 41}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    extends WWParameters;
    Real I_SP, M_SP, I_COD, M_COD, I_TSS, M_TSS, I_NH, M_NH, I_NT, M_NT;
  equation
    S.I = T_data.y[1];
    S.S = T_data.y[2];
    X.I = T_data.y[3];
    X.S = T_data.y[4];
    X.H = T_data.y[5];
    X.A = T_data.y[6];
    X.P = T_data.y[7];
    S.O2 = T_data.y[8];
    S.NO = T_data.y[9];
    S.NH = T_data.y[10];
    S.ND = T_data.y[11];
    X.ND = T_data.y[12];
    S.ALK = T_data.y[13];
    S.P = T_data.y[14];
    S.Fe2 = 0.0;
    S.Fe3 = 0.0;
    S.Al3 = 0.0;
    X.FeP = 0.0;
    X.AlP = 0.0;
    X.FeOH = 0.0;
    X.AlOH = 0.0;
    Out1.S = S;
    Out1.X = X;
    Out1.Q = -abs(T_data.y[15]);
    Out1.T = T;
    der(I_SP) = S.P;
    M_SP = I_SP / time;
    der(I_COD) = COD;
    M_COD = I_COD / time;
    der(I_TSS) = TSS;
    M_TSS = I_TSS / time;
    der(I_NH) = S.NH;
    M_NH = I_NH / time;
    der(I_NT) = NT;
    M_NT = I_NT / time;
  
    OutSensor1.Value.O2 = Out1.S.O2;
    OutSensor1.Value.NH = Out1.S.NH;
    OutSensor1.Value.NO = Out1.S.NO;
    OutSensor1.Value.S = Out1.S.S;
    OutSensor1.Value.P = Out1.S.P;
    OutSensor1.Value.Q = abs(T_data.y[15]);
    OutSensor1.Value.TSS = 0.75 * (Out1.X.I + Out1.X.S + Out1.X.H + Out1.X.A + Out1.X.P);
    annotation(
      defaultComponentName = "Inflow",
      Documentation(info = "<html>
            <H1 style=\"font-size:20px\">Inflow (Data from file) </H1>
            <p style=\"font-size:20px\">
            This Model reads inflow parameters of the ASM1 from a file. The format of
            the file has to fullfill the requirements of the
            Modelica.Blocks.Sources.CombiTimeTable.
            </p>
        </html>"));
  end Inflow;

  model Effluent
    extends OpenWasteWater.Icons.EffluentSink;
    InPipe In1 annotation(
      Placement(visible = true, transformation(origin = {-100, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends WWParameters;
    Real Q, T;
    Real I_SP, M_SP, I_COD, M_COD, I_TSS, M_TSS, I_NH, M_NH, I_NT, M_NT;
  initial equation
    I_SP = 0;
    I_COD = 0;
    I_TSS = 0;
    I_NH = 0;
    I_NT = 0;
  equation
    Q = In1.Q;
    T = In1.T;
    S = In1.S;
    X = In1.X;
    der(I_SP) = S.P;
    M_SP = I_SP / time;
    der(I_COD) = COD;
    M_COD = I_COD / time;
    der(I_TSS) = TSS;
    M_TSS = I_TSS / time;
    der(I_NH) = S.NH;
    M_NH = I_NH / time;
    der(I_NT) = NT;
    M_NT = I_NT / time;
    annotation(
      defaultComponentName = "SinkWater",
      Documentation(info = "<html>
                  <H1 style=\"font-size:20px\">SinkWater (for Water and Sludge) </H1>
                  <p style=\"font-size:20px\">
                  This Model recieves as inflow the values of an
                  ASMConnector to fullfill the condition that
                  the sum of all flows is equal zero<br />
                  The use of the inStream function is not required here,
                  but to have a look the the effluent, is helpfull.<br />
                  All units are g/m3 except temperature in °C and S_ALK im mol/m3
                  </p>
          </html>"));
  end Effluent;

  model WasteSludge
    extends OpenWasteWater.Icons.SludgeSink;
    InPipe In1 annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends WWParameters;
    Real Q, T;
    Real I_SP, M_SP, I_COD, M_COD, I_TSS, M_TSS, I_NH, M_NH, I_NT, M_NT;
  initial equation
    I_SP = 0;
    I_COD = 0;
    I_TSS = 0;
    I_NH = 0;
    I_NT = 0;
  equation
    Q = In1.Q;
    T = In1.T;
    S = In1.S;
    X = In1.X;
    der(I_SP) = S.P;
    M_SP = I_SP / time;
    der(I_COD) = COD;
    M_COD = I_COD / time;
    der(I_TSS) = TSS;
    M_TSS = I_TSS / time;
    der(I_NH) = S.NH;
    M_NH = I_NH / time;
    der(I_NT) = NT;
    M_NT = I_NT / time;
    annotation(
      defaultComponentName = "SinkWater",
      Documentation(info = "<html>
                  <H1 style=\"font-size:20px\">SinkWater (for Water and Sludge) </H1>
                  <p style=\"font-size:20px\">
                  This Model recieves as inflow the values of an
                  ASMConnector to fullfill the condition that
                  the sum of all flows is equal zero<br />
                  The use of the inStream function is not required here,
                  but to have a look the the effluent, is helpfull.<br />
                  All units are g/m3 except temperature in °C and S_ALK im mol/m3
                  </p>
          </html>"));
  end WasteSludge;

  model NitrificationTank
    extends OpenWasteWater.Icons.nitri;
    parameter Real H = 6 "m";
    extends ASM1_CSTR;
    inQ InQair annotation(
      Placement(visible = true, transformation(origin = {0, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    kLa = fkLa(Q = InQair.Q, V_R = V_R, H = H);
  end NitrificationTank;

  model DenitrificationTank
    extends OpenWasteWater.Icons.deni;
    parameter Real kkLa = 2.0 "1/d  (for deni tanks 2 - 4 1/d)";
    extends ASM1_CSTR;
  equation
    kLa = kkLa;
    annotation(
      Diagram);
  end DenitrificationTank;

  partial model ASM1_CSTR
    OpenWasteWater.ASM1P.InPipe In1 annotation(
      Placement(transformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0)));
    OpenWasteWater.ASM1P.OutPipe Out1 annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}}, rotation = 0)));
    parameter Real V_R = 1000 "m3";
    parameter OpenWasteWater.ASM1P.Soluble Sini(I = 30, S = 1.15, O2 = 4.0, NH = 0.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0);
    parameter OpenWasteWater.ASM1P.Particulate Xini(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0);
    extends OpenWasteWater.ASM1P.kinetic;
    extends OpenWasteWater.ASM1P.WWParameters;
    Real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, rA;
    Real kLa;
    OpenWasteWater.ASM1P.TechUnits.outWWSensor OutSensor1 annotation(
      Placement(visible = true, transformation(origin = {52, 92}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {53, 53}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
  initial equation
    S = Sini;
    X = Xini;
  equation
// Temperature adaption
    TT = In1.T;
    r1 = muemH * S.S / (K_S + S.S) * S.O2 / (K_O2H + S.O2) * S.P / (K_P + S.P) * X.H;
    r2 = muemH * S.S / (K_S + S.S) * S.P / (K_P + S.P) * S.NO / (K_NO + S.NO) * K_O2H / (K_O2H + S.O2) * eta_D * X.H;
    r3 = muemA * S.NH / (K_NH + S.NH) * S.O2 / (K_O2A + S.O2) * S.P / (K_P + S.P) * X.A;
    r4 = b_H * X.H;
    r5 = b_A * X.A;
    r6 = k_a * S.ND * X.H;
    r7 = k_h * X.S / X.H / (K_X + X.S / X.H) * (S.O2 / (K_O2H + S.O2) + eta_H * K_O2H / (K_O2H + S.O2) * S.NO / (K_NO + S.NO)) * X.H;
    r8 = r7 * X.ND / X.S;
    r9 = k_Feox * S.Fe2 * S.O2;
    r10 = k_PRE * X.FeOH * S.P;
    r11 = k_PRE * X.AlOH * S.P;
    r12 = k_RED * X.FeP * S.ALK / (K_ALK + S.ALK);
    r13 = k_RED * X.AlP * S.ALK / (K_ALK + S.ALK);
    r14 = k_POH * S.Fe3 * S.ALK / (K_ALK + S.ALK);
    r15 = k_POH * S.Al3 * S.ALK / (K_ALK + S.ALK);
    rA = kLa * (S_O2sat - S.O2);
    der(S.I) = In1.Q / V_R * (In1.S.I - S.I);
    der(S.S) = In1.Q / V_R * (In1.S.S - S.S) - 1 / Y_H * r1 - 1 / Y_H * r2 + r7;
    der(X.I) = In1.Q / V_R * (In1.X.I - X.I);
    der(X.S) = In1.Q / V_R * (In1.X.S - X.S) - r7 + (1 - fp) * (r4 + r5);
    der(X.H) = In1.Q / V_R * (In1.X.H - X.H) + r1 + r2 - r4;
    der(X.A) = In1.Q / V_R * (In1.X.A - X.A) + r3 - r5;
    der(X.P) = In1.Q / V_R * (In1.X.P - X.P) + fp * (r4 + r5);
    der(S.O2) = In1.Q / V_R * (In1.S.O2 - S.O2) + (Y_H - 1) / Y_H * r1 + (Y_A - 4.57) / Y_A * r3 - 0.143 * r9 + rA;
    der(S.NO) = In1.Q / V_R * (In1.S.NO - S.NO) - (1 - Y_H) / (2.86 * Y_H) * r2 + 1 / Y_A * r3;
    der(S.NH) = In1.Q / V_R * (In1.S.NH - S.NH) - iXB * (r1 + r2) - (1 / Y_A + iXB) * r3 + r6;
    der(S.ND) = In1.Q / V_R * (In1.S.ND - S.ND) - r6 + r8;
    der(X.ND) = In1.Q / V_R * (In1.X.ND - X.ND) + (iXB - iXP * fp) * (r4 + r5) - r8;
    der(S.ALK) = In1.Q / V_R * (In1.S.ALK - S.ALK) - iXB / 14 * (r1 + r2 + r3) + (1 - Y_H) / (14 * 2.86 * Y_H) * r2 - 1 / 7 / Y_A * r3 + 1 / 14 * r6 - 3 / 55.845 * r14 - 3 / 26.982 * r15;
    der(S.P) = In1.Q / V_R * (In1.S.P - S.P) - r10 - r11 - iPBM * (r1 + r2 + r3) + iPb * (1 - fp) * (r4 + r5) + 0.002 * r7 + r12 / 4.869 + r13 / 3.927;
    der(S.Fe2) = In1.Q / V_R * (In1.S.Fe2 - S.Fe2) - r9;
    der(S.Fe3) = In1.Q / V_R * (In1.S.Fe3 - S.Fe3) - r14 + r9;
    der(S.Al3) = In1.Q / V_R * (In1.S.Al3 - S.Al3) - r15;
    der(X.FeP) = In1.Q / V_R * (In1.X.FeP - X.FeP) + 4.869 * r10 - r12;
    der(X.AlP) = In1.Q / V_R * (In1.X.AlP - X.AlP) + 3.927 * r11 - r13;
    der(X.FeOH) = In1.Q / V_R * (In1.X.FeOH - X.FeOH) + 1.914 * r14 + 0.709 * r12 - 3.45 * r10;
    der(X.AlOH) = In1.Q / V_R * (In1.X.AlOH - X.AlOH) + 2.93 * r15 + 0.754 * r13 - 2.54 * r11;
    Out1.S = S;
    Out1.X = X;
    0 = Out1.Q + In1.Q;
    Out1.T = In1.T;
    OutSensor1.Value.O2 = S.O2;
    OutSensor1.Value.NH = S.NH;
    OutSensor1.Value.NO = S.NO;
    OutSensor1.Value.S = S.S;
    OutSensor1.Value.P = S.P;
    OutSensor1.Value.Q = In1.Q;
    OutSensor1.Value.TSS = TSS;
    annotation(
      defaultComponentName = "AMS1_CSTR",
      Documentation(info = "<html>
                                                        <H1 style=\"font-size:20px\">ASM1_CSTR </H1>
                                                        <p style=\"font-size:20px\">
                                                        This Model is a model for an CSTR hat work like described for the ASM no. 1
                                                        Input parameters are the oxygen mass transfer coefficient kLa and the
                                                        reactor volume V_R. If the model works under nitrification or denitrification
                                                        conditions can be control with the kLa-value.
                                                        </p>
                        </html>"));
  end ASM1_CSTR;

  package TechUnits
    record MeasuredValues
      Real O2 "gO2/m3 dissolved oxygen";
      Real S "gCOD/m3 readily biodegradable organic substances";
      Real P "gP/m3 Po4-P";
      Real NH "gN/m3 ammonium + ammonia N";
      Real NO "gN/m3 nitrite + nitrate N";
      Real Q "m3/d flow rate ";
      Real TSS "g/m3 total suspended solids";
    end MeasuredValues;

    connector inWWSensor
      input MeasuredValues Value;
      annotation(
        defaultComponentName = "InWWSensor",
        Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {148, 0, 204}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {0, 0, 15}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
        Documentation(info = "<html>
        <H1 style=\"font-size:20px\">Connector (measured values in WWTP)  </H1>
        <p style=\"font-size:20px\">
        Used for control of flowrates to enhance the performance of the WWTP and/or operate economic.
        </p>
      </html>"));
    end inWWSensor;

    connector outWWSensor
      output MeasuredValues Value;
      annotation(
        defaultComponentName = "OutWWSensor",
        Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {255, 182, 203}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {155, 155, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
        Documentation(info = "<html>
        <H1 style=\"font-size:20px\">Connector (for a Flow only (Air and Water)  </H1>
        <p style=\"font-size:20px\">
        Used for control of flowrates to enhance the performance of the WWTP and/or operate economic.
        </p>
      </html>"));
    end outWWSensor;

    connector inV
      input Real V;
      annotation(
        defaultComponentName = "InV",
        Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {148, 0, 204}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {0, 0, 15}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
        Documentation(info = "<html>
        <H1 style=\"font-size:20px\">Connector (Volume (level) in Tank)  </H1>
        <p style=\"font-size:20px\">
        Used for control of flowrates to enhance the performance of the WWTP and/or operate economic.
        </p>
      </html>"));
    end inV;

    connector outV
      output Real V;
      annotation(
        defaultComponentName = "OutV",
        Window(x = 0.45, y = 0.01, width = 0.35, height = 0.49),
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {78, 154, 6}, fillColor = {255, 182, 203}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {155, 155, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-88, 92}, {88, -94}}, textString = "%name")}),
        Documentation(info = "<html>
        <H1 style=\"font-size:20px\">Connector (Volume (level) in Tank)  </H1>
        <p style=\"font-size:20px\">
        Used for control of flowrates to enhance the performance of the WWTP and/or operate economic.
        </p>
      </html>"));
    end outV;

    model rElement
      OpenWasteWater.ASM1P.TechUnits.inWWSensor In1;
      OpenWasteWater.ASM1P.TechUnits.outWWSensor Out1;
      parameter Real t_D = 0.001 "delay of signal in d";
      Real t_R;
    initial equation
      Out1.Value.S = 1.0;
      Out1.Value.P = 1.0;
      Out1.Value.NH = 1.0;
      Out1.Value.NO = 1.0;
      Out1.Value.O2 = 1.0;
      Out1.Value.TSS = 1.0;
    equation
      if t_D <= 0.000001 then
        t_R = 0.00001;
      else
        t_R = t_D;
      end if;
      der(Out1.Value.S) = (In1.Value.S - Out1.Value.S) / t_R;
      der(Out1.Value.P) = (In1.Value.P - Out1.Value.P) / t_R;
      der(Out1.Value.NH) = (In1.Value.NH - Out1.Value.NH) / t_R;
      der(Out1.Value.NO) = (In1.Value.NO - Out1.Value.NO) / t_R;
      der(Out1.Value.O2) = (In1.Value.O2 - Out1.Value.O2) / t_R;
      der(Out1.Value.TSS) = (In1.Value.TSS - Out1.Value.TSS) / t_R;
      Out1.Value.Q = In1.Value.Q;
    end rElement;

    model Sensor
      parameter Real t_D = 0.001 "delay of signal in d";
      OpenWasteWater.ASM1P.TechUnits.rElement T1(t_D = t_D / 6), T2(t_D = t_D / 6), T3(t_D = t_D / 6), T4(t_D = t_D / 6), T5(t_D = t_D / 6), T6(t_D = t_D / 6);
      OpenWasteWater.ASM1P.TechUnits.inWWSensor S_In1 annotation(
        Placement(visible = true, transformation(origin = {0, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.outWWSensor S_Out1 annotation(
        Placement(visible = true, transformation(origin = {2, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {2, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(S_In1, T1.In1);
      connect(T1.Out1, T2.In1);
      connect(T2.Out1, T3.In1);
      connect(T3.Out1, T4.In1);
      connect(T4.Out1, T5.In1);
      connect(T5.Out1, T6.In1);
      connect(T6.Out1, S_Out1);
      annotation(
        Icon(graphics = {Rectangle(origin = {-5, 9}, lineThickness = 1.5, extent = {{-35, 69}, {47, -85}}), Ellipse(origin = {1, 1}, lineThickness = 1.5, extent = {{-33, 33}, {33, -33}}, endAngle = 360), Line(origin = {9, 7}, points = {{-9, -9}, {7, 7}}, color = {32, 74, 135}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 6, smooth = Smooth.Bezier), Line(origin = {-21, 18}, points = {{-3, 2}, {3, -2}}, color = {164, 0, 0}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {0, 29}, points = {{0, 5}, {0, -5}}, color = {164, 0, 0}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {-27.3081, -0.345964}, points = {{-4, 0}, {4, 0}, {4, 0}}, color = {164, 0, 0}, thickness = 1), Line(origin = {29.051, -0.0329492}, points = {{-4, 0}, {4, 0}, {4, 0}}, color = {164, 0, 0}, thickness = 1), Line(origin = {21, 19}, points = {{3, 3}, {-3, -3}}, color = {164, 0, 0}, thickness = 1)}));
    end Sensor;

    model Pump
      extends OpenWasteWater.Icons.pump;
      parameter Real Qmax = 70000 "maximum performance, m3/d";
      parameter Real Qmin = 0 "minimum performance, m3/d";
      parameter Real Qini = 50.0 "initial flow rate m3/d";
      parameter Real k_RT = 5 "response time constant 1/d";
      Real Qactual "m3/d";
      InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-98, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {98, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inQ Qset annotation(
        Placement(visible = true, transformation(origin = {-100, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      Qactual = Qini;
    equation
      der(Qactual) = k_RT * ((if Qset.Q > Qmax then Qmax else if Qset.Q < Qmin then Qmin else Qset.Q) - Qactual);
      In1.Q = Qactual;
      In1.Q + Out1.Q = 0;
      In1.S = Out1.S;
      In1.X = Out1.X;
      In1.T = Out1.T;
      annotation(
        defaultComponentName = "pump",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">pumping unit</H1>
      <p style=\"font-size:20px\">
      The pump contains in and outflow connectors and is placed in between the components of the model.
      The outflow is equal to the input value Q given by the controller.

      The concentration of the compounds in the stream remains the same.

      The maximum and minimum pumping performance can not be exceeded by the pump.</p>
      </html>"));
    end Pump;

    model Blower
      extends OpenWasteWater.Icons.blower;
      parameter Real Qmax = 15361 "max. performance, m3/d";
      parameter Real Qmin = 0 "min. performance, m3/d";
      parameter Real Qini = 5000.0 "initial air flow rate m3/d";
      parameter Real k_RT = 50 "response time constant 1/d";
      inQ Qset annotation(
        Placement(visible = true, transformation(origin = {102, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      outQ Qair annotation(
        Placement(visible = true, transformation(origin = {-8, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-8, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      Qair.Q = Qini;
    equation
      der(Qair.Q) = k_RT * ((if Qset.Q > Qmax then Qmax else if Qset.Q < Qmin then Qmin else Qset.Q) - Qair.Q);
      annotation(
        defaultComponentName = "blower",
        Documentation(info = "<html>
  <H1 style=\"font-size:20px\">aeration unit</H1>
  <p style=\"font-size:20px\">
  The blower contains an outflow connectors and is placed at each of the activated sludge tanks.
  The outflow is equal to the input value Q given by the controller.

  The air flow rate will only be used to calculate the kLa in the tanks without any impact on the compound concentrations in the water.

  The maximum and minimum blowing performance was chosen to be the performance of reactor 3 and 4 and can not be exceeded by the blower.
  </p>
  </html>"));
    end Blower;

    model SetQ
      parameter Real Q = 100.0;
      OpenWasteWater.ASM1P.outQ OutQ1 annotation(
        Placement(visible = true, transformation(origin = {-102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-96, 1}, extent = {{-16, -17}, {16, 17}}, rotation = 0)));
    equation
      OutQ1.Q = Q;
      annotation(
        Icon(graphics = {Rectangle(lineThickness = 3, extent = {{-100, 40}, {100, -40}})}));
    end SetQ;

    model SludgeControl
      parameter Real SP_TSS = 60000.0 "g/m3 TSS for digestion";
      parameter Real SP_Q_min = 1 "lower limit SP_Q, m3/d";
      parameter Real SP_Q_max = 500 "upper limit SP_Q, m3/d";
      parameter Real k_P_Q = 2 "k-value proportional part";
      parameter Real k_I_Q = 10 "k-value integral part";
      parameter Real k_D_Q = 0 "k-value of the differential part";
      Real P_Q "proportional controller part";
      Real I_Q "integral controller part";
      Real D_Q "differential controler part";
      Real e_Q(start = 1.0) ", e_O2 control derivation, g/m3";
      Real e_Q_i "integral of control derivation";
      inWWSensor In1WWSensor annotation(
        Placement(visible = true, transformation(origin = {0, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      outQ Out1Q annotation(
        Placement(visible = true, transformation(origin = {0, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    algorithm

    initial equation
      e_Q_i = 0.0;
    equation
      e_Q = (-SP_TSS) + In1WWSensor.Value.TSS;
      der(e_Q_i) = e_Q;
      P_Q = k_P_Q * e_Q;
      I_Q = k_I_Q * e_Q_i;
      D_Q = k_D_Q * der(e_Q);
      Out1Q.Q = min(max(SP_Q_min, 20.0 + P_Q + I_Q + D_Q), SP_Q_max);
      annotation(
        Icon(graphics = {Text(origin = {-1, 0}, lineThickness = 3, extent = {{-37, 18}, {37, -18}}, textString = "Q-Control"), Rectangle(origin = {-3, 0}, lineThickness = 3, extent = {{-51, 38}, {51, -38}})}));
    end SludgeControl;

    model VolumeControl
      parameter Real SP_V = 50.0 "m3 Sludge for feed for digestion";
      parameter Real SP_Q_min = 10 "lower limit SP_Q, m3/d";
      parameter Real SP_Q_max = 200 "upper limit SP_Q, m3/d";
      parameter Real k_P_Q = 2 "k-value proportional part";
      parameter Real k_I_Q = 10 "k-value integral part";
      parameter Real k_D_Q = 0 "k-value of the differential part";
      Real P_Q "proportional controller part";
      Real I_Q "integral controller part";
      Real D_Q "differential controler part";
      Real e_Q ", e_O2 control derivation, g/m3";
      Real e_Q_i "integral of control derivation";
      inV In1V annotation(
        Placement(visible = true, transformation(origin = {0, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      outQ Out1Q annotation(
        Placement(visible = true, transformation(origin = {0, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    algorithm

    initial equation
      e_Q_i = 0;
    equation
      e_Q = (-SP_V) + In1V.V;
      der(e_Q_i) = e_Q;
      P_Q = k_P_Q * e_Q;
      I_Q = k_I_Q * e_Q_i;
      D_Q = k_D_Q * der(e_Q);
      Out1Q.Q = min(max(SP_Q_min, 60.0 + P_Q + I_Q + D_Q), SP_Q_max);
      annotation(
        Icon(graphics = {Text(origin = {-1, 0}, lineThickness = 3, extent = {{-37, 18}, {37, -18}}, textString = "Q-Control"), Rectangle(origin = {-3, 0}, lineThickness = 3, extent = {{-51, 38}, {51, -38}})}));
    end VolumeControl;

    model controller
      parameter Real Q_RS = 18841;
      parameter Real Q_WS = 385;
      parameter Real Q_REC = 55338;
      parameter Real Q_air_N1 = 15361 "constant Air flow to N1";
      parameter Real Q_air_N2 = 15361 "constant Air flow to N2";
      parameter Real Q_air_N3 = 10000 "mean Air flow to N3";
      parameter Real SP_NO = 1.0 "setpoint S.NO in DN2, gN/m3";
      parameter Real SP_NH = 3.0 "setpoint S.O2 in N2, gN/m3";
      Real SP_O2 "setpoint O2, gO2/m3";
      parameter Real SP_O2_min = 1 "lower limit SP_O2, gO2/m3";
      parameter Real SP_O2_max = 4 "upper limit SP_O2, gO2/m3";
      Real P_NO "proportional controller part";
      Real I_NO "integral controller part";
      Real D_NO "differential controler part";
      parameter Real k_P_NO = 2000 "k-value proportional part";
      parameter Real k_I_NO = 1000 "k-value integral part";
      parameter Real k_D_NO = 0 "k-value of the differential part";
      Real P_O2 "proportional controller part";
      Real I_O2 "integral controller part";
      Real D_O2 "differential controler part";
      parameter Real k_P_O2 = 2000 "k-value proportional part";
      parameter Real k_I_O2 = 1000 "k-value integral part";
      parameter Real k_D_O2 = 0 "k-value of the differential part";
      Real e_NO ", e_O2 control derivation, g/m3";
      Real e_NO_i "integral of control derivation";
      Real e_NH ", e_NH control derivation, g/m3";
      Real e_O2 ", e_O2 control derivation, g/m3";
      Real e_O2_i "integral of control derivation";
      Real f, h_O2, h_NO "help variable";
      MeasuredValues Inf, DN2, N2, N3;
      OpenWasteWater.ASM1P.outQ C_RS annotation(
        Placement(visible = true, transformation(origin = {-82, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.outQ C_WS annotation(
        Placement(visible = true, transformation(origin = {-50, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.outQ C_REC annotation(
        Placement(visible = true, transformation(origin = {-14, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.outQ CairR3 annotation(
        Placement(visible = true, transformation(origin = {14, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.outQ CairR4 annotation(
        Placement(visible = true, transformation(origin = {50, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.outQ CairR5 annotation(
        Placement(visible = true, transformation(origin = {84, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.inWWSensor InSensorInf annotation(
        Placement(visible = true, transformation(origin = {-72, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.inWWSensor InSensorDN2 annotation(
        Placement(visible = true, transformation(origin = {-24, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-28, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.inWWSensor InSensorN2 annotation(
        Placement(visible = true, transformation(origin = {24, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {24, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.inWWSensor InSensorN3 annotation(
        Placement(visible = true, transformation(origin = {68, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {68, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      e_NO_i = 1.0;
      e_O2_i = 1.0;
    equation
      Inf = InSensorInf.Value;
      DN2 = InSensorDN2.Value;
      N2 = InSensorN2.Value;
      N3 = InSensorN3.Value;
// controlled Return Sludge RS = 1
      C_RS.Q = InSensorInf.Value.Q;
//constant pumps
//C_RS.Q = Q_RS;
      C_WS.Q = Q_WS;
// control of recirculation
      e_NO = SP_NO - DN2.NO;
      e_NO = der(e_NO_i);
      P_NO = e_NO * k_P_NO;
      I_NO = e_NO_i * k_I_NO;
      D_NO = der(e_NO) * k_D_NO;
      C_REC.Q = max(10, Q_REC + P_NO + I_NO + D_NO);
      h_NO = C_REC.Q / Q_REC;
// Aeration
//constant blowers
      CairR3.Q = Q_air_N1;
      CairR4.Q = Q_air_N2;
// conrolles aeration of N3
// O2 setpoint with respect to NH concentration
      e_NH = SP_NH - N3.NH;
      f = 2 ^ (-e_NH);
      SP_O2 = if f < SP_O2_min then SP_O2_min else if f > SP_O2_max then SP_O2_max else f;
      e_O2 = SP_O2 - N3.O2;
      der(e_O2_i) = e_O2;
      P_O2 = e_O2 * k_P_O2;
      I_O2 = e_O2_i * k_I_O2;
      D_O2 = der(e_O2) * k_D_O2;
      CairR5.Q = max(100, Q_air_N3 + P_O2 + I_O2 + D_O2);
      h_O2 = CairR5.Q / 6000;
      annotation(
        Icon(graphics = {Text(origin = {0, 78}, lineThickness = 3, extent = {{-56, 12}, {56, -12}}, textString = "Controler"), Text(origin = {-36, 1}, lineThickness = 1.5, extent = {{-32, 15}, {32, -15}}, textString = "Waster", fontSize = 9), Text(origin = {44, 0}, lineThickness = 1.5, extent = {{-44, 6}, {44, -6}}, textString = "Air", fontSize = 9), Rectangle(origin = {1, 1}, lineThickness = 1.5, extent = {{-99, 99}, {99, -99}}), Line(origin = {1, 59}, points = {{-99, 1}, {99, -1}}, thickness = 1.5, smooth = Smooth.Bezier), Line(origin = {1, -66.4545}, points = {{-99, 1}, {99, -1}}, thickness = 1.5, smooth = Smooth.Bezier), Line(origin = {0, -4}, points = {{0, 62}, {0, -62}}, thickness = 1.5, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
        defaultComponentName = "Controller",
        Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Controller for models including ADM</H1>
  <p style=\"font-size:20px\">
  This controller includes the control for the dewatering unit and is therfore defined for the models including ADM.

  The controller contains input connectors from all five activated sludge tanks, the dewatering and the thickener unit.

  The controller gives output values for all pumps and blowers in the system.
  The fix values were taken from the BSM1 and BSM2.
    'John B. Copp. The cost simulation benchmark - description and simulator manual. COST.'
    'U. Jeppsson, M.-N. Pons, I. Nopens, J. Alex, J.B. Copp, K.V. Gernaey, C. Rosen, J.-P.
  Steyer, and P.A. Vanrolleghem. Benchmark simulation model no 2: general protocol and
  exploratory case studies. Water Science & Technology, 56(8):67, oct 2007.'


  The desired sludge concentration [gTSS/m3] can be adjusted.

  The internal recycle flow can be chosen as fix flow (reference) or as controlled flow. Watch the different initial conditions for the WWTP_ADM and PC_WWTP_ADM.
  Chosen can be by commenting out the not desired flows (e.g. for controlled WWTP_ADM, comment out the Q_int for reference and controlled PC_WWTP_ADM)

  The aeration of the fifth reactor can be chosen as fix air flow (reference) or as controlled flow.

  The calculation as PI-controller is based on
   'Rainer Froriep Heinz Mann, Horst Schielgen. Einführung in die Regelungstechnik. Hanser
  Verlag, 2009.'

  The k-values k_I and k_P as well as the conversion factor f were chosen by testing
  </p>
  </html>"));
    end controller;

    model controller_WWTP
      output Real Q_int, Q_r, Q_w "water/sludge flow, m3/d";
      output Real Q_air_R1, Q_air_R2, Q_air_R3, Q_air_R4, Q_air_R5 "air flow, m3/d";
      ASM.InPipe In1, In2, In3, In4, In5;
      //controlling aeration
      parameter Real SP_NH = 3 "setpoint S.NH, gN/m3";
      Real SP_O2 "setpoint O2, gO2/m3";
      Real SP_O2_min = 1 "lower limit SP_O2, gO2/m3";
      Real SP_O2_max = 4 "upper limit SP_O2, gO2/m3";
      Real P_O2 "proportional controller part";
      Real I_O2 "integral controller part";
      parameter Real k_P_O2 = 0 "k-value proportional part";
      parameter Real k_I_O2 = 5000 "k-value integral part";
      Real e_NH, e_O2 "control derivation, g/m3";
      Real e_O2_i "integral of control derivation";
      Real f "help variable";
      //controlling internal recycle flow
      parameter Real SP_NO = 3 "setpoint S.NO, gN/m3";
      Real P_NO "proportional controller part";
      Real I_NO "integral controller part";
      parameter Real k_P_NO = 1000 "k-value proportional part";
      parameter Real k_I_NO = 1000 "k-value integral part";
      Real e_NO "control deviation, gN/m3";
      Real e_NO_i "Integral of control deviation";
      Real speed_int, speed_bR5 "changing speed";
    initial equation
      e_NO_i = 0;
      e_O2_i = 0;
    equation
//constant pumps
      Q_w = 385;
      Q_r = 18446;
//controlled pumps
////internal recycle flow
      e_NO = SP_NO - In2.S.NO;
      P_NO = e_NO * k_P_NO;
      e_NO = der(e_NO_i);
      I_NO = e_NO_i * k_I_NO;
      Q_int = 32000 + P_NO + I_NO;
//controlled
//Q_int = 55338; //reference
//constant blowers
      Q_air_R1 = 0.01;
      Q_air_R2 = 0.01;
      Q_air_R3 = 15361;
      Q_air_R4 = 15361;
//controlled blower
////blower reactor 5
      e_NH = SP_NH - In5.S.NH;
      f = 2 ^ (-e_NH);
      SP_O2 = if f < SP_O2_min then SP_O2_min else if f > SP_O2_max then SP_O2_max else f;
      e_O2 = SP_O2 - In5.S.O2;
      e_O2 = der(e_O2_i);
      P_O2 = e_O2 * k_P_O2;
      I_O2 = e_O2_i * k_I_O2;
      Q_air_R5 = 3500 + P_O2 + I_O2;
//controlled
//Q_air_R5 = 4781; //reference
//changing speed
      der(Q_air_R5) = speed_bR5;
      der(Q_int) = speed_int;
      annotation(
        defaultComponentName = "Controller",
        Documentation(info = "<html>
  <H1 style=\"font-size:20px\">Controller for models without ADM</H1>
  <p style=\"font-size:20px\">
  This controller does !not! include the control for the dewatering unit and is therfore not compatible with the ADM.

  The controller contains input connectors from all five activated sludge tanks.

  The controller gives output values for all pumps and blowers in the WWTP system.
  The fix values were taken from the BSM1.
    'John B. Copp. The cost simulation benchmark - description and simulator manual. COST.'

  The internal recycle flow can be chosen as fix flow (reference) or as controlled flow.
  Chosen can be by commenting out the not desired flow (e.g. for reference, comment out the Q_int for controlled)

  The aeration of the fifth reactor can be chosen as fix air flow (reference) or as controlled flow.

  The calculation as PI-controller is based on
   'Rainer Froriep Heinz Mann, Horst Schielgen. Einführung in die Regelungstechnik. Hanser
  Verlag, 2009.'

  The k-values k_I and k_P as well as the conversion factor f were chosen by testing
  </p>
  </html>"));
    end controller_WWTP;
  end TechUnits;

  package SedTank
    model Tank
      parameter Real V = 6000;
      parameter Soluble Sini(I = 30, S = 1.15, O2 = 4.0, NH = 0.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 0.5, Fe2 = 0.5, Fe3 = 0.5, Al3 = 0.5);
      parameter Particulate Xini(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 2, FeOH = 2, AlP = 2, AlOH = 2);
      OpenWasteWater.ASM1P.InPipe In1;
      OpenWasteWater.ASM1P.OutPipe Out1;
      OpenWasteWater.ASM1P.Soluble S;
      OpenWasteWater.ASM1P.Particulate X;
    initial equation
      S = Sini;
      X = Xini;
    equation
      der(S.I) = In1.Q / V * (In1.S.I - S.I);
      der(S.S) = In1.Q / V * (In1.S.S - S.S);
      der(X.I) = In1.Q / V * (In1.X.I - X.I);
      der(X.S) = In1.Q / V * (In1.X.S - X.S);
      der(X.H) = In1.Q / V * (In1.X.H - X.H);
      der(X.A) = In1.Q / V * (In1.X.A - X.A);
      der(X.P) = In1.Q / V * (In1.X.P - X.P);
      der(S.O2) = In1.Q / V * (In1.S.O2 - S.O2);
      der(S.NO) = In1.Q / V * (In1.S.NO - S.NO);
      der(S.NH) = In1.Q / V * (In1.S.NH - S.NH);
      der(S.ND) = In1.Q / V * (In1.S.ND - S.ND);
      der(X.ND) = In1.Q / V * (In1.X.ND - X.ND);
      der(S.ALK) = In1.Q / V * (In1.S.ALK - S.ALK);
      der(S.P) = In1.Q / V * (In1.S.P - S.P);
      der(S.Fe2) = In1.Q / V * (In1.S.Fe2 - S.Fe2);
      der(S.Fe3) = In1.Q / V * (In1.S.Fe3 - S.Fe3);
      der(S.Al3) = In1.Q / V * (In1.S.Al3 - S.Al3);
      der(X.FeP) = In1.Q / V * (In1.X.FeP - X.FeP);
      der(X.AlP) = In1.Q / V * (In1.X.AlP - X.AlP);
      der(X.FeOH) = In1.Q / V * (In1.X.FeOH - X.FeOH);
      der(X.AlOH) = In1.Q / V * (In1.X.AlOH - X.AlOH);
      0 = Out1.Q + In1.Q;
      Out1.T = In1.T;
      Out1.S = S;
      Out1.X = X;
      annotation(
        defaultComponentName = "Tank",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Tank (no reactions) </H1>
      <p style=\"font-size:20px\">
      This Model is a model for the retention of the components of ASM no. 1
      for a given volume V. It can be used as equilibrium tank or as a module
      for a simplified secoundary clarifier model.
      </p>
    </html>"));
    end Tank;

    model StorageTank
      extends OpenWasteWater.Icons.storage;
      parameter Real Vini = 100 "m3";
      parameter Soluble Sini(I = 30, S = 1.15, O2 = 4.0, NH = 0.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 0.5, Fe2 = 0.5, Fe3 = 0.5, Al3 = 0.5);
      parameter Particulate Xini(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 2, FeOH = 2, AlP = 2, AlOH = 2);
      output Real V;
      OpenWasteWater.ASM1P.Soluble S;
      OpenWasteWater.ASM1P.Particulate X;
      OpenWasteWater.ASM1P.InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-104, 0}, extent = {{-14, -14}, {14, 14}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {101, -1}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.outV OutV_StorateTank annotation(
        Placement(visible = true, transformation(origin = {79, 93}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {76, 96}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      V = Vini;
      S = Sini;
      X = Xini;
    equation
      der(V) = In1.Q + Out1.Q;
      der(S.I) = In1.Q / V * In1.S.I + (Out1.Q - der(V)) / V * S.I;
      der(S.S) = In1.Q / V * In1.S.S + (Out1.Q - der(V)) / V * S.S;
      der(X.I) = In1.Q / V * In1.X.I + (Out1.Q - der(V)) / V * X.I;
      der(X.S) = In1.Q / V * In1.X.S + (Out1.Q - der(V)) / V * X.S;
      der(X.H) = In1.Q / V * In1.X.H + (Out1.Q - der(V)) / V * X.H;
      der(X.A) = In1.Q / V * In1.X.A + (Out1.Q - der(V)) / V * X.A;
      der(X.P) = In1.Q / V * In1.X.P + (Out1.Q - der(V)) / V * X.P;
      der(S.O2) = In1.Q / V * In1.S.O2 + (Out1.Q - der(V)) / V * S.O2;
      der(S.NO) = In1.Q / V * In1.S.NO + (Out1.Q - der(V)) / V * S.NO;
      der(S.NH) = In1.Q / V * In1.S.NH + (Out1.Q - der(V)) / V * S.NH;
      der(S.ND) = In1.Q / V * In1.S.ND + (Out1.Q - der(V)) / V * S.ND;
      der(X.ND) = In1.Q / V * In1.X.ND + (Out1.Q - der(V)) / V * X.ND;
      der(S.ALK) = In1.Q / V * In1.S.ALK + (Out1.Q - der(V)) / V * S.ALK;
      der(S.P) = In1.Q / V * In1.S.P + (Out1.Q - der(V)) / V * S.P;
      der(S.Fe2) = In1.Q / V * In1.S.Fe2 + (Out1.Q - der(V)) / V * S.Fe2;
      der(S.Fe3) = In1.Q / V * In1.S.Fe3 + (Out1.Q - der(V)) / V * S.Fe3;
      der(S.Al3) = In1.Q / V * In1.S.Al3 + (Out1.Q - der(V)) / V * S.Al3;
      der(X.FeP) = In1.Q / V * In1.X.FeP + (Out1.Q - der(V)) / V * X.FeP;
      der(X.AlP) = In1.Q / V * In1.X.AlP + (Out1.Q - der(V)) / V * X.AlP;
      der(X.FeOH) = In1.Q / V * In1.X.FeOH + (Out1.Q - der(V)) / V * X.FeOH;
      der(X.AlOH) = In1.Q / V * In1.X.AlOH + (Out1.Q - der(V)) / V * X.AlOH;
      Out1.T = In1.T;
      Out1.S = S;
      Out1.X = X;
      OutV_StorateTank.V = V;
      annotation(
        defaultComponentName = "STorageTank",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Tank (no reactions) </H1>
      <p style=\"font-size:20px\">
      This Model is a model for the retention of the components of ASM no. 1
      for variable volume V. It can be used as equilibrium tank.
      </p>
    </html>"));
    end StorageTank;

    model Separator
      parameter Real fns = 0.00228 "fraction of not sedimentable soldis";
      OpenWasteWater.ASM1P.InPipe In1;
      OpenWasteWater.ASM1P.OutPipe Out1, Out2;
    equation
      0 = Out1.Q + Out2.Q + In1.Q;
      Out1.T = In1.T;
      Out2.T = In1.T;
      Out1.S = In1.S;
      Out1.X.I = fns * In1.X.I;
      Out1.X.S = fns * In1.X.S;
      Out1.X.H = fns * In1.X.H;
      Out1.X.A = fns * In1.X.A;
      Out1.X.P = fns * In1.X.P;
      Out1.X.ND = fns * In1.X.ND;
      Out1.X.FeP = fns * In1.FeP;
      Out1.X.AlP = fns * In1.AlP;
      Out1.X.FeOH = fns * In1.FeOH;
      Out1.X.AlOH = fns * In1.AlOH;
      Out2.S = In1.S;
      Out2.Q * Out2.X.I + In1.Q * In1.X.I + Out1.Q * Out1.X.I = 0;
      Out2.Q * Out2.X.S + In1.Q * In1.X.S + Out1.Q * Out1.X.S = 0;
      Out2.Q * Out2.X.H + In1.Q * In1.X.H + Out1.Q * Out1.X.H = 0;
      Out2.Q * Out2.X.A + In1.Q * In1.X.A + Out1.Q * Out1.X.A = 0;
      Out2.Q * Out2.X.P + In1.Q * In1.X.P + Out1.Q * Out1.X.P = 0;
      Out2.Q * Out2.X.ND + In1.Q * In1.X.ND + Out1.Q * Out1.X.ND = 0;
      Out2.Q * Out2.X.FeP + In1.Q * In1.X.FeP + Out1.Q * Out1.X.FeP = 0;
      Out2.Q * Out2.X.AlP + In1.Q * In1.X.AlP + Out1.Q * Out1.X.AlP = 0;
      Out2.Q * Out2.X.FeOH + In1.Q * In1.X.FeOH + Out1.Q * Out1.X.FeOH = 0;
      Out2.Q * Out2.X.AlOH + In1.Q * In1.X.AlOH + Out1.Q * Out1.X.AlOH = 0;
      annotation(
        defaultComponentName = "Separator",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Separator </H1>
      <p style=\"font-size:20px\">
      This Model is a model of a separator of the components of ASM no. 1.
      The dissolved components \"S\" are splited; the particulat components
      \"X\" are treated in a way that the fraction of non sedimentable solids
      \"fns\" is treated like dissoled and the rest is separated and concentrated
      in the sludge connetor. The clarified fraction is discharged via the
      clear connector.
      It can be used as a module for a simplified secoundary clarifier model.
      </p>
      <p style=\"font-size:20px\">
      <strong>Example:</strong> <br />
      separator Cl(Q2=100, fns=0.002) <br />
      connect(inlet, <strong>Cl.in1</strong>); <br />
      connect(<strong>Cl.clear</strong>, sink1.inlet); <br />
      connect(<strong>Cl.sludge</strong>, sink2.inlet); <br />
      </p>
    </html>"));
    end Separator;

    model Centrifuge
      extends OpenWasteWater.Icons.centrifuge;
      parameter Real fns = 0.02 "fraction of not sedimentable soldis";
      OpenWasteWater.ASM1P.InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-103, -19}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {-100, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {-102, -64}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {-100, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out2 annotation(
        Placement(visible = true, transformation(origin = {77, -93}, extent = {{-11, -11}, {11, 11}}, rotation = 0), iconTransformation(origin = {78, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.outWWSensor Out2WWSensor annotation(
        Placement(visible = true, transformation(origin = {96, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {96, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      0 = Out1.Q + Out2.Q + In1.Q;
      Out1.T = In1.T;
      Out2.T = In1.T;
      Out1.S = In1.S;
      Out1.X.I = fns * In1.X.I;
      Out1.X.S = fns * In1.X.S;
      Out1.X.H = fns * In1.X.H;
      Out1.X.A = fns * In1.X.A;
      Out1.X.P = fns * In1.X.P;
      Out1.X.ND = fns * In1.X.ND;
      Out1.X.FeP = fns * In1.X.FeP;
      Out1.X.AlP = fns * In1.X.AlP;
      Out1.X.FeOH = fns * In1.X.FeOH;
      Out1.X.AlOH = fns * In1.X.AlOH;
      Out2.S = In1.S;
      Out2.Q * Out2.X.I + In1.Q * In1.X.I + Out1.Q * Out1.X.I = 0;
      Out2.Q * Out2.X.S + In1.Q * In1.X.S + Out1.Q * Out1.X.S = 0;
      Out2.Q * Out2.X.H + In1.Q * In1.X.H + Out1.Q * Out1.X.H = 0;
      Out2.Q * Out2.X.A + In1.Q * In1.X.A + Out1.Q * Out1.X.A = 0;
      Out2.Q * Out2.X.P + In1.Q * In1.X.P + Out1.Q * Out1.X.P = 0;
      Out2.Q * Out2.X.ND + In1.Q * In1.X.ND + Out1.Q * Out1.X.ND = 0;
      Out2.Q * Out2.X.FeP + In1.Q * In1.X.FeP + Out1.Q * Out1.X.FeP = 0;
      Out2.Q * Out2.X.AlP + In1.Q * In1.X.AlP + Out1.Q * Out1.X.AlP = 0;
      Out2.Q * Out2.X.FeOH + In1.Q * In1.X.FeOH + Out1.Q * Out1.X.FeOH = 0;
      Out2.Q * Out2.X.AlOH + In1.Q * In1.X.AlOH + Out1.Q * Out1.X.AlOH = 0;
      Out2WWSensor.Value.TSS = (Out2.X.I + Out2.X.S + Out2.X.H + Out2.X.A + Out2.X.P) * 0.75 + Out2.X.FeP + Out2.X.AlP + Out2.X.FeOH + Out2.X.AlOH;
      Out2WWSensor.Value.O2 = Out2.S.O2;
      Out2WWSensor.Value.S = Out2.S.S;
      Out2WWSensor.Value.P = Out2.S.P;
      Out2WWSensor.Value.NH = Out2.S.NH;
      Out2WWSensor.Value.NO = Out2.S.NO;
      Out2WWSensor.Value.Q = In1.Q;
      annotation(
        defaultComponentName = "Separator",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Separator </H1>
      <p style=\"font-size:20px\">
      This Model is a model of a separator of the components of ASM no. 1.
      The dissolved components \"S\" are splited; the particulat components
      \"X\" are treated in a way that the fraction of non sedimentable solids
      \"fns\" is treated like dissoled and the rest is separated and concentrated
      in the sludge connetor. The clarified fraction is discharged via the
      clear connector.
      It can be used as a module for a simplified secoundary clarifier model.
      </p>
      <p style=\"font-size:20px\">
      <strong>Example:</strong> <br />
      separator Cl(Q2=100, fns=0.002) <br />
      connect(inlet, <strong>Cl.in1</strong>); <br />
      connect(<strong>Cl.clear</strong>, sink1.inlet); <br />
      connect(<strong>Cl.sludge</strong>, sink2.inlet); <br />
      </p>
    </html>"));
    end Centrifuge;

    model PreThick
      parameter Real V = 600;
      extends WWParameters;
      OpenWasteWater.ASM1P.InPipe In1;
      OpenWasteWater.ASM1P.OutPipe Out1, Out2;
      OpenWasteWater.ASM1P.Soluble S;
      OpenWasteWater.ASM1P.Particulate X;
      OpenWasteWater.ASM1P.TechUnits.outWWSensor Out2WWSensor;
      Real HRT_h, n_COD, n_X, H;
    equation
      S = In1.S;
      X = In1.X;
      0 = Out1.Q + Out2.Q + In1.Q;
      HRT_h = V / In1.Q * 24.0;
      n_COD = 2.7 * (2 * log(HRT_h) + 9) / 100;
      H = n_COD * COD / TSS * 0.75;
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;
      Out1.T = In1.T;
      Out2.T = In1.T;
      Out1.S = In1.S;
      Out1.X.I = (1 - n_X) * In1.X.I;
      Out1.X.S = (1 - n_X) * In1.X.S;
      Out1.X.H = (1 - n_X) * In1.X.H;
      Out1.X.A = (1 - n_X) * In1.X.A;
      Out1.X.P = (1 - n_X) * In1.X.P;
      Out1.X.ND = (1 - n_X) * In1.X.ND;
      Out1.X.FeP = (1 - n_X) * In1.X.FeP;
      Out1.X.AlP = (1 - n_X) * In1.X.AlP;
      Out1.X.FeOH = (1 - n_X) * In1.X.FeOH;
      Out1.X.AlOH = (1 - n_X) * In1.X.AlOH;
      Out2.S = In1.S;
      Out2.Q * Out2.X.I + In1.Q * In1.X.I + Out1.Q * Out1.X.I = 0;
      Out2.Q * Out2.X.S + In1.Q * In1.X.S + Out1.Q * Out1.X.S = 0;
      Out2.Q * Out2.X.H + In1.Q * In1.X.H + Out1.Q * Out1.X.H = 0;
      Out2.Q * Out2.X.A + In1.Q * In1.X.A + Out1.Q * Out1.X.A = 0;
      Out2.Q * Out2.X.P + In1.Q * In1.X.P + Out1.Q * Out1.X.P = 0;
      Out2.Q * Out2.X.ND + In1.Q * In1.X.ND + Out1.Q * Out1.X.ND = 0;
      Out2.Q * Out2.X.FeP + In1.Q * In1.X.FeP + Out1.Q * Out1.X.FeP = 0;
      Out2.Q * Out2.X.AlP + In1.Q * In1.X.AlP + Out1.Q * Out1.X.AlP = 0;
      Out2.Q * Out2.X.FeOH + In1.Q * In1.X.FeOH + Out1.Q * Out1.X.FeOH = 0;
      Out2.Q * Out2.X.AlOH + In1.Q * In1.X.AlOH + Out1.Q * Out1.X.AlOH = 0;
      Out2WWSensor.Value.O2 = Out2.S.O2;
      Out2WWSensor.Value.S = Out2.S.S;
      Out2WWSensor.Value.NH = Out2.S.NH;
      Out2WWSensor.Value.NO = Out2.S.NO;
      Out2WWSensor.Value.P = Out2.S.P;
      Out2WWSensor.Value.Q = abs(Out2.Q);
      Out2WWSensor.Value.TSS = 0.75 * (Out2.X.I + Out2.X.S + Out2.X.H + Out2.X.A + Out2.X.P);
      annotation(
        defaultComponentName = "PreThick",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">Separator </H1>
      <p style=\"font-size:20px\">
      This Model is a model of a separator of the components of ASM no. 1.
      The dissolved components \"S\" are splited; the particulat components
      \"X\" are treated in a way that the fraction of non sedimentable solids
      \"fns\" is treated like dissoled and the rest is separated and concentrated
      in the sludge connetor. The clarified fraction is discharged via the
      clear connector.
      It can be used as a module for a simplified secoundary clarifier model.
      </p>
      <p style=\"font-size:20px\">
      <strong>Example:</strong> <br />
      separator Cl(Q2=100, fns=0.002) <br />
      connect(inlet, <strong>Cl.in1</strong>); <br />
      connect(<strong>Cl.clear</strong>, sink1.inlet); <br />
      connect(<strong>Cl.sludge</strong>, sink2.inlet); <br />
      </p>
    </html>"));
    end PreThick;

    model SCT
      extends OpenWasteWater.Icons.SecClarSimple;
      parameter Real V = 6000 "Volume of the SCT";
      parameter Real fns = 0.00228 "fraction of not sedimentable soldis";
      parameter Integer N = 6;
      Tank[N] T(each V = V/N);
      Separator S(fns = fns);
      OpenWasteWater.ASM1P.InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-100, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-97, 14}, extent = {{-13, -14}, {13, 14}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {98, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 60}, extent = {{-12, -14}, {12, 14}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out2 annotation(
        Placement(visible = true, transformation(origin = {-2, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-1, -98}, extent = {{-13, -12}, {13, 12}}, rotation = 0)));
    equation
      connect(In1, T[1].In1);
    
      for i in 1:N - 1 loop
        connect(T[i].Out1, T[i + 1].In1);
      end for;
    
      connect(T[N].Out1, S.In1);
      connect(S.Out1, Out1);
      connect(S.Out2, Out2);
      annotation(
        defaultComponentName = "SCT",
        Documentation(info = "<html>
            <H1 style=\"font-size:20px\">SCT (Secondary Clarifier Tank) </H1>
            <p style=\"font-size:20px\">
            This Model is a model of a secondary clarifier for the components of ASM no. 1.
            The dissolved components \"S\" are splited; the particulat components
            \"X\" are treated in a way that the fraction of non sedimentable solids
            \"fns\" is treated like dissoled and the rest is separated and concentrated
            in the sludge connector. The clarified fraction is discharged via the
            clear connector.
            </p>
            <p style=\"font-size:20px\">
            <strong>Example:</strong> <br />
            SCT Cl(V=6000, Q2=1000, fns=0.002) <br />
            connect(inlet, <strong>Cl.in1</strong>); <br />
            connect(<strong>Cl.clear</strong>, sink.inlet); <br />
            connect(<strong>Cl.sludge</strong>, sink.inlet); <br />
            </p>
        </html>"));
    end SCT;

    model PreClar
      extends OpenWasteWater.Icons.preclar1;
      parameter Real V = 600 "Volume of the SCT";
      constant Integer N = 6;
      Tank T1(V = V / N), T2(V = V / N), T3(V = V / N), T4(V = V / N), T5(V = V / N), T6(V = V / N);
      PreThick PST(V = V);
      OpenWasteWater.ASM1P.InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-100, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-99, 4}, extent = {{-13, -14}, {13, 14}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {96, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 2}, extent = {{-12, -14}, {12, 14}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out2 annotation(
        Placement(visible = true, transformation(origin = {-2, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-1, -98}, extent = {{-13, -12}, {13, 12}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.outWWSensor OutTSSWWSensor annotation(
        Placement(visible = true, transformation(origin = {62, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(In1, T1.In1) annotation(
        Line);
      connect(T1.Out1, T2.In1);
      connect(T2.Out1, T3.In1);
      connect(T3.Out1, T4.In1);
      connect(T4.Out1, T5.In1);
      connect(T5.Out1, T6.In1);
      connect(T6.Out1, PST.In1);
      connect(PST.Out1, Out1);
      connect(PST.Out2, Out2);
      connect(PST.Out2WWSensor, OutTSSWWSensor);
      annotation(
        defaultComponentName = "SCT",
        Documentation(info = "<html>
            <H1 style=\"font-size:20px\">SCT (Secondary Clarifier Tank) </H1>
            <p style=\"font-size:20px\">
            This Model is a model of a secondary clarifier for the components of ASM no. 1.
            The dissolved components \"S\" are splited; the particulat components
            \"X\" are treated in a way that the fraction of non sedimentable solids
            \"fns\" is treated like dissoled and the rest is separated and concentrated
            in the sludge connector. The clarified fraction is discharged via the
            clear connector.
            </p>
            <p style=\"font-size:20px\">
            <strong>Example:</strong> <br />
            SCT Cl(V=6000, Q2=1000, fns=0.002) <br />
            connect(inlet, <strong>Cl.in1</strong>); <br />
            connect(<strong>Cl.clear</strong>, sink.inlet); <br />
            connect(<strong>Cl.sludge</strong>, sink.inlet); <br />
            </p>
        </html>"));
    end PreClar;
  end SedTank;

  package SecClar
    connector InClarConnector
      Real T;
      flow Real Q;
      input Soluble S;
      input Real Xn;
      input Real Xs;
    end InClarConnector;

    connector OutClarConnector
      Real T;
      flow Real Q;
      output Soluble S;
      output Real Xn;
      output Real Xs;
    end OutClarConnector;

    connector inXconnect
      input Real Xs;
    end inXconnect;

    connector outXconnect
      output Real Xs;
    end outXconnect;

    connector inFrac
      input Real relation[10];
    end inFrac;

    connector outFrac
      output Real relation[10];
    end outFrac;

    partial model ClarZone
      parameter Real A = 1500 "m2";
      parameter Real z = 0.5 "m";
      parameter Soluble Sini(I = 30, S = 1, O2 = 2, NH = 0.72, NO = 5, ND = 1, ALK = 3.61, P = 0.5, Fe2 = 0.0, Fe3 = 0.0, Al3 = 0.0);
      Real Xs;
      Real Xn;
      Soluble S;
      Real V;
      Real TSS;
    equation
      V = A * z;
      TSS = 0.75 * (Xs + Xn);
    end ClarZone;

    model InputPart
      parameter Real fns = 0.02 "-";
      Real Xtot, Xsum;
      Real frac[10];
      InPipe In1;
      OutClarConnector Out1, Out2;
      inXconnect iXup, iXdown;
      outFrac oB, oT;
      extends ClarZone;
    initial equation
      S = Sini;
      Xs = 600;
      Xn = 6.0;
    equation
      0 = Out1.Q + Out2.Q + In1.Q;
      Out1.T = In1.T;
      Out2.T = In1.T;
      Xsum = In1.X.I + In1.X.S + In1.X.H + In1.X.A + In1.X.P + In1.X.ND + In1.X.FeP + In1.X.AlP + In1.X.FeOH + In1.X.AlOH;
      Xtot = In1.X.I + In1.X.S + In1.X.H + In1.X.A + In1.X.P + In1.X.ND + (In1.X.FeP + In1.X.AlP + In1.X.FeOH + In1.X.AlOH) / 0.75;
      frac[1] = In1.X.I / Xsum;
      frac[2] = In1.X.S / Xsum;
      frac[3] = In1.X.H / Xsum;
      frac[4] = In1.X.A / Xsum;
      frac[5] = In1.X.P / Xsum;
      frac[6] = In1.X.ND / Xsum;
      frac[7] = In1.X.FeP / Xsum;
      frac[8] = In1.X.AlP / Xsum;
      frac[9] = In1.X.FeOH / Xsum;
      frac[10] = In1.X.AlOH / Xsum;
      der(Xs) = In1.Q / V * ((1 - fns) * Xtot - Xs) + (min(fJS(Xs), fJS(iXup.Xs)) - min(fJS(Xs), fJS(iXdown.Xs))) / z;
      der(Xn) = In1.Q / V * (fns * Xtot - Xn);
      der(S.I) = In1.Q / V * (In1.S.I - S.I);
      der(S.S) = In1.Q / V * (In1.S.S - S.S);
      der(S.O2) = In1.Q / V * (In1.S.O2 - S.O2);
      der(S.NO) = In1.Q / V * (In1.S.NO - S.NO);
      der(S.NH) = In1.Q / V * (In1.S.NH - S.NH);
      der(S.ND) = In1.Q / V * (In1.S.ND - S.ND);
      der(S.ALK) = In1.Q / V * (In1.S.ALK - S.ALK);
      der(S.P) = In1.Q / V * (In1.S.P - S.P);
      der(S.Fe2) = In1.Q / V * (In1.S.Fe2 - S.Fe2);
      der(S.Fe3) = In1.Q / V * (In1.S.Fe3 - S.Fe3);
      der(S.Al3) = In1.Q / V * (In1.S.Al3 - S.Al3);
      Out1.S = S;
      Out1.Xn = Xn;
      Out1.Xs = Xs;
      oB.relation = frac;
      Out2.S = S;
      Out2.Xn = Xn;
      Out2.Xs = Xs;
      oT.relation = frac;
    end InputPart;

    model TopPart
      parameter Real A = 1500 "m2";
      parameter Real z = 0.5 "m";
      Real Xtot;
      Real frac[10];
      OutPipe Out1;
      InClarConnector In1;
      outXconnect oXdown;
      inFrac iT;
      extends ClarZone;
    initial equation
      S = Sini;
      Xs = 0;
      Xn = 4.6;
    equation
      0 = Out1.Q + In1.Q;
      Out1.T = In1.T;
      frac = iT.relation;
      der(Xs) = In1.Q / V * (In1.Xs - Xs) - min(fJS(In1.Xs), fJS(Xs)) / z;
      der(Xn) = In1.Q / V * (In1.Xn - Xn);
      der(S.I) = In1.Q / V * (In1.S.I - S.I);
      der(S.S) = In1.Q / V * (In1.S.S - S.S);
      der(S.O2) = In1.Q / V * (In1.S.O2 - S.O2);
      der(S.NO) = In1.Q / V * (In1.S.NO - S.NO);
      der(S.NH) = In1.Q / V * (In1.S.NH - S.NH);
      der(S.ND) = In1.Q / V * (In1.S.ND - S.ND);
      der(S.ALK) = In1.Q / V * (In1.S.ALK - S.ALK);
      der(S.P) = In1.Q / V * (In1.S.P - S.P);
      der(S.Fe2) = In1.Q / V * (In1.S.Fe2 - S.Fe2);
      der(S.Fe3) = In1.Q / V * (In1.S.Fe3 - S.Fe3);
      der(S.Al3) = In1.Q / V * (In1.S.Al3 - S.Al3);
      Xtot = Xs + Xn;
      Out1.S = S;
      Out1.X.I = Xtot * frac[1];
      Out1.X.S = Xtot * frac[2];
      Out1.X.H = Xtot * frac[3];
      Out1.X.A = Xtot * frac[4];
      Out1.X.P = Xtot * frac[5];
      Out1.X.ND = Xtot * frac[6];
      Out1.X.FeP = Xtot * frac[7] * 0.75;
      Out1.X.AlP = Xtot * frac[8] * 0.75;
      Out1.X.FeOH = Xtot * frac[9] * 0.75;
      Out1.X.AlOH = Xtot * frac[10] * 0.75;
      oXdown.Xs = Xs;
    end TopPart;

    model UpPart
      OutClarConnector Out1;
      InClarConnector In1;
      outXconnect oXdown;
      inXconnect iXup;
      extends ClarZone;
    initial equation
      S = Sini;
      Xs = 50;
      Xn = 4.6;
    equation
      0 = In1.Q + Out1.Q;
      Out1.T = In1.T;
      der(Xs) = In1.Q / V * (In1.Xs - Xs) + (min(fJS(iXup.Xs), fJS(Xs)) - min(fJS(Xs), fJS(In1.Xs))) / z;
      der(Xn) = In1.Q / V * (In1.Xn - Xn);
      der(S.I) = In1.Q / V * (In1.S.I - S.I);
      der(S.S) = In1.Q / V * (In1.S.S - S.S);
      der(S.O2) = In1.Q / V * (In1.S.O2 - S.O2);
      der(S.NO) = In1.Q / V * (In1.S.NO - S.NO);
      der(S.NH) = In1.Q / V * (In1.S.NH - S.NH);
      der(S.ND) = In1.Q / V * (In1.S.ND - S.ND);
      der(S.ALK) = In1.Q / V * (In1.S.ALK - S.ALK);
      der(S.P) = In1.Q / V * (In1.S.P - S.P);
      der(S.Fe2) = In1.Q / V * (In1.S.Fe2 - S.Fe2);
      der(S.Fe3) = In1.Q / V * (In1.S.Fe3 - S.Fe3);
      der(S.Al3) = In1.Q / V * (In1.S.Al3 - S.Al3);
      Out1.S = S;
      Out1.Xn = Xn;
      Out1.Xs = Xs;
      oXdown.Xs = Xs;
    end UpPart;

    model DownPart
      InClarConnector In1;
      OutClarConnector Out1;
      outXconnect oXup;
      inXconnect iXdown;
      extends ClarZone;
    initial equation
      S = Sini;
      Xs = 5300;
      Xn = 4.6;
    equation
      0 = In1.Q + Out1.Q;
      In1.T = Out1.T;
      der(Xs) = In1.Q / V * (In1.Xs - Xs) + (min(fJS(In1.Xs), fJS(Xs)) - min(fJS(Xs), fJS(iXdown.Xs))) / z;
      der(Xn) = In1.Q / V * (In1.Xn - Xn);
      der(S.I) = In1.Q / V * (In1.S.I - S.I);
      der(S.S) = In1.Q / V * (In1.S.S - S.S);
      der(S.O2) = In1.Q / V * (In1.S.O2 - S.O2);
      der(S.NO) = In1.Q / V * (In1.S.NO - S.NO);
      der(S.NH) = In1.Q / V * (In1.S.NH - S.NH);
      der(S.ND) = In1.Q / V * (In1.S.ND - S.ND);
      der(S.ALK) = In1.Q / V * (In1.S.ALK - S.ALK);
      der(S.P) = In1.Q / V * (In1.S.P - S.P);
      der(S.Fe2) = In1.Q / V * (In1.S.Fe2 - S.Fe2);
      der(S.Fe3) = In1.Q / V * (In1.S.Fe3 - S.Fe3);
      der(S.Al3) = In1.Q / V * (In1.S.Al3 - S.Al3);
      Out1.S = S;
      Out1.Xn = Xn;
      Out1.Xs = Xs;
      oXup.Xs = Xs;
    end DownPart;

    model BottomPart
      Real frac[10];
      Real Xtot;
      OutPipe Out1;
      InClarConnector In1;
      outXconnect oXup;
      inFrac iB;
      extends ClarZone;
    initial equation
      S = Sini;
      Xs = 8000;
      Xn = 10.6;
    equation
      0 = In1.Q + Out1.Q;
      In1.T = Out1.T;
      frac = iB.relation;
      der(Xs) = In1.Q / V * (In1.Xs - Xs) + min(fJS(In1.Xs), fJS(Xs)) / z;
      der(Xn) = In1.Q / V * (In1.Xn - Xn);
      der(S.I) = In1.Q / V * (In1.S.I - S.I);
      der(S.S) = In1.Q / V * (In1.S.S - S.S);
      der(S.O2) = In1.Q / V * (In1.S.O2 - S.O2);
      der(S.NO) = In1.Q / V * (In1.S.NO - S.NO);
      der(S.NH) = In1.Q / V * (In1.S.NH - S.NH);
      der(S.ND) = In1.Q / V * (In1.S.ND - S.ND);
      der(S.ALK) = In1.Q / V * (In1.S.ALK - S.ALK);
      der(S.P) = In1.Q / V * (In1.S.P - S.P);
      der(S.Fe2) = In1.Q / V * (In1.S.Fe2 - S.Fe2);
      der(S.Fe3) = In1.Q / V * (In1.S.Fe3 - S.Fe3);
      der(S.Al3) = In1.Q / V * (In1.S.Al3 - S.Al3);
      Xtot = Xs + Xn;
      Out1.S = S;
      Out1.X.I = Xtot * frac[1];
      Out1.X.S = Xtot * frac[2];
      Out1.X.H = Xtot * frac[3];
      Out1.X.A = Xtot * frac[4];
      Out1.X.P = Xtot * frac[5];
      Out1.X.ND = Xtot * frac[6];
      Out1.X.FeP = Xtot * frac[7] * 0.75;
      Out1.X.AlP = Xtot * frac[8] * 0.75;
      Out1.X.FeOH = Xtot * frac[9] * 0.75;
      Out1.X.AlOH = Xtot * frac[10] * 0.75;
      oXup.Xs = Xs;
    end BottomPart;

    model SCL
      // parameter Real Q2 = 18881 "Return sludge flow + waste sludge";
      parameter Real fns = 0.00228 "fraction of not sedimentable soldis";
      parameter Real A = 1500.0 "m2";
      parameter Real z = 0.4 "m";
      extends OpenWasteWater.Icons.SecClar;
      InPipe In1 annotation(
        Placement(transformation(extent = {{-110, 4}, {-90, 24}}, rotation = 0)));
      OutPipe Out1 annotation(
        Placement(transformation(extent = {{92, 47}, {112, 67}}, rotation = 0)));
      TopPart Top(A = A, z = z);
      UpPart Up1(A = A, z = 0.4), Up2(A = A, z = z), Up3(A = A, z = z);
      InputPart InC(A = A, z = z, fns = fns);
      DownPart Down1(A = A, z = z), Down2(A = A, z = z), Down3(A = A, z = z), Down4(A = A, z = z);
      BottomPart Bottom(A = A, z = z);
      OutPipe Out2 annotation(
        Placement(visible = true, transformation(origin = {0, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(In1, InC.In1);
      connect(Top.Out1, Out1);
      connect(Bottom.Out1, Out2);
//////////////////////////////////////////////////////
      connect(InC.oT, Top.iT);
      connect(InC.oB, Bottom.iB);
//////////////////////////////////////////////////////
      connect(InC.Out1, Up1.In1);
      connect(InC.Out2, Down1.In1);
      connect(Down1.oXup, InC.iXdown);
      connect(Up1.oXdown, InC.iXup);
//////////////////////////////////////////////////////
      connect(Up1.Out1, Up2.In1);
      connect(Up2.oXdown, Up1.iXup);
      connect(Up2.Out1, Up3.In1);
      connect(Up3.oXdown, Up2.iXup);
      connect(Up3.Out1, Top.In1);
      connect(Top.oXdown, Up3.iXup);
//////////////////////////////////////////////////////
      connect(Down1.Out1, Down2.In1);
      connect(Down2.oXup, Down1.iXdown);
      connect(Down2.Out1, Down3.In1);
      connect(Down3.oXup, Down2.iXdown);
      connect(Down3.Out1, Down4.In1);
      connect(Down4.oXup, Down3.iXdown);
      connect(Down4.Out1, Bottom.In1);
      connect(Bottom.oXup, Down4.iXdown);
      annotation(
        defaultComponentName = "SCL",
        Documentation(info = "<html>
      <H1 style=\"font-size:20px\">SCT (Secondary Clarifier Takacs) </H1>
      <p style=\"font-size:20px\">
      This Model is a model of a secondary clarifier for the components of ASM no. 1.
      The dissolved components \"S\" are splited; the particulat components
      \"X\" are treated in a way that the fraction of non sedimentable solids
      \"fns\" is treated like dissoled and the rest is separated and concentrated
      in the sludge connector. The models base on an 10 zone sedementation tank like
      described in the literature. The clarified fraction is discharged via the
      clear connector.
      </p>
      <p style=\"font-size:20px\">
      <strong>Example:</strong> <br />
      SCT Cl(V=6000, Q2=1000, fns=0.002, A=1000, z=0.4) <br />
      connect(inlet, <strong>Cl.in1</strong>); <br />
      connect(<strong>Cl.clear</strong>, sink.inlet); <br />
      connect(<strong>Cl.sludge</strong>, sink.inlet); <br />
      </p>
    </html>"),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end SCL;
  end SecClar;

  package WWTP_Examples
    model WWTP_COST
      OpenWasteWater.ASM1P.Inflow Inf(T = 20) annotation(
        Placement(visible = true, transformation(origin = {-90, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Mixer3 M1 annotation(
        Placement(visible = true, transformation(origin = {-64, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN1(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 30.2, NO = 8.9, ND = 0.9, ALK = 4.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 1500, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-38, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN2(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 30.2, NO = 1.9, ND = 0.9, ALK = 5.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 1500, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-14, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N1(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 10.2, NO = 8.9, ND = 0.9, ALK = 4.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 2000, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {14, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N2(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 5.2, NO = 8.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 2000, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {40, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N3(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 0.2, NO = 10.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 2000, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {64, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B1 annotation(
        Placement(visible = true, transformation(origin = {14, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B2 annotation(
        Placement(visible = true, transformation(origin = {40, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B3 annotation(
        Placement(visible = true, transformation(origin = {66, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.controller PLS(Q_REC = 20000, Q_RS = 18446, Q_WS = 350, Q_air_N3 = 10000, k_D_O2 = 0, k_I_NO = 5000, k_I_O2 = 5000, k_P_NO = 2000, k_P_O2 = 2500) annotation(
        Placement(visible = true, transformation(origin = {36, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump REC_Pump(Qmax = 92239) annotation(
        Placement(visible = true, transformation(origin = {-60, 12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump RS_Pump annotation(
        Placement(visible = true, transformation(origin = {-38, -84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump WS_Pump annotation(
        Placement(visible = true, transformation(origin = {42, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.SecClar.SCL SC annotation(
        Placement(visible = true, transformation(origin = {2, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Effluent Clear annotation(
        Placement(visible = true, transformation(origin = {82, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.WasteSludge Sludge annotation(
        Placement(visible = true, transformation(origin = {82, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S1 annotation(
        Placement(visible = true, transformation(origin = {-24, 8}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S2 annotation(
        Placement(visible = true, transformation(origin = {2, -68}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.TechUnits.Sensor SensorInf(t_D = 0.01) annotation(
        Placement(visible = true, transformation(origin = {-42, 92}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.Mixer2 M2 annotation(
        Placement(visible = true, transformation(origin = {-58, 44}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Dosing Dosing_Fe(InAl3 = 7500, InFe2 = 20000, InFe3 = 30000, betaP = 5) annotation(
        Placement(visible = true, transformation(origin = {-32, 36}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(SensorInf.S_Out1, PLS.InSensorInf) annotation(
        Line(points = {{-32, 92}, {98, 92}, {98, -38}, {30, -38}, {30, 4}, {29, 4}, {29, -19}}, color = {78, 154, 6}));
      connect(Inf.OutSensor1, SensorInf.S_In1) annotation(
        Line(points = {{-82, 92}, {-50, 92}, {-50, 92}, {-52, 92}}, color = {78, 154, 6}));
      connect(PLS.CairR5, B3.Qset) annotation(
        Line(points = {{45, -14}, {79, -14}, {79, 31}, {76, 31}}, color = {78, 154, 6}));
      connect(PLS.CairR4, B2.Qset) annotation(
        Line(points = {{45, -10}, {53, -10}, {53, 41}, {50, 41}}, color = {78, 154, 6}));
      connect(PLS.CairR3, B1.Qset) annotation(
        Line(points = {{45, -6}, {47, -6}, {47, 31.8}, {24, 31.8}, {24, 31}}, color = {78, 154, 6}));
      connect(PLS.C_WS, WS_Pump.Qset) annotation(
        Line(points = {{27, -14}, {21.2, -14}, {21.2, -81}, {32, -81}}, color = {78, 154, 6}));
      connect(PLS.C_RS, RS_Pump.Qset) annotation(
        Line(points = {{27, -10}, {-18.8, -10}, {-18.8, -81}, {-28, -81}}, color = {78, 154, 6}));
      connect(PLS.C_REC, REC_Pump.Qset) annotation(
        Line(points = {{27, -6}, {7.2, -6}, {7.2, 15.6}, {-12.8, 15.6}, {-12.8, 15}, {-50, 15}}, color = {78, 154, 6}));
      connect(DN2.OutSensor1, PLS.InSensorDN2) annotation(
        Line(points = {{-9, 69}, {-9, 89.3}, {91.3, 89.3}, {91.3, -30.7}, {33, -30.7}, {33, -19}}, color = {78, 154, 6}));
      connect(N2.OutSensor1, PLS.InSensorN2) annotation(
        Line(points = {{45, 69}, {45, 87.3}, {87.3, 87.3}, {87.3, -6.7}, {38, -6.7}, {38, -19}}, color = {78, 154, 6}));
      connect(N3.OutSensor1, PLS.InSensorN3) annotation(
        Line(points = {{69, 69}, {69, 85.3}, {85.3, 85.3}, {85.3, -2.7}, {43, -2.7}, {43, -19}}, color = {78, 154, 6}));
      connect(SC.Out1, Clear.In1) annotation(
        Line(points = {{12, -34}, {72, -34}}));
      connect(RS_Pump.Out1, M1.In2) annotation(
        Line(points = {{-48, -81}, {-88, -81}, {-88, 78}, {-74, 78}}));
      connect(S1.Out1, SC.In1) annotation(
        Line(points = {{-33, 5}, {-47, 5}, {-47, -39}, {-8, -39}}));
      connect(N3.Out1, S1.In1) annotation(
        Line(points = {{74, 64}, {84, 64}, {84, 4}, {-4, 4}, {-4, 8}, {-14, 8}}));
      connect(REC_Pump.Out1, M1.In3) annotation(
        Line(points = {{-70, 15}, {-82, 15}, {-82, 74}, {-74, 74}}));
      connect(S1.Out2, REC_Pump.In1) annotation(
        Line(points = {{-33, 9.5}, {-43, 9.5}, {-43, 9}, {-50, 9}}));
      connect(S2.Out2, RS_Pump.In1) annotation(
        Line(points = {{0.5, -77}, {0.5, -87}, {-28, -87}}));
      connect(S2.Out1, WS_Pump.In1) annotation(
        Line(points = {{5, -77}, {6.6, -77}, {6.6, -87}, {32, -87}}));
      connect(SC.Out2, S2.In1) annotation(
        Line(points = {{2, -50}, {2, -58}}));
      connect(WS_Pump.Out1, Sludge.In1) annotation(
        Line(points = {{52, -81}, {72, -81}, {72, -82}}));
      connect(Inf.Out1, M1.In1) annotation(
        Line(points = {{-80, 82}, {-74, 82}, {-74, 82}, {-74, 82}}));
      connect(N2.Out1, N3.In1) annotation(
        Line(points = {{50, 64}, {54, 64}}));
      connect(N1.Out1, N2.In1) annotation(
        Line(points = {{24, 64}, {30, 64}}));
      connect(DN2.Out1, N1.In1) annotation(
        Line(points = {{-4, 64}, {4, 64}}));
      connect(DN1.Out1, DN2.In1) annotation(
        Line(points = {{-28, 64}, {-24, 64}}));
      connect(B3.Qair, N3.InQair) annotation(
        Line(points = {{65, 44}, {65, 64.8}, {64, 64.8}, {64, 54}}, color = {78, 154, 6}));
      connect(B2.Qair, N2.InQair) annotation(
        Line(points = {{39, 54}, {39, 64.8}, {40, 64.8}, {40, 54}}, color = {78, 154, 6}));
      connect(B1.Qair, N1.InQair) annotation(
        Line(points = {{13, 44}, {13, 64.8}, {14, 64.8}, {14, 54}}, color = {78, 154, 6}));
      connect(M2.Out1, DN1.In1) annotation(
        Line(points = {{-48, 43.5}, {-38, 43.5}, {-38, 52}, {-50, 52}, {-50, 64}, {-48, 64}}));
      connect(Inf.OutSensor1, Dosing_Fe.SensorInfDos) annotation(
        Line(points = {{-82, 92}, {-80, 92}, {-80, 24}, {-23, 24}, {-23, 37}}, color = {78, 154, 6}));
      connect(Dosing_Fe.Out1, M2.In2) annotation(
        Line(points = {{-42, 30}, {-74, 30}, {-74, 42}, {-68, 42}}));
      connect(M1.Out1, M2.In1) annotation(
        Line(points = {{-54, 78}, {-48, 78}, {-48, 72}, {-52, 72}, {-52, 58}, {-74, 58}, {-74, 46}, {-68, 46}}));
      annotation(
        experiment(StartTime = 0, StopTime = 14, Tolerance = 1e-06, Interval = 0.02));
    end WWTP_COST;

    model WWTP_BWWT
      OpenWasteWater.ASM1P.Inflow Inf(Inf_File = "/home/awwjb/Git/OpenWasteWater/Resources/ASM1P/Inf_strm_P.txt", T = 20) annotation(
        Placement(visible = true, transformation(origin = {-90, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Mixer3 M1 annotation(
        Placement(visible = true, transformation(origin = {-64, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN1(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 30.2, NO = 8.9, ND = 0.9, ALK = 4.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 3250, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-38, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN2(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 30.2, NO = 1.9, ND = 0.9, ALK = 5.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 3250, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-14, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N1(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 10.2, NO = 8.9, ND = 0.9, ALK = 4.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 3450, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {14, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N2(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 5.2, NO = 8.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 3450, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {40, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N3(Sini = Soluble(I = 30, S = 1.15, O2 = 4.0, NH = 0.2, NO = 10.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.0, Fe3 = 0.1, Al3 = 0.0), V_R = 3450, Xini = Particulate(I = 1168, S = 28, H = 2088.0, A = 136.0, P = 564.0, ND = 2.55, FeP = 600.0, FeOH = 50.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {64, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B1 annotation(
        Placement(visible = true, transformation(origin = {14, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B2 annotation(
        Placement(visible = true, transformation(origin = {40, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B3 annotation(
        Placement(visible = true, transformation(origin = {66, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.controller PLS(Q_REC = 20000, Q_RS = 18446, Q_WS = 350, Q_air_N3 = 10000, k_D_O2 = 0, k_I_NO = 5000, k_I_O2 = 5000, k_P_NO = 2000, k_P_O2 = 2500) annotation(
        Placement(visible = true, transformation(origin = {36, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump REC_Pump(Qmax = 92239) annotation(
        Placement(visible = true, transformation(origin = {-60, 12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump RS_Pump annotation(
        Placement(visible = true, transformation(origin = {-38, -84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump WS_Pump annotation(
        Placement(visible = true, transformation(origin = {42, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.SecClar.SCL SC(A = 2000)  annotation(
        Placement(visible = true, transformation(origin = {2, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Effluent Clear annotation(
        Placement(visible = true, transformation(origin = {82, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.WasteSludge Sludge annotation(
        Placement(visible = true, transformation(origin = {82, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S1 annotation(
        Placement(visible = true, transformation(origin = {-24, 8}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S2 annotation(
        Placement(visible = true, transformation(origin = {2, -68}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.TechUnits.Sensor SensorInf(t_D = 0.01) annotation(
        Placement(visible = true, transformation(origin = {-42, 92}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.Mixer2 M2 annotation(
        Placement(visible = true, transformation(origin = {-58, 44}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Dosing Dosing_Fe(InAl3 = 7500, InFe2 = 20000, InFe3 = 30000) annotation(
        Placement(visible = true, transformation(origin = {-32, 36}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(SensorInf.S_Out1, PLS.InSensorInf) annotation(
        Line(points = {{-32, 92}, {98, 92}, {98, -38}, {30, -38}, {30, 4}, {29, 4}, {29, -19}}, color = {78, 154, 6}));
      connect(Inf.OutSensor1, SensorInf.S_In1) annotation(
        Line(points = {{-82, 92}, {-50, 92}, {-50, 92}, {-52, 92}}, color = {78, 154, 6}));
      connect(PLS.CairR5, B3.Qset) annotation(
        Line(points = {{45, -14}, {79, -14}, {79, 31}, {76, 31}}, color = {78, 154, 6}));
      connect(PLS.CairR4, B2.Qset) annotation(
        Line(points = {{45, -10}, {53, -10}, {53, 41}, {50, 41}}, color = {78, 154, 6}));
      connect(PLS.CairR3, B1.Qset) annotation(
        Line(points = {{45, -6}, {47, -6}, {47, 31.8}, {24, 31.8}, {24, 31}}, color = {78, 154, 6}));
      connect(PLS.C_WS, WS_Pump.Qset) annotation(
        Line(points = {{27, -14}, {21.2, -14}, {21.2, -81}, {32, -81}}, color = {78, 154, 6}));
      connect(PLS.C_RS, RS_Pump.Qset) annotation(
        Line(points = {{27, -10}, {-18.8, -10}, {-18.8, -81}, {-28, -81}}, color = {78, 154, 6}));
      connect(PLS.C_REC, REC_Pump.Qset) annotation(
        Line(points = {{27, -6}, {7.2, -6}, {7.2, 15.6}, {-12.8, 15.6}, {-12.8, 15}, {-50, 15}}, color = {78, 154, 6}));
      connect(DN2.OutSensor1, PLS.InSensorDN2) annotation(
        Line(points = {{-9, 69}, {-9, 89.3}, {91.3, 89.3}, {91.3, -30.7}, {33, -30.7}, {33, -19}}, color = {78, 154, 6}));
      connect(N2.OutSensor1, PLS.InSensorN2) annotation(
        Line(points = {{45, 69}, {45, 87.3}, {87.3, 87.3}, {87.3, -6.7}, {38, -6.7}, {38, -19}}, color = {78, 154, 6}));
      connect(N3.OutSensor1, PLS.InSensorN3) annotation(
        Line(points = {{69, 69}, {69, 85.3}, {85.3, 85.3}, {85.3, -2.7}, {43, -2.7}, {43, -19}}, color = {78, 154, 6}));
      connect(SC.Out1, Clear.In1) annotation(
        Line(points = {{12, -34}, {72, -34}}));
      connect(RS_Pump.Out1, M1.In2) annotation(
        Line(points = {{-48, -81}, {-88, -81}, {-88, 78}, {-74, 78}}));
      connect(S1.Out1, SC.In1) annotation(
        Line(points = {{-33, 5}, {-47, 5}, {-47, -39}, {-8, -39}}));
      connect(N3.Out1, S1.In1) annotation(
        Line(points = {{74, 64}, {84, 64}, {84, 4}, {-4, 4}, {-4, 8}, {-14, 8}}));
      connect(REC_Pump.Out1, M1.In3) annotation(
        Line(points = {{-70, 15}, {-82, 15}, {-82, 74}, {-74, 74}}));
      connect(S1.Out2, REC_Pump.In1) annotation(
        Line(points = {{-33, 9.5}, {-43, 9.5}, {-43, 9}, {-50, 9}}));
      connect(S2.Out2, RS_Pump.In1) annotation(
        Line(points = {{0.5, -77}, {0.5, -87}, {-28, -87}}));
      connect(S2.Out1, WS_Pump.In1) annotation(
        Line(points = {{5, -77}, {6.6, -77}, {6.6, -87}, {32, -87}}));
      connect(SC.Out2, S2.In1) annotation(
        Line(points = {{2, -50}, {2, -58}}));
      connect(WS_Pump.Out1, Sludge.In1) annotation(
        Line(points = {{52, -81}, {72, -81}, {72, -82}}));
      connect(Inf.Out1, M1.In1) annotation(
        Line(points = {{-80, 82}, {-74, 82}, {-74, 82}, {-74, 82}}));
      connect(N2.Out1, N3.In1) annotation(
        Line(points = {{50, 64}, {54, 64}}));
      connect(N1.Out1, N2.In1) annotation(
        Line(points = {{24, 64}, {30, 64}}));
      connect(DN2.Out1, N1.In1) annotation(
        Line(points = {{-4, 64}, {4, 64}}));
      connect(DN1.Out1, DN2.In1) annotation(
        Line(points = {{-28, 64}, {-24, 64}}));
      connect(B3.Qair, N3.InQair) annotation(
        Line(points = {{65, 44}, {65, 64.8}, {64, 64.8}, {64, 54}}, color = {78, 154, 6}));
      connect(B2.Qair, N2.InQair) annotation(
        Line(points = {{39, 54}, {39, 64.8}, {40, 64.8}, {40, 54}}, color = {78, 154, 6}));
      connect(B1.Qair, N1.InQair) annotation(
        Line(points = {{13, 44}, {13, 64.8}, {14, 64.8}, {14, 54}}, color = {78, 154, 6}));
      connect(M2.Out1, DN1.In1) annotation(
        Line(points = {{-48, 43.5}, {-38, 43.5}, {-38, 52}, {-50, 52}, {-50, 64}, {-48, 64}}));
      connect(Inf.OutSensor1, Dosing_Fe.SensorInfDos) annotation(
        Line(points = {{-82, 92}, {-80, 92}, {-80, 24}, {-23, 24}, {-23, 37}}, color = {78, 154, 6}));
      connect(Dosing_Fe.Out1, M2.In2) annotation(
        Line(points = {{-42, 30}, {-74, 30}, {-74, 42}, {-68, 42}}));
      connect(M1.Out1, M2.In1) annotation(
        Line(points = {{-54, 78}, {-48, 78}, {-48, 72}, {-52, 72}, {-52, 58}, {-74, 58}, {-74, 46}, {-68, 46}}));
      annotation(
        experiment(StartTime = 0, StopTime = 14, Tolerance = 1e-06, Interval = 0.02));
    end WWTP_BWWT;

    model WWTP_ADM
      OpenWasteWater.ASM1P.WasteSludge PS_and_WS annotation(
        Placement(visible = true, transformation(origin = {82, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Effluent Eff annotation(
        Placement(visible = true, transformation(origin = {90, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Inflow Inflow1(Inf_File = "/home/awwjb/Git/OpenWasteWater/Resources/ASM1P/Inf_Raw_dry_P.txt") annotation(
        Placement(visible = true, transformation(origin = {-82, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.SedTank.PreClar Vorklaerung annotation(
        Placement(visible = true, transformation(origin = {-46, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump PS_Pump(Qmax = 700, k_RT = 50) annotation(
        Placement(visible = true, transformation(origin = {-80, 14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.Sub_WWTP_Modells.WWTP_Part wWTP_Part1 annotation(
        Placement(visible = true, transformation(origin = {44, 66}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      OpenWasteWater.ASM1P.SedTank.Centrifuge WasteSludgeCentrifuge annotation(
        Placement(visible = true, transformation(origin = {0, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump WS_Pump(Qmax = 100, k_RT = 50) annotation(
        Placement(visible = true, transformation(origin = {12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.Mixer2 M2 annotation(
        Placement(visible = true, transformation(origin = {-76, -24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.TechUnits.SludgeControl sludgeControl_ES(SP_Q_max = 100, SP_Q_min = 2, k_I_Q = 0.2, k_P_Q = 0.2) annotation(
        Placement(visible = true, transformation(origin = {40, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.SedTank.StorageTank StorageTank1(Vini = 10) annotation(
        Placement(visible = true, transformation(origin = {-56, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump Pump_Dig(Qini = 100.0, k_RT = 50) annotation(
        Placement(visible = true, transformation(origin = {-2, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.VolumeControl volumeControl1(SP_Q_max = 80, SP_Q_min = 50, SP_V = 20, k_D_Q = 0.5, k_I_Q = 2, k_P_Q = 1) annotation(
        Placement(visible = true, transformation(origin = {-22, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ADM_P.inflow_ASM inflow_ADM1 annotation(
        Placement(visible = true, transformation(origin = {24, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ADM_P.outflow_ASM outflow_ASM1 annotation(
        Placement(visible = true, transformation(origin = {-30, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ADM_P.adm_models.ADM1 ADM1(Gini = OpenWasteWater.ADM_P.Gaseous_ADM(CO2 = 14.6, H2 = 0.0037, NH3 = 0.005, CH4 = 1646.0), Sini = OpenWasteWater.ADM_P.Solubles_ADM(su = 11.5, aa = 5.4, fa = 100.0, va = 10.5, bu = 14.3, pro = 16.85, ac = 44.0, H2 = 0.00024, ch4 = 130, IC = 55.5, IN = 660.0, I = 5700.0, P = 54.5, Fe2 = 0, Fe3 = 0, Al3 = 0), V_D = 1200, V_H = 120, Xini = OpenWasteWater.ADM_P.Particulates_ADM(c = 5500, ch = 53, pr = 54, li = 67, su = 835, aa = 630, fa = 560, c4 = 273, pro = 129, ac = 820, H2 = 390, I = 34000, FeP = 550, FeOH = 220, AlP = 0, AlOH = 0)) annotation(
        Placement(visible = true, transformation(origin = {58, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Sub_WWTP_Modells.SludgeThickening sludgeThickening1(k_I_Q = 1, k_P_Q = 1) annotation(
        Placement(visible = true, transformation(origin = {28, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Mixer3 M1 annotation(
        Placement(visible = true, transformation(origin = {-6, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  TechUnits.SetQ setQPS(Q = 30)  annotation(
        Placement(visible = true, transformation(origin = {-56, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Inflow1.Out1, Vorklaerung.In1) annotation(
        Line(points = {{-72, 75}, {-47, 75}, {-47, 74}, {-56, 74}}));
      connect(Vorklaerung.Out2, PS_Pump.In1) annotation(
        Line(points = {{-46, 64}, {-46, 56}, {-83, 56}, {-83, 24}}));
      connect(sludgeThickening1.Out2, PS_and_WS.In1) annotation(
        Line(points = {{37, -89}, {60, -89}, {60, -80}, {72, -80}}));
      connect(M1.Out1, wWTP_Part1.In1) annotation(
        Line(points = {{4, 74}, {10, 74}, {10, 82}, {24, 82}}));
      connect(WasteSludgeCentrifuge.Out1, M1.In3) annotation(
        Line(points = {{-10, 22}, {-12, 22}, {-12, 16}, {-30, 16}, {-30, 69}, {-16, 69}}));
      connect(outflow_ASM1.Out1, sludgeThickening1.In1) annotation(
        Line(points = {{-24, -86}, {-4, -86}, {-4, -79}, {19, -79}}));
      connect(ADM1.Out1, outflow_ASM1.In1) annotation(
        Line(points = {{68, -50}, {82, -50}, {82, -66}, {-50, -66}, {-50, -90}, {-36, -90}, {-36, -86}}, color = {233, 185, 110}));
      connect(inflow_ADM1.Out1, ADM1.In1) annotation(
        Line(points = {{30, -50}, {48, -50}, {48, -50}, {50, -50}}, color = {233, 185, 110}));
      connect(Pump_Dig.Out1, inflow_ADM1.In1) annotation(
        Line(points = {{8, -50}, {18, -50}, {18, -50}, {18, -50}}));
      connect(volumeControl1.Out1Q, Pump_Dig.Qset) annotation(
        Line(points = {{-22, -40}, {-22, -51}, {-12, -51}}, color = {78, 154, 6}));
      connect(StorageTank1.Out1, Pump_Dig.In1) annotation(
        Line(points = {{-46, -54}, {-25, -54}, {-25, -57}, {-12, -57}}));
      connect(StorageTank1.OutV_StorateTank, volumeControl1.In1V) annotation(
        Line(points = {{-48, -44}, {-48, -26}, {-22, -26}, {-22, -32}}, color = {78, 154, 6}));
      connect(M2.Out1, StorageTank1.In1) annotation(
        Line(points = {{-76, -34}, {-76, -54}, {-66, -54}}));
      connect(sludgeControl_ES.Out1Q, WS_Pump.Qset) annotation(
        Line(points = {{40, 24}, {40, 24}, {40, 16}, {14, 16}, {14, 10}, {14, 10}}, color = {78, 154, 6}));
      connect(WasteSludgeCentrifuge.Out2WWSensor, sludgeControl_ES.In1WWSensor) annotation(
        Line(points = {{10, 32}, {8, 32}, {8, 40}, {40, 40}, {40, 32}, {40, 32}}, color = {78, 154, 6}));
      connect(WasteSludgeCentrifuge.Out2, WS_Pump.In1) annotation(
        Line(points = {{8, 19}, {8, 10}}));
      connect(wWTP_Part1.Out2, WasteSludgeCentrifuge.In1) annotation(
        Line(points = {{64, 54}, {70, 54}, {70, 44}, {-24, 44}, {-24, 26}, {-10, 26}}));
      connect(WS_Pump.Out1, M2.In1) annotation(
        Line(points = {{16, -10}, {16, -10}, {16, -20}, {-48, -20}, {-48, -4}, {-74, -4}, {-74, -14}, {-74, -14}}));
      connect(PS_Pump.Out1, M2.In2) annotation(
        Line(points = {{-76, 4}, {-78, 4}, {-78, -14}, {-78, -14}}));
      connect(Inflow1.OutSensor1, wWTP_Part1.InfWWSensor) annotation(
        Line(points = {{-74, 86}, {24, 86}, {24, 84}}, color = {78, 154, 6}));
      connect(wWTP_Part1.Out1, Eff.In1) annotation(
        Line(points = {{64, 62}, {76, 62}, {76, 26}, {80, 26}, {80, 26}}));
      connect(ADM1.Out1, outflow_ASM1.In1) annotation(
        Line(points = {{68, -50}, {80, -50}, {80, -72}, {-52, -72}, {-52, -86}, {-36, -86}}, color = {233, 185, 110}));
      connect(outflow_ASM1.Out1, sludgeThickening1.In1) annotation(
        Line(points = {{-24, -86}, {6, -86}, {6, -78}, {18, -78}}));
      connect(Vorklaerung.Out1, M1.In1) annotation(
        Line(points = {{-36, 74}, {-32, 74}, {-32, 82}, {-20, 82}, {-20, 78}, {-16, 78}}));
      connect(sludgeThickening1.Out1, M1.In2) annotation(
        Line(points = {{18, -88}, {4, -88}, {4, -90}, {-92, -90}, {-92, 60}, {-34, 60}, {-34, 72}, {-26, 72}, {-26, 74}, {-16, 74}}));
  connect(setQPS.OutQ1, PS_Pump.Qset) annotation(
        Line(points = {{-66, 36}, {-78, 36}, {-78, 24}}, color = {78, 154, 6}));
      annotation(
        experiment(StartTime = 0, StopTime = 14, Tolerance = 1e-06, Interval = 0.0104167));
    end WWTP_ADM;

    model Test
      OpenWasteWater.ASM1P.Inflow Inf1(Inf_File = "/home/awwjb/Git/OpenWasteWater/Resources/ASM1P/Inf_strm_P.txt") annotation(
        Placement(visible = true, transformation(origin = {-70, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Effluent Effluent1 annotation(
        Placement(visible = true, transformation(origin = {84, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN1(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {0, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN2(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {34, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N1(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 1.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), V_R = 1333, Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-40, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N2(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {4, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N3(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), V_R = 1333, Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {50, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B1 annotation(
        Placement(visible = true, transformation(origin = {-40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B2 annotation(
        Placement(visible = true, transformation(origin = {6, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B3 annotation(
        Placement(visible = true, transformation(origin = {52, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SetQ setQ1(Q = 15000) annotation(
        Placement(visible = true, transformation(origin = {-18, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SetQ setQ2(Q = 15000) annotation(
        Placement(visible = true, transformation(origin = {30, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SetQ setQ3(Q = 5000) annotation(
        Placement(visible = true, transformation(origin = {76, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      OpenWasteWater.ASM1P.SecClar.SCL SCL1 annotation(
        Placement(visible = true, transformation(origin = {20, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.WasteSludge Sludge1 annotation(
        Placement(visible = true, transformation(origin = {84, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Mixer3 M3 annotation(
        Placement(visible = true, transformation(origin = {-40, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S1 annotation(
        Placement(visible = true, transformation(origin = {80, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S2 annotation(
        Placement(visible = true, transformation(origin = {20, -74}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.TechUnits.Pump P_Rec annotation(
        Placement(visible = true, transformation(origin = {82, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      OpenWasteWater.ASM1P.TechUnits.Pump P_RS annotation(
        Placement(visible = true, transformation(origin = {-46, -68}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump P_WS annotation(
        Placement(visible = true, transformation(origin = {56, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SetQ setQ4(Q = 18000) annotation(
        Placement(visible = true, transformation(origin = {-20, -64}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SetQ setQ5(Q = 365) annotation(
        Placement(visible = true, transformation(origin = {36, -62}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SetQ setQ6(Q = 55000) annotation(
        Placement(visible = true, transformation(origin = {64, 42}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
    equation
      connect(B1.Qair, N1.InQair) annotation(
        Line(points = {{-40, 0}, {6, 0}}, color = {78, 154, 6}));
      connect(B2.Qair, N2.InQair) annotation(
        Line(points = {{6, 0}, {4, 0}, {4, 10}}, color = {78, 154, 6}));
      connect(B3.Qair, N3.InQair) annotation(
        Line(points = {{52, -2}, {50, -2}, {50, 10}}, color = {78, 154, 6}));
      connect(setQ1.OutQ1, B1.Qset) annotation(
        Line(points = {{-24, -14}, {-30, -14}, {-30, -12}}, color = {78, 154, 6}));
      connect(setQ2.OutQ1, B2.Qset) annotation(
        Line(points = {{24, -14}, {16, -14}, {16, -12}}, color = {78, 154, 6}));
      connect(setQ3.OutQ1, B3.Qset) annotation(
        Line(points = {{70, -14}, {62, -14}}, color = {78, 154, 6}));
      connect(M3.Out1, DN1.In1) annotation(
        Line(points = {{-30, 64}, {-10, 64}}));
      connect(DN1.Out1, DN2.In1) annotation(
        Line(points = {{10, 64}, {24, 64}}));
      connect(DN2.Out1, N1.In1) annotation(
        Line(points = {{44, 64}, {56, 64}, {56, 34}, {-62, 34}, {-62, 20}, {-50, 20}}));
      connect(N1.Out1, N2.In1) annotation(
        Line(points = {{-30, 20}, {-6, 20}}));
      connect(N2.Out1, N3.In1) annotation(
        Line(points = {{14, 20}, {40, 20}}));
      connect(N3.Out1, S1.In1) annotation(
        Line(points = {{60, 20}, {70, 20}}));
      connect(S1.Out1, P_Rec.In1) annotation(
        Line(points = {{90, 22}, {94, 22}, {94, 40}, {86, 40}, {86, 56}}));
      connect(S1.Out2, SCL1.In1) annotation(
        Line(points = {{90, 18}, {94, 18}, {94, -24}, {-10, -24}, {-10, -46}, {10, -46}}));
      connect(SCL1.Out2, S2.In1) annotation(
        Line(points = {{20, -58}, {20, -64}}));
      connect(S2.Out1, P_WS.In1) annotation(
        Line(points = {{22, -82}, {24, -82}, {24, -94}, {38, -94}, {38, -78}, {46, -78}}));
      connect(P_WS.Out1, Sludge1.In1) annotation(
        Line(points = {{66, -70}, {68, -70}, {68, -74}, {74, -74}}));
      connect(SCL1.Out1, Effluent1.In1) annotation(
        Line(points = {{30, -42}, {74, -42}, {74, -44}}));
      connect(S2.Out2, P_RS.In1) annotation(
        Line(points = {{18, -82}, {18, -92}, {-12, -92}, {-12, -72}, {-36, -72}}));
      connect(setQ4.OutQ1, P_RS.Qset) annotation(
        Line(points = {{-26, -64}, {-31, -64}, {-31, -66}, {-36, -66}}, color = {78, 154, 6}));
      connect(setQ5.OutQ1, P_WS.Qset) annotation(
        Line(points = {{42, -62}, {46, -62}, {46, -72}}, color = {78, 154, 6}));
      connect(setQ6.OutQ1, P_Rec.Qset) annotation(
        Line(points = {{70, 42}, {80, 42}, {80, 56}}, color = {78, 154, 6}));
      connect(Inf1.Out1, M3.In1) annotation(
        Line(points = {{-60, 70}, {-50, 70}, {-50, 68}}));
      connect(P_Rec.Out1, M3.In2) annotation(
        Line(points = {{80, 76}, {78, 76}, {78, 94}, {-94, 94}, {-94, 64}, {-50, 64}}));
      connect(P_RS.Out1, M3.In3) annotation(
        Line(points = {{-56, -64}, {-82, -64}, {-82, 60}, {-50, 60}}));
      annotation(
        Icon(coordinateSystem(initialScale = 0.5)),
        experiment(StartTime = 0, StopTime = 14, Tolerance = 1e-06, Interval = 0.02));
    end Test;

    model Test1
      OpenWasteWater.ASM1P.Effluent Effluent1 annotation(
        Placement(visible = true, transformation(origin = {76, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WasteSludge sinkWater annotation(
        Placement(visible = true, transformation(origin = {76, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Sub_WWTP_Modells.WWTP_Part wWTP_Part annotation(
        Placement(visible = true, transformation(origin = {-16, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Inflow inflow(Inf_File = "/home/awwjb/Git/OpenWasteWater/Resources/ASM1P/Inf_dry_P.txt") annotation(
        Placement(visible = true, transformation(origin = {-70, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  OpenWasteWater.ADM_P.inflow_ASM inflow_ADM1 annotation(
        Placement(visible = true, transformation(origin = {-42, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ADM_P.adm_models.ADM1 adm1 annotation(
        Placement(visible = true, transformation(origin = {-10, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ADM_P.outflow_ASM outflow_ASM annotation(
        Placement(visible = true, transformation(origin = {28, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(wWTP_Part.Out1, Effluent1.In1) annotation(
        Line(points = {{-6, 74}, {66, 74}, {66, 72}}));
      connect(inflow.OutSensor1, wWTP_Part.InfWWSensor) annotation(
        Line(points = {{-62, 88}, {-26, 88}, {-26, 85}}, color = {78, 154, 6}));
  connect(wWTP_Part.Out2, inflow_ADM1.In1) annotation(
        Line(points = {{-6, 70}, {6, 70}, {6, 54}, {-70, 54}, {-70, 30}, {-48, 30}, {-48, 28}}));
  connect(inflow_ADM1.Out1, adm1.In1) annotation(
        Line(points = {{-36, 28}, {-18, 28}}, color = {233, 185, 110}));
  connect(adm1.Out1, outflow_ASM.In1) annotation(
        Line(points = {{0, 28}, {22, 28}}, color = {233, 185, 110}));
  connect(outflow_ASM.Out1, sinkWater.In1) annotation(
        Line(points = {{34, 28}, {66, 28}, {66, 26}}));
  connect(inflow.Out1, wWTP_Part.In1) annotation(
        Line(points = {{-60, 78}, {-36, 78}, {-36, 84}, {-26, 84}}));
      annotation(
        experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-6, Interval = 0.0102639));
    end Test1;
  end WWTP_Examples;

  package Sub_WWTP_Modells
    model WWTP_Part
      OpenWasteWater.ASM1P.Mixer3 M1 annotation(
        Placement(visible = true, transformation(origin = {-68, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN1(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-38, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.DenitrificationTank DN2(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {-12, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N1(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 1.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), V_R = 1333, Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {16, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N2(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {42, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.NitrificationTank N3(Sini = Soluble(I = 30, S = 5.15, O2 = 4.0, NH = 15.2, NO = 16.9, ND = 0.9, ALK = 3.54, P = 1.0, Fe2 = 0.1, Fe3 = 0.1, Al3 = 0.0), V_R = 1333, Xini = Particulate(I = 78.0, S = 37.4, H = 1216.0, A = 136.0, P = 208.0, ND = 2.55, FeP = 40.0, FeOH = 34.0, AlP = 0.0, AlOH = 0.0)) annotation(
        Placement(visible = true, transformation(origin = {66, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B1 annotation(
        Placement(visible = true, transformation(origin = {16, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B2 annotation(
        Placement(visible = true, transformation(origin = {42, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Blower B3 annotation(
        Placement(visible = true, transformation(origin = {66, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.controller PLS(Q_REC = 20000, Q_RS = 18446, Q_WS = 350, Q_air_N3 = 10000, k_D_O2 = 0, k_I_NO = 5000, k_I_O2 = 5000, k_P_NO = 2000, k_P_O2 = 2500) annotation(
        Placement(visible = true, transformation(origin = {36, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump REC_Pump(Qini = 55000.0, Qmax = 92239, k_RT = 50) annotation(
        Placement(visible = true, transformation(origin = {-58, 28}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump RS_Pump(Qini = 18000.0, k_RT = 50) annotation(
        Placement(visible = true, transformation(origin = {-36, -74}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump WS_Pump(Qini = 385.0, k_RT = 50) annotation(
        Placement(visible = true, transformation(origin = {44, -76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.SecClar.SCL SC annotation(
        Placement(visible = true, transformation(origin = {4, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S1 annotation(
        Placement(visible = true, transformation(origin = {-24, 22}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Divider2 S2 annotation(
        Placement(visible = true, transformation(origin = {4, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {96, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out2 annotation(
        Placement(visible = true, transformation(origin = {98, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {98, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Sensor InflowSensor annotation(
        Placement(visible = true, transformation(origin = {-40, 94}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      OpenWasteWater.ASM1P.InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      TechUnits.inWWSensor InfWWSensor annotation(
        Placement(visible = true, transformation(origin = {-100, 94}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Mixer2 M2 annotation(
        Placement(visible = true, transformation(origin = {-58, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.Dosing DosP annotation(
        Placement(visible = true, transformation(origin = {-30, 52}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(In1, M1.In1) annotation(
        Line(points = {{-102, 0}, {-96, 0}, {-96, 82}, {-78, 82}}, color = {0, 0, 255}));
      connect(InfWWSensor, InflowSensor.S_In1) annotation(
        Line(points = {{-100, 94}, {-50, 94}, {-50, 94}, {-50, 94}}, color = {78, 154, 6}));
      connect(InflowSensor.S_Out1, PLS.InSensorInf) annotation(
        Line(points = {{-30, 94}, {-14, 94}, {-14, 98}, {98, 98}, {98, -36}, {28, -36}, {28, 6}}, color = {78, 154, 6}));
      connect(WS_Pump.Out1, Out2) annotation(
        Line(points = {{53.8, -73}, {97.8, -73}}));
      connect(SC.Out1, Out1) annotation(
        Line(points = {{14.2, -28.3}, {96.2, -28.3}}));
      connect(PLS.CairR3, B1.Qset) annotation(
        Line(points = {{45, 18}, {56, 18}, {56, 30}, {30, 30}, {30, 48}, {26, 48}}, color = {78, 154, 6}));
      connect(PLS.CairR4, B2.Qset) annotation(
        Line(points = {{45, 14}, {60, 14}, {60, 34}, {56, 34}, {56, 48}, {52, 48}}, color = {78, 154, 6}));
      connect(PLS.CairR5, B3.Qset) annotation(
        Line(points = {{45, 10}, {80, 10}, {80, 46}, {76, 46}, {76, 48}}, color = {78, 154, 6}));
      connect(PLS.C_REC, REC_Pump.Qset) annotation(
        Line(points = {{27, 18}, {2, 18}, {2, 36}, {-48, 36}, {-48, 31}}, color = {78, 154, 6}));
      connect(PLS.C_RS, RS_Pump.Qset) annotation(
        Line(points = {{27, 14}, {-16, 14}, {-16, -71}, {-26, -71}}, color = {78, 154, 6}));
      connect(PLS.C_WS, WS_Pump.Qset) annotation(
        Line(points = {{27, 10}, {24, 10}, {24, -73}, {34, -73}}, color = {78, 154, 6}));
      connect(DN2.OutSensor1, PLS.InSensorDN2) annotation(
        Line(points = {{-6, 82}, {-6, 88}, {94, 88}, {94, -32}, {33, -32}, {33, 5}}, color = {78, 154, 6}));
      connect(N2.OutSensor1, PLS.InSensorN2) annotation(
        Line(points = {{48, 82}, {48, 86}, {90, 86}, {90, -8}, {38, -8}, {38, 5}}, color = {78, 154, 6}));
      connect(N3.OutSensor1, PLS.InSensorN3) annotation(
        Line(points = {{72, 82}, {72, 84}, {88, 84}, {88, -4}, {43, -4}, {43, 5}}, color = {78, 154, 6}));
      connect(RS_Pump.Out1, M1.In2) annotation(
        Line(points = {{-46, -71}, {-88, -71}, {-88, 78}, {-78, 78}}));
      connect(S1.Out1, SC.In1) annotation(
        Line(points = {{-33, 19}, {-44, 19}, {-44, -33}, {-6, -33}}));
      connect(N3.Out1, S1.In1) annotation(
        Line(points = {{76, 76}, {86, 76}, {86, 2}, {-2, 2}, {-2, 22}, {-14, 22}}));
      connect(REC_Pump.Out1, M1.In3) annotation(
        Line(points = {{-68, 31}, {-82, 31}, {-82, 73}, {-78, 73}}));
      connect(S1.Out2, REC_Pump.In1) annotation(
        Line(points = {{-33, 23.5}, {-40, 23.5}, {-40, 25}, {-48, 25}}));
      connect(S2.Out2, RS_Pump.In1) annotation(
        Line(points = {{2.5, -67}, {2.5, -67}, {2.5, -79}, {-25.5, -79}, {-25.5, -79}}));
      connect(S2.Out1, WS_Pump.In1) annotation(
        Line(points = {{6.6, -67}, {8.6, -67}, {8.6, -81}, {34.6, -81}, {34.6, -81}}));
      connect(SC.Out2, S2.In1) annotation(
        Line(points = {{4, -43.8}, {4, -47.8}}));
      connect(N2.Out1, N3.In1) annotation(
        Line(points = {{52, 76}, {56, 76}, {56, 76}, {56, 76}}));
      connect(N1.Out1, N2.In1) annotation(
        Line(points = {{26, 76}, {32, 76}, {32, 76}, {32, 76}}));
      connect(DN2.Out1, N1.In1) annotation(
        Line(points = {{-2, 76}, {6, 76}, {6, 76}, {6, 76}}));
      connect(DN1.Out1, DN2.In1) annotation(
        Line(points = {{-28, 76}, {-22, 76}, {-22, 76}, {-22, 76}}));
      connect(B3.Qair, N3.InQair) annotation(
        Line(points = {{66, 60}, {66, 60}, {66, 66}, {66, 66}}, color = {78, 154, 6}));
      connect(B2.Qair, N2.InQair) annotation(
        Line(points = {{42, 60}, {42, 60}, {42, 66}, {42, 66}}, color = {78, 154, 6}));
      connect(B1.Qair, N1.InQair) annotation(
        Line(points = {{16, 60}, {16, 60}, {16, 66}, {16, 66}}, color = {78, 154, 6}));
      connect(DosP.Out1, M2.In1) annotation(
        Line(points = {{-40, 46}, {-76, 46}, {-76, 54}, {-68, 54}}));
      connect(M1.Out1, M2.In2) annotation(
        Line(points = {{-58, 78}, {-54, 78}, {-54, 68}, {-78, 68}, {-78, 58}, {-68, 58}}));
      connect(M2.Out1, DN1.In1) annotation(
        Line(points = {{-48, 56}, {-46, 56}, {-46, 62}, {-52, 62}, {-52, 76}, {-48, 76}}));
      connect(InflowSensor.S_Out1, DosP.SensorInfDos) annotation(
        Line(points = {{-30, 94}, {-26, 94}, {-26, 64}, {-10, 64}, {-10, 54}, {-20, 54}}, color = {78, 154, 6}));
      annotation(
        experiment(StartTime = 0, StopTime = 14, Tolerance = 1e-6, Interval = 0.02),
        Diagram(coordinateSystem(initialScale = 0.1)));
    end WWTP_Part;

    model SludgeThickening
      parameter Real SP_TSS = 220000.0 "g/m3 TSS";
      parameter Real fns = 0.02;
      parameter Real SP_Q_min = 5.0 "m3/d";
      parameter Real SP_Q_max = 15 "m3/d";
      parameter Real k_P_Q = 0.2;
      parameter Real k_I_Q = 0.2;
      parameter Real k_D_Q = 0.0;
      OpenWasteWater.ASM1P.InPipe In1 annotation(
        Placement(visible = true, transformation(origin = {-92, 54}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {-94, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out1 annotation(
        Placement(visible = true, transformation(origin = {-94, -52}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {-94, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.OutPipe Out2 annotation(
        Placement(visible = true, transformation(origin = {93, -53}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {91, -51}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.SludgeControl sludgeControl1(SP_TSS = SP_TSS, SP_Q_min = SP_Q_min, SP_Q_max = SP_Q_max, k_P_Q = k_P_Q, k_I_Q = k_I_Q, k_D_Q = k_D_Q) annotation(
        Placement(visible = true, transformation(origin = {28, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SedTank.Centrifuge Centrifuge(fns = fns) annotation(
        Placement(visible = true, transformation(origin = {-6, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      OpenWasteWater.ASM1P.TechUnits.Pump pump(k_RT = 5) annotation(
        Placement(visible = true, transformation(origin = {24, -12}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
      connect(Centrifuge.Out1, Out1) annotation(
        Line(points = {{-16, 20}, {-16, 20}, {-16, 0}, {-52, 0}, {-52, -52}, {-94, -52}, {-94, -52}}));
      connect(In1, Centrifuge.In1) annotation(
        Line(points = {{-92, 54}, {-50, 54}, {-50, 26}, {-16, 26}, {-16, 24}, {-16, 24}}, color = {0, 0, 255}));
      connect(pump.Out1, Out2) annotation(
        Line(points = {{28, -22}, {30, -22}, {30, -54}, {94, -54}, {94, -52}}));
      connect(Centrifuge.Out2, pump.In1) annotation(
        Line(points = {{2, 16}, {2, 16}, {2, 4}, {20, 4}, {20, -2}, {20, -2}}));
      connect(sludgeControl1.Out1Q, pump.Qset) annotation(
        Line(points = {{28, 32}, {26, 32}, {26, -2}, {26, -2}}, color = {78, 154, 6}));
      connect(Centrifuge.Out2WWSensor, sludgeControl1.In1WWSensor) annotation(
        Line(points = {{4, 30}, {4, 30}, {4, 52}, {28, 52}, {28, 40}, {28, 40}}, color = {78, 154, 6}));
      annotation(
        Icon(graphics = {Text(origin = {-3, 50}, lineColor = {204, 0, 0}, lineThickness = 1.5, extent = {{-39, 18}, {35, -12}}, textString = "Sludge Tickening", textStyle = {TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold}), Rectangle(lineThickness = 0.75, extent = {{-74, 72}, {74, -72}}), Line(points = {{74, 72}, {-74, -72}}, thickness = 0.75)}, coordinateSystem(initialScale = 0.1)));
    end SludgeThickening;
  end Sub_WWTP_Modells;
end ASM1P;