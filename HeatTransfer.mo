within ;
package HeatTransfer
  connector FlowInlet
    import Modelica.Units.SI;
  SI.Pressure p;
  flow SI.VolumeFlowRate q;
  SI.Temperature T;
    annotation (Icon(graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={0,140,72},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0})}));
  end FlowInlet;

  connector WallConnector
    import Modelica.Units.SI;

    SI.VolumeFlowRate q;
    SI.Temperature T1, T2, T;
    flow SI.HeatFlowRate Phi;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None)}),                           Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end WallConnector;

  record DuctParameters
    import Modelica.Units.SI;
    SI.Volume V;
    SI.Area Cv;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end DuctParameters;

  model BasicLiquid
    import Modelica.Units.SI;
    SI.Pressure p;
    SI.Temperature T;
    SI.Density rho = 1e3;
    SI.SpecificHeatCapacity c=1e3;
  end BasicLiquid;

  model Redeclare_Test
    import Modelica.Units.SI;
    model BasicWaterModel_1 = BasicLiquid (rho=1000, c=4180);
    model BasicWaterModel_2 = BasicLiquid (redeclare parameter Real rho=1000,
          redeclare parameter SI.SpecificHeatCapacity c=4180);

    BasicWaterModel_1 basic_1;
    BasicWaterModel_2 basic_2;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Redeclare_Test;

  model HEXDuctSection
    import Modelica.Units.SI;
    FlowInlet Port1
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    FlowInlet Port2
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    WallConnector Wall
      annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));

    parameter DuctParameters Pars;

    model Virtual_Liquid = BasicLiquid;

    Virtual_Liquid L1 "The liquid at Port1.";
    Virtual_Liquid L2 "The liquid at Port2.";

    SI.Temperature T;
    SI.Density rho;
    SI.SpecificHeatCapacity c;
    SI.Enthalpy H;


  equation
    // Hydraulics
    Port1.q + Port2.q = 0;
    Port1.p - Port2.p = 1e-3; //rho/Pars.Cv^2*abs(Port1.q)*Port1.q;

    //Thermodynamics
    der(H) = Wall.Phi +
    L1.c*L1.rho*Port1.q*Port1.T +
    L2.c*L2.rho*Port2.q*Port2.T;


    //H = rho*Pars.V*c*T;
    T = if Port1.q>0 then Port2.T else Port1.T;
    rho = if Port1.q>0 then L2.rho else L1.rho;
    c = if Port1.q>0 then L2.c else L1.c;

    // Communication
    Wall.q = Port1.q;
    Wall.T = T;
    Wall.T1 = Port1.T;
    Wall.T2 = Port2.T;

    L1.p = Port1.p;
    L2.p = Port2.p;
    L1.T = Port1.T;
    L2.T = Port2.T;




    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end HEXDuctSection;

  record ThermalSurfaceParameters
    import Modelica.Units.SI;
    SI.SurfaceCoefficientOfHeatTransfer h0;
    SI.VolumeFlowRate q0;
    Real n;
    SI.LinearTemperatureCoefficient ah;
    SI.Temperature T0 "nomial mean";

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end ThermalSurfaceParameters;

  record WallParameters
    import Modelica.Units.SI;
    SI.Area A;
    SI.Thickness d;
    SI.ThermalConductivity lambda;
    SI.Area SA, SB;
    Real Y;
    SI.ThermalResistance Rf "Fouling";

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end WallParameters;

  model HEXWallSection
    import Modelica.Units.SI;

    WallConnector WA
      annotation (Placement(transformation(extent={{-8,88},{12,108}})));
    WallConnector WB
      annotation (Placement(transformation(extent={{-8,-108},{12,-88}})));
    parameter WallParameters Pars_Wall;
    parameter ThermalSurfaceParameters Pars_Thermal;
    SI.ThermalResistance R, RA, RB, Rw;

  protected
    SI.Temperature DTlm, DT1, DT2;
  equation
    //Heat transfer
    WA.Phi + WB.Phi = 0;
    WB.Phi = DTlm/R;
    // Log-mean temperature difference

    DT1 = WA.T1 - WB.T1;
    DT2 = WA.T2 - WB.T2;
    DTlm =
      if (abs(DT1-DT2)) > 0.05*max(abs(DT1), abs(DT2))
        then (DT1-DT2)/ln(DT1/DT2)
      else if DT1*DT2 == 0 then 0.5*(DT1+DT2)
      else 0.5*(DT1+DT2)*
        (1-sqrt(DT1-DT2)/(DT1*DT2)*
        (1+sqrt(DT1-DT2)/(DT1*DT2)/2)/12);

    RA = 1/(Pars_Wall.SA*Pars_Thermal.h0*abs(WA.q/Pars_Wall.SA*Pars_Thermal.q0)^Pars_Wall.SA*Pars_Thermal.n*
    (1+Pars_Wall.SA*Pars_Thermal.ah*(WA.T-Pars_Wall.SA*Pars_Thermal.T0))*Pars_Wall.A);
    RB = 1/(Pars_Wall.SB*Pars_Thermal.h0*abs(WB.q/Pars_Wall.SB*Pars_Thermal.q0)^Pars_Wall.SB*Pars_Thermal.n*
    (1+Pars_Wall.SB*Pars_Thermal.ah*(WB.T-Pars_Wall.SB*Pars_Thermal.T0))*Pars_Wall.A);
    Rw = Pars_Wall.d/(Pars_Wall.lambda*Pars_Wall.Y*Pars_Wall.A);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end HEXWallSection;
end HeatTransfer;
