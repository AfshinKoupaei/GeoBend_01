// Grasshopper Script Instance
#region Usings
using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
#endregion

//   ██████  ██████  ██████  ███████     ██████  ██    ██                                                                                     
//  ██      ██    ██ ██   ██ ██          ██   ██  ██  ██                                                                                      
//  ██      ██    ██ ██   ██ █████       ██████    ████                                                                                       
//  ██      ██    ██ ██   ██ ██          ██   ██    ██                                                                                        
//   ██████  ██████  ██████  ███████     ██████     ██                                                                                        
//                                                                                                                                            
//                                                                                                                                            
//   █████  ███████ ███████ ██   ██ ██ ███    ██     ██   ██  ██████  ██    ██ ██████   █████  ███████ ██                                     
//  ██   ██ ██      ██      ██   ██ ██ ████   ██     ██  ██  ██    ██ ██    ██ ██   ██ ██   ██ ██      ██                                     
//  ███████ █████   ███████ ███████ ██ ██ ██  ██     █████   ██    ██ ██    ██ ██████  ███████ █████   ██                                     
//  ██   ██ ██           ██ ██   ██ ██ ██  ██ ██     ██  ██  ██    ██ ██    ██ ██      ██   ██ ██      ██                                     
//  ██   ██ ██      ███████ ██   ██ ██ ██   ████     ██   ██  ██████   ██████  ██      ██   ██ ███████ ██                                     
//                                                                                                                                            
//                                                                                                                                            
//  ██   ██ ████████ ████████ ██████         ██     ██  █████  ██   ██  █████  ██████  ██       █████  ██████      ██████  ██████  ███    ███ 
//  ██   ██    ██       ██    ██   ██ ██    ██     ██  ██   ██ ██  ██  ██   ██ ██   ██ ██      ██   ██ ██   ██    ██      ██    ██ ████  ████ 
//  ███████    ██       ██    ██████       ██     ██   ███████ █████   ███████ ██   ██ ██      ███████ ██████     ██      ██    ██ ██ ████ ██ 
//  ██   ██    ██       ██    ██      ██  ██     ██    ██   ██ ██  ██  ██   ██ ██   ██ ██      ██   ██ ██   ██    ██      ██    ██ ██  ██  ██ 
//  ██   ██    ██       ██    ██         ██     ██     ██   ██ ██   ██ ██   ██ ██████  ███████ ██   ██ ██████  ██  ██████  ██████  ██      ██ 
//                                                                                                                                            
//  All Rights Reserved, The use of parts or all of the code is permitted with the condition of referencing to the Author  : 
//  Afshin Koupaei, 2022 , http://akadlab.com


public class Script_Instance : GH_ScriptInstance
{
    #region Notes
    /* 
      Members:
        RhinoDoc RhinoDocument
        GH_Document GrasshopperDocument
        IGH_Component Component
        int Iteration

      Methods (Virtual & overridable):
        Print(string text)
        Print(string format, params object[] args)
        Reflect(object obj)
        Reflect(object obj, string method_name)
    */
    #endregion

    private void RunScript(
	Polyline StartSection,
	double KerfWidth,
	double LeftOverLayer,
	double FullPieceLength,
	double StartR,
	double StartL,
	List<double> DistListR,
	List<double> DistListL,
	bool blnRotateSection,
	ref object ExtensionSection,
	ref object KERFLines,
	ref object CrvPolygon,
	ref object BendingInformation)
    {
    //Declarations
    List<Curve> crvHrLines = new List<Curve>();
    List<Curve> crvCutPoly = new List<Curve>();
    Polyline crvEndSection = new Polyline();


    
    //Function
    crvEndSection = GeoBend(KerfWidth, LeftOverLayer, StartSection, blnRotateSection, FullPieceLength, StartR, StartL, DistListR, DistListL, ref crvHrLines, ref crvCutPoly);

    // OUTPUT
    ExtensionSection = crvEndSection;
    KERFLines = crvHrLines;
    BendingInformation = bndInfo;
    CrvPolygon = crvCutPoly;
      
    }

//   ██████  ███████  ██████      ██████  ███████ ███    ██ ██████      ███████ ██ ███    ███ ██████  ██      ███████ 
//  ██       ██      ██    ██     ██   ██ ██      ████   ██ ██   ██     ██      ██ ████  ████ ██   ██ ██      ██      
//  ██   ███ █████   ██    ██     ██████  █████   ██ ██  ██ ██   ██     ███████ ██ ██ ████ ██ ██████  ██      █████   
//  ██    ██ ██      ██    ██     ██   ██ ██      ██  ██ ██ ██   ██          ██ ██ ██  ██  ██ ██      ██      ██      
//   ██████  ███████  ██████      ██████  ███████ ██   ████ ██████      ███████ ██ ██      ██ ██      ███████ ███████ 
//                                                                                                                    
//                                                                                                                    

     public Polyline GeoBend(
    double kw, double LOL, Polyline crvSection, bool rotateSection, double pieceLength,
    double stR, double stL,
    List<double> distR, List<double> distL,
    ref List<Curve> crvHrLines, ref List<Curve> crvCutPoly)
  {
    ///<summary> The Function Bends the geometry according to parameters and returns the Section Curve for Continuation</summary>
    double docTolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;

    Curve crvBase = new PolylineCurve();
    Polyline crvEndSection = new Polyline();

    // <<<<<<<<<<<<<<<<<<<< SectionCurve Analysis and Base Curve production >>>>>>>>>>>>>>>>>>>>
    // Adjusting the short and long edge
    crvSection = plLengthAdjust(crvSection);
    Info("sectionCurve point count", crvSection.Count);
    // Defines the side of Operations
    if (rotateSection) plShiftRight(ref crvSection);
    // Drawing the intermediate curve as the basis for bending
    crvBase = plInterim(crvSection, LOL).ToPolylineCurve();
    //Defining the width of the material
    Info("MaterialWidth", (crvSection[1] - crvSection[0]).Length);
    //Defining the materialThickness
    double materialThickness = (crvSection[2] - crvSection[1]).Length;
    Info("MaterialThickness", materialThickness);
    Info("Lefover Layer", LOL);
    //KD retrieving KerfDepth
    double kd = materialThickness - LOL;
    Info("Kerf Depth", kd);
    //Defining Rotation Angle
    double rotAngle = Math.Asin((kw / 2) / Math.Sqrt(kd * kd + (kw / 2) * (kw / 2))) * 2;
    Info("Rotation Angle RAD", rotAngle);
    Info("Rotation Angle Degrees", RhinoMath.ToDegrees(rotAngle));

    //Defining the Section Plane // Should be after Curve Adjustment
    Plane plnSection = new Plane(crvSection[0], crvSection[1] - crvSection[0], crvSection[3] - crvSection[0]);
    //Retrieving main piece Direction according to the Sectopm Plane
    Vector3d vecDirection = plnSection.ZAxis;
    vecDirection.Unitize();

    //Start point and End Point of Base Curve
    Point3d ptR = crvBase.PointAtStart;
    Point3d ptL = crvBase.PointAtEnd;

    //workPlane  - Surface plane, to be used for polygons and for Milling Plane
    Plane plnSurface = new Plane(ptR, ptL - ptR, vecDirection);


    //vecExtension = new Vector3d (vecDirection);
    //crvEndSection cc ;

    crvEndSection = crvSection.Duplicate();
    crvEndSection.Transform(Rhino.Geometry.Transform.Translation(vecDirection * pieceLength));

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Defining the OFFSET VECTOR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Vector3d vecOffsetDir = new Vector3d (crvBase.PointAtStart - crvBase.PointAtEnd);
    vecOffsetDir.Rotate(Math.PI / 2, vecDirection);
    vecOffsetDir.Unitize();


    // defining the start point of KERFs on side Right and Side Left
    double stRLength = 0, stLLength = 0;
    if(stR > 0 && stR <= 1.0) {stRLength = pieceLength * stR;} else{stRLength = stR;}
    if(stL > 0 && stL <= 1.0) {stLLength = pieceLength * stL;} else{stLLength = stL;}


    //Curve List Definition
    List<Curve> _crvHrLines = new List<Curve>();


    //<<<<<<<<<<<<<<<<<<<<<<<<<<<   Generating initial HrLines, the kerf lines along the sheet direction  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    double distAddR=0, distAddL = 0;

    _crvHrLines.Add(new LineCurve(ptR, ptL));       // First Line on the ground
    for ( int i = 0 ; i < distR.Count + 1 ; i++ )   // Kerf Lines
    {
      _crvHrLines.Add(new LineCurve(vecDirection * (stRLength + distAddR) + ptR, vecDirection * (stLLength + distAddL) + ptL));
      if (i < distR.Count) { distAddR += distR[i];  distAddL += distL[i];}
    }
    _crvHrLines.Add(new LineCurve(vecDirection * pieceLength + ptR, vecDirection * pieceLength + ptL)); //Last Line
    Info("crvHrLines Count ", _crvHrLines.Count);


    //<<<<<<<<<<<<<<<<<<<<<<<<<<<   Defining the Cut polygons >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    List<Curve> _crvCutPoly = new List<Curve>();
    _crvCutPoly = polyFromHrLines(_crvHrLines, vecDirection, plnSurface, kw);



    //<<<<<<<<<<<<<<<<<<<<<<<<<<<  Main Curve Bending
    for (int i = 1 ; i < _crvHrLines.Count - 1 ; i++)
    {
      Vector3d rotAxis = new Vector3d(_crvHrLines[i].PointAtEnd - _crvHrLines[i].PointAtStart);
      for (int j = i ; j < _crvHrLines.Count - 1 ; j++)
      {

        _crvHrLines[j + 1].Rotate(rotAngle, rotAxis, _crvHrLines[i].PointAtStart);
        _crvCutPoly[j].Rotate(rotAngle, rotAxis, _crvHrLines[i].PointAtStart);
      }
      crvEndSection.Transform(Rhino.Geometry.Transform.Rotation(rotAngle, rotAxis, _crvHrLines[i].PointAtStart));
    }



    //<<<<<<<<<<<<<<<<<<<<<<<<<<< Returning values and setting by ref parameters

    crvCutPoly = _crvCutPoly;
    crvHrLines = _crvHrLines;

    return crvEndSection;

  }

//  ██████  ███████  ██████ ██       █████  ██████   █████  ████████ ██  ██████  ███    ██ ███████ 
//  ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██ ██      
//  ██   ██ █████   ██      ██      ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██ ███████ 
//  ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██      ██ 
//  ██████  ███████  ██████ ███████ ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████ ███████ 
//                                                                                                 
//                                                                                                 
//                                                                                                                                                       

  public PolylineCurve triFromPoints(Point3d pt0, Point3d pt1, Point3d pt2)
  {
    List<Point3d> polyPoints = new List<Point3d>();
    polyPoints.Add(pt0);
    polyPoints.Add(pt1);
    polyPoints.Add(pt2);
    polyPoints.Add(pt0);
    return new PolylineCurve(polyPoints);

  }

  public List<Curve> polyFromHrLines(List<Curve> _crvHrLines, Vector3d vecDirection, Plane plnSurface, double kw)
  {
    List<Curve> _crvCutPoly = new List<Curve>();
    for(int i = 0 ; i < _crvHrLines.Count - 1 ; i++)
    {
      List<Point3d> polyPoints = new List<Point3d>();

      double kerfMovDn = kw / 2 / Math.Sin(Vector3d.VectorAngle(vecDirection, _crvHrLines[i].PointAtStart - _crvHrLines[i].PointAtEnd, plnSurface));
      double kerfMovUp = kw / 2 / Math.Sin(Vector3d.VectorAngle(vecDirection, _crvHrLines[i + 1].PointAtStart - _crvHrLines[i + 1].PointAtEnd, plnSurface));

      polyPoints.Add(_crvHrLines[i].PointAtStart + vecDirection * (i == 0 ? 0 : kerfMovDn ));             // Point at Right , down
      polyPoints.Add(_crvHrLines[i].PointAtEnd + vecDirection * (i == 0 ? 0 : kerfMovDn ));               //Point at Left Down
      polyPoints.Add(_crvHrLines[i + 1].PointAtEnd + vecDirection * (i + 1 == _crvHrLines.Count - 1 ? 0 : kerfMovUp ) * -1);      //Point at Left Up
      polyPoints.Add(_crvHrLines[i + 1].PointAtStart + vecDirection * (i + 1 == _crvHrLines.Count - 1 ? 0 : kerfMovUp ) * -1);    //Point at Right Up
      polyPoints.Add(_crvHrLines[i].PointAtStart + vecDirection * (i == 0 ? 0 : kerfMovDn ));             // Point at Right down to close the curve
      _crvCutPoly.Add(new PolylineCurve(polyPoints));                                //making the poly line from the points
    }
    return _crvCutPoly;
  }

  public Brep ruledLoft(List<Curve> crvIn)
  {
    List<Brep> brepOut = new List<Brep>();

    for( int i = 0 ; i < crvIn.Count - 1; i++  )
    {
      brepOut.Add(NurbsSurface.CreateRuledSurface(crvIn[i], crvIn[i + 1]).ToBrep());
    }
    return Brep.JoinBreps(brepOut, 0.001)[0];

  }

  public List<Brep> ruledLoftList(List<Curve> crvIn)
  {
    List<Brep> brepOut = new List<Brep>();

    for( int i = 0 ; i < crvIn.Count - 1; i++  )
    {
      brepOut.Add(NurbsSurface.CreateRuledSurface(crvIn[i], crvIn[i + 1]).ToBrep());
    }
    return brepOut;

  }
  
  public enum sideDef {
    UpperSide,
  BottomSide}

  public double iDist(double vA, double vB, double rotationAngle, double kerfMove, double patchDistance)
  {
    //Rhino.Geometry.Triangle3d
    double vC = Math.Sqrt(vA * vA + vB * vB - 2 * vA * vB * Math.Cos(Math.PI - rotationAngle));
    double sinBeta = (Math.Sin(Math.PI - rotationAngle) * vB) / vC;

    return Math.Tan(Math.Asin(sinBeta)) * (vA - kerfMove) + patchDistance;
  }

  public double iDist(Triangle3d tri, sideDef side, double kerfMove)
  {
    if (side == sideDef.UpperSide){
      return Math.Tan(tri.AngleA) * (tri.AB.Length - kerfMove);
    }
    else{
      return Math.Tan(tri.AngleC) * (tri.BC.Length - kerfMove);
    }


  }

  public Polyline plLengthAdjust(Polyline src)
  { // test if the first segment of the rectangle is the longest segment
    if((src[1] - src[0]).Length < (src[2] - src[1]).Length)
    {/**/
      Polyline plTmp = new Polyline(5);
      plTmp.Add(src[1]);
      plTmp.Add(src[2]);
      plTmp.Add(src[3]);
      plTmp.Add(src[0]);
      plTmp.Add(src[1]);
      src = plTmp;
      /**/
      //src.Reverse();

    }
    return src;
  }

  public void plShiftRight(ref Polyline src)
  { // rotates the rectangle segments so the segments 0 and 2 are changed
    Polyline plTmp = new Polyline(5);
    plTmp.Add(src[2]);
    plTmp.Add(src[3]);
    plTmp.Add(src[0]);
    plTmp.Add(src[1]);
    plTmp.Add(src[2]);
    src = plTmp;

  }
  public Polyline plInterim (Polyline src, double dist)
  {
    Point3d ptR = new Point3d(), ptL = new Point3d();
    Vector3d movVec = new Vector3d(src[2] - src[1]);
    movVec.Unitize();
    ptR = src[1] + dist * movVec;
    ptL = src[0] + dist * movVec;
    return new Polyline(new List<Point3d>(){ptR,ptL});
  }


//  ██ ███    ██ ███████  ██████      ██████  ███████  ██████ 
//  ██ ████   ██ ██      ██    ██     ██   ██ ██      ██      
//  ██ ██ ██  ██ █████   ██    ██     ██   ██ █████   ██      
//  ██ ██  ██ ██ ██      ██    ██     ██   ██ ██      ██      
//  ██ ██   ████ ██       ██████      ██████  ███████  ██████ 
//                                                            
//                                                            



  public void pp(string txt, int val)
  {
    Print(txt + " : " + val.ToString());
  }
  public void pp(string txt, double val)
  {
    Print(txt + " : " + val.ToString());
  }
  public void pp(string txt)
  {
    Print(txt);
  }
  public void pp(string txt, List<String> itm)
  {
    String txt2 = "";
    foreach(String items in itm)
    {
      txt2 += items.ToString() + " , ";
    }
    Print(txt + " : " + txt2);
  }
  public void pp(string txt, List<int> itm)
  {
    String txt2 = "";
    foreach(int items in itm)
    {
      txt2 += items.ToString() + " , ";
    }
    Print(txt + " : " + txt2);
  }
  public void pp(string txt, List<double> itm)
  {
    String txt2 = "";
    foreach(double items in itm)
    {
      txt2 += items.ToString() + " , ";
    }
    Print(txt + " : " + txt2);
  }
  public void Info(string information)
  {
    bndInfo += information + "\n";

  }
  public void Info(int ln, string information)
  {
    bndInfo += "Code Line : " + ln.ToString() + " // " + information + "\n";
  }
  public void Info(string descriptionA, string descriptionB)
  {
    bndInfo += descriptionA + " // " + descriptionB + "\n";

  }

  public void Info(string descriptionA, double val)
  {
    bndInfo += descriptionA + " // " + val.ToString() + "\n";
  }
  string bndInfo = "Bending results and details : \n";

}
