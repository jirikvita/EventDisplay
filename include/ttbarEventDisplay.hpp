/* File ttbarEventDisplay.hpp
 *
 * Created       : Fri Jan 11 20:30:05 CST 2008
 * Author        : kvita using Dag's caf_mc_util::MCEventDisplay
 * Purpose       : 
 * Last modified : $Id: ttbarEventDisplay.hpp,v 1.3 2008/01/18 00:31:47 kvita Exp $
 * Comments      : 
 * Usage         : See TopEventDisplay.C
 *
 *   particles can be added in several ways using
 *   ttbarEventDisplay::AddParticle(TLorentzVector *vec, int pdgid);
 *   ttbarEventDisplay::AddParticle(TLorentzVector &vec, int pdgid);
 *   ttbarEventDisplay::AddParticlePxPyPzE(double Px, double Py, double Pz, double E, int pdgid);
 * See also
 *   display -> SetEtaMax(4.7);
 *   display -> SetLegoPhiTheta(-0.7*TMath::Pi()/4., 0.);
 *
 */

#ifndef ttbarEventDisplay_HPP_
#define ttbarEventDisplay_HPP_

#include <string>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
// #include "TPostScript.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLine.h"
//#include "TH2D.h"
//#include "TH3D.h"
//#include "TGraph2D.h"
#include "TString.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TCurlyLine.h"
//#include "TMath.h"

//#include "TGeoBBox.h"

using namespace std;

class ttbarEventDisplay {
  
 public:

 // Constructor, destructor: 
  ttbarEventDisplay(int runno = 0, int evtno = 0, bool doOneCavas = false);
  ~ttbarEventDisplay();
  
  // dag + jiri's modifications
  Int_t GetPdgColor( int id );
  TString GetPdgName( int id );
  TEllipse* DrawDot( float x, float y, float r, int pConservation );
  TCurlyLine* DrawParticleLine( float x, float y, float x2, float y2, int pdgid);
  TCurlyLine* DrawAxisLine( float x, float y, float x2, float y2, int lstyle = 2, int lwidth = 1, int lcol = 1);
  TLatex* DrawLatex( float x, float y, TString textstring, int col, float tsize );

  // jiri:
  void Clear();
  void Reset(int runno, int evtno);
  void AddParticle(TLorentzVector *vec, int pdgid);
  void AddParticle(TLorentzVector &vec, int pdgid);
  void AddParticlePxPyPzE(double Px, double Py, double Pz, double E, int pdgid);
  double GetMaxPt();
  double GetMaxE();
  double GetMaxP();

  void Translate3Dto2D(double x, double y, double z, double &xx, double &yy);
  void SetLegoPhiTheta(float phi, float theta = 0.);
  void SetEtaMax(float maxeta);

  void DrawRZ();
  void DrawXY();
  void DrawLego();
  void DrawAll();
  void PrintEps();
  void Print(TString suffix = "eps");

 private:
  
  int _runno;
  int _evtno;
  bool _doOneCavas;

  // pad 2D ranges:
  float _x0;
  float _y0;
  float _offx;
  float _offy;
  float _xmin, _xmax, _ymin, _ymax;

  // ranges for 3D lego:
  float _phi;
  float _theta;
  float _etamin, _etamax, _phimin, _phimax, _ptmin, _ptmax;
    

  int _screenReso;

  TCanvas *_globalCan;
  TCanvas *_cRZ;
  TCanvas *_cXY;
  TCanvas *_cLego;

  //  TH3D *_legoHist;

  vector< pair<TLorentzVector*, int> > _particles;
  vector<TCurlyLine*> _lines;
  vector<TEllipse*> _ellipses;
  vector<TLatex*> _latex;
  //  vector<TGraph2D*> _graphs2d;
};

#endif
