/* File ttbarEventDisplay.hpp
 *
 * Created       : Fri Jan 11 20:30:05 CST 2008
 * Author        : kvita using Dag's caf_mc_util::MCEventDisplay
 * Purpose       : 
 * Last modified : $Id: ttbarEventDisplay.cpp,v 1.4 2008/07/01 00:46:46 kvita Exp $
 * Comments      : 
 */

#include <iostream>

 
// __________________________________________________________________


// jk, reminder (30.6.2008) from january:
// for compiling the package:
#include "ttbarEventDisplay.hpp"
// for compiling the script DataTopEventDisplay.C:
// #include "top_ptspectrum_cafe/top_ptspectrum_cafe/ttbarEventDisplay.hpp"

//#if !(defined(__CINT__) || defined(__MAKECINT__)) 
//# ifndef __CINT__
//# ifndef ROOT_VERSION
//#   include "top_ptspectrum_cafe/ttbarEventDisplay.hpp"
//# endif
// #if defined__CINT__ || __MAKECINT__
//# ifdef __CINT__
//# ifdef ROOT_VERSION
//#   include "top_ptspectrum_cafe/top_ptspectrum_cafe/ttbarEventDisplay.hpp"
//# endif



// __________________________________________________________________

// Constructor, destructor: 
ttbarEventDisplay::ttbarEventDisplay(int runno, int evtno, bool doOneCavas) : _cXY(0), _cRZ(0), 
									      _cLego(0), _globalCan(0)
{

  _doOneCavas = doOneCavas;

  SetEtaMax(4.7);
  _theta = 0.;
  //  _phi = 3*TMath::Pi()/8.;
  //  _phi = -TMath::Pi()/4.;
  //  _phi = -3*TMath::Pi()/8.;
  //  _phi = 0.;
  _phi = -0.7*TMath::Pi()/4.;

  Reset(runno, evtno);

}


ttbarEventDisplay::~ttbarEventDisplay()
{

}
// __________________________________________________________________

void ttbarEventDisplay::Clear()

{

  for (int i = 0; i < _particles.size(); ++i)
    delete _particles[i].first;
  _particles.clear();
  for (int i = 0; i < _lines.size(); ++i)
    delete _lines[i];
  _lines.clear();
  for (int i = 0; i < _ellipses.size(); ++i)
    delete _ellipses[i];
  _ellipses.clear();
  for (int i = 0; i < _latex.size(); ++i)
    delete _latex[i];
  _latex.clear();
  
  if (_doOneCavas) {
    if (_globalCan)
      delete _globalCan; 
  } else {
    if (_cXY)
      delete _cXY;
    if (_cRZ)
      delete _cRZ;
    if (_cLego)
      delete _cLego;
  }
  
  _globalCan = 0;
  _cXY = 0;
  _cRZ = 0;
  _cLego = 0;
  

}


// __________________________________________________________________

void ttbarEventDisplay::Reset(int runno, int evtno)
{

  _runno = runno;
  _evtno = evtno;


  _x0 = 0.;
  _y0 = 0.;
  _screenReso = 1200;
  
  Clear();

  int width = _screenReso / 3;
  int height = width;
  int offset = 15;
  if (_doOneCavas) {
    TString name = Form("GlobalDisplay_%i_%i", _runno, _evtno);
    _globalCan = new TCanvas(name, name, 0, 0, 3*width, height);
    _globalCan -> Divide(3, 1);
    //    _globalCan = new TCanvas(name, name, 0, 0, 2*width, 2*height);
    //    _globalCan -> Divide(2, 2);
    _globalCan -> cd(1);
    gPad -> SetName("ttbarEventDisplayXY");
    _globalCan -> cd(2);
    gPad -> SetName("ttbarEventDisplayRZ");
    _globalCan -> cd(3);
    gPad -> SetName("ttbarEventDisplayLego");
  } else {
    _cXY = new TCanvas(Form("ttbarEventDisplayXY_%i_%i", _runno, _evtno), 
		       Form("ttbarEventDisplayXY_%i_%i", _runno, _evtno), 0, 0, width, height);
    _cRZ = new TCanvas(Form("ttbarEventDisplayRZ_%i_%i", _runno, _evtno), 
		       Form("ttbarEventDisplayRZ_%i_%i", _runno, _evtno), width + offset, 0, width, height);
    _cLego = new TCanvas(Form("ttbarEventDisplayLego_%i_%i", _runno, _evtno), 
			 Form("ttbarEventDisplayLego_%i_%i", _runno, _evtno), width + offset, height + 2*offset, width, height);
  }

  
}

// __________________________________________________________________

void ttbarEventDisplay::AddParticle(TLorentzVector &vec, int pdgid)
{
  TLorentzVector *particle = new TLorentzVector(vec);
  pair<TLorentzVector*,int> ptcl;
  ptcl.first = particle;
  ptcl.second = pdgid;
  _particles.push_back(ptcl);
  // _particles.push_back(make_pair<TLorentzVector*,int>(particle, pdgid));
}

// __________________________________________________________________

void ttbarEventDisplay::AddParticle(TLorentzVector *vec, int pdgid)
{
  pair<TLorentzVector*,int> ptcl;
  ptcl.first = vec;
  ptcl.second = pdgid;
  _particles.push_back(ptcl);
  // _particles.push_back(make_pair<TLorentzVector*,int>(vec, pdgid));
}


// __________________________________________________________________

void ttbarEventDisplay::AddParticlePxPyPzE(double Px, double Py, double Pz, double E, int pdgid)
{
  
  TLorentzVector *vec = new TLorentzVector();
  vec -> SetPxPyPzE(Px, Py, Pz, E);
  pair<TLorentzVector*,int> ptcl;
  ptcl.first = vec;
  ptcl.second = pdgid;
  _particles.push_back(ptcl);
  // _particles.push_back(make_pair<TLorentzVector*,int>(vec, pdgid));
}


// __________________________________________________________________

double ttbarEventDisplay::GetMaxE()
{
  double maxpt = -1.;
  for (int i = 0; i < _particles.size(); ++i)
    if (_particles[i].first -> E() > maxpt)
      maxpt = _particles[i].first -> E();
  
  return maxpt;
}

// __________________________________________________________________

double ttbarEventDisplay::GetMaxP()
{
  double maxpt = -1.;
  for (int i = 0; i < _particles.size(); ++i)
    if (_particles[i].first -> P() > maxpt)
      maxpt = _particles[i].first -> P();
  
  return maxpt;
}

// __________________________________________________________________

double ttbarEventDisplay::GetMaxPt()
{
  double maxpt = -1.;
  for (int i = 0; i < _particles.size(); ++i)
    if (_particles[i].first -> Pt() > maxpt)
      maxpt = _particles[i].first -> Pt();
  
  return maxpt;
}

// __________________________________________________________________

void ttbarEventDisplay::PrintEps()
{
  Print("eps");
}


// __________________________________________________________________

void ttbarEventDisplay::Print(TString suffix)
{

  if (_doOneCavas) {
    _globalCan -> Print( Form("%s.%s", _globalCan -> GetName(), suffix.Data() ));
  } else {
    _cRZ -> cd();
    _cRZ -> Print( Form("%s.%s", _cRZ -> GetName(), suffix.Data() ));
    _cXY -> cd();
    _cXY -> Print( Form("%s.%s", _cXY -> GetName(), suffix.Data() ));
    _cLego -> cd();
    _cLego -> Print( Form("%s.%s", _cLego -> GetName(), suffix.Data() ));
  }
}

// __________________________________________________________________

void ttbarEventDisplay::DrawRZ()
{

  //  double sizeE = GetMaxPt();
  double sizeE = GetMaxP();
  double sizePt = GetMaxP();


  float x, y, x2, y2;
  float textsize = 0.05;

  if (_doOneCavas) 
    _globalCan -> cd(2);
  else
    _cRZ -> cd();

  //  gPad->Range(-SF*sizeE,-SF*sizePt,SF*sizeE,SF*sizePt);

  double SFx = 1.25;
  double SFy = 1.50;
  gPad->Range(-SFx*sizeE,-SFy*sizeE,SFx*sizeE,SFy*sizeE);

  // draw box
  _lines.push_back(DrawAxisLine(_x0 - sizeE, _y0 + sizePt, _x0 + sizeE, _y0 + sizePt));
  _lines.push_back(DrawAxisLine(_x0 - sizeE, _y0 - sizePt, _x0 + sizeE, _y0 - sizePt));
  _lines.push_back(DrawAxisLine(_x0 - sizeE, _y0 - sizePt, _x0 - sizeE, _y0 + sizePt));
  _lines.push_back(DrawAxisLine(_x0 + sizeE, _y0 - sizePt, _x0 + sizeE, _y0 + sizePt));
  _lines.push_back(DrawAxisLine(_x0 - sizeE, _y0, _x0 + sizeE, _y0));
  _lines.push_back(DrawAxisLine(_x0, _y0 - sizePt, _x0, _y0 + sizePt));

  TLatex *PtText = new TLatex(_x0 - 0.67*sizeE, _y0 - 1.13*sizePt,Form("E Scale: %3.1f GeV", sizePt));
  PtText->SetTextAlign(22);
  PtText->SetTextColor(1);
  PtText->SetTextSize(textsize);
  PtText->Draw();
  _latex.push_back(PtText);

  _latex.push_back(DrawLatex(_x0 - 0.40*sizeE, _y0 + 1.1*sizeE,Form("Run: %i, Event: %i", _runno, _evtno),
			     1, textsize));

  double phioff = 0.;
  //  double phioff = TMath::Pi();

  // draw particles:
  for (int i = 0; i < _particles.size(); ++i) {
    x = _x0;
    y = _y0;
    float signum = sin(_particles[i].first -> Phi()  - phioff) > 0 ? 1. : -1.;
    x2 = _particles[i].first -> Pz();
    y2 = _particles[i].first -> Pt()*signum;
    float SF = 0.10;
    float x2_long = (SF*sizeE + _particles[i].first -> Pz());
    float y2_long = (SF*sizeE + _particles[i].first -> Pt())*signum;

    _lines.push_back(DrawParticleLine(x, y, x2, y2, _particles[i].second));
    _latex.push_back(DrawLatex(x2_long, y2_long, 
			       GetPdgName(_particles[i].second),
			       GetPdgColor(_particles[i].second), textsize));

  }
}

// __________________________________________________________________

void ttbarEventDisplay::DrawXY()
{

  double size = GetMaxPt();
  float SF = 1.25;

  float x, y, x2, y2;
  float textsize = 0.05;

  if (_doOneCavas) 
    _globalCan -> cd(1);
  else
    _cXY -> cd();

  gPad -> Range(-SF*size,-SF*size,SF*size,SF*size);

  // draw axes:
  TEllipse *radius = new TEllipse(_x0, _y0, size, size);
  radius -> SetLineColor(1);
  radius -> SetLineWidth(1);
  radius -> SetLineStyle(2);
  radius -> Draw();
  _ellipses.push_back(radius);

  _lines.push_back(DrawAxisLine(_x0 - size, _y0, _x0 + size, _y0));
  _lines.push_back(DrawAxisLine(_x0, _y0 - size, _x0, _y0 + size));

  _latex.push_back(DrawLatex(_x0 - 0.67*size, _y0 - 1.13*size,Form("p_{T} Scale: %3.1f GeV", size),
			     1, textsize));
  _latex.push_back(DrawLatex(_x0 - 0.40*size, _y0 + 1.1*size,Form("Run: %i, Event: %i", _runno, _evtno),
			     1, textsize));

  double phioff = 0.;
  //  double phioff = - TMath::Pi() / 2.;

  float SFcap = 1.1;
  float axistextsize = 0.04;
  float testphi =  - phioff;
  _latex.push_back(DrawLatex(SFcap*size*cos(testphi),
			     SFcap*size*sin(testphi),
			     Form("x"), 1, axistextsize));
  //			     Form("#phi = 0"), 1, axistextsize));


  testphi = TMath::Pi() / 2. - phioff;
  _latex.push_back(DrawLatex(SFcap*size*cos(testphi),
			     SFcap*size*sin(testphi),
			     Form("y"), 1, axistextsize));
  //			     Form("#phi = 0"), 1, axistextsize));


  testphi = TMath::Pi() - phioff;
  //  _latex.push_back(DrawLatex(SFcap*size*cos(testphi),SFcap*size*sin(testphi),Form("#phi = #pi"), 1, axistextsize));
  testphi = TMath::Pi()/2 - phioff;
  //  _latex.push_back(DrawLatex(SFcap*size*cos(testphi),SFcap*size*sin(testphi),Form("#phi = #pi/2"), 1, axistextsize));
  testphi = -TMath::Pi()/2 - phioff;
  //  _latex.push_back(DrawLatex(SFcap*size*cos(testphi),SFcap*size*sin(testphi),Form("#phi = -#pi/2"), 1, axistextsize));
  

  // draw particles:
  for (int i = 0; i < _particles.size(); ++i) {
    x = _x0;
    y = _y0;
    float xphi = _particles[i].first -> Phi()- phioff;
    
    x2 = _particles[i].first -> Pt() * cos(xphi);
    y2 = _particles[i].first -> Pt() * sin(xphi);

    float SF = 0.10;
    float x2_long = (SF*size + _particles[i].first -> Pt()) * cos(xphi);
    float y2_long = (SF*size + _particles[i].first -> Pt()) * sin(xphi);

    _lines.push_back(DrawParticleLine(x, y, x2, y2, _particles[i].second));
    _latex.push_back(DrawLatex(x2_long, y2_long, 
			       GetPdgName(_particles[i].second),
			       GetPdgColor(_particles[i].second), textsize));

  }

}

// __________________________________________________________________

void ttbarEventDisplay::Translate3Dto2D(double x, double y, double z, 
					double &xx, double &yy)
{

  xx = _x0;
  yy = _y0;


  // phi:
  xx += x*sin(_phi) / (_phimax - _phimin) * sqrt( pow(_ymax - _ymin - _offy, 2) + pow(_xmax - _xmin - _offx, 2) );
  yy += x*cos(_phi) / (_phimax - _phimin) * sqrt( pow(_ymax - _ymin - _offy, 2) + pow(_xmax - _xmin - _offx, 2) );

  xx += y*cos(_phi)/(_etamax - _etamin) * (_xmax - _xmin - _offx);
  yy -= y*sin(_phi)/(_etamax - _etamin) * (_ymax - _ymin - _offy);


  // theta:
  xx += z*sin(_theta) / (_ptmax - _ptmin) * (_ymax - _ymin - _offy);
  yy += z*cos(_theta) / (_ptmax - _ptmin) * (_ymax - _ymin - _offy);

  xx += y*cos(_theta)/(_etamax - _etamin) * (_xmax - _xmin - _offx);
  yy -= y*sin(_theta)/(_etamax - _etamin) * (_xmax - _xmin - _offx);


}


// __________________________________________________________________

void ttbarEventDisplay::SetLegoPhiTheta(float phi, float theta)
{

  _phi = phi;
  _theta = theta;

}

// __________________________________________________________________

void ttbarEventDisplay::SetEtaMax(float etamax)
{
  _etamax = etamax;
  _etamin = -_etamax;
}
// __________________________________________________________________

void ttbarEventDisplay::DrawLego()
{

  //  gSystem -> Load("libGeom");

  _x0 = 0.25;
  _y0 = -0.60;

  int nbins = 100;
  //  float etamax = 4.7;

  //  _phimin = -TMath::Pi();
  //  _phimax = TMath::Pi();

  _phimin = 0.;
  _phimax = TMath::TwoPi();

  _ptmax = 1.2*GetMaxPt();
  _ptmin = 0.;
  //  _ptmin = -0.2*_ptmax;
  _xmin = -1.; _xmax = 1.; _ymin = -1; _ymax = 1.;

  if (_doOneCavas) 
    _globalCan -> cd(3);
  else
    _cLego->cd();
  gPad -> Range(_xmin, _ymin, _xmax, _ymax);

  float textsize = 0.05;
  _latex.push_back(DrawLatex(_xmin + 0.4*(_xmax - _xmin), 0.9*_ymax,Form("Run: %i, Event: %i", _runno, _evtno),
			     1, textsize));
  _latex.push_back(DrawLatex(_xmin + 0.75*(_xmax - _xmin), 0.9*_ymin,Form("p_{T} Scale: %3.1f GeV", GetMaxPt()),
			     1, textsize));
  textsize = 0.04;

  double SFy = 3.5;
  double SFx = 3.5;
  _xmin /= SFx;
  _xmax /= SFx;
  _ymin /= SFy;
  _ymax /= SFy;

  //  _offx = 0.2;
  //  _offy = 0.5;
  //  _offx = 0.3;
  _offy = 0.;
  _offx = 0.;

  /*
    float maxZ = GetMaxPt();
    _legoHist = new TH3D(TString(Form("LegoHist_%i_%i", _runno, _evtno)), 
    TString("Lego View"),
    nbins, 0., TMath::TwoPi(), 
    nbins, -etamax, etamax,
    nbins, 0., maxZ);
    _legoHist -> SetStats(0);
    _legoHist -> Draw();
  */

  //  double origin[] = {0., -etamax, 0.};
  //  TGeoBBox *box = new TGeoBBox("LegoFrame", TMath::TwoPi(), etamax, maxZ, origin);
  //  box -> Draw();

  //  int ip = 0;
  //  TGraph2D *graph = new TGraph2D();

  // draw axes:
    double x1, x2, y1, y2;
    Translate3Dto2D(_phimin, _etamin, 0., x1, y1);
    Translate3Dto2D(_phimax, _etamin, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));
    //    Translate3Dto2D(_phimin, _etamin, _ptmax, x2, y2);
    //    _lines.push_back(DrawAxisLine(x1, y1, x2, y2));
    Translate3Dto2D(_phimin, _etamax, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2,1 ));

    Translate3Dto2D(_phimax, _etamax, 0., x1, y1);
    Translate3Dto2D(_phimin, _etamax, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));
    Translate3Dto2D(_phimax, _etamin, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));
    // a
    Translate3Dto2D(_phimax, _etamax, _ptmax, x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));    

    //    Translate3Dto2D(_phimin, _etamax, _ptmax, x1, y1);
    //    Translate3Dto2D(_phimin, _etamin, _ptmax, x2, y2);
    //    _lines.push_back(DrawAxisLine(x1, y1, x2, y2));


    // a
    int Nsub = 10;
    int lstyle = 2;
    for (int k = 1; k <= Nsub; ++k) {
      Translate3Dto2D(_phimin, _etamax, k*_ptmax / Nsub, x1, y1);
      Translate3Dto2D(_phimax, _etamax, k*_ptmax / Nsub, x2, y2);
      if (k == Nsub)
	lstyle = 1;
      _lines.push_back(DrawAxisLine(x1, y1, x2, y2, lstyle));
    }
    Translate3Dto2D(_phimin, _etamax, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1)); 

    // a
    //    int Nsub = 10;
    lstyle = 2;
    for (int k = 1; k <= Nsub; ++k) {
      Translate3Dto2D(_phimax, _etamax, k*_ptmax / Nsub, x1, y1);
      Translate3Dto2D(_phimax, _etamin, k*_ptmax / Nsub, x2, y2);
      if (k == Nsub)
	lstyle = 1;
      _lines.push_back(DrawAxisLine(x1, y1, x2, y2, lstyle));
    }

    Translate3Dto2D(_phimax, _etamin, _ptmax, x1, y1);
    Translate3Dto2D(_phimax, _etamin, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));

    Translate3Dto2D(_phimax, _etamin, _ptmax, x1, y1);
    //    Translate3Dto2D(_phimin, _etamin, _ptmax, x2, y2);
    //    _lines.push_back(DrawAxisLine(x1, y1, x2, y2));

 
    // Cross:
    Translate3Dto2D(_phimin, 0., 0., x1, y1);
    Translate3Dto2D(_phimax, 0., 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));
    Translate3Dto2D(_phimin + 0.5*(_phimax - _phimin), _etamin, 0., x1, y1);
    Translate3Dto2D(_phimin + 0.5*(_phimax - _phimin), _etamax, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2, 1));

    // Subcross:
    /*
      float feta = _etamin + 0.25*(_etamax - _etamin)
      Translate3Dto2D(_phimin, feta, 0., x1, y1);
      Translate3Dto2D(_phimax, feta, 0., x2, y2);
      _lines.push_back(DrawAxisLine(x1, y1, x2, y2));
      feta = _etamin + 0.75*(_etamax - _etamin);
      Translate3Dto2D(_phimin, feta, 0., x1, y1);
      Translate3Dto2D(_phimax, feta, 0., x2, y2);
      _lines.push_back(DrawAxisLine(x1, y1, x2, y2));
    */

    for (int k = 1; k < 5; ++k) {
      float feta = 1.*k;
      Translate3Dto2D(_phimin, feta, 0., x1, y1);
      Translate3Dto2D(_phimax, feta, 0., x2, y2);
      _lines.push_back(DrawAxisLine(x1, y1, x2, y2));
      feta = -1.*k;
      Translate3Dto2D(_phimin, feta, 0., x1, y1);
      Translate3Dto2D(_phimax, feta, 0., x2, y2);
      _lines.push_back(DrawAxisLine(x1, y1, x2, y2));

    }

    Translate3Dto2D(_phimin + 0.25*(_phimax - _phimin), _etamin, 0., x1, y1);
    Translate3Dto2D(_phimin + 0.25*(_phimax - _phimin), _etamax, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2));

    Translate3Dto2D(_phimin + 0.75*(_phimax - _phimin), _etamin, 0., x1, y1);
    Translate3Dto2D(_phimin + 0.75*(_phimax - _phimin), _etamax, 0., x2, y2);
    _lines.push_back(DrawAxisLine(x1, y1, x2, y2));

    Translate3Dto2D(_phimin, _etamin, 0., x1, y1);
    if ( TMath::Abs(_phimin + TMath::Pi()) < 1.e-3)
      _latex.push_back(DrawLatex(x1 - 0.25*(_xmax - _xmin), y1, TString("#phi = -#pi"), 1, textsize));
    else     if ( TMath::Abs(_phimin) < 1.e-3)
      _latex.push_back(DrawLatex(x1 - 0.25*(_xmax - _xmin), y1, TString("#phi = 0"), 1, textsize));

    Translate3Dto2D(_phimax, _etamin, 0., x1, y1);
    if ( TMath::Abs(_phimax - TMath::Pi()) < 1.e-3)
      _latex.push_back(DrawLatex(x1 - 0.25*(_xmax - _xmin), y1, TString("#phi = #pi"), 1, textsize));
    else     if ( TMath::Abs(_phimax - TMath::TwoPi()) < 1.e-3)
      _latex.push_back(DrawLatex(x1 - 0.25*(_xmax - _xmin), y1, TString("#phi = 2#pi"), 1, textsize));
    // half:
    double val = 0.5*(_phimax+_phimin);
    Translate3Dto2D(val, _etamin, 0., x1, y1);
    if ( TMath::Abs(val - TMath::Pi()) < 1.e-3)
      _latex.push_back(DrawLatex(x1 - 0.25*(_xmax - _xmin), y1, TString("#phi = #pi"), 1, textsize));
    else     if ( TMath::Abs(val) < 1.e-3)
      _latex.push_back(DrawLatex(x1 - 0.25*(_xmax - _xmin), y1, TString("#phi = 0"), 1, textsize));

    Translate3Dto2D(_phimin, _etamin, 0., x1, y1);
    _latex.push_back(DrawLatex(x1, y1 - 0.2*(_ymax - _ymin), Form("#eta = %3.1f", _etamin), 1, textsize));
    Translate3Dto2D(_phimin, _etamax, 0., x1, y1);
    _latex.push_back(DrawLatex(x1, y1 - 0.2*(_ymax - _ymin), Form("#eta = %3.1f", _etamax), 1, textsize));
     Translate3Dto2D(_phimin, 0., 0., x1, y1);
    _latex.push_back(DrawLatex(x1 + 0.05*(_xmax - _xmin), y1 - 0.2*(_ymax - _ymin), "#eta = 0", 1, textsize));
   

    for (int i = 0; i < _particles.size(); ++i) {
    
    /*
      double delta = 0.01;
      TH2D *_ptcHist = new TH2D(Form("Particle%i_run%i_evt%i", i, _runno, _evtno),
      Form("Particle%i_run%i_evt%i", i, _runno, _evtno),
      1,_particles[i].first -> Phi() - delta, _particles[i].first -> Phi() + delta, 
      1,_particles[i].first -> Eta() - delta, _particles[i].first -> Eta() + delta
      );
      //			      nbins, 0., TMath::TwoPi(), 
      //			      nbins, -etamax, etamax);
      _ptcHist -> Sumw2();
      _ptcHist -> SetStats(0);
      _ptcHist -> Fill( _particles[i].first -> Phi(), 
      _particles[i].first -> Eta(),
      _particles[i].first -> Pt());
      int col = GetPdgColor(_particles[i].second);
      _ptcHist -> SetLineColor(col);    
      _ptcHist -> SetFillColor(col);
      _ptcHist -> SetMarkerColor(col);
      _ptcHist -> Draw("SameLego");
    */
    /*
      graph -> SetPoint(ip++,
      _particles[i].first -> Phi(), 
      _particles[i].first -> Eta(),
      _particles[i].first -> Pt());
    */

      double phioff = 0.;
      //  double phioff = TMath::Pi();

      double xphi = _particles[i].first -> Phi() - phioff;
      if (xphi > _phimax)
	xphi = _phimin + xphi - _phimax;
      else
	if (xphi < _phimin)
	  xphi = _phimax - (_phimin - xphi);
      
    Translate3Dto2D(xphi, _particles[i].first -> Eta(), 0.,
		    x1, y1);
    Translate3Dto2D(xphi, _particles[i].first -> Eta(), _particles[i].first -> Pt(),
		    x2, y2);

    _lines.push_back(DrawParticleLine(x1, y1, x2, y2, _particles[i].second));

  }
  //  graph -> Draw("P");
  //  _graphs2d.push_back(graph);
  
  

  _x0 = 0.;
  _y0 = 0.;
  

}

// __________________________________________________________________

void ttbarEventDisplay::DrawAll()
{
  DrawXY();
  DrawRZ();
  DrawLego();
}


// __________________________________________________________________

TCurlyLine* ttbarEventDisplay::DrawAxisLine(float x, float y, float x2, float y2, int lstyle, int lwidth, int lcol)
{
  TCurlyLine *l = new TCurlyLine(x,y,x2,y2);
  l->SetAmplitude(0.0);
  // l->SetWavy();
  l->SetLineWidth(lwidth);
  l->SetLineStyle(lstyle);
  l->SetLineColor(lcol);
  l->Draw();
  return l;
}

// __________________________________________________________________

Int_t ttbarEventDisplay::GetPdgColor( int id )
{
  if      (     id  == 21 ) return 14;  // g
  else if ( abs(id) ==  0 ) return 1; // ttbar
  else if ( abs(id) ==  6 ) return 2;   // t
  else if ( abs(id) == 24 ) return 108;   // W
  else if ( abs(id) == 11 ) return 102;   // e
  else if ( abs(id) == 12 ) return 102;   // e
  else if ( abs(id) == 22 ) return 102;   // gamma
  else if ( abs(id) == 13 ) return 103; // mu
  else if ( abs(id) == 14 ) return 103; // mu
  else if ( abs(id) == 15 ) return 105; // tau
  else if ( abs(id) == 16 ) return 105; // tau
  else if ( abs(id) ==  1 ) return 4;   // d
  else if ( abs(id) ==  2 ) return 104; // u
  else if ( abs(id) ==  3 ) return 6;   // s
  else if ( abs(id) ==  4 ) return 106; // c
  else if ( abs(id) ==  5 ) return 142; // b
  else 
    return 1;
}

TEllipse* ttbarEventDisplay::DrawDot( float x, float y, float r, int pConservation )
{
  TEllipse *c = new TEllipse(x,y,r,r);
  if ( pConservation == 1 )
    c->SetFillColor(103);
  else if ( pConservation == 0 )
    c->SetFillColor(2);
  else if ( pConservation == -1 )
    c->SetFillColor(6);
  c->Draw();

  return c;
}

TCurlyLine* ttbarEventDisplay::DrawParticleLine( float x, float y, float x2, float y2, int pdgid)
{
  TCurlyLine *l = new TCurlyLine(x,y,x2,y2);
  l->SetLineWidth(2);
  int abspdgid = TMath::Abs(pdgid);
  if ( abspdgid == 22 || abspdgid == 23 || abspdgid == 24) {
    l->SetAmplitude(0.005);
    l->SetWavy();
  }
  else if ( abspdgid == 21 ) {
    l->SetAmplitude(0.005);
    l->SetCurly();
    l->SetLineWidth(1);
  }
  else
    l->SetAmplitude(0);
  if ( abs(pdgid) == 5)
    //    l->SetLineWidth(4);
    l->SetLineWidth(2);
  if ( abspdgid == 12 || abspdgid == 14 || abspdgid == 16 )
    l->SetLineStyle(2);
  l->SetLineColor(GetPdgColor(pdgid));
  l->Draw();

  return l;
}

    
TLatex* ttbarEventDisplay::DrawLatex( float x, float y, TString textstring, int col, float tsize )
{
  TLatex *text = new TLatex(x,y,textstring);
  text->SetTextAlign(22);
  text->SetTextColor(col);
  text->SetTextSize(tsize);
  text->Draw();
  return text;
}

TString ttbarEventDisplay::GetPdgName( int id )
{

  // jiri:
  if ( id ==   0 ) return "t#bar{t}";

  // dag:
  if ( id ==   1 ) return "d";
  if ( id ==  -1 ) return "#bar{d}";
  if ( id ==   2 ) return "u";
  if ( id ==  -2 ) return "#bar{u}";
  if ( id ==   3 ) return "s";
  if ( id ==  -3 ) return "#bar{s}";
  if ( id ==   4 ) return "c";
  if ( id ==  -4 ) return "#bar{c}";
  if ( id ==   5 ) return "b";
  if ( id ==  -5 ) return "#bar{b}";
  if ( id ==   6 ) return "t";
  if ( id ==  -6 ) return "#bar{t}";
  if ( id ==  11 ) return "e^{-}";
  if ( id == -11 ) return "e^{+}";
  if ( id ==  12 ) return "#nu_{e}";
  if ( id == -12 ) return "#bar{#nu}_{e}";
  if ( id ==  13 ) return "#mu^{-}";
  if ( id == -13 ) return "#mu^{+}";
  if ( id ==  14 ) return "#nu_{#mu}";
  if ( id == -14 ) return "#bar{#nu}_{#mu}";
  if ( id ==  15 ) return "#tau^{-}";
  if ( id == -15 ) return "#tau^{+}";
  if ( id ==  16 ) return "#nu_{#tau}";
  if ( id == -16 ) return "#bar{#nu}_{#tau}";
  if ( id ==  21 ) return "g";
  if ( id ==  23 ) return "Z";
  if ( id ==  24 ) return "W^{+}";
  if ( id == -24 ) return "W^{-}";
  if ( id ==  22 ) return "#gamma";

  // mesons
  if ( id ==  111 ) return "#pi^{0}";
  if ( id ==  211 ) return "#pi^{+}";
  if ( id == -211 ) return "#pi^{-}";
  if ( id ==  113 ) return "#rho^{0}";
  if ( id ==  213 ) return "#rho^{+}";
  if ( id == -213 ) return "#rho^{-}";
  if ( id ==  130 ) return "K^{0}_{L}";
  if ( id ==  310 ) return "K^{0}_{S}";
  if ( abs(id) ==  311 ) return "K^{0}";
  if ( id ==  321 ) return "K^{+}";
  if ( id == -321 ) return "K^{-}";
  if ( abs(id) ==  313 ) return "K*^{0}";
  if ( id ==  323 ) return "K*^{+}";
  if ( id == -323 ) return "K*^{-}";
  if ( id == -511 ) return "#bar{B}^{0}";
  if ( id ==  511 ) return "B^{0}";
  if ( id ==  513 ) return "#bar{B}*^{0}";
  if ( id ==  513 ) return "B*^{0}";
  if ( id ==  521 ) return "B^{+}";
  if ( id == -521 ) return "B^{-}";
  if ( id ==  513 ) return "B*^{0}";
  if ( id ==  523 ) return "B*^{+}";
  if ( id == -523 ) return "B*^{-}";
  if ( abs(id) ==  533 ) return "B*^{0}_{s}";
  if ( id ==  221 ) return "#eta^{0}";
  if ( id ==  331 ) return "#eta'^{0}";
  if ( id ==  333 ) return "#phi^{0}";
  if ( id ==  223 ) return "#omega^{0}";
  // di-quarks
  if ( id ==  2101 ) return "ud";
  if ( id == -2101 ) return "#bar{ud}";
  if ( id ==  2103 ) return "ud";
  if ( id == -2103 ) return "#bar{ud}";
  if ( id ==  1103 ) return "dd";
  if ( id == -1103 ) return "#bar{dd}";
  if ( id ==  2203 ) return "uu";
  if ( id == -2203 ) return "#bar{uu}";
  // baryons
  if ( id ==  2212 ) return "p";
  if ( id == -2212 ) return "#bar{p}";
  if ( id ==  2112 ) return "n";
  if ( id == -2112 ) return "#bar{n}";
  if ( id ==  2224 ) return "#Delta^{++}";
  if ( id ==  2114 ) return "#Delta^{0}";
  if ( id ==  1114 ) return "#Delta^{-}";
  if ( id ==  2214 ) return "#Delta^{+}";
  if ( id == -2224 ) return "#bar{#Delta}^{--}";
  if ( id == -2114 ) return "#bar{#Delta}^{0}";
  if ( id == -1114 ) return "#bar{#Delta}^{+}";
  if ( id == -2214 ) return "#bar{#Delta}^{-}";
  if ( id ==  3122 ) return "#Lambda";
  if ( id == -3122 ) return "#bar{#Lambda}";
  if ( id ==  3212 ) return "#Sigma^{0}";
  if ( id == -3212 ) return "#bar{#Sigma}^{0}";
  if ( id ==  3222 ) return "#Sigma^{+}";
  if ( id == -3222 ) return "#bar{#Sigma}^{-}";
  if ( id ==  3112 ) return "#Sigma^{-}";
  if ( id == -3112 ) return "#bar{#Sigma}^{+}";
  if ( id ==  3214 ) return "#Sigma*^{0}";
  if ( id == -3214 ) return "#bar{#Sigma*}^{0}";
  if ( id ==  3224 ) return "#Sigma*^{+}";
  if ( id == -3224 ) return "#bar{#Sigma*}^{-}";
  if ( id ==  3114 ) return "#Sigma*^{-}";
  if ( id == -3114 ) return "#bar{#Sigma*}^{+}";

  char sbuf[256];
  sprintf(sbuf,"%d",id);
  //printf(" ------ %6d\n",id);
  return TString(sbuf);
}


// __________________________________________________________________
// __________________________________________________________________
// __________________________________________________________________
