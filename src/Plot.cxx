
#include "TApplication.h"
#include "ttbarEventDisplay.hpp"
//#include "TFile.h"

// jk 24.8.2018

int main()
{

  // TFile *outfile = new TFile("outfile.root", "recreate");
  
  ttbarEventDisplay *ttevt = new ttbarEventDisplay(0, 0, true);

  ttevt -> AddParticlePxPyPzE(10., 20., 100., 500., 1);
  ttevt -> DrawRZ();
  ttevt -> DrawXY();
  ttevt -> DrawLego();

  ttevt -> Print("png");
  ttevt -> Print("C");
  // ttevt -> Print("pdf");
  
  // gApplication -> Run();

  return 0;
  
}
