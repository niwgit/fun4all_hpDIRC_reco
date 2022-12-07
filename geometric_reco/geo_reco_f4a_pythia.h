#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TArc.h"
#include "PrtLutNode.h"
#include "TVector3.h"
#include <set>
#include <map>
#include <vector>

TFile *fFile;
TTree *fTree;
TChain *fChain;

int fVerbose;
int nevents;
int fMethod;
TString fInputFile;
int fp1;
int fp2;

TF1 *fFit;
TSpectrum *fSpect;
TH1F *hthetac[5];
TH1F *hthetacd[5];
TH1F *hnph[5];
double fAngle[5];
TF1 *fFunc[5];
TH1F *fLnDiff[5];
double fCriticalAngle;
TString fCorrFile;

void FindPeak(double (&cherenkovreco)[5], double (&spr)[5]);
int FindPdg(double mom, double cangle);
void FitRing(double& x0, double& y0, double& theta);
//double FindStartTime(PrtEvent *e);
int fDetectorID;
double fBboxNum,fPipehAngle,fDphi,fBarPhi;
double fSigma[5];

TClonesArray *fLut[10];
void drawTheoryLines(double mom=6);
