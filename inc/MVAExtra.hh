#include <algorithm>
#include <cstdlib>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <stdio.h>

using namespace std;



class MVAExtra
{

public:

  MVAExtra(string name, TH1* mvaS, TH1* mvaB);

  TH1D* GetROC(double rbin);

  double GetSeparation();
  
private:

  double GetEffForRoot(double theCut);
  double Root(double refValue);
  
  TH1* GetCumulativeDist(TH1* h);

  TH1 *fmvaS, *fmvaB;
  TH1 *fmvaScumul, *fmvaBcumul;
  
  //TSpline *fSplmvaCumS, *fSplmvaCumB;
  
  int fCutOrientation;
  int fMaxIter;
  int fNbins;
  
  double fAbsTol;
  double fXmax;
  double fXmin;

  string mvaeName;
 
};

MVAExtra::MVAExtra(string name, TH1* mvaS, TH1* mvaB):
  fMaxIter(100),
  fAbsTol(0.0)

{

  mvaeName=name;
  
  fmvaS =  mvaS; fmvaS->SetTitle("MVA Signal");
  fmvaB =  mvaB; fmvaB->SetTitle("MVA Backgr");
  fNbins = fmvaS->GetNbinsX();
  fXmax = fmvaS->GetXaxis()->GetXmax();
  fXmin = fmvaS->GetXaxis()->GetXmin(); 
  fCutOrientation = (fmvaS->GetMean() > fmvaB->GetMean()) ? +1 : -1;
}

double MVAExtra::GetEffForRoot(double theCut){
  
  double retVal=0;
 
  // retrieve the class object
  retVal = fmvaScumul->GetBinContent( fmvaScumul->FindBin( theCut ) );
    
  // caution: here we take some "forbidden" action to hide a problem:
  // in some cases, in particular for likelihood, the binned efficiency distributions
  // do not equal 1, at xmin, and 0 at xmax; of course, in principle we have the
  // unbinned information available in the trees, but the unbinned minimization is
  // too slow, and we don't need to do a precision measurement here. Hence, we force
  // this property.
  double eps = 1.0e-5;
  if      (theCut-fXmin < eps) retVal = (fCutOrientation > 0) ? 1.0 : 0.0;
  else if (fXmax-theCut < eps) retVal = (fCutOrientation > 0) ? 0.0 : 1.0;
 
 
  return retVal;
}

double MVAExtra::Root(double refValue){
  
  double a  = fXmin, b = fXmax;
  double fa = GetEffForRoot( a ) - refValue;
  double fb = GetEffForRoot( b ) - refValue;
  if (fb*fa > 0) {
    cout << "<Root> initial interval w/o root: "
	 << "(a=" << a << ", b=" << b << "),"
	 << " (Eff_a=" << GetEffForRoot( a ) 
	 << ", Eff_b=" << GetEffForRoot( b ) << "), "
	 << "(fa=" << fa << ", fb=" << fb << "), "
	 << "refValue = " << refValue << endl;
    return 1;
  }

  bool   ac_equal(kFALSE);
  double fc = fb;
  double c  = 0, d = 0, e = 0;
  for (int iter= 0; iter <= fMaxIter; iter++) {
    if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {

      // Rename a,b,c and adjust bounding interval d
      ac_equal = kTRUE;
      c  = a; fc = fa;
      d  = b - a; e  = b - a;
    }
  
    if (TMath::Abs(fc) < TMath::Abs(fb)) {
      ac_equal = kTRUE;
      a  = b;  b  = c;  c  = a;
      fa = fb; fb = fc; fc = fa;
    }

    double tol = 0.5 * 2.2204460492503131e-16 * TMath::Abs(b);
    double m   = 0.5 * (c - b);
    if (fb == 0 || TMath::Abs(m) <= tol || TMath::Abs(fb) < fAbsTol) return b;
  
    // Bounds decreasing too slowly: use bisection
    if (TMath::Abs (e) < tol || TMath::Abs (fa) <= TMath::Abs (fb)) { d = m; e = m; }      
    else {
      // Attempt inverse cubic interpolation
      double p, q, r;
      double s = fb / fa;
      
      if (ac_equal) { p = 2 * m * s; q = 1 - s; }
      else {
	q = fa / fc; r = fb / fc;
	p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
	q = (q - 1) * (r - 1) * (s - 1);
      }
      // Check whether we are in bounds
      if (p > 0) q = -q;
      else       p = -p;
      
      double min1 = 3 * m * q - TMath::Abs (tol * q);
      double min2 = TMath::Abs (e * q);
      if (2 * p < (min1 < min2 ? min1 : min2)) {
	// Accept the interpolation
	e = d;        d = p / q;
      }
      else { d = m; e = m; } // Interpolation failed: use bisection.
    }
    // Move last best guess to a
    a  = b; fa = fb;
    // Evaluate new trial root
    if (TMath::Abs(d) > tol) b += d;
    else                     b += (m > 0 ? +tol : -tol);

    fb = GetEffForRoot( b ) - refValue;

  }

  // Return our best guess if we run out of iterations
  cout << "<Root> maximum iterations (" << fMaxIter 
       << ") reached before convergence" << endl;

  return b;
}

TH1* MVAExtra::GetCumulativeDist(TH1* h){
  
  // get the cumulative distribution of a histogram
  TH1* cumulativeDist= (TH1*) h->Clone(Form("%sCumul",h->GetTitle()));
  //cumulativeDist->Smooth(5); // with this, I get less beautiful ROC curves, hence out!
 
  float partialSum = 0;
  float inverseSum = 0.;
 
  float val;
  for (int ibinEnd=1, ibin=cumulativeDist->GetNbinsX(); ibin >=ibinEnd ; ibin--){
    val = cumulativeDist->GetBinContent(ibin);
    if (val>0) inverseSum += val;
  }
  inverseSum = 1/inverseSum; // as I learned multiplications are much faster than division, and later I need one per bin. Well, not that it would really matter here I guess :)
 
  for (int ibinEnd=1, ibin=cumulativeDist->GetNbinsX(); ibin >=ibinEnd ; ibin--){
    val = cumulativeDist->GetBinContent(ibin);
    if (val>0) partialSum += val;
    cumulativeDist->SetBinContent(ibin,partialSum*inverseSum);
  }
  return cumulativeDist;
}
    
TH1D* MVAExtra::GetROC(double rbin=2){

  int fCutOrientation = (fmvaS->GetMean() > fmvaB->GetMean()) ? +1 : -1;
  
  fmvaScumul = GetCumulativeDist(fmvaS);
  fmvaBcumul = GetCumulativeDist(fmvaB);
  fmvaScumul->Scale(1.0/TMath::Max(numeric_limits<double>::epsilon(),fmvaScumul->GetMaximum()));
  fmvaBcumul->Scale(1.0/TMath::Max(numeric_limits<double>::epsilon(),fmvaBcumul->GetMaximum()));
  fmvaScumul->SetMinimum(0);
  fmvaBcumul->SetMinimum(0);
  //   fmvaScumul->Draw("hist");
  //   fmvaBcumul->Draw("histsame");

  // background rejection (=1-eff.) versus signal efficiency
  TH1D* rejBvsS = new TH1D(mvaeName.c_str(), "", fNbins, 0, 1);
  rejBvsS->SetXTitle("Signal eff");
  rejBvsS->SetYTitle("Backgr rejection (1-eff)");

  //fSplmvaCumS  = new TSpline1( "spline2_signal",     new TGraph( fmvaScumul ) );
  //fSplmvaCumB  = new TSpline1( "spline2_background", new TGraph( fmvaBcumul ) );
  // verify spline sanity
  //gTools().CheckSplines( fmvaScumul, fSplmvaCumS );
  //gTools().CheckSplines( fmvaBcumul, fSplmvaCumB );

  double effB = 0;
  for (int bini=1; bini<=fNbins; bini++) {

    // find cut value corresponding to a given signal efficiency
    double effS = rejBvsS->GetBinCenter( bini );
    double cut  = Root(effS);

    // retrieve background efficiency for given cut
    effB = fmvaBcumul->GetBinContent(fmvaBcumul->FindBin(cut));
    //effB = fSplmvaCumB->Eval(cut);
    
    // and fill histogram
    rejBvsS->SetBinContent(bini, 1.0-effB );
  }

  rejBvsS->Rebin(rbin); rejBvsS->Scale(1./rbin);
  
  return rejBvsS;
}


double MVAExtra::GetSeparation(){

  double s,b,sum=0;
  double dx = fXmax-(fXmin)/(double)fNbins;
  for (int i=0;i<fNbins;i++){
    s=fmvaS->GetBinContent(i+1)/(fmvaS->GetSumOfWeights()*dx);
    b=fmvaB->GetBinContent(i+1)/(fmvaB->GetSumOfWeights()*dx);
    if (s+b>0)
      sum+= (0.5*(s-b)*(s-b)/(s+b))*dx;
  }

  return sum;

}
