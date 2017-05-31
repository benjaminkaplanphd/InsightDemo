#include <vector>
#include <map>

#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include "myDSCBPdf.h"

using namespace std;
using namespace RooFit;

TString vA[]={"3p7","3p8","4","4p5","5","6","7","8","9","10","11","12","15","20"};

struct fitResults_t{
  vector<double> a;
  map<TString,vector<double> > y,e;
};
fitResults_t myFits;

void PrepareWorkspace(RooWorkspace *wks);
void BuildModel(RooWorkspace *wks, TString model);
void BuildGenericModel(RooWorkspace *wks);
void BuildSimplifiedModel(RooWorkspace *wks);
void LoadData(RooWorkspace *wks);
void FitModel(RooWorkspace *wks, TString model);
void FitGenericModel(RooWorkspace *wks);
void FitSimplifiedModel(RooWorkspace *wks);
void PlotFitResults(fitResults_t fitres);
void PlotFitResults(RooWorkspace *wks);

   
void resModel( bool redoWks=true, TString model="simplified" ){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  //All data and mdoels will be stored in a RooWorkspace
  RooWorkspace *wks = 0;
  if( redoWks ){
    wks = new RooWorkspace("wks");
    PrepareWorkspace(wks);

    //The first step is to study the generic double-sided crystall ball in the simulated data
    BuildModel(wks,model);

    //Need to load the simulated data
    LoadData(wks);
    wks->Print();
    TFile * fileOut = new TFile("workspace_"+model+"Model.root","RECREATE");
    wks->Write();
    fileOut->Close();
    delete wks;
  }
  TFile *fileIn = new TFile("workspace_"+model+"Model.root");
  wks = (RooWorkspace*)fileIn->FindObjectAny("wks");  
  
  FitModel(wks,model);
  if( model=="generic" ) PlotFitResults(myFits);
  else PlotFitResults(wks);
}

void PrepareWorkspace( RooWorkspace *wks ){
  wks->factory("Mmumu[0,100]");  //The is the primary discriminant
  wks->factory("PassPreSelection[0,1]"); //Pre-defined variable to filter out data entries that might be interesting
  wks->factory("Passfilter[0,1]"); //Pre-defined variable to filter out data entries that might be interesting
  wks->factory("weight0[0,2]"); //Data is weighted
  return;
}

void BuildModel( RooWorkspace *wks, TString model ){

  if( model=="generic" ) BuildGenericModel(wks);
  else if( model=="simplified" ) BuildSimplifiedModel(wks);
  else{
    cout << "model " << model << " note recognized.  exiting..." << endl;
    exit(0);
  }

}

void BuildGenericModel( RooWorkspace *wks ){

  //Define variables needed for generic model.
  wks->factory("aL[1,0,10]");
  wks->factory("aH[1,0,10]");
  wks->factory("nL[5,0,50]");
  wks->factory("nH[5,0,50]");
  wks->factory("mu[5,0,100]");
  wks->factory("sigma[1,0,10]");
  
  RooArgSet params;
  params.add(*wks->arg("aL"));
  params.add(*wks->arg("aH"));
  params.add(*wks->arg("nL"));
  params.add(*wks->arg("nH"));
  params.add(*wks->arg("mu"));
  params.add(*wks->arg("sigma"));
  
  wks->defineSet("paramsGen",params);
    
  //Construct the double-sided Crystal Ball, add it to the workspace (creates a copy) and remove the local instance from memory
  myDSCBPdf * myPdf = new myDSCBPdf("dsCBGen","",*(RooAbsReal*)wks->var("Mmumu"),*(RooAbsReal*)wks->var("mu"),*(RooAbsReal*)wks->var("sigma"),*(RooAbsReal*)wks->var("aL"),*(RooAbsReal*)wks->var("aH"),*(RooAbsReal*)wks->var("nL"),*(RooAbsReal*)wks->var("nH"));
  wks->import(*myPdf);
  delete myPdf;

  return;
}

void BuildSimplifiedModel( RooWorkspace *wks ){

  //Define variables needed for generic model.
  wks->factory("alpha_CB[1,0,10]");
  wks->factory("nL[2.5]");
  wks->factory("nH[10]");
  wks->factory("a_mu[1,0,2]");
  wks->factory("a_sigma[0.02,0,1]");
  
  RooArgSet params;
  params.add(*wks->arg("alpha_CB"));
  params.add(*wks->arg("a_mu"));
  params.add(*wks->arg("a_sigma"));
  
  wks->defineSet("paramsSim",params);
  
  char c='a';  //Character counter
  TString catList="",pdfList="";
  int NA = sizeof(vA)/sizeof(TString);
  for( int iA=0;iA<NA;iA++ ){
    TString sA = vA[iA];
    TString fA = sA; fA.ReplaceAll("p",".");
    //Now, mu and sigma need to be defined separately for each mass hypothesis
    wks->factory("A"+sA+"["+fA+"]");
    wks->factory("expr::mu"+sA+"('a_mu*A"+sA+"',a_mu,A"+sA+")");
    wks->factory("expr::sigma"+sA+"('a_sigma*A"+sA+"',a_sigma,A"+sA+")");
    
    //Construct the double-sided Crystal Ball for each mass, add it to the workspace (creates a copy) and remove the local instance from memory
    myDSCBPdf * myPdf = new myDSCBPdf("dsCBSim"+sA,"",*(RooAbsReal*)wks->var("Mmumu"),*(RooAbsReal*)wks->function("mu"+sA),*(RooAbsReal*)wks->function("sigma"+sA),*(RooAbsReal*)wks->var("alpha_CB"),*(RooAbsReal*)wks->var("alpha_CB"),*(RooAbsReal*)wks->var("nL"),*(RooAbsReal*)wks->var("nH"));
    wks->import(*myPdf);
    delete myPdf;
    if( iA!=0 ){
      catList += ",";
      pdfList += ",";
    }
    catList += TString(c)+"="; catList += iA;
    pdfList += TString(c)+"=dsCBSim"+sA;
    c++;
  }
  wks->factory("catSim["+catList+"]");
  wks->factory("SIMUL::dsCBSim(catSim,"+pdfList+")");
  
  return;
}

void LoadData( RooWorkspace *wks ){

  int NA = sizeof(vA)/sizeof(TString);
  
  for( int iA=0;iA<NA;iA++ ){
    TString sA = vA[iA];
    //Open the file and retriese the data container
    TFile *fileS = new TFile("/local/data/kaplan/results/DiMuon/DiMuonAnalysis.MC12C.H125_aa"+sA+"_2mu2tauv3.root");
    TTree *treeS = (TTree*)fileS->FindObjectAny("data");

    RooArgSet args;
    args.add(*wks->arg("Mmumu"));
    args.add(*wks->arg("PassPreSelection"));
    args.add(*wks->arg("Passfilter"));
    args.add(*wks->arg("weight0"));

    //Create a RooDataSet from the data
    RooDataSet * data = new RooDataSet( "data"+sA, "data"+sA, args, Import(*treeS),WeightVar("weight0"),Cut("PassPreSelection==1&&Passfilter==1") );
    wks->import(*data);
    delete data;
  }
  return;
}

void FitModel( RooWorkspace *wks, TString model ){

  if( model=="generic" ) FitGenericModel(wks);
  else if( model=="simplified" ) FitSimplifiedModel(wks);
  else{
    cout << "model " << model << " note recognized.  exiting..." << endl;
    exit(0);
  }

}

void FitGenericModel( RooWorkspace *wks ){

  int NA = sizeof(vA)/sizeof(TString);
  
  for( int iA=0;iA<NA;iA++ ){
    TString sA = vA[iA];
    TString fA = sA; fA.ReplaceAll("p",".");
    double a = fA.Atof();
      
    RooAbsData *data = wks->data("data"+sA);
    wks->pdf("dsCBGen")->fitTo(*data,Range(0.5*a,1.5*a));

    myFits.a.push_back(a);
    TIterator *itr = wks->set("paramsGen")->createIterator();
    while( RooRealVar *var = (RooRealVar*)itr->Next() ){
      myFits.y[var->GetName()].push_back(var->getVal());
      myFits.e[var->GetName()].push_back(var->getError());
    }     
  }
}

void FitSimplifiedModel( RooWorkspace *wks ){

  int NA = sizeof(vA)/sizeof(TString);
  
  map<string,RooDataSet*> data_map;
  char c='a';
  for( int iA=0;iA<NA;iA++ ){
    TString sA = vA[iA];
      
    data_map[string(1,c)] = (RooDataSet*)wks->data("data"+sA);
    c++;
  }
  RooDataSet *dataSim = new RooDataSet("dataSim","",RooArgSet(*wks->arg("Mmumu"),*wks->arg("catSim")),Index(*wks->cat("catSim")),Import(data_map));
  RooFitResult *frSim = wks->pdf("dsCBSim")->fitTo(*dataSim,Save());
  frSim->SetName("frSim");
  wks->import(*frSim);
  delete frSim;
}

void PlotFitResults( fitResults_t fitres ){
  
  TGraphErrors *gNum=0;
  TGraphErrors *gDen=0;
  int iC=0;
  map<TString,vector<double> >::iterator itr = fitres.y.begin();
  for(;itr!=fitres.y.end();++itr){
    TString sC="";sC+=iC;
    TString name = (*itr).first;
    vector<double> x = fitres.a;
    vector<double> y = (*itr).second;
    vector<double> e = fitres.e[name];
    TGraphErrors *g = new TGraphErrors(x.size(),&x[0],&y[0],0,&e[0]);
    g->SetMarkerStyle(20);
    TCanvas *c = new TCanvas("c"+sC,"Fit: "+name,800,600);
    g->Draw("ap");
    g->GetXaxis()->SetTitle(name);
    cout << name << endl;
    if( name=="mu"||name=="sigma" ) g->Fit("pol1");
    else g->Fit("pol0");
    iC++;
    if( name=="aL" ) gNum = (TGraphErrors*)g->Clone("gNum");
    if( name=="aH" ) gDen = (TGraphErrors*)g->Clone("gDen");
  }
  if( gNum ){
    TString sC="";sC+=iC;
    for( int i=0;i<gNum->GetN();i++ ){
       double x(0),n(0),d(0),eN(0),eD(0);
       gNum->GetPoint(i,x,n);
       gDen->GetPoint(i,x,d);
       eN = gNum->GetErrorY(i);
       eD = gDen->GetErrorY(i);
       double y = n/d;
       double e = y*sqrt(pow(eN/n,2)+pow(eD/d,2));
       gNum->SetPoint(i,x,y);
       gNum->SetPointError(i,0,e);
    }
    TCanvas *c = new TCanvas("c"+sC,"Fit: aL/aH",800,600);
    gNum->Draw("ap");
    gNum->GetXaxis()->SetTitle("aL/aH");
    cout << "aL/aH" << endl;
    gNum->Fit("pol0");
  }
}

void PlotFitResults( RooWorkspace *wks ){


  int NA = sizeof(vA)/sizeof(TString);
  
  int iC=0;
  for( int iA=0;iA<NA;iA++ ){
    TString sA = vA[iA];
    TString fA = sA; fA.ReplaceAll("p",".");
    double a = fA.Atof();
    
    RooDataSet *data = (RooDataSet*)wks->data("data"+sA);
    RooAbsPdf  *pdf  = wks->pdf("dsCBSim"+sA);
    
    RooPlot *plot = wks->var("Mmumu")->frame(Title("Fit Projection, mA = "+sA),Range(0.5*a,1.5*a));
    data->plotOn(plot,DataError(RooAbsData::SumW2));
    pdf->plotOn(plot,VisualizeError(*(RooFitResult*)wks->genobj("frSim"),1),FillColor(kOrange+2));
    pdf->plotOn(plot);
    TString sC="";sC+=iC;
    TCanvas *c = new TCanvas("c"+sC,"",1200,400);
    c->Divide(2,1);
    c->cd(1);
    plot->Draw();
    c->cd(2);
    plot->Draw();
    plot->SetMinimum(0.01);
    c->cd(2)->SetLogy();
    iC++;
  }

}
