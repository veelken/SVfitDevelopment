
#include <TFile.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TArrayF.h>
#include <TF1.h>

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <assert.h>

enum { kNormByQuantiles, kNormByNegLogMax, kNormByValue };

void showHistogram2d(TH2* histogram, 
		     const std::string& xAxisTitle, const std::string& yAxisTitle,
		     int zAxisNormOption, double zMin, double zMax, 
		     Float_t* genX, Float_t* genY, 
		     const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetTopMargin(0.10);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.14);
  canvas->SetBottomMargin(0.12);

  histogram->SetTitle("");
  histogram->SetStats(false);
  int numBinsX = histogram->GetNbinsX();
  int numBinsY = histogram->GetNbinsY();
  if ( zAxisNormOption == kNormByQuantiles ) {
    std::vector<double> binContents;
    for ( int iBinX = 1; iBinX <= numBinsX; ++iBinX ) {
      for ( int iBinY = 1; iBinY <= numBinsY; ++iBinY ) {
	binContents.push_back(histogram->GetBinContent(iBinX, iBinY));
      }
    }
    std::sort(binContents.begin(), binContents.end());
    histogram->SetMinimum(binContents[TMath::Nint(0.05*binContents.size())]);
    histogram->SetMaximum(binContents[TMath::Nint(0.95*binContents.size())]);
  } else if ( zAxisNormOption == kNormByNegLogMax ) {
    double maxBinContent = 0.;
    for ( int iBinX = 1; iBinX <= numBinsX; ++iBinX ) {
      for ( int iBinY = 1; iBinY <= numBinsY; ++iBinY ) {
	double binContent = histogram->GetBinContent(iBinX, iBinY);
	if ( binContent > maxBinContent ) {
	  std::cout << "binX = " << iBinX << " (x = " << histogram->GetXaxis()->GetBinCenter(iBinX) << ")," 
		    << " binY = " << iBinY << " (y = " << histogram->GetYaxis()->GetBinCenter(iBinY) << "): maxBinContent = " << maxBinContent << std::endl;
	  maxBinContent = binContent;
	}
      }
    }
    double logMaxBinContent = TMath::Log(maxBinContent);
    for ( int iBinX = 1; iBinX <= numBinsX; ++iBinX ) {
      for ( int iBinY = 1; iBinY <= numBinsY; ++iBinY ) {
	double binContent = histogram->GetBinContent(iBinX, iBinY);
	if ( binContent > 0. ) {
	  histogram->SetBinContent(iBinX, iBinY, -TMath::Log(binContent) + logMaxBinContent);
	} else {
	  histogram->SetBinContent(iBinX, iBinY, -1.);
	}
      }
    }
    histogram->SetMinimum(0.);
    histogram->SetMaximum(zMax);
  } else if ( zAxisNormOption == kNormByValue ) {
    histogram->SetMinimum(zMin);
    histogram->SetMaximum(zMax);
  } else assert(0);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(1.15);

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(1.30);

  gStyle->SetPalette(1,0);
  histogram->Draw("COLZ");

  TMarker* genMarker = 0;
  if ( genX && genY ) {
    genMarker = new TMarker(*genX, *genY, 34);
    genMarker->SetMarkerColor(1);
    genMarker->SetMarkerSize(2);
    genMarker->Draw();
  }

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete genMarker;
  delete canvas;  
}

void showHistogram1d(TH1* histogram, 
		     const std::string& xAxisTitle,
		     Float_t* genX, 
		     const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetTopMargin(0.10);
  canvas->SetLeftMargin(0.16);
  canvas->SetRightMargin(0.14);
  canvas->SetBottomMargin(0.12);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(1.15);

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle("Sampling Points");
  yAxis->SetTitleOffset(1.60);

  histogram->SetLineColor(1);
  histogram->SetLineWidth(2);
  histogram->SetMarkerColor(1);
  histogram->SetMarkerStyle(20);
  histogram->Draw("e1p");

  TMarker* genMarker = 0;
  if ( genX ) {
    genMarker = new TMarker(*genX, 0.10, 34);
    genMarker->SetMarkerColor(1);
    genMarker->SetMarkerSize(2);
    genMarker->Draw();
  }

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete genMarker;
  delete canvas;  
}

void showGraph(TGraph* graph, 
	       const std::string& xAxisTitle,
	       Float_t* genX, 
	       double yMin, double yMax, const std::string& yAxisTitle,
	       const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetTopMargin(0.10);
  canvas->SetLeftMargin(0.16);
  canvas->SetRightMargin(0.14);
  canvas->SetBottomMargin(0.12);

  int numPoints = graph->GetN();
  double xMin, xMax, yDummy;
  graph->GetPoint(0, xMin, yDummy);
  graph->GetPoint(numPoints - 1, xMax, yDummy);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numPoints/100, xMin, xMax);
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(1.15);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(1.15);

  dummyHistogram->Draw("axis");

  TGraph* genMarker = 0;
  if ( genX ) {
    genMarker = new TGraph(2);
    genMarker->SetPoint(0, xMin, *genX);
    genMarker->SetPoint(1, xMax, *genX);
    genMarker->SetLineColor(8);
    genMarker->SetLineWidth(1);
    genMarker->SetMarkerColor(8);
    genMarker->SetMarkerStyle(20);
    genMarker->SetMarkerSize(1);
    genMarker->Draw("L");
  }

  graph->SetLineColor(1);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(1);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->Draw("L");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete dummyHistogram;
  delete genMarker;
  delete canvas;  
}

int convertToInt(const TObjArray& run_ls_event, int idx)
{
  TObjString* entry = dynamic_cast<TObjString*>(run_ls_event.At(idx));
  assert(entry);
  return atoi(entry->GetString().Data());
}

void getGen(const std::string& genNtupleFileName, int run, int ls, int event, Float_t& genX1, Float_t& genX2, Float_t& genPt, Float_t& genEta, Float_t& genPhi, Float_t& genMass)
{
  TFile* genNtupleFile = new TFile(genNtupleFileName.data());
  std::string genNtupleName = "neuralMtautauNtupleProducer/neuralMtautauNtuple";
  TTree* genNtuple = dynamic_cast<TTree*>(genNtupleFile->Get(genNtupleName.data()));
  if ( !genNtuple ) {
    std::cerr << "Failed to find tree with name = " << genNtupleName << " in file = " << genNtupleFileName << " --> aborting !!" << std::endl;
    assert(0);
  }

  Float_t genNtuple_X1, genNtuple_X2, genNtuple_pt, genNtuple_eta, genNtuple_phi, genNtuple_mass;
  genNtuple->SetBranchAddress("genX1", &genNtuple_X1);
  genNtuple->SetBranchAddress("genX2", &genNtuple_X2);
  genNtuple->SetBranchAddress("genDiTauPt", &genNtuple_pt);
  genNtuple->SetBranchAddress("genDiTauEta", &genNtuple_eta);
  genNtuple->SetBranchAddress("genDiTauPhi", &genNtuple_phi);
  genNtuple->SetBranchAddress("genMtautau", &genNtuple_mass);
  
  ULong64_t genNtuple_run, genNtuple_ls, genNtuple_event;
  genNtuple->SetBranchAddress("run", &genNtuple_run);
  genNtuple->SetBranchAddress("lumisection", &genNtuple_ls);
  genNtuple->SetBranchAddress("event", &genNtuple_event);
    
  bool isFound = false;

  int numEvents = genNtuple->GetEntries();
  std::cout << "numEvents = " << numEvents << std::endl;
  for ( int iEvent = 0; iEvent < numEvents; ++iEvent ) {
    genNtuple->GetEvent(iEvent);
    std::cout << "run " << genNtuple_run << ", ls = " << genNtuple_ls << ", event " << genNtuple_event << std::endl;
    if ( genNtuple_run == run && genNtuple_ls == ls && genNtuple_event == event ) {
      genX1   = genNtuple_X1;
      genX2   = genNtuple_X2;
      genPt   = genNtuple_pt;
      genEta  = genNtuple_eta;
      genPhi  = genNtuple_phi;
      genMass = genNtuple_mass;
      isFound = true;
    }
  }

  if ( !isFound ) {
    std::cerr << "Failed to find run# " << run << ", #ls " << ls << ", #event " << event << " in file = " << genNtupleFileName << " --> aborting !!" << std::endl;
    assert(0);
  }

  delete genNtupleFile;
}

TH1* bookMassHistogram(const std::string& histogramName)
{
  double xMin = 1.e+1;
  double xMax = 1.e+3;
  double logBinWidth = 1.025;
  int numBins = TMath::Log(xMax/xMin)/TMath::Log(logBinWidth);
  TArrayF binning(numBins + 1);
  double x = xMin;  
  for ( int iBin = 0; iBin <= numBins; ++iBin ) {
    binning[iBin] = x;
    //std::cout << "binning[" << iBin << "] = " << binning[iBin] << std::endl;
    x *= logBinWidth;
  }  
  TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, binning.GetArray());
  return histogram;
}

double square(double x)
{
  return x*x;
}

void divideByBinWidth(TH1* histogram)
{
  int numBins = histogram->GetNbinsX();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double binError = histogram->GetBinError(iBin);
    double binWidth = histogram->GetBinWidth(iBin);
    histogram->SetBinContent(iBin, binContent/binWidth);
    histogram->SetBinError(iBin, binError/binWidth);
  }
}

void printMaximum(const TH1* histogram)
{
  double xMax = -1.;
  double yMax = -1.;
  int numBins = histogram->GetNbinsX();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double binCenter = histogram->GetBinCenter(iBin);
    if ( binContent > yMax ) {
      xMax = binCenter;
      yMax = binContent;
    }
  }
  std::cout << "yMax = " << yMax << " @ xMax = " << xMax << std::endl;
}

double compDeltaFunctionJacobiFactor(double X)
{
  double jacobiFactor = 1./(X*X*X);
  return jacobiFactor;
}

void plotMarkovChainMonitorTree()
{
  std::string inputFilePath = "/data1/veelken/tmp/svFitStudies/";

  std::vector<std::string> eventsToProcess;
  eventsToProcess.push_back("1:102093:40804761");
  //eventsToProcess.push_back("1:3369:1346645");
  //eventsToProcess.push_back("1:89621:35819785");
  //eventsToProcess.push_back("1:79157:31638382");

  std::string genNtupleFileName = "/data1/veelken/tmp/svFitStudies/testNSVfitTrackLikelihoods3_ntuple_Higgs_mutau_125_2012Nov30.root";

  gROOT->SetBatch(true);

  for ( std::vector<std::string>::const_iterator eventToProcess = eventsToProcess.begin();
	eventToProcess != eventsToProcess.end(); ++eventToProcess ) {
    std::cout << "processing " << (*eventToProcess) << std::endl;
    TObjArray* run_ls_event = TString(eventToProcess->data()).Tokenize(":");
    assert(run_ls_event);
    if ( run_ls_event->GetEntries() != 3 ) {
      std::cerr << "Failed to parse " << (*eventToProcess) << " !!" << std::endl;
      assert(0);
    }

    int run = convertToInt(*run_ls_event, 0);
    int ls = convertToInt(*run_ls_event, 1);
    int event = convertToInt(*run_ls_event, 2);

    std::string inputFileName = Form(
      "%s/nSVfitProducerByIntegration2W%sMaxGenVertex_run%i_ls%i_ev%i_mc.root", 
      inputFilePath.data(), "decayKinePlusMEt", run, ls, event);
    TFile* inputFile = new TFile(inputFileName.data());
    assert(inputFile);

    std::string treeName = "monitorTree";
    TTree* monitorTree = dynamic_cast<TTree*>(inputFile->Get(treeName.data()));
    assert(monitorTree);

    int numEntries = monitorTree->GetEntries();

    TH2* histogramDensity_vs_X1_and_X2           = new TH2D("histogramDensity_vs_X1_and_X2",           "histogramDensity_vs_X1_and_X2",            100, 0., 1., 100, 0., 1.);
    TH2* histogramMass_vs_X1_and_X2              = new TH2D("histogramMass_vs_X1_and_X2",              "histogramMass_vs_X1_and_X2",               100, 0., 1., 100, 0., 1.);

    TH2* histogramDensity_vs_decayDist1_and_phi1 = new TH2D("histogramDensity_vs_decayDist1_and_phi1", "histogramDensity_vs_decayDist1_and_phi1",  200, -0.2, +0.2, 180, -TMath::Pi(), +TMath::Pi());
    TH2* histogramDensity_vs_X1_and_phi1         = new TH2D("histogramDensity_vs_X1_and_phi1",         "histogramDensity_vs_X1_and_phi1",          100,  0.,   1.,  180, -TMath::Pi(), +TMath::Pi());
    TH2* histogramDensity_vs_decayDist2_and_phi2 = new TH2D("histogramDensity_vs_decayDist2_and_phi2", "histogramDensity_vs_decayDist2_and_phi2",  200, -0.2, +0.2, 180, -TMath::Pi(), +TMath::Pi());
    TH2* histogramDensity_vs_X2_and_phi2         = new TH2D("histogramDensity_vs_X2_and_phi2",         "histogramDensity_vs_X2_and_phi2",          100,  0.,   1.,  180, -TMath::Pi(), +TMath::Pi());
        
    TH1* histogramDiTauPtDistribution            = new TH1D("histogramDiTauPtDistribution",            "histogramDiTauPtDistribution",             250, 0., 250.);
    TH1* histogramDiTauEtaDistribution           = new TH1D("histogramDiTauEtaDistribution",           "histogramDiTauEtaDistribution",            198, -9.9, +9.9);
    TH1* histogramDiTauPhiDistribution           = new TH1D("histogramDiTauPhiDistribution",           "histogramDiTauPhiDistribution",            360, -TMath::Pi(), +TMath::Pi());
    TH1* histogramDiTauMassDistribution          = new TH1D("histogramDiTauMassDistribution",          "histogramDiTauMassDistribution",           250, 0., 250.);
    
    TH1* histogramSVfitMass_unweighted                      = bookMassHistogram("histogramSVfitMass_unweighted");
    TH1* histogramSVfitMass_gjAngle_labGt0_00025_unweighted = bookMassHistogram("histogramSVfitMass_gjAngle_labGt0_00025_unweighted");
    TH1* histogramSVfitMass_gjAngle_labLt0_00025_unweighted = bookMassHistogram("histogramSVfitMass_gjAngle_labLt0_00025_unweighted");
    TH1* histogramSVfitMass_weighted                        = bookMassHistogram("histogramSVfitMass_weighted");
    TH1* histogramSVfitMass_gjAngle_labGt0_00025_weighted   = bookMassHistogram("histogramSVfitMass_gjAngle_labGt0_00025_weighted");
    TH1* histogramSVfitMass_gjAngle_labLt0_00025_weighted   = bookMassHistogram("histogramSVfitMass_gjAngle_labLt0_00025_weighted");

    TH1* histogramLeg1PtDistribution             = new TH1D("histogramLeg1PtDistribution",             "histogramLeg1PtDistribution",              250, 0., 250.);
    TH1* histogramLeg1EtaDistribution            = new TH1D("histogramLeg1EtaDistribution",            "histogramLeg1EtaDistribution",             198, -9.9, +9.9);
    TH1* histogramLeg1PhiDistribution            = new TH1D("histogramLeg1PhiDistribution",            "histogramLeg1PhiDistribution",             360, -TMath::Pi(), +TMath::Pi());
    TH1* histogramLeg1dEtaDistribution           = new TH1D("histogramLeg1dEtaDistribution",           "histogramLeg1dPhiDistribution",            400, -0.05, +0.05);
    TH1* histogramLeg1dPhiDistribution           = new TH1D("histogramLeg1dPhiDistribution",           "histogramLeg1dEtaDistribution",            400, -0.05, +0.05);
    TH2* histogramLeg1dPhi_vs_dEtaDistribution   = new TH2D("histogramLeg1dPhi_vs_dEtaDistribution",   "histogramLeg1dPhi_vs_dEtaDistribution",    400, -0.05, +0.05, 400, -0.05, +0.05);

    TH1* histogramLeg2PtDistribution             = new TH1D("histogramLeg2PtDistribution",             "histogramLeg1PtDistribution",              250, 0., 250.);
    TH1* histogramLeg2EtaDistribution            = new TH1D("histogramLeg2EtaDistribution",            "histogramLeg2EtaDistribution",             198, -9.9, +9.9);
    TH1* histogramLeg2PhiDistribution            = new TH1D("histogramLeg2PhiDistribution",            "histogramLeg2PhiDistribution",             360, -TMath::Pi(), +TMath::Pi());
    TH1* histogramLeg2dEtaDistribution           = new TH1D("histogramLeg2dEtaDistribution",           "histogramLeg2dPhiDistribution",            400, -0.05, +0.05);
    TH1* histogramLeg2dPhiDistribution           = new TH1D("histogramLeg2dPhiDistribution",           "histogramLeg2dEtaDistribution",            400, -0.05, +0.05);
    TH2* histogramLeg2dPhi_vs_dEtaDistribution   = new TH2D("histogramLeg2dPhi_vs_dEtaDistribution",   "histogramLeg2dPhi_vs_dEtaDistribution",    400, -0.05, +0.05, 400, -0.05, +0.05);

    TH1* histogramMoveDistribution               = new TH1D("histogramMoveDistribution",               "histogramMoveDistribution",               1000, 0., monitorTree->GetEntries());

    TGraph* graphEvolution_X1                    = new TGraph(numEntries);
    TGraph* graphEvolution_phiLab1               = new TGraph(numEntries);
    TGraph* graphEvolution_decayDist1            = new TGraph(numEntries);
    TGraph* graphEvolution_X2                    = new TGraph(numEntries);
    TGraph* graphEvolution_phiLab2               = new TGraph(numEntries);
    TGraph* graphEvolution_decayDist2            = new TGraph(numEntries);

    Float_t X1, phiLab1, decayDistance1, X2, phiLab2, decayDistance2;
    monitorTree->SetBranchAddress("x0", &X1);
    monitorTree->SetBranchAddress("x1", &phiLab1);
    //monitorTree->SetBranchAddress("x2", &decayDistance1);
    //monitorTree->SetBranchAddress("x3", &X2);
    //monitorTree->SetBranchAddress("x4", &phiLab2);
    //monitorTree->SetBranchAddress("x5", &decayDistance2);
    //monitorTree->SetBranchAddress("x0", &X1);
    //monitorTree->SetBranchAddress("x1", &phiLab1);
    monitorTree->SetBranchAddress("x2", &X2);
    monitorTree->SetBranchAddress("x3", &phiLab2);
    //monitorTree->SetBranchAddress("x1", &X2);
    //monitorTree->SetBranchAddress("x0", &phiLab1);
    //monitorTree->SetBranchAddress("x1", &phiLab2);

    Float_t diTauPt, diTauEta, diTauPhi, diTauMass;
    monitorTree->SetBranchAddress("APt", &diTauPt);
    monitorTree->SetBranchAddress("AEta", &diTauEta);
    monitorTree->SetBranchAddress("APhi", &diTauPhi);
    monitorTree->SetBranchAddress("AMass", &diTauMass);

    Float_t leg1Pt, leg1Eta, leg1Phi, leg1VisPt, leg1VisEta, leg1VisPhi;
    monitorTree->SetBranchAddress("A_leg1Pt", &leg1Pt);
    monitorTree->SetBranchAddress("A_leg1Eta", &leg1Eta);
    monitorTree->SetBranchAddress("A_leg1Phi", &leg1Phi);
    monitorTree->SetBranchAddress("A_leg1VisPt", &leg1VisPt);
    monitorTree->SetBranchAddress("A_leg1VisEta", &leg1VisEta);
    monitorTree->SetBranchAddress("A_leg1VisPhi", &leg1VisPhi);

    Float_t leg2Pt, leg2Eta, leg2Phi, leg2VisPt, leg2VisEta, leg2VisPhi;
    monitorTree->SetBranchAddress("A_leg2Pt", &leg2Pt);
    monitorTree->SetBranchAddress("A_leg2Eta", &leg2Eta);
    monitorTree->SetBranchAddress("A_leg2Phi", &leg2Phi);
    monitorTree->SetBranchAddress("A_leg2VisPt", &leg2VisPt);
    monitorTree->SetBranchAddress("A_leg2VisEta", &leg2VisEta);
    monitorTree->SetBranchAddress("A_leg2VisPhi", &leg2VisPhi);
    
    Int_t move;
    monitorTree->SetBranchAddress("move", &move);

    TF1* pdfWeight_vs_eta_m90 = new TF1("pdfWeight_vs_eta", "[0] + x*([1] + x*([2] + x*([3] + x*([4] + x*([5]  + x*[6])))))", -2., +2.); 
    pdfWeight_vs_eta_m90->SetParameter(0, 0.173744);
    pdfWeight_vs_eta_m90->SetParameter(1, 1.00115e-18);
    pdfWeight_vs_eta_m90->SetParameter(2, -0.0370055);
    pdfWeight_vs_eta_m90->SetParameter(3, 3.40875e-17);
    pdfWeight_vs_eta_m90->SetParameter(4, -0.014232);
    pdfWeight_vs_eta_m90->SetParameter(5, -1.07345e-17);
    pdfWeight_vs_eta_m90->SetParameter(6, 0.00320732);

    for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
      monitorTree->GetEntry(iEntry);
      
      if ( (iEntry % 100000) == 0 ) 
      //if ( TMath::Abs(leg1Eta - leg1VisEta) < 0.000250 && TMath::Abs(leg1Phi - leg1VisPhi) < 0.000250 )
	std::cout << "entry #" << iEntry << ":" 
		  << " X1 = " << X1 << ", phiLab1 = " << phiLab1 << ", decayDistance1 = " << decayDistance1 << "," 
		  << " X2 = " << X2 << ", phiLab2 = " << phiLab2 << ", decayDistance2 = " << decayDistance2 << ", Mass = " << diTauMass << std::endl;
      
      //if ( !(leg2Eta > 0.0485 && leg2Eta < 0.0495) ) continue;
      //if ( !(phiLab1 > 1. && phiLab1 < 2.) ) continue;
      //if ( !(phiLab1 > -2. && phiLab1 < -1.) ) continue;

      histogramDensity_vs_X1_and_X2->Fill(X1, X2);
      histogramMass_vs_X1_and_X2->Fill(X1, X2, diTauMass);

      histogramDensity_vs_decayDist1_and_phi1->Fill(decayDistance1, phiLab1);
      histogramDensity_vs_X1_and_phi1->Fill(X1, phiLab1);
      histogramDensity_vs_decayDist2_and_phi2->Fill(decayDistance2, phiLab2);
      histogramDensity_vs_X2_and_phi2->Fill(X2, phiLab2);
    
      histogramDiTauPtDistribution->Fill(diTauPt);
      histogramDiTauEtaDistribution->Fill(diTauEta);
      histogramDiTauPhiDistribution->Fill(diTauPhi);
      histogramDiTauMassDistribution->Fill(diTauMass);

      double gjAngle_lab = TMath::Sqrt(square(leg1Eta - leg1VisEta) + square(leg1Phi - leg1VisPhi));
      double jacobiFactor = 2.*TMath::Pi()*gjAngle_lab;
      //double weight = jacobiFactor;
      //double weight = pdfWeight_vs_eta_m90->Eval(diTauEta);
      double weight = compDeltaFunctionJacobiFactor(X1)*compDeltaFunctionJacobiFactor(X2);
      //double weight = compDeltaFunctionJacobiFactor(X1)*compDeltaFunctionJacobiFactor(X2)*(X1*X2);      
      std::cout << "eta = " << diTauEta << ": weight = " << weight << std::endl;
      histogramSVfitMass_unweighted->Fill(diTauMass);
      histogramSVfitMass_weighted->Fill(diTauMass, weight);
      if ( gjAngle_lab > 0.00025 ) {
	histogramSVfitMass_gjAngle_labGt0_00025_unweighted->Fill(diTauMass);
	histogramSVfitMass_gjAngle_labGt0_00025_weighted->Fill(diTauMass, weight);
      } else {
	histogramSVfitMass_gjAngle_labLt0_00025_unweighted->Fill(diTauMass);
	histogramSVfitMass_gjAngle_labLt0_00025_weighted->Fill(diTauMass, weight);
      }

      histogramLeg1PtDistribution->Fill(leg1Pt);
      histogramLeg1EtaDistribution->Fill(leg1Eta);
      histogramLeg1PhiDistribution->Fill(leg1Phi);
      histogramLeg1dEtaDistribution->Fill(leg1Eta - leg1VisEta);
      histogramLeg1dPhiDistribution->Fill(leg1Phi - leg1VisPhi);
      histogramLeg1dPhi_vs_dEtaDistribution->Fill(leg1Eta - leg1VisEta, leg1Phi - leg1VisPhi);
      //std::cout << "leg1Eta = " << leg1Eta << ", leg1VisEta = " << leg1VisEta << " --> dEta = " << (leg1Eta - leg1VisEta) << ";"
      //	  << " leg1Phi = " << leg1Phi << ", leg1VisPhi = " << leg1VisPhi << " --> dPhi = " << (leg1Phi - leg1VisPhi) << std::endl;

      histogramLeg2PtDistribution->Fill(leg2Pt);
      histogramLeg2EtaDistribution->Fill(leg2Eta);
      histogramLeg2PhiDistribution->Fill(leg2Phi);
      histogramLeg2dEtaDistribution->Fill(leg2Eta - leg2VisEta);
      histogramLeg2dPhiDistribution->Fill(leg2Phi - leg2VisPhi);
      histogramLeg2dPhi_vs_dEtaDistribution->Fill(leg2Eta - leg2VisEta, leg2Phi - leg2VisPhi);
      //std::cout << "leg2Eta = " << leg2Eta << ", leg2VisEta = " << leg2VisEta << " --> dEta = " << (leg2Eta - leg2VisEta) << ";"
      //	  << " leg2Phi = " << leg2Phi << ", leg2VisPhi = " << leg2VisPhi << " --> dPhi = " << (leg2Phi - leg2VisPhi) << std::endl;

      histogramMoveDistribution->Fill(move);

      graphEvolution_X1->SetPoint(iEntry, move, X1);
      graphEvolution_phiLab1->SetPoint(iEntry, move, phiLab1);
      graphEvolution_decayDist1->SetPoint(iEntry, move, decayDistance1);
      graphEvolution_X2->SetPoint(iEntry, move, X2);
      graphEvolution_phiLab2->SetPoint(iEntry, move, phiLab2);
      graphEvolution_decayDist2->SetPoint(iEntry, move, decayDistance2);
    }
    
    histogramMass_vs_X1_and_X2->Divide(histogramDensity_vs_X1_and_X2);
    
    Float_t genX1, genX2;
    Float_t genDiTauPt, genDiTauEta, genDiTauPhi, genDiTauMass;
    //getGen(genNtupleFileName, run, ls, event, genX1, genX2, genDiTauPt, genDiTauEta, genDiTauPhi, genDiTauMass);
/*
    showHistogram2d(
      histogramDensity_vs_X1_and_X2, "X_{1}", "X_{2}", kNormByNegLogMax, 0.,   2., &genX1, &genX2, 
      Form("plots/plotMarkovChainMonitorTree_density_vs_X1andX2_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram2d(
      histogramMass_vs_X1_and_X2, "X_{1}", "X_{2}", kNormByValue, 0., 250., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_mass_vs_X1andX2_run%i_ls%i_ev%i.root", run, ls, event));

    showHistogram2d(
      histogramDensity_vs_decayDist1_and_phi1, "D^{decay}_{1}", "#phi^{lab}_{1}", kNormByNegLogMax, 0., 2., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_density_vs_decayDist1andPhi1_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram2d(
      histogramDensity_vs_X1_and_phi1, "X_{1}", "#phi^{lab}_{1}", kNormByNegLogMax, 0., 2., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_density_vs_X1andPhi1_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram2d(
      histogramDensity_vs_decayDist2_and_phi2, "D^{decay}_{2}", "#phi^{lab}_{2}", kNormByNegLogMax, 0., 2., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_density_vs_decayDist2andPhi2_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram2d(
      histogramDensity_vs_X2_and_phi2, "X_{2}", "#phi^{lab}_{2}", kNormByNegLogMax, 0., 2., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_density_vs_X2andPhi2_run%i_ls%i_ev%i.root", run, ls, event));
    
    showHistogram1d(
      histogramDiTauPtDistribution, "P_{T}^{#tau#tau}", &genDiTauPt, 
      Form("plots/plotMarkovChainMonitorTree_diTauPtDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramDiTauEtaDistribution, "#eta_{#tau#tau}", &genDiTauEta, 
      Form("plots/plotMarkovChainMonitorTree_diTauEtaDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramDiTauPhiDistribution, "#phi_{#tau#tau}", &genDiTauPhi, 
      Form("plots/plotMarkovChainMonitorTree_diTauPhiDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramDiTauMassDistribution, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_diTauMassDistribution_run%i_ls%i_ev%i.root", run, ls, event));
 */
    divideByBinWidth(histogramSVfitMass_unweighted);
    std::cout << "SVfit mass (unweighted):" << std::endl;
    printMaximum(histogramSVfitMass_unweighted);
    showHistogram1d(
      histogramSVfitMass_unweighted, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_svFitMass_unweighted_run%i_ls%i_ev%i.root", run, ls, event));
    divideByBinWidth(histogramSVfitMass_gjAngle_labGt0_00025_unweighted);
    showHistogram1d(
      histogramSVfitMass_gjAngle_labGt0_00025_unweighted, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_svFitMass_gjAngle_labGt0_00025_unweighted_run%i_ls%i_ev%i.root", run, ls, event));
    divideByBinWidth(histogramSVfitMass_gjAngle_labLt0_00025_unweighted);    
    showHistogram1d(
      histogramSVfitMass_gjAngle_labLt0_00025_unweighted, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_svFitMass_gjAngle_labLt0_00025_unweighted_run%i_ls%i_ev%i.root", run, ls, event));
    divideByBinWidth(histogramSVfitMass_weighted);
    std::cout << "SVfit mass (weighted by Jacobi factor):" << std::endl;
    printMaximum(histogramSVfitMass_weighted);
    showHistogram1d(
      histogramSVfitMass_weighted, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_svFitMass_weighted_run%i_ls%i_ev%i.root", run, ls, event));
    divideByBinWidth(histogramSVfitMass_gjAngle_labGt0_00025_weighted);
    showHistogram1d(
      histogramSVfitMass_gjAngle_labGt0_00025_weighted, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_svFitMass_gjAngle_labGt0_00025_weighted_run%i_ls%i_ev%i.root", run, ls, event));
    divideByBinWidth(histogramSVfitMass_gjAngle_labLt0_00025_weighted);    
    showHistogram1d(
      histogramSVfitMass_gjAngle_labLt0_00025_weighted, "M_{#tau#tau}", &genDiTauMass, 
      Form("plots/plotMarkovChainMonitorTree_svFitMass_gjAngle_labLt0_00025_weighted_run%i_ls%i_ev%i.root", run, ls, event));
/*
    showHistogram1d(
      histogramLeg1PtDistribution, "P_{T}^{1}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg1PtDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg1EtaDistribution, "#eta_{1}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg1EtaDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg1PhiDistribution, "#phi_{1}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg1PhiDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg1dEtaDistribution, "#eta_{1} - #eta_{vis}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg1dEta_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg1dPhiDistribution, "#phi_{1} - #phi_{vis}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg1dPhi_run%i_ls%i_ev%i.root", run, ls, event));
 */
    showHistogram2d(
      histogramLeg1dPhi_vs_dEtaDistribution, "#eta_{1} - #eta_{vis}", "#phi_{1} - #phi_{vis}", kNormByNegLogMax, 0., 2., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_leg1dPhi_vs_dEta_run%i_ls%i_ev%i.root", run, ls, event));
/*
    showHistogram1d(
      histogramLeg2PtDistribution, "P_{T}^{2}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg2PtDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg2EtaDistribution, "#eta_{2}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg2EtaDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg2PhiDistribution, "#phi_{2}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg2PhiDistribution_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg2dEtaDistribution, "#eta_{2} - #eta_{vis}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg2dEta_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram1d(
      histogramLeg2dPhiDistribution, "#phi_{2} - #phi_{vis}", 0, 
      Form("plots/plotMarkovChainMonitorTree_leg2dPhi_run%i_ls%i_ev%i.root", run, ls, event));
    showHistogram2d(
      histogramLeg2dPhi_vs_dEtaDistribution, "#eta_{2} - #eta_{vis}", "#phi_{2} - #phi_{vis}", kNormByNegLogMax, 0., 2., 0, 0, 
      Form("plots/plotMarkovChainMonitorTree_leg2dPhi_vs_dEta_run%i_ls%i_ev%i.root", run, ls, event));

    showHistogram1d(
      histogramMoveDistribution, "Move", 0, 
      Form("plots/plotMarkovChainMonitorTree_moveDistribution_run%i_ls%i_ev%i.root", run, ls, event));
 */
/*    
    showGraph(
      graphEvolution_X1, "Move", &genX1, 0., 1., "X_{1}", 
      Form("plots/plotMarkovChainMonitorTree_evolution_X1_run%i_ls%i_ev%i.root", run, ls, event));
    showGraph(
      graphEvolution_phiLab1, "Move", 0, -TMath::Pi(), +TMath::Pi(), "#phi^{lab}_{1}", 
      Form("plots/plotMarkovChainMonitorTree_evolution_phiLab1_run%i_ls%i_ev%i.root", run, ls, event));
    showGraph(
      graphEvolution_decayDist1, "Move", 0, -1., +1., "D^{decay}_{1}",
      Form("plots/plotMarkovChainMonitorTree_evolution_decayDist1_run%i_ls%i_ev%i.root", run, ls, event));
    showGraph(
      graphEvolution_X2, "Move", &genX2, 0., 1., "X_{2}",   
      Form("plots/plotMarkovChainMonitorTree_evolution_X2_run%i_ls%i_ev%i.root", run, ls, event));
    showGraph(
      graphEvolution_phiLab2, "Move", 0, -TMath::Pi(), +TMath::Pi(), "#phi^{lab}_{2}",   
      Form("plots/plotMarkovChainMonitorTree_evolution_phiLab2_run%i_ls%i_ev%i.root", run, ls, event));
    showGraph(
      graphEvolution_decayDist2, "Move", 0, -1., +1., "D^{decay}_{2}", 
      Form("plots/plotMarkovChainMonitorTree_evolution_decayDist2_run%i_ls%i_ev%i.root", run, ls, event));
 */
    delete histogramDensity_vs_X1_and_X2;
    delete histogramMass_vs_X1_and_X2;

    delete histogramDensity_vs_decayDist1_and_phi1;
    delete histogramDensity_vs_X1_and_phi1;
    delete histogramDensity_vs_decayDist2_and_phi2;
    delete histogramDensity_vs_X2_and_phi2;

    delete histogramDiTauPtDistribution;
    delete histogramDiTauEtaDistribution;
    delete histogramDiTauPhiDistribution;
    delete histogramDiTauMassDistribution;

    delete histogramSVfitMass_unweighted;
    delete histogramSVfitMass_gjAngle_labGt0_00025_unweighted;
    delete histogramSVfitMass_gjAngle_labLt0_00025_unweighted;
    delete histogramSVfitMass_weighted;
    delete histogramSVfitMass_gjAngle_labGt0_00025_weighted;
    delete histogramSVfitMass_gjAngle_labLt0_00025_weighted;

    delete histogramLeg1PtDistribution;
    delete histogramLeg1EtaDistribution;
    delete histogramLeg1PhiDistribution;
    delete histogramLeg1dPhi_vs_dEtaDistribution;

    delete histogramLeg2PtDistribution;
    delete histogramLeg2EtaDistribution;
    delete histogramLeg2PhiDistribution;
    delete histogramLeg2dPhi_vs_dEtaDistribution;

    delete histogramMoveDistribution;

    delete graphEvolution_X1;
    delete graphEvolution_phiLab1;
    delete graphEvolution_decayDist1;
    delete graphEvolution_X2;
    delete graphEvolution_phiLab2;
    delete graphEvolution_decayDist2;

    delete inputFile;
  }
}
