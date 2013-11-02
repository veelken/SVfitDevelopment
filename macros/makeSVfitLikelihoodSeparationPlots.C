
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

TH1* getHistogram(TFile* inputFile, const TString& directory, const TString& histogramName)
{  
  TString histogramName_full = TString("DQMData").Append("/").Append(directory);
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(histogramName);

  TH1* histogram = (TH1*)inputFile->Get(histogramName_full.Data());

  if ( histogram && !histogram->GetSumw2N() ) histogram->Sumw2();
  else if ( !histogram) 
    std::cerr << "Failed to load histogram = " << histogramName_full << " from file = " << inputFile->GetName() << " !!" << std::endl;

  return histogram;
}

double getIntegral(const TH1* histogram, double xMin, double xMax)
{
  double integral = 0.;
  int numBins = histogram->GetNbinsX();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double binCenter = histogram->GetBinCenter(iBin);
    if ( binCenter >= xMin && binCenter < xMax ) integral += histogram->GetBinContent(iBin);
  }
  return integral;
}

TGraph* makeSeparationCurve(const TH1* histogram_Z, const TH1* histogram_Higgs)
{
  const double xMin =   0.;
  const double xMax = 200.;

  double normalization_Z     = getIntegral(histogram_Z,     xMin, xMax);
  double normalization_Higgs = getIntegral(histogram_Higgs, xMin, xMax);

  assert(histogram_Z->GetNbinsX() == histogram_Higgs->GetNbinsX());

  const double epsilon = 1.e-3;

  double runningSum_Z     = 0.;
  double runningSum_Higgs = 0.;

  int numBins = histogram_Z->GetNbinsX();

  TGraph* graph = new TGraph(numBins);

  for ( int iBin = numBins; iBin >= 1; --iBin ) {
    assert(TMath::Abs(histogram_Z->GetBinCenter(iBin) - histogram_Higgs->GetBinCenter(iBin)) < epsilon);
    double binCenter = histogram_Z->GetBinCenter(iBin);
    if ( binCenter >= xMin && binCenter < xMax ) {
      runningSum_Z     += histogram_Z->GetBinContent(iBin);
      runningSum_Higgs += histogram_Higgs->GetBinContent(iBin);
    }
    graph->SetPoint(iBin - 1, runningSum_Higgs/normalization_Higgs, 1. - runningSum_Z/normalization_Z);
  }
  
  return graph;
}

TGraph* makeIntGraph(const TH1* histogram)
{
  const double xMin =   0.;
  const double xMax = 200.;

  double normalization = getIntegral(histogram, xMin, xMax);

  double runningSum = 0.;

  int numBins = histogram->GetNbinsX();

  TGraph* graph = new TGraph(numBins);

  for ( int iBin = numBins; iBin >= 1; --iBin ) {
    double binCenter = histogram->GetBinCenter(iBin);
    if ( binCenter >= xMin && binCenter < xMax ) {
      runningSum += histogram->GetBinContent(iBin);
    }
    graph->SetPoint(iBin - 1, binCenter, runningSum/normalization);
  }
  
  return graph;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		const std::vector<std::string>& entriesToPlot,
		std::map<std::string, TGraph*>& graphs,
		std::map<std::string, std::string>& legendEntries,
		int colors[], int lineWidths[], int lineStyles[],
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 100, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);

  size_t numEntriesToPlot = entriesToPlot.size();
  for ( size_t iEntryToPlot = 0; iEntryToPlot < numEntriesToPlot; ++iEntryToPlot ) {
    const std::string& entryToPlot = entriesToPlot[iEntryToPlot];
    TGraph* graph = graphs[entryToPlot];
    assert(graph);
    graph->SetLineColor(colors[iEntryToPlot]);
    graph->SetLineWidth(lineWidths[iEntryToPlot]);
    graph->SetLineStyle(lineStyles[iEntryToPlot]);
    graph->SetMarkerColor(colors[iEntryToPlot]);
    graph->Draw("L");

    legend->AddEntry(graph, legendEntries[entryToPlot].data(), "l");
  }

  legend->Draw();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete dummyHistogram;
  delete legend;
  delete canvas;  
}

void makeSVfitLikelihoodSeparationPlots()
{
  TString inputFileName_Z = "/data1/veelken/tmp/svFitStudies/testNSVfitTrackLikelihoods4_plots_Z_tautau_90_2013Aug06.root";
  TFile* inputFile_Z = new TFile(inputFileName_Z.Data());

  TString inputFileName_Higgs = "/data1/veelken/tmp/svFitStudies/testNSVfitTrackLikelihoods4_plots_Higgs_tautau_125_2013Aug06.root";
  TFile* inputFile_Higgs = new TFile(inputFileName_Higgs.Data());

  gROOT->SetBatch(true);

  std::vector<std::string> likelihoodOptions;
  //likelihoodOptions.push_back("const");
  //likelihoodOptions.push_back("decayKine");
  //likelihoodOptions.push_back("met");
  //likelihoodOptions.push_back("track");
  //likelihoodOptions.push_back("decayKinePlusMEt");
  //likelihoodOptions.push_back("trackPlusDecayKine");
  //likelihoodOptions.push_back("trackPlusMEt");
  //likelihoodOptions.push_back("trackPlusDecayKinePlusMEt");
  likelihoodOptions.push_back("decayKinePlusMEtVEGAS_max");
  likelihoodOptions.push_back("decayKinePlusMEtVEGAS_meanGt80%*max");
  likelihoodOptions.push_back("decayKinePlusMEtMarkovChain_max");
  likelihoodOptions.push_back("decayKinePlusMEtMarkovChain_meanGt80%*max");
  likelihoodOptions.push_back("trackPlusDecayKinePlusMEtVEGAS_max");
  likelihoodOptions.push_back("trackPlusDecayKinePlusMEtVEGAS_meanGt80%*max");
  likelihoodOptions.push_back("trackPlusDecayKinePlusMEtMarkovChain_max");
  likelihoodOptions.push_back("trackPlusDecayKinePlusMEtMarkovChain_meanGt80%*max");

  std::map<std::string, std::string> dqmDirectories; // key = likelihoodOption
  //dqmDirectories["const"]                     = "nSVfitProducerByIntegration2W";
  //dqmDirectories["decayKine"]                 = "nSVfitProducerByIntegration2WdecayKine";
  //dqmDirectories["met"]                       = "nSVfitProducerByIntegration2Wmet";
  //dqmDirectories["track"]                     = "nSVfitProducerByIntegration2Wtrack";
  //dqmDirectories["decayKinePlusMEt"]          = "nSVfitProducerByIntegration2WdecayKinePlusMEt";
  //dqmDirectories["trackPlusDecayKine"]        = "nSVfitProducerByIntegration2WtrackPlusDecayKine";
  //dqmDirectories["trackPlusMEt"]              = "nSVfitProducerByIntegration2WtrackPlusMEt";
  //dqmDirectories["trackPlusDecayKinePlusMEt"] = "nSVfitProducerByIntegration2WtrackPlusDecayKinePlusMEt";
  dqmDirectories["decayKinePlusMEtVEGAS_max"]                          = "nSVfitProducerByIntegration2WdecayKinePlusMEtMax";
  dqmDirectories["decayKinePlusMEtVEGAS_meanGt80%*max"]                = "nSVfitProducerByIntegration2WdecayKinePlusMEtMeanGt80pc";
  dqmDirectories["decayKinePlusMEtMarkovChain_max"]                    = "nSVfitProducerByIntegration2WtrackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEtMax";
  dqmDirectories["decayKinePlusMEtMarkovChain_meanGt80%*max"]          = "nSVfitProducerByIntegration2WtrackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEtMeanGt80pc";
  dqmDirectories["trackPlusDecayKinePlusMEtVEGAS_max"]                 = "nSVfitProducerByIntegrationWdecayKinePlusMEtMax";
  dqmDirectories["trackPlusDecayKinePlusMEtVEGAS_meanGt80%*max"]       = "nSVfitProducerByIntegrationWdecayKinePlusMEtMeanGt80pc";
  dqmDirectories["trackPlusDecayKinePlusMEtMarkovChain_max"]           = "nSVfitProducerByIntegrationWtrackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEtMax";
  dqmDirectories["trackPlusDecayKinePlusMEtMarkovChain_meanGt80%*max"] = "nSVfitProducerByIntegrationWtrackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEtMeanGt80pc";
  
  std::map<std::string, std::string> legendEntries; // key = likelihoodOption
  //legendEntries["const"]                      = "1";
  //legendEntries["decayKine"]                  = "Kine";
  //legendEntries["met"]                        = "MET";
  //legendEntries["track"]                      = "Trk";
  //legendEntries["decayKinePlusMEt"]           = "Kine+MET";
  //legendEntries["trackPlusDecayKine"]         = "Kine+Trk";
  //legendEntries["trackPlusMEt"]               = "MET+Trk";
  //legendEntries["trackPlusDecayKinePlusMEt"]  = "Kine+MET+Trk";
  legendEntries["decayKinePlusMEtVEGAS_max"]                           = "Kine+MET VEGAS_{max}";
  legendEntries["decayKinePlusMEtVEGAS_meanGt80%*max"]                 = "Kine+MET VEGAS_{mean}";
  legendEntries["decayKinePlusMEtMarkovChain_max"]                     = "Kine+MET MC_{max}";
  legendEntries["decayKinePlusMEtMarkovChain_meanGt80%*max"]           = "Kine+MET MC_{mean}";
  legendEntries["trackPlusDecayKinePlusMEtVEGAS_max"]                  = "Kine+MET+Trk VEGAS_{max}";
  legendEntries["trackPlusDecayKinePlusMEtVEGAS_meanGt80%*max"]        = "Kine+MET+Trk VEGAS_{mean}";
  legendEntries["trackPlusDecayKinePlusMEtMarkovChain_max"]            = "Kine+MET+Trk MC_{max}";
  legendEntries["trackPlusDecayKinePlusMEtMarkovChain_meanGt80%*max"]  = "Kine+MET+Trk MC_{mean}";

  std::string histogramName = "plotEntryType/svFitMass";

  std::map<std::string, const TH1*> histograms_Z;     // key = likelihoodOption
  std::map<std::string, const TH1*> histograms_Higgs; // key = likelihoodOption

  std::map<std::string, TGraph*> separationCurves; // key = likelihoodOption

  for ( std::vector<std::string>::const_iterator likelihoodOption = likelihoodOptions.begin();
	likelihoodOption != likelihoodOptions.end(); ++likelihoodOption ) {
    TH1* histogram_Z = getHistogram(inputFile_Z, dqmDirectories[*likelihoodOption], histogramName);
    histograms_Z[*likelihoodOption] = histogram_Z;

    TH1* histogram_Higgs = getHistogram(inputFile_Higgs, dqmDirectories[*likelihoodOption], histogramName);
    histograms_Higgs[*likelihoodOption] = histogram_Higgs;

    TGraph* separationCurve = makeSeparationCurve(histogram_Z, histogram_Higgs);
    separationCurves[*likelihoodOption] = separationCurve;
  }

  int colors[8]     = {  1,  2,  8,  4, 28,  6,  7, 14 };
  int lineWidths[8] = {  1,  1,  1,  1,  2,  2,  2,  1 }; 
  int lineStyles[8] = {  1,  1,  1,  1,  7,  7,  7,  1 }; 

  showGraphs(800, 600, 
	     likelihoodOptions,
	     separationCurves,
	     legendEntries,
	     colors, lineWidths, lineStyles,
	     0.050, 0.165, 0.165, 0.40, 0.40,
	     0., 1., "Int_{H}^{M > X}", 1.2,
	     0., 1., "1 - Int_{Z}^{M > X}", 1.2,
	     "makeSVfitLikelihoodSeparationPlots.pdf");

  TGraph* intGraph_Z_opt1     = makeIntGraph(histograms_Z["decayKinePlusMEtVEGAS_max"]);
  TGraph* intGraph_Higgs_opt1 = makeIntGraph(histograms_Higgs["decayKinePlusMEtVEGAS_meanGt80%*max"]);
  TGraph* intGraph_Z_opt2     = makeIntGraph(histograms_Z["trackPlusDecayKinePlusMEtVEGAS_max"]);
  TGraph* intGraph_Higgs_opt2 = makeIntGraph(histograms_Higgs["decayKinePlusMEtMarkovChain_meanGt80%*max"]);
  std::map<std::string, TGraph*> intGraphs;
  intGraphs["Z_opt1"]     = intGraph_Z_opt1;
  intGraphs["Higgs_opt1"] = intGraph_Higgs_opt1;
  intGraphs["Z_opt2"]     = intGraph_Z_opt2;
  intGraphs["Higgs_opt2"] = intGraph_Higgs_opt2;
  
  std::map<std::string, std::string> intLegendEntries;
  intLegendEntries["Z_opt1"]     = "Z/#gamma^{*} #rightarrow #tau#tau (Kine+MET)";
  intLegendEntries["Higgs_opt1"] = "Higgs #rightarrow #tau#tau (Kine+MET)";
  intLegendEntries["Z_opt2"]     = "Z/#gamma^{*} #rightarrow #tau#tau (Kine+MET+Trk)";
  intLegendEntries["Higgs_opt2"] = "Higgs #rightarrow #tau#tau (Kine+MET+Trk)";
  
  std::vector<std::string> intEntries;
  intEntries.push_back("Z_opt1");
  intEntries.push_back("Higgs_opt1");
  intEntries.push_back("Z_opt2");
  intEntries.push_back("Higgs_opt2");
  
  showGraphs(800, 600, 
	     intEntries,
	     intGraphs,
	     intLegendEntries,
	     colors, lineWidths, lineStyles,
	     0.035, 0.48, 0.725, 0.39, 0.16,
	     0., 250., "X / GeV", 1.2,
	     0., 1., "Int^{M > X}", 1.2,
	     "makeSVfitLikelihoodSeparationPlots_intGraphs.pdf");

  delete inputFile_Z;
  delete inputFile_Higgs;
}
