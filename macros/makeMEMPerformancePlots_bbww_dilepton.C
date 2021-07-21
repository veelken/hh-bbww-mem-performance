
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

enum { kUndefined, kSignal_lo, kSignal_nlo, kBackground_lo, kBackground_nlo };

enum { kProbSignal, kProbBackground, kLR };

TH1* loadHistogram(TFile* inputFile, const std::string& directory_part1, const std::string& directory_part2, int signal_or_background, const std::string& histogramName)
{  
  TString histogramName_full = directory_part1.data();
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(directory_part2.data());
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  if      ( signal_or_background == kSignal_lo      ) histogramName_full.Append("signal_lo");
  else if ( signal_or_background == kSignal_nlo     ) histogramName_full.Append("signal_nlo");
  else if ( signal_or_background == kBackground_lo  ) histogramName_full.Append("background_lo");
  else if ( signal_or_background == kBackground_nlo ) histogramName_full.Append("background_nlo");
  else assert(0);
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(histogramName.data());
  TH1* histogram = (TH1*)inputFile->Get(histogramName_full.Data());
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName_full.Data() << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  double integral = histogram->Integral(1, histogram->GetNbinsX()); // CV: exclude underflow and overflow bins
  if ( integral > 0. ) histogram->Scale(1./integral);
  return histogram;
}

TH1* addHistograms(const std::string& histogramSumName, const TH1* histogram1, const TH1* histogram2, const TH1* histogram3)
{
  TH1* histogramSum = (TH1*)histogram1->Clone(histogramSumName.data());
  histogramSum->Reset();
  if ( !histogramSum->GetSumw2N() ) histogramSum->Sumw2();
  histogramSum->Add(histogram1);
  histogramSum->Add(histogram2);
  if ( histogram3 ) histogramSum->Add(histogram2);
  double integral = histogramSum->Integral(1, histogramSum->GetNbinsX()); // CV: exclude underflow and overflow bins
  if ( integral > 0. ) histogramSum->Scale(1./integral);
  return histogramSum;
}

TH1* rebinHistogram(const TH1* histogram, int numBinsX)
{
  TH1* histogram_rebinned = (TH1*)histogram->Clone(Form("%s_rebinned", histogram->GetName()));
  if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();
  if ( numBinsX < histogram->GetNbinsX() && (histogram->GetNbinsX() % numBinsX) == 0 ) {
    histogram_rebinned->Rebin(histogram->GetNbinsX() / numBinsX);
  }
  return histogram_rebinned;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    TH1* histogram3, const std::string& legendEntry3,
		    TH1* histogram4, const std::string& legendEntry4,
		    int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
		    const std::string& labelText, double labelTextSize,
		    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		    int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = ( outputFileName.find("_prob") != std::string::npos ) ? 0.162 : 0.150;
  const double margin_bottom   = ( outputFileName.find("_prob") != std::string::npos ) ? 0.140 : 0.125;
  const double margin_right    = 0.015;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(margin_top);
  canvas->SetLeftMargin(margin_left);
  canvas->SetBottomMargin(margin_bottom);
  canvas->SetRightMargin(margin_right); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();

  assert(histogram1);
  double numBinsX_forRebinning;
  if ( xMax > xMin ) {
    const TAxis* xAxis = histogram1->GetXaxis();
    assert(xAxis->GetXmax() > xAxis->GetXmin());
    numBinsX_forRebinning = TMath::Nint(numBinsX*(xAxis->GetXmax() - xAxis->GetXmin())/(xMax - xMin));
    //std::cout << "numBinsX = " << numBinsX << ": numBinsX_forRebinning = " << numBinsX_forRebinning << std::endl;
  } else {
    numBinsX_forRebinning = numBinsX;
  }
  
  assert(histogram1);
  TH1* histogram1_rebinned = rebinHistogram(histogram1, numBinsX_forRebinning);
  histogram1_rebinned->SetFillColor(0);
  histogram1_rebinned->SetFillStyle(0);
  histogram1_rebinned->SetLineColor(colors[0]);
  histogram1_rebinned->SetLineStyle(lineStyles[0]);
  histogram1_rebinned->SetLineWidth(lineWidths[0]);
  histogram1_rebinned->SetMarkerColor(colors[0]);
  histogram1_rebinned->SetMarkerStyle(markerStyles[0]);
  histogram1_rebinned->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_rebinned = rebinHistogram(histogram2, numBinsX_forRebinning);
  histogram2_rebinned->SetFillColor(0);
  histogram2_rebinned->SetFillStyle(0);
  histogram2_rebinned->SetLineColor(colors[1]);
  histogram2_rebinned->SetLineStyle(lineStyles[1]);
  histogram2_rebinned->SetLineWidth(lineWidths[1]);
  histogram2_rebinned->SetMarkerColor(colors[1]);
  histogram2_rebinned->SetMarkerStyle(markerStyles[1]);
  histogram2_rebinned->SetMarkerSize(markerSizes[1]);

  TH1* histogram3_rebinned = nullptr;     
  if ( histogram3 ) {
    histogram3_rebinned = rebinHistogram(histogram3, numBinsX_forRebinning);
    histogram3_rebinned->SetFillColor(0);
    histogram3_rebinned->SetFillStyle(0);
    histogram3_rebinned->SetLineColor(colors[2]);
    histogram3_rebinned->SetLineStyle(lineStyles[2]);
    histogram3_rebinned->SetLineWidth(lineWidths[2]);
    histogram3_rebinned->SetMarkerColor(colors[2]);
    histogram3_rebinned->SetMarkerStyle(markerStyles[2]);
    histogram3_rebinned->SetMarkerSize(markerSizes[2]);
  }

  TH1* histogram4_rebinned = nullptr; 
  if ( histogram4 ) {
    histogram4_rebinned = rebinHistogram(histogram4, numBinsX_forRebinning);
    histogram4_rebinned->SetFillColor(0);
    histogram4_rebinned->SetFillStyle(0);
    histogram4_rebinned->SetLineColor(colors[3]);
    histogram4_rebinned->SetLineStyle(lineStyles[3]);
    histogram4_rebinned->SetLineWidth(lineWidths[3]);
    histogram4_rebinned->SetMarkerColor(colors[3]);
    histogram4_rebinned->SetMarkerStyle(markerStyles[3]);
    histogram4_rebinned->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis = histogram1_rebinned->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  if ( xMax > xMin ) {
    xAxis->SetRangeUser(xMin, xMax);
  }
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram1_rebinned->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(60);
  yAxis->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis->SetRangeUser(yMin, yMax);
  }
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.055);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);

  histogram1_rebinned->SetTitle("");
  histogram1_rebinned->SetStats(false);

  histogram1_rebinned->Draw(Form("%ssame", drawOptions[0].data()));
  histogram2_rebinned->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3_rebinned ) histogram3_rebinned->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4_rebinned ) histogram4_rebinned->Draw(Form("%ssame", drawOptions[3].data()));
  histogram1_rebinned->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntry1 != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogram1_rebinned, legendEntry1.data(), legendOptions[0].data());
    legend->AddEntry(histogram2_rebinned, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3 ) legend->AddEntry(histogram3_rebinned, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4 ) legend->AddEntry(histogram4_rebinned, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete histogram1_rebinned;
  delete histogram2_rebinned;
  delete histogram3_rebinned;
  delete histogram4_rebinned;
  delete canvas;
}

double square(double x)
{
  return x*x;
}

TH1* compRatioHistogram(const TH1* histogram_numerator, const TH1* histogram_denominator)
{
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());

  std::string histogramRatioName = Form("%s_div_%s", histogram_numerator->GetName(), histogram_denominator->GetName());
  TH1* histogramRatio = (TH1*)histogram_numerator->Clone(histogramRatioName.data());
  histogramRatio->Reset();
  if ( !histogramRatio->GetSumw2N() ) histogramRatio->Sumw2();
  histogramRatio->SetTitle("");
  histogramRatio->SetStats(false);

  int numBinsX = histogram_numerator->GetNbinsX();
  for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin ) {
    double binContent_numerator   = histogram_numerator->GetBinContent(idxBin);
    double binError_numerator     = histogram_numerator->GetBinError(idxBin);
    double binContent_denominator = histogram_denominator->GetBinContent(idxBin);
    double binError_denominator   = histogram_denominator->GetBinError(idxBin);

    if ( binContent_denominator > 0. ) {
      //double binContent_ratio = binContent_numerator/binContent_denominator - 1.;
      double binContent_ratio = binContent_numerator/binContent_denominator;
      double binErr2_ratio = 0.;
      binErr2_ratio += square(binError_numerator/binContent_denominator);
      binErr2_ratio += square(binError_denominator*binContent_numerator/square(binContent_denominator));
      double binErr_ratio = TMath::Sqrt(binErr2_ratio);
      histogramRatio->SetBinContent(idxBin, binContent_ratio);
      histogramRatio->SetBinError(idxBin, binErr_ratio);
    }
  }
  
  return histogramRatio;
}

void copyHistogramStyle(const TH1* histogram_source, TH1* histogram_target)
{
  histogram_target->SetMarkerColor(histogram_source->GetMarkerColor());
  histogram_target->SetMarkerSize(histogram_source->GetMarkerSize());
  histogram_target->SetMarkerStyle(histogram_source->GetMarkerStyle());
  histogram_target->SetLineColor(histogram_source->GetLineColor());
  histogram_target->SetLineWidth(histogram_source->GetLineWidth());
  histogram_target->SetLineStyle(histogram_source->GetLineStyle());
  histogram_target->SetFillColor(histogram_source->GetFillColor());
  histogram_target->SetFillStyle(histogram_source->GetFillStyle());
}

void showHistograms_wRatio(double canvasSizeX, double canvasSizeY,
			   TH1* histogramRef, const std::string& legendEntryRef,
			   TH1* histogram2, const std::string& legendEntry2,
			   TH1* histogram3, const std::string& legendEntry3,
			   TH1* histogram4, const std::string& legendEntry4,
			   int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
			   double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
			   const std::string& labelText, double labelTextSize,
			   double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
			   int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
			   bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
			   const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.155;
  const double margin_right    = 0.015;
  
  const double topPad_sizeY    = 0.69;
  const double clearance_sizeY = 0.02;
  const double bottomPad_sizeY = 0.29;
  const double sum_sizeY       = topPad_sizeY + clearance_sizeY + bottomPad_sizeY;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.);
  canvas->SetLeftMargin(0.);
  canvas->SetBottomMargin(0.);
  canvas->SetRightMargin(0.); 
  canvas->Draw();
  canvas->cd();

  const double topPad_x0 = 0.;
  const double topPad_x1 = 1.;
  const double topPad_y1 = 1.;
  const double topPad_y0 = topPad_y1 - topPad_sizeY/sum_sizeY;
  assert(topPad_x1 > topPad_x0 && topPad_y1 > topPad_y0);

  TPad* topPad = new TPad("topPad", "topPad", topPad_x0, topPad_y0, topPad_x1, topPad_y1);
  topPad->SetFillColor(10);
  topPad->SetBorderSize(0);
  topPad->SetTopMargin(margin_top);
  topPad->SetLeftMargin(margin_left);
  topPad->SetBottomMargin(0.000);
  topPad->SetRightMargin(margin_right); 
  topPad->SetLogx(false);
  topPad->SetLogy(useLogScale);
  topPad->Draw();
  topPad->cd();

  assert(histogramRef);
  double numBinsX_forRebinning;
  if ( xMax > xMin ) {
    const TAxis* xAxis = histogramRef->GetXaxis();
    assert(xAxis->GetXmax() > xAxis->GetXmin());
    numBinsX_forRebinning = TMath::Nint(numBinsX*(xAxis->GetXmax() - xAxis->GetXmin())/(xMax - xMin));
    //std::cout << "numBinsX = " << numBinsX << ": numBinsX_forRebinning = " << numBinsX_forRebinning << std::endl;
  } else {
    numBinsX_forRebinning = numBinsX;
  }
  
  assert(histogramRef);
  TH1* histogramRef_rebinned = rebinHistogram(histogramRef, numBinsX_forRebinning);
  histogramRef_rebinned->SetFillColor(0);
  histogramRef_rebinned->SetFillStyle(0);
  histogramRef_rebinned->SetLineColor(colors[0]);
  histogramRef_rebinned->SetLineStyle(lineStyles[0]);
  histogramRef_rebinned->SetLineWidth(lineWidths[0]);
  histogramRef_rebinned->SetMarkerColor(colors[0]);
  histogramRef_rebinned->SetMarkerStyle(markerStyles[0]);
  histogramRef_rebinned->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_rebinned = rebinHistogram(histogram2, numBinsX_forRebinning);
  histogram2_rebinned->SetFillColor(0);
  histogram2_rebinned->SetFillStyle(0);
  histogram2_rebinned->SetLineColor(colors[1]);
  histogram2_rebinned->SetLineStyle(lineStyles[1]);
  histogram2_rebinned->SetLineWidth(lineWidths[1]);
  histogram2_rebinned->SetMarkerColor(colors[1]);
  histogram2_rebinned->SetMarkerStyle(markerStyles[1]);
  histogram2_rebinned->SetMarkerSize(markerSizes[1]);

  TH1* histogram3_rebinned = nullptr;     
  if ( histogram3 ) {
    histogram3_rebinned = rebinHistogram(histogram3, numBinsX_forRebinning);
    histogram3_rebinned->SetFillColor(0);
    histogram3_rebinned->SetFillStyle(0);
    histogram3_rebinned->SetLineColor(colors[2]);
    histogram3_rebinned->SetLineStyle(lineStyles[2]);
    histogram3_rebinned->SetLineWidth(lineWidths[2]);
    histogram3_rebinned->SetMarkerColor(colors[2]);
    histogram3_rebinned->SetMarkerStyle(markerStyles[2]);
    histogram3_rebinned->SetMarkerSize(markerSizes[2]);
  }

  TH1* histogram4_rebinned = nullptr; 
  if ( histogram4 ) {
    histogram4_rebinned = rebinHistogram(histogram4, numBinsX_forRebinning);
    histogram4_rebinned->SetFillColor(0);
    histogram4_rebinned->SetFillStyle(0);
    histogram4_rebinned->SetLineColor(colors[3]);
    histogram4_rebinned->SetLineStyle(lineStyles[3]);
    histogram4_rebinned->SetLineWidth(lineWidths[3]);
    histogram4_rebinned->SetMarkerColor(colors[3]);
    histogram4_rebinned->SetMarkerStyle(markerStyles[3]);
    histogram4_rebinned->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis_top = histogramRef_rebinned->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitleSize(65);
  xAxis_top->SetTitleFont(43);
  //xAxis_top->SetLabelOffset(-0.01);
  xAxis_top->SetLabelSize(0.050);
  if ( xMax > xMin ) {
    xAxis_top->SetRangeUser(xMin, xMax);
  }
  xAxis_top->SetLabelFont(42);
  xAxis_top->SetTickLength(0.040);
  xAxis_top->SetNdivisions(505);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogramRef_rebinned->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(65);
  yAxis_top->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis_top->SetRangeUser(yMin, yMax);
  }
  //yAxis_top->SetLabelOffset(0.010);
  yAxis_top->SetLabelSize(0.055);
  yAxis_top->SetLabelFont(42);
  yAxis_top->SetTickLength(0.040);  
  yAxis_top->SetNdivisions(505);

  histogramRef_rebinned->SetTitle("");
  histogramRef_rebinned->SetStats(false);

  histogramRef_rebinned->Draw(drawOptions[0].data());
  histogram2_rebinned->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3_rebinned ) histogram3_rebinned->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4_rebinned ) histogram4_rebinned->Draw(Form("%ssame", drawOptions[3].data()));
  histogramRef_rebinned->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntryRef != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogramRef_rebinned, legendEntryRef.data(), legendOptions[0].data());
    legend->AddEntry(histogram2_rebinned, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3 ) legend->AddEntry(histogram3_rebinned, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4 ) legend->AddEntry(histogram4_rebinned, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->cd();

  const double bottomPad_x0 = 0.;
  const double bottomPad_x1 = 1.;
  const double bottomPad_y0 = 0.;
  const double bottomPad_y1 = bottomPad_y0 + bottomPad_sizeY/sum_sizeY;
  assert(bottomPad_x1 > bottomPad_x0 && bottomPad_y1 > bottomPad_y0);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", bottomPad_x0, bottomPad_y0, bottomPad_x1, bottomPad_y1);
  bottomPad->SetFillColor(10);
  bottomPad->SetBorderSize(0);
  bottomPad->SetTopMargin(0.000);
  bottomPad->SetLeftMargin(margin_left);
  bottomPad->SetBottomMargin(margin_bottom*topPad_sizeY/bottomPad_sizeY);
  bottomPad->SetRightMargin(margin_right);
  bottomPad->SetLogx(false);
  bottomPad->SetLogy(false);
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogramRatio2_rebinned = compRatioHistogram(histogram2_rebinned, histogramRef_rebinned);
  copyHistogramStyle(histogram2_rebinned, histogramRatio2_rebinned);

  TH1* histogramRatio3_rebinned = nullptr;
  if ( histogram3_rebinned ) {
    histogramRatio3_rebinned = compRatioHistogram(histogram3_rebinned, histogramRef_rebinned);
    copyHistogramStyle(histogram3_rebinned, histogramRatio3_rebinned);
  }

  TH1* histogramRatio4_rebinned = nullptr;
  if ( histogram4_rebinned ) {
    histogramRatio4_rebinned = compRatioHistogram(histogram4_rebinned, histogramRef_rebinned);
    copyHistogramStyle(histogram4_rebinned, histogramRatio4_rebinned);
  }

  TAxis* xAxis_bottom = histogramRatio2_rebinned->GetXaxis();
  xAxis_bottom->SetTitle(xAxisTitle.data());
  xAxis_bottom->SetTitleOffset(xAxisOffset);
  xAxis_bottom->SetTitleSize(65);
  xAxis_bottom->SetTitleFont(43);
  //xAxis_bottom->SetLabelOffset(-0.01);
  xAxis_bottom->SetLabelSize(0.132);
  if ( xMax > xMin ) {
    xAxis_bottom->SetRangeUser(xMin, xMax);
  }
  xAxis_bottom->SetLabelFont(42);
  xAxis_bottom->SetTickLength(0.105);
  xAxis_bottom->SetNdivisions(505);

  TAxis* yAxis_bottom = histogramRatio2_rebinned->GetYaxis();
  yAxis_bottom->SetTitle("Ratio");
  yAxis_bottom->SetTitleOffset(yAxisOffset);
  yAxis_bottom->SetTitleSize(65);
  yAxis_bottom->SetTitleFont(43);
  //yAxis_bottom->SetLabelOffset(0.010);
  yAxis_bottom->SetLabelSize(0.132);
  yAxis_bottom->SetLabelFont(42);
  yAxis_bottom->SetTickLength(0.060);
  yAxis_bottom->SetNdivisions(505);

  histogramRatio2_rebinned->SetMinimum(yMin_ratio);
  histogramRatio2_rebinned->SetMaximum(yMax_ratio);

  histogramRatio2_rebinned->Draw(drawOptions[1].data());
  
  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_bottom->GetXmin(), 1.);
  graph_line->SetPoint(1, xAxis_bottom->GetXmax(), 1.);
  graph_line->SetLineColor(colors[0]);
  graph_line->SetLineStyle(lineStyles[0]);
  graph_line->SetLineWidth(lineWidths[0]);
  graph_line->Draw("L");

  histogramRatio2_rebinned->Draw("same");

  if ( histogramRatio3_rebinned ) histogramRatio3_rebinned->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogramRatio4_rebinned ) histogramRatio4_rebinned->Draw(Form("%ssame", drawOptions[3].data()));
  histogramRatio2_rebinned->Draw("axissame");

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete histogramRef_rebinned;
  delete histogram2_rebinned;
  delete histogramRatio2_rebinned;
  delete histogram3_rebinned;
  delete histogramRatio3_rebinned;
  delete histogram4_rebinned;
  delete histogramRatio4_rebinned;
  delete graph_line;
  delete topPad;
  delete bottomPad;
  delete canvas;
}

TGraph* compGraphEfficiency(const std::string& graphName, const TH1* histogram)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numPoints = xAxis->GetNbins();
  TGraph* graphEfficiency = new TGraph(numPoints + 2);
  graphEfficiency->SetName(graphName.data());
  graphEfficiency->SetPoint(0, xAxis->GetXmin(), 1.0);
  double integral = histogram->Integral(1, histogram->GetNbinsX());
  double sum = 0.;
  for ( int idxPoint = 1; idxPoint <= numPoints; ++idxPoint ) {
    int idxBin = idxPoint;
    double binCenter = xAxis->GetBinCenter(idxBin);
    double binContent = histogram->GetBinContent(idxBin);
    sum += binContent;
    graphEfficiency->SetPoint(idxPoint, binCenter, 1.0 - (sum/integral));
  }
  graphEfficiency->SetPoint(numPoints + 1, xAxis->GetXmax(), 0.0);
  return graphEfficiency;
} 

TGraph* compGraphROC(const std::string& graphName, const TGraph* graphEfficiency_signal, const TGraph* graphEfficiency_background, bool useLogScale)
{
  //std::cout << "<compGraphROC>:" << std::endl;
  //std::cout << " graphName = " << graphName << std::endl;
  assert(graphEfficiency_signal->GetN() == graphEfficiency_background->GetN());
  int numPoints = graphEfficiency_signal->GetN();
  TGraph* graphROC = new TGraph(numPoints);
  graphROC->SetName(graphName.data());
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, efficiency_signal;
    graphEfficiency_signal->GetPoint(idxPoint, x, efficiency_signal);
    double efficiency_background = graphEfficiency_background->Eval(x);
    double xROC = efficiency_signal;
    double yROC;
    if ( useLogScale ) yROC = efficiency_background;
    else yROC = 1.0 - efficiency_background;
    //std::cout << "point #" << idxPoint << ": x = " << xROC << ", y = " << yROC << std::endl;
    graphROC->SetPoint(idxPoint, xROC, yROC);
  }
  return graphROC;
}

TGraph* compGraphROC(const std::string& graphName, const TH1* histogram_signal, const TH1* histogram_background, bool useLogScale)
{
  //std::cout << "<compGraphROC>:" << std::endl;
  //std::cout << " graphName = " << graphName << std::endl;
  assert(histogram_signal->GetNbinsX() == histogram_background->GetNbinsX());
  TGraph* graphEfficiency_signal = compGraphEfficiency(Form("%s_signal", graphName.data()), histogram_signal);
  TGraph* graphEfficiency_background = compGraphEfficiency(Form("%s_background", graphName.data()), histogram_background);
  TGraph* graphROC = compGraphROC(graphName, graphEfficiency_signal, graphEfficiency_background, useLogScale);
  delete graphEfficiency_signal;
  delete graphEfficiency_background;
  return graphROC;
}

struct graphPoint
{
  graphPoint(double x, double y)
    : x_(x)
    , y_(y)
  {}
  ~graphPoint() {}
  double x_;
  double y_;
};

TGraph* sparsifyGraph(TGraph* graph, double minDeltaX = 0.025)
{
  std::vector<graphPoint> graphPoints_sparsified;
  double x_last = -1.e+3;
  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    if ( x < 0.01 ) continue; // CV: prevent point @ zero signal efficiency and zero background rate from being drawn
    //if ( x > 0.99 ) continue; // CV: prevent point @ 100% signal efficiency and 100% background rate from being drawn
    if ( idxPoint == 0 || TMath::Abs(x - x_last) > minDeltaX || idxPoint == (numPoints - 1) ) {
      graphPoints_sparsified.push_back(graphPoint(x, y));
      x_last = x;
    }
  }
  int numPoints_sparsified = graphPoints_sparsified.size();
  TGraph* graph_sparsified = new TGraph(numPoints_sparsified);
  graph_sparsified->SetName(Form("%s_sparsified", graph->GetName()));
  for ( int idxPoint = 0; idxPoint < numPoints_sparsified; ++idxPoint ) {
    const graphPoint& graphPoint_sparsified = graphPoints_sparsified[idxPoint];
    graph_sparsified->SetPoint(idxPoint, graphPoint_sparsified.x_, graphPoint_sparsified.y_);
  }
  return graph_sparsified;
}

void setLineStyle(TLine* line)
{
  line->SetLineColor(14);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions,
		const std::string& labelText, double labelTextSize,
		double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.145;
  const double margin_right    = 0.015;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(margin_top);
  canvas->SetLeftMargin(margin_left);
  canvas->SetBottomMargin(margin_bottom);
  canvas->SetRightMargin(margin_right); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();

  assert(graph1);
  TGraph* graph1_sparsified = sparsifyGraph(graph1);
  graph1_sparsified->SetLineColor(colors[0]);
  graph1_sparsified->SetLineStyle(lineStyles[0]);
  graph1_sparsified->SetLineWidth(lineWidths[0]);
  graph1_sparsified->SetMarkerColor(colors[0]);
  graph1_sparsified->SetMarkerStyle(markerStyles[0]);
  graph1_sparsified->SetMarkerSize(markerSizes[0]);

  TGraph* graph2_sparsified = nullptr;
  if ( graph2 ) {    
    graph2_sparsified = sparsifyGraph(graph2);
    graph2_sparsified->SetLineColor(colors[1]);
    graph2_sparsified->SetLineStyle(lineStyles[1]);
    graph2_sparsified->SetLineWidth(lineWidths[1]);
    graph2_sparsified->SetMarkerColor(colors[1]);
    graph2_sparsified->SetMarkerStyle(markerStyles[1]);
    graph2_sparsified->SetMarkerSize(markerSizes[1]);
  }
  
  TGraph* graph3_sparsified = nullptr;
  if ( graph3 ) {    
    graph3_sparsified = sparsifyGraph(graph3);
    graph3_sparsified->SetLineColor(colors[2]);
    graph3_sparsified->SetLineStyle(lineStyles[2]);
    graph3_sparsified->SetLineWidth(lineWidths[2]);
    graph3_sparsified->SetMarkerColor(colors[2]);
    graph3_sparsified->SetMarkerStyle(markerStyles[2]);
    graph3_sparsified->SetMarkerSize(markerSizes[2]);
  }
  
  TGraph* graph4_sparsified = nullptr;
  if ( graph4 ) {
    graph4_sparsified = sparsifyGraph(graph4);
    graph4_sparsified->SetLineColor(colors[3]);
    graph4_sparsified->SetLineStyle(lineStyles[3]);
    graph4_sparsified->SetLineWidth(lineWidths[3]);
    graph4_sparsified->SetMarkerColor(colors[3]);
    graph4_sparsified->SetMarkerStyle(markerStyles[3]);
    graph4_sparsified->SetMarkerSize(markerSizes[3]);
  }

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numBinsX, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  assert(yMax > yMin);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(55);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(55);
  yAxis->SetTitleFont(43);
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.050);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);
  
  dummyHistogram->Draw("axis");
  graph1_sparsified->Draw(drawOptions[0].data());
  if ( graph2_sparsified ) graph2_sparsified->Draw(drawOptions[1].data());
  if ( graph3_sparsified ) graph3_sparsified->Draw(drawOptions[2].data());
  if ( graph4_sparsified ) graph4_sparsified->Draw(drawOptions[3].data());

  std::vector<TLine*> lines;
  if ( outputFileName == "hh_bbwwMEM_dilepton_effectOfFakes_2graphs_ROC.pdf" ) {
    double x = 0.35;
    assert(graph1);
    double y1 = graph1->Eval(x);
    TLine* line1_horizontal = new TLine(0., y1, x, y1);
    setLineStyle(line1_horizontal);
    line1_horizontal->Draw();
    lines.push_back(line1_horizontal);
    assert(graph2);
    double y2 = graph2->Eval(x);
    TLine* line2_horizontal = new TLine(0., y2, x, y2);
    setLineStyle(line2_horizontal);
    line2_horizontal->Draw();
    lines.push_back(line2_horizontal);
    TLine* line_vertical = new TLine(x, yMin, x, TMath::Max(y1, y2));
    setLineStyle(line_vertical);
    line_vertical->Draw();
    lines.push_back(line_vertical);
  } else if ( outputFileName == "hh_bbwwMEM_dilepton_effectOfFakes_ROC_missingBJet.pdf" ) {
    double x = 0.35;
    assert(graph1);
    double y = graph1->Eval(x);
    TLine* line_horizontal = new TLine(0., y, x, y);
    setLineStyle(line_horizontal);
    line_horizontal->Draw();
    lines.push_back(line_horizontal);
    TLine* line_vertical = new TLine(x, yMin, x, y);
    setLineStyle(line_vertical);
    line_vertical->Draw();
    lines.push_back(line_vertical);
  }

  dummyHistogram->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( graph2_sparsified || graph3_sparsified || graph4_sparsified ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(graph1_sparsified, legendEntry1.data(), legendOptions[0].data());
    if ( graph2_sparsified ) legend->AddEntry(graph2_sparsified, legendEntry2.data(), legendOptions[1].data());
    if ( graph3_sparsified ) legend->AddEntry(graph3_sparsified, legendEntry3.data(), legendOptions[2].data());
    if ( graph4_sparsified ) legend->AddEntry(graph4_sparsified, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete graph1_sparsified;
  delete graph2_sparsified;
  delete graph3_sparsified;
  delete graph4_sparsified;
  for ( std::vector<TLine*>::iterator line = lines.begin();
        line != lines.end(); ++line ) {
    delete (*line);
  }
  delete canvas;
}

TGraph* compRatioGraph(const TGraph* graph_numerator, const TGraph* graph_denominator)
{
  assert(graph_numerator->GetN() == graph_denominator->GetN());

  std::string graphRatioName = Form("%s_div_%s", graph_numerator->GetName(), graph_denominator->GetName());
  int numPoints = graph_numerator->GetN();
  TGraph* graphRatio = new TGraph(numPoints);

  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y_numerator;
    graph_numerator->GetPoint(idxPoint, x, y_numerator);

    double y_denominator = graph_denominator->Eval(x);

    if ( y_denominator > 0. ) {
      //double y_ratio = y_numerator/y_denominator - 1.;
      double y_ratio = y_numerator/y_denominator;
      graphRatio->SetPoint(idxPoint, x, y_ratio);
    }
  }
  
  return graphRatio;
}

void copyGraphStyle(const TGraph* graph_source, TGraph* graph_target)
{
  graph_target->SetMarkerColor(graph_source->GetMarkerColor());
  graph_target->SetMarkerSize(graph_source->GetMarkerSize());
  graph_target->SetMarkerStyle(graph_source->GetMarkerStyle());
  graph_target->SetLineColor(graph_source->GetLineColor());
  graph_target->SetLineWidth(graph_source->GetLineWidth());
  graph_target->SetLineStyle(graph_source->GetLineStyle());
}

void showGraphs_wRatio(double canvasSizeX, double canvasSizeY,
		       TGraph* graphRef, const std::string& legendEntryRef,
		       TGraph* graph2, const std::string& legendEntry2,
		       TGraph* graph3, const std::string& legendEntry3,
		       TGraph* graph4, const std::string& legendEntry4,
		       int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		       double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
		       const std::string& labelText, double labelTextSize,
		       double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		       int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		       bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
		       const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.145;
  const double margin_right    = 0.015;
  
  const double topPad_sizeY    = 0.69;
  const double clearance_sizeY = 0.02;
  const double bottomPad_sizeY = 0.29;
  const double sum_sizeY       = topPad_sizeY + clearance_sizeY + bottomPad_sizeY;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.);
  canvas->SetLeftMargin(0.);
  canvas->SetBottomMargin(0.);
  canvas->SetRightMargin(0.); 
  canvas->Draw();
  canvas->cd();

  const double topPad_x0 = 0.;
  const double topPad_x1 = 1.;
  const double topPad_y1 = 1.;
  const double topPad_y0 = topPad_y1 - topPad_sizeY/sum_sizeY;
  assert(topPad_x1 > topPad_x0 && topPad_y1 > topPad_y0);

  TPad* topPad = new TPad("topPad", "topPad", topPad_x0, topPad_y0, topPad_x1, topPad_y1);
  topPad->SetFillColor(10);
  topPad->SetBorderSize(0);
  topPad->SetTopMargin(margin_top);
  topPad->SetLeftMargin(margin_left);
  topPad->SetBottomMargin(0.000);
  topPad->SetRightMargin(margin_right); 
  topPad->SetLogx(false);
  topPad->SetLogy(useLogScale);
  topPad->Draw();
  topPad->cd();

  assert(graphRef);
  TGraph* graphRef_sparsified = sparsifyGraph(graphRef);
  graphRef_sparsified->SetLineColor(colors[0]);
  graphRef_sparsified->SetLineStyle(lineStyles[0]);
  graphRef_sparsified->SetLineWidth(lineWidths[0]);
  graphRef_sparsified->SetMarkerColor(colors[0]);
  graphRef_sparsified->SetMarkerStyle(markerStyles[0]);
  graphRef_sparsified->SetMarkerSize(markerSizes[0]);

  assert(graph2);
  TGraph* graph2_sparsified = sparsifyGraph(graph2);
  graph2_sparsified->SetLineColor(colors[1]);
  graph2_sparsified->SetLineStyle(lineStyles[1]);
  graph2_sparsified->SetLineWidth(lineWidths[1]);
  graph2_sparsified->SetMarkerColor(colors[1]);
  graph2_sparsified->SetMarkerStyle(markerStyles[1]);
  graph2_sparsified->SetMarkerSize(markerSizes[1]);
  
  TGraph* graph3_sparsified = nullptr;
  if ( graph3 ) {    
    graph3_sparsified = sparsifyGraph(graph3);
    graph3_sparsified->SetLineColor(colors[2]);
    graph3_sparsified->SetLineStyle(lineStyles[2]);
    graph3_sparsified->SetLineWidth(lineWidths[2]);
    graph3_sparsified->SetMarkerColor(colors[2]);
    graph3_sparsified->SetMarkerStyle(markerStyles[2]);
    graph3_sparsified->SetMarkerSize(markerSizes[2]);
  }
  
  TGraph* graph4_sparsified = nullptr;
  if ( graph4 ) {
    graph4_sparsified = sparsifyGraph(graph4);
    graph4_sparsified->SetLineColor(colors[3]);
    graph4_sparsified->SetLineStyle(lineStyles[3]);
    graph4_sparsified->SetLineWidth(lineWidths[3]);
    graph4_sparsified->SetMarkerColor(colors[3]);
    graph4_sparsified->SetMarkerStyle(markerStyles[3]);
    graph4_sparsified->SetMarkerSize(markerSizes[3]);
  }

  TH1* dummyHistogram_top = new TH1D("dummyHistogram_top", "dummyHistogram_top", numBinsX, xMin, xMax);
  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  assert(yMax > yMin);
  dummyHistogram_top->SetMinimum(yMin);
  dummyHistogram_top->SetMaximum(yMax);

  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitleSize(60);
  xAxis_top->SetTitleFont(43);
  //xAxis_top->SetLabelOffset(-0.01);
  xAxis_top->SetLabelSize(0.050);
  xAxis_top->SetLabelFont(42);
  xAxis_top->SetTickLength(0.040);
  xAxis_top->SetNdivisions(505);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(60);
  yAxis_top->SetTitleFont(43);
  //yAxis_top->SetLabelOffset(0.010);
  yAxis_top->SetLabelSize(0.055);
  yAxis_top->SetLabelFont(42);
  yAxis_top->SetTickLength(0.040);  
  yAxis_top->SetNdivisions(505);
  
  dummyHistogram_top->Draw("axis");
  graphRef_sparsified->Draw(drawOptions[0].data());
  graph2_sparsified->Draw(drawOptions[1].data());
  if ( graph3_sparsified ) graph3_sparsified->Draw(drawOptions[2].data());
  if ( graph4_sparsified ) graph4_sparsified->Draw(drawOptions[3].data());
  dummyHistogram_top->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntryRef != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(graphRef_sparsified, legendEntryRef.data(), legendOptions[0].data());
    legend->AddEntry(graph2_sparsified, legendEntry2.data(), legendOptions[1].data());
    if ( graph3_sparsified ) legend->AddEntry(graph3_sparsified, legendEntry3.data(), legendOptions[2].data());
    if ( graph4_sparsified ) legend->AddEntry(graph4_sparsified, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->cd();

  const double bottomPad_x0 = 0.;
  const double bottomPad_x1 = 1.;
  const double bottomPad_y0 = 0.;
  const double bottomPad_y1 = bottomPad_y0 + bottomPad_sizeY/sum_sizeY;
  assert(bottomPad_x1 > bottomPad_x0 && bottomPad_y1 > bottomPad_y0);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", bottomPad_x0, bottomPad_y0, bottomPad_x1, bottomPad_y1);
  bottomPad->SetFillColor(10);
  bottomPad->SetBorderSize(0);
  bottomPad->SetTopMargin(0.000);
  bottomPad->SetLeftMargin(margin_left);
  bottomPad->SetBottomMargin(margin_bottom*topPad_sizeY/bottomPad_sizeY);
  bottomPad->SetRightMargin(margin_right);
  bottomPad->SetLogx(false);
  bottomPad->SetLogy(false);
  bottomPad->Draw();
  bottomPad->cd();

  TH1* dummyHistogram_bottom = new TH1D("dummyHistogram_bottom", "dummyHistogram_bottom", numBinsX, xMin, xMax);
  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);

  TGraph* graphRatio2_sparsified = compRatioGraph(graph2_sparsified, graphRef_sparsified);
  copyGraphStyle(graph2_sparsified, graphRatio2_sparsified);

  TGraph* graphRatio3_sparsified = nullptr;
  if ( graph3_sparsified ) {
    graphRatio3_sparsified = compRatioGraph(graph3_sparsified, graphRef_sparsified);
    copyGraphStyle(graph3_sparsified, graphRatio3_sparsified);
  }

  TGraph* graphRatio4_sparsified = nullptr;
  if ( graph4_sparsified ) {
    graphRatio4_sparsified = compRatioGraph(graph4_sparsified, graphRef_sparsified);
    copyGraphStyle(graph4_sparsified, graphRatio4_sparsified);
  }

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle(xAxisTitle.data());
  xAxis_bottom->SetTitleOffset(xAxisOffset);
  xAxis_bottom->SetTitleSize(60);
  xAxis_bottom->SetTitleFont(43);
  //xAxis_bottom->SetLabelOffset(-0.01);
  xAxis_bottom->SetLabelSize(0.132);
  xAxis_bottom->SetLabelFont(42);
  xAxis_bottom->SetTickLength(0.105);
  xAxis_bottom->SetNdivisions(505);

  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("Ratio");
  yAxis_bottom->SetTitleOffset(yAxisOffset);
  yAxis_bottom->SetTitleSize(60);
  yAxis_bottom->SetTitleFont(43);
  //yAxis_bottom->SetLabelOffset(0.010);
  yAxis_bottom->SetLabelSize(0.132);
  yAxis_bottom->SetLabelFont(42);
  yAxis_bottom->SetTickLength(0.060);
  yAxis_bottom->SetNdivisions(505);

  dummyHistogram_bottom->SetMinimum(yMin_ratio);
  dummyHistogram_bottom->SetMaximum(yMax_ratio);

  dummyHistogram_bottom->Draw("axis");

  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_bottom->GetXmin(), 1.);
  graph_line->SetPoint(1, xAxis_bottom->GetXmax(), 1.);
  graph_line->SetLineColor(colors[0]);
  graph_line->SetLineStyle(lineStyles[0]);
  graph_line->SetLineWidth(lineWidths[0]);
  graph_line->Draw("L");

  graphRatio2_sparsified->Draw(drawOptions[1].data());
  if ( graphRatio3_sparsified ) graphRatio3_sparsified->Draw(drawOptions[2].data());
  if ( graphRatio4_sparsified ) graphRatio4_sparsified->Draw(drawOptions[3].data());
  dummyHistogram_bottom->Draw("axissame");

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete dummyHistogram_top;
  delete dummyHistogram_bottom;
  delete graphRef_sparsified;
  delete graph2_sparsified;
  delete graphRatio2_sparsified;
  delete graph3_sparsified;
  delete graphRatio3_sparsified;
  delete graph4_sparsified;
  delete graphRatio4_sparsified;
  delete graph_line;
  delete canvas;
}

void makeMEMPerformancePlots_bbww_dilepton()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  bool makePlots_signal_vs_background = true;
  bool makePlots_effectOfFakes = true;
  bool makePlots_effectOfSmearing = true;
  bool makePlots_effectOfHigherOrders = true;

  std::string inputFilePath = "/hdfs/local/veelken/hhAnalysis/2016/2021May17/histograms/hh_bbwwMEM_dilepton/";
  std::string inputFileName = "histograms_harvested_stage2_hh_bbwwMEM_dilepton.root";
  TString inputFileName_full = inputFilePath.data();
  if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
  inputFileName_full.Append(inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.Data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName_full.Data() << " !!" << std::endl;
    assert(0);
  }

  std::map<bool, std::map<bool, std::string>> directories_part1; // key = apply_jetSmearing, apply_metSmearing
  directories_part1[false][false]      = "hh_bbwwMEM_dilepton_jetSmearingDisabled_metSmearingDisabled";
  directories_part1[false][true]       = "hh_bbwwMEM_dilepton_jetSmearingDisabled_metSmearingEnabled";
  directories_part1[true][false]       = "hh_bbwwMEM_dilepton_jetSmearingEnabled_metSmearingDisabled";
  directories_part1[true][true]        = "hh_bbwwMEM_dilepton_jetSmearingEnabled_metSmearingEnabled";
  
  std::map<bool, std::map<bool, std::string>> directories_part2; // key = selGenBJet_lead_isFake, selGenBJet_sublead_isFake
  directories_part2[false][false]      = "sel/mem_genuineLeadingBJet_genuineSubleadingBJet";
  directories_part2[false][true]       = "sel/mem_genuineLeadingBJet_fakeSubleadingBJet";
  directories_part2[true][false]       = "sel/mem_fakeLeadingBJet_genuineSubleadingBJet";
  directories_part2[true][true]        = "sel/mem_fakeLeadingBJet_fakeSubleadingBJet";

  std::map<bool, std::string> directories_part2_missingBJet; // key = selGenBJet_isFake
  directories_part2_missingBJet[false] = "sel/mem_missingBJet_genuineBJet";
  directories_part2_missingBJet[true]  = "sel/mem_missingBJet_fakeBJet";

  std::map<int, std::string> histogramNames; // key = { kProbSignal, kProbBackground, kLR }
  histogramNames[kProbSignal]        = "log_memProb_signal";
  histogramNames[kProbBackground]    = "log_memProb_background";
  histogramNames[kLR]                = "memLR";

  std::string labelText_signal = "HH #rightarrow b#bar{b} WW^{*} #rightarrow b#bar{b} l^{+}#nu l^{-}#bar{#nu}";
  std::string labelText_background = "t#bar{t} #rightarrow bW #bar{b}W #rightarrow b l^{+}#nu #bar{b} l^{-}#bar{#nu}";
  //std::string labelText_signal_vs_background = Form("%s vs %s", labelText_signal.data(), labelText_background.data());
  std::string labelText_signal_vs_background = "";

  int showHistograms_canvasSizeX = 1050;
  int showHistograms_canvasSizeY =  950;
  int showHistograms_canvasSizeY_wRatio = 1150;
  double showHistograms_xAxisOffset = 0.96;
  double showHistograms_xAxisOffset_wRatio = 2.90;
  double showHistograms_yAxisOffset = 1.21;
  double showHistograms_yAxisOffset_wRatio = 1.45;
  int showHistograms_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showHistograms_markerStyles[4] = { 20, 24, 21, 25 };
  int showHistograms_markerSizes[4]  = { 2, 2, 2, 2 };
  int showHistograms_lineStyles[4]   = { 1, 1, 1, 7 };
  int showHistograms_lineWidths[4]   = { 3, 2, 2, 3 };
  std::vector<std::string> showHistograms_drawOptions = { "hist", "ep", "ep", "hist" };
  std::vector<std::string> showHistograms_legendOptions = { "l", "p", "p", "l" };

  int showHistograms_signal_vs_background_colors[2]       = { kBlack, kRed };
  int showHistograms_signal_vs_background_markerStyles[2] = { 20, 24 };
  int showHistograms_signal_vs_background_markerSizes[2]  = { 2, 2 };
  int showHistograms_signal_vs_background_lineStyles[2]   = { 1, 1 };
  int showHistograms_signal_vs_background_lineWidths[2]   = { 2, 2 };
  std::vector<std::string> showHistograms_signal_vs_background_drawOptions = { "ep", "ep" };
  std::vector<std::string> showHistograms_signal_vs_background_legendOptions = { "p", "p" };

  int numBinsX_memLR = 40;
  double xMin_memLR = 0.;
  double xMax_memLR = 1.;
  double yMin_memLR = 1.1e-5;
  double yMin_memLR_wRatio = 1.1e-5;
  double yMax_memLR = 1.9e0;

  int numBinsX_probS = 40;
  double xMin_probS = -70.;
  double xMax_probS = +10.;
  double yMin_probS = 1.1e-5;
  double yMin_probS_wRatio = 1.1e-5;
  double yMax_probS = 1.9e0;

  int numBinsX_probB = 40;
  double xMin_probB = -70.;
  double xMax_probB = +10.;
  double yMin_probB = 3.1e-5;
  double yMin_probB_wRatio = 3.1e-5;
  double yMax_probB = 1.9e0;

  int showGraphs_canvasSizeX = 1050;
  int showGraphs_canvasSizeY =  950;
  int showGraphs_canvasSizeY_wRatio = 1150;
  double showGraphs_xAxisOffset = 1.18;
  double showGraphs_xAxisOffset_wRatio = 3.00;
  double showGraphs_yAxisOffset = 1.44;
  double showGraphs_yAxisOffset_wRatio = 1.65;
  int showGraphs_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showGraphs_markerStyles[4] = { 20, 24, 21, 25 };
  int showGraphs_markerSizes[4]  = { 2, 2, 2, 2 };
  int showGraphs_lineStyles[4]   = { 1, 1, 1, 1 };
  int showGraphs_lineWidths[4]   = { 2, 2, 2, 2 };
  std::vector<std::string> showGraphs_drawOptions = { "Lp", "Lp", "Lp", "Lp" };
  std::vector<std::string> showGraphs_legendOptions = { "lp", "lp", "lp", "lp" };

  if ( makePlots_signal_vs_background ) {
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_noSmearing_2genuineBJets_signal, "Signal",
      histogram_probS_noSmearing_2genuineBJets_background, "Background",
      nullptr, "",
      nullptr, "",
      showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
      showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
      "", 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}", showHistograms_yAxisOffset, 
      "hh_bbwwMEM_dilepton_signal_vs_background_probS_unsmeared.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_noSmearing_2genuineBJets_signal, "Signal",
      histogram_probB_noSmearing_2genuineBJets_background, "Background",
      nullptr, "",
      nullptr, "",
      showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
      showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
      "", 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}", showHistograms_xAxisOffset,
      true, yMin_probB, yMax_probB, "dN/dlog w_{1}", showHistograms_yAxisOffset, 
      "hh_bbwwMEM_dilepton_signal_vs_background_probB_unsmeared.pdf");

    TH1* histogram_memLR_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_noSmearing_2genuineBJets_signal, "Signal",
      histogram_memLR_noSmearing_2genuineBJets_background, "Background",
      nullptr, "",
      nullptr, "",
      showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
      showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
      "", 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_signal_vs_background_memLR_unsmeared.pdf");

    TGraph* graph_ROC_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_noSmearing_2genuineBJets",
      histogram_memLR_noSmearing_2genuineBJets_signal, histogram_memLR_noSmearing_2genuineBJets_background, true);

    showGraphs(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY,
      graph_ROC_noSmearing_2genuineBJets_logScale, "",
      nullptr, "",
      nullptr, "",
      nullptr, "",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      0.055, 0.23, 0.86, 0.33, 0.08, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
      true, 1.e-3, 1.e0, "Background Rate", showGraphs_yAxisOffset, 
      "hh_bbwwMEM_dilepton_ROC_unsmeared.pdf");
    //-------------------------------------------------------------------------------------------------
  }

  if ( makePlots_effectOfFakes ) {
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_fakeLeadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_fakeSubleadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_2fakeBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kSignal_lo, histogramNames[kProbSignal]);
  
    TH1* histogram_probS_noSmearing_1fakeBJet_signal = addHistograms("histogram_probS_noSmearing_1fakeBJet_signal",
      histogram_probS_noSmearing_fakeLeadBJet_signal, 
      histogram_probS_noSmearing_fakeSubleadBJet_signal, 
      histogram_probS_noSmearing_2fakeBJets_signal);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
      histogram_probS_noSmearing_1fakeBJet_signal, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}", showHistograms_yAxisOffset, 
      "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probS_signal.pdf");
    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
      histogram_probS_noSmearing_fakeLeadBJet_signal, "fake lead. b-jet",
      histogram_probS_noSmearing_fakeSubleadBJet_signal, "fake sublead. b-jet",
      histogram_probS_noSmearing_2fakeBJets_signal, "2 fake b-jets",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}", showHistograms_yAxisOffset, 
      "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probS_signal.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_noSmearing_fakeBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kSignal_lo, histogramNames[kProbSignal]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_missingBJet_noSmearing_genuineBJet_signal, "genuine b-jet",
      histogram_probS_missingBJet_noSmearing_fakeBJet_signal, "fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}^{m}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}^{m}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_probS_missingBJet_signal.pdf");
    
    TH1* histogram_probB_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeLeadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeSubleadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_2fakeBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kSignal_lo, histogramNames[kProbBackground]);
  
    TH1* histogram_probB_noSmearing_1fakeBJet_signal = addHistograms("histogram_probB_noSmearing_1fakeBJet_signal",
      histogram_probB_noSmearing_fakeLeadBJet_signal, 
      histogram_probB_noSmearing_fakeSubleadBJet_signal, 
      histogram_probB_noSmearing_2fakeBJets_signal);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
      histogram_probB_noSmearing_1fakeBJet_signal, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}", showHistograms_xAxisOffset,
      true, yMin_probB, yMax_probB, "dN/dlog w_{1}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probB_signal.pdf");
    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
      histogram_probB_noSmearing_fakeLeadBJet_signal, "fake lead. b-jet",
      histogram_probB_noSmearing_fakeSubleadBJet_signal, "fake sublead. b-jet",
      histogram_probB_noSmearing_2fakeBJets_signal, "2 fake b-jets",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.045, 0.25, 0.69, 0.28, 0.23, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}",showHistograms_xAxisOffset,
      true, 1.e-4, 1.e0, "dN/dlog w_{1}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probB_signal.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_noSmearing_fakeBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kSignal_lo, histogramNames[kProbBackground]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_missingBJet_noSmearing_genuineBJet_signal, "genuine b-jet",
      histogram_probB_missingBJet_noSmearing_fakeBJet_signal, "fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}^{m}",showHistograms_xAxisOffset,
      true, yMin_probB, yMax_probB, "dN/dlog w_{1}^{m}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_probB_missingBJet_signal.pdf");
    //-------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_fakeLeadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_fakeSubleadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_noSmearing_2fakeBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kBackground_lo, histogramNames[kProbSignal]);
  
    TH1* histogram_probS_noSmearing_1fakeBJet_background = addHistograms("histogram_probS_noSmearing_1fakeBJet_background",
      histogram_probS_noSmearing_fakeLeadBJet_background, 
      histogram_probS_noSmearing_fakeSubleadBJet_background, 
      histogram_probS_noSmearing_2fakeBJets_background);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_noSmearing_2genuineBJets_background, "2 genuine b-jets",
      histogram_probS_noSmearing_1fakeBJet_background, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probS_background.pdf");
    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_noSmearing_2genuineBJets_background, "2 genuine b-jets",
      histogram_probS_noSmearing_fakeLeadBJet_background, "fake lead. b-jet",
      histogram_probS_noSmearing_fakeSubleadBJet_background, "fake sublead. b-jet",
      histogram_probS_noSmearing_2fakeBJets_background, "2 fake b-jets",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.045, 0.25, 0.69, 0.28, 0.23, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probS_background.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_noSmearing_fakeBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kBackground_lo, histogramNames[kProbSignal]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probS_missingBJet_noSmearing_genuineBJet_background, "genuine b-jet",
      histogram_probS_missingBJet_noSmearing_fakeBJet_background, "fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}^{m}", showHistograms_xAxisOffset,
      true, yMin_probS, yMax_probS, "dN/dlog w_{0}^{m}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_probS_missingBJet_background.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeLeadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_fakeSubleadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_noSmearing_2fakeBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kBackground_lo, histogramNames[kProbBackground]);
  
    TH1* histogram_probB_noSmearing_1fakeBJet_background = addHistograms("histogram_probB_noSmearing_1fakeBJet_background",
      histogram_probB_noSmearing_fakeLeadBJet_background, 
      histogram_probB_noSmearing_fakeSubleadBJet_background, 
      histogram_probB_noSmearing_2fakeBJets_background);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_noSmearing_2genuineBJets_background, "2 genuine b-jets",
      histogram_probB_noSmearing_1fakeBJet_background, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}",showHistograms_xAxisOffset,
      true, yMin_probB, yMax_probB, "dN/dlog w_{1}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_probB_background.pdf");
    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_noSmearing_2genuineBJets_background, "2 genuine b-jets",
      histogram_probB_noSmearing_fakeLeadBJet_background, "fake lead. b-jet",
      histogram_probB_noSmearing_fakeSubleadBJet_background, "fake sublead. b-jet",
      histogram_probB_noSmearing_2fakeBJets_background, "2 fake b-jets",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.045, 0.25, 0.69, 0.28, 0.23, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}",showHistograms_xAxisOffset,
      true, yMin_probB, yMax_probB, "dN/dlog w_{1}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_probB_background.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_noSmearing_fakeBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kBackground_lo, histogramNames[kProbBackground]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_probB_missingBJet_noSmearing_genuineBJet_background, "genuine b-jet",
      histogram_probB_missingBJet_noSmearing_fakeBJet_background, "fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}^{m}",showHistograms_xAxisOffset,
      true, yMin_probB, yMax_probB, "dN/dlog w_{1}^{m}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_probB_missingBJet_background.pdf");
    //-------------------------------------------------------------------------------------------------
  
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_memLR_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeLeadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeSubleadBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_2fakeBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kSignal_lo, histogramNames[kLR]);
  
    TH1* histogram_memLR_noSmearing_1fakeBJet_signal = addHistograms("histogram_memLR_noSmearing_1fakeBJet_signal",
      histogram_memLR_noSmearing_fakeLeadBJet_signal, 
      histogram_memLR_noSmearing_fakeSubleadBJet_signal, 
      histogram_memLR_noSmearing_2fakeBJets_signal);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
      histogram_memLR_noSmearing_1fakeBJet_signal, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.40, 0.72, 0.33, 0.15, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_memLR_signal.pdf");
    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_noSmearing_2genuineBJets_signal, "2 genuine b-jets",
      histogram_memLR_noSmearing_fakeLeadBJet_signal, "fake lead. b-jet",
      histogram_memLR_noSmearing_fakeSubleadBJet_signal, "fake sublead. b-jet",
      histogram_memLR_noSmearing_2fakeBJets_signal, "2 fake b-jets",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.045, 0.40, 0.69, 0.28, 0.23, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_memLR_signal.pdf");

    TH1* histogram_memLR_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeLeadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_fakeSubleadBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][true], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_noSmearing_2fakeBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[true][true], kBackground_lo, histogramNames[kLR]);
  
    TH1* histogram_memLR_noSmearing_1fakeBJet_background = addHistograms("histogram_memLR_noSmearing_1fakeBJet_background",
      histogram_memLR_noSmearing_fakeLeadBJet_background, 
      histogram_memLR_noSmearing_fakeSubleadBJet_background, 
      histogram_memLR_noSmearing_2fakeBJets_background);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_noSmearing_2genuineBJets_background, "2 genuine b-jets",
      histogram_memLR_noSmearing_1fakeBJet_background, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.40, 0.72, 0.33, 0.15, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_2histograms_memLR_background.pdf");
    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_noSmearing_2genuineBJets_background, "2 genuine b-jets",
      histogram_memLR_noSmearing_fakeLeadBJet_background, "fake lead. b-jet",
      histogram_memLR_noSmearing_fakeSubleadBJet_background, "fake sublead. b-jet",
      histogram_memLR_noSmearing_2fakeBJets_background, "2 fake b-jets",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.045, 0.40, 0.69, 0.28, 0.23, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_4histograms_memLR_background.pdf");

    TGraph* graph_ROC_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_noSmearing_2genuineBJets",
      histogram_memLR_noSmearing_2genuineBJets_signal, histogram_memLR_noSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_noSmearing_fakeLeadBJet_logScale = compGraphROC("graph_ROC_noSmearing_fakeLeadBJet",
      histogram_memLR_noSmearing_fakeLeadBJet_signal, histogram_memLR_noSmearing_fakeLeadBJet_background, true);
    TGraph* graph_ROC_noSmearing_fakeSubleadBJet_logScale = compGraphROC("graph_ROC_noSmearing_fakeSubleadBJet",
      histogram_memLR_noSmearing_fakeSubleadBJet_signal, histogram_memLR_noSmearing_fakeSubleadBJet_background, true);
    TGraph* graph_ROC_noSmearing_2fakeBJets_logScale = compGraphROC("graph_ROC_noSmearing_2fakeBJets",
      histogram_memLR_noSmearing_2fakeBJets_signal, histogram_memLR_noSmearing_2fakeBJets_background, true);

    TGraph* graph_ROC_noSmearing_1fakeBJet_logScale = compGraphROC("graph_ROC_noSmearing_",
      histogram_memLR_noSmearing_1fakeBJet_signal, histogram_memLR_noSmearing_1fakeBJet_background, true);

    showGraphs(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY,
      graph_ROC_noSmearing_2genuineBJets_logScale, "2 genuine b-jets",
      graph_ROC_noSmearing_1fakeBJet_logScale, "#geq 1 fake b-jet",
      nullptr, "",
      nullptr, "",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      0.055, 0.23, 0.72, 0.33, 0.15, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
      true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
      "hh_bbwwMEM_dilepton_effectOfFakes_2graphs_ROC.pdf");
    showGraphs(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY,
      graph_ROC_noSmearing_2genuineBJets_logScale, "2 genuine b-jets",
      graph_ROC_noSmearing_fakeLeadBJet_logScale, "fake lead. b-jet",
      graph_ROC_noSmearing_fakeSubleadBJet_logScale, "fake sublead. b-jet",
      graph_ROC_noSmearing_2fakeBJets_logScale, "2 fake b-jets",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      0.045, 0.23, 0.66, 0.33, 0.28, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
      true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
      "hh_bbwwMEM_dilepton_effectOfFakes_4graphs_ROC.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_noSmearing_fakeBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kSignal_lo, histogramNames[kLR]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_missingBJet_noSmearing_genuineBJet_signal, "genuine b-jet",
      histogram_memLR_missingBJet_noSmearing_fakeBJet_signal, "fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.40, 0.72, 0.33, 0.15, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P_{m}", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP_{m}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_memLR_missingBJet_signal.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_noSmearing_fakeBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[true], kBackground_lo, histogramNames[kLR]);

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_missingBJet_noSmearing_genuineBJet_background, "genuine b-jet",
      histogram_memLR_missingBJet_noSmearing_fakeBJet_background, "fake b-jet",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.055, 0.40, 0.72, 0.33, 0.15, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P_{m}", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP_{m}", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_memLR_missingBJet_background.pdf");

    showHistograms(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY,
      histogram_memLR_missingBJet_noSmearing_genuineBJet_signal, "Signal",
      histogram_memLR_missingBJet_noSmearing_genuineBJet_background, "Background",
      nullptr, "",
      nullptr, "",
      showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
      showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
      0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
      "", 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset,
      true, yMin_memLR, yMax_memLR, "dN/dP", showHistograms_yAxisOffset,
      "hh_bbwwMEM_dilepton_effectOfFakes_memLR_missingBJet.pdf");

    TGraph* graph_ROC_missingBJet_noSmearing_genuineBJet_logScale = compGraphROC("graph_ROC_missingBJet_noSmearing_genuineBJet",
      histogram_memLR_missingBJet_noSmearing_genuineBJet_signal, histogram_memLR_missingBJet_noSmearing_genuineBJet_background, true);
    TGraph* graph_ROC_missingBJet_noSmearing_fakeBJet_logScale = compGraphROC("graph_ROC_missingBJet_noSmearing_fakeBJet",
      histogram_memLR_missingBJet_noSmearing_fakeBJet_signal, histogram_memLR_missingBJet_noSmearing_fakeBJet_background, true);

    showGraphs(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY,
      //graph_ROC_missingBJet_noSmearing_genuineBJet_logScale, "genuine b-jet",
      //graph_ROC_missingBJet_noSmearing_fakeBJet_logScale, "fake b-jet",
      graph_ROC_missingBJet_noSmearing_genuineBJet_logScale, "",
      nullptr, "",
      nullptr, "",
      nullptr, "",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      0.055, 0.23, 0.79, 0.33, 0.15, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
      true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
      "hh_bbwwMEM_dilepton_effectOfFakes_ROC_missingBJet.pdf");
    //-------------------------------------------------------------------------------------------------
  }

  if ( makePlots_effectOfSmearing ) {
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probS_noSmearing_2genuineBJets_signal, "MC truth",
      histogram_probS_jetSmearing_2genuineBJets_signal, "E_{b} smearing",
      histogram_probS_metSmearing_2genuineBJets_signal, "#rho smearing",
      //histogram_probS_jet_and_metSmearing_2genuineBJets_signal, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset_wRatio,
      true, yMin_probS_wRatio, yMax_probS, 1. - 0.29, 1. + 0.29, "dN/dlog w_{0}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probS_signal.pdf");

    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jetSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbSignal]);
    
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probS_missingBJet_noSmearing_genuineBJet_signal, "MC truth",
      histogram_probS_missingBJet_jetSmearing_genuineBJet_signal, "E_{b} smearing",
      histogram_probS_missingBJet_metSmearing_genuineBJet_signal, "#rho smearing",
      //histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_signal, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}^{m}", showHistograms_xAxisOffset_wRatio,
      true, yMin_probS_wRatio, yMax_probS, 1. - 0.29, 1. + 0.29, "dN/dlog w_{0}^{m}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probS_missingBJet_signal.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probB_noSmearing_2genuineBJets_signal, "MC truth",
      histogram_probB_jetSmearing_2genuineBJets_signal, "E_{b} smearing",
      histogram_probB_metSmearing_2genuineBJets_signal, "#rho smearing",
      //histogram_probB_jet_and_metSmearing_2genuineBJets_signal, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}",showHistograms_xAxisOffset_wRatio,
      true, yMin_probB_wRatio, yMax_probB, 1. - 0.29, 1. + 0.29, "dN/dlog w_{1}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probB_signal.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jetSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kProbBackground]);
    
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probB_missingBJet_noSmearing_genuineBJet_signal, "MC truth",
      histogram_probB_missingBJet_jetSmearing_genuineBJet_signal, "E_{b} smearing",
      histogram_probB_missingBJet_metSmearing_genuineBJet_signal, "#rho smearing",
      //histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_signal, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}^{m}",showHistograms_xAxisOffset_wRatio,
      true, yMin_probB_wRatio, yMax_probB, 1. - 0.29, 1. + 0.29, "dN/dlog w_{1}^{m}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probB_missingBJet_signal.pdf");
    //-------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probS_noSmearing_2genuineBJets_background, "MC truth",
      histogram_probS_jetSmearing_2genuineBJets_background, "E_{b} smearing",
      histogram_probS_metSmearing_2genuineBJets_background, "#rho smearing",
      //histogram_probS_jet_and_metSmearing_2genuineBJets_background, "E_{b} + #rho smearing",
      nullptr, "", 
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}", showHistograms_xAxisOffset_wRatio,
      true, yMin_probS_wRatio, yMax_probS, 1. - 0.29, 1. + 0.29, "dN/dlog w_{0}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probS_background.pdf");
    
    TH1* histogram_probS_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jetSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbSignal]);
    
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probS_missingBJet_noSmearing_genuineBJet_background, "MC truth",
      histogram_probS_missingBJet_jetSmearing_genuineBJet_background, "E_{b} smearing",
      histogram_probS_missingBJet_metSmearing_genuineBJet_background, "#rho smearing",
      //histogram_probS_missingBJet_jet_and_metSmearing_genuineBJet_background, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}^{m}", showHistograms_xAxisOffset_wRatio,
      true, yMin_probS_wRatio, yMax_probS, 1. - 0.29, 1. + 0.29, "dN/dlog w_{0}^{m}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probS_missingBJet_background.pdf");

    TH1* histogram_probB_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probB_noSmearing_2genuineBJets_background, "MC truth",
      histogram_probB_jetSmearing_2genuineBJets_background, "E_{b} smearing",
      histogram_probB_metSmearing_2genuineBJets_background, "#rho smearing",
      //histogram_probB_jet_and_metSmearing_2genuineBJets_background, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}", showHistograms_xAxisOffset_wRatio,
      true, yMin_probB_wRatio, yMax_probB, 1. - 0.29, 1. + 0.29, "dN/dlog w_{1}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probB_background.pdf");

    TH1* histogram_probB_missingBJet_noSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jetSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kProbBackground]);
    
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probB_missingBJet_noSmearing_genuineBJet_background, "MC truth",
      histogram_probB_missingBJet_jetSmearing_genuineBJet_background, "E_{b} smearing",
      histogram_probB_missingBJet_metSmearing_genuineBJet_background, "#rho smearing",
      //histogram_probB_missingBJet_jet_and_metSmearing_genuineBJet_background, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.23, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}^{m}", showHistograms_xAxisOffset_wRatio,
      true, yMin_probB_wRatio, yMax_probB, 1. - 0.29, 1. + 0.29, "dN/dlog w_{1}^{m}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_probB_missingBJet_background.pdf");
    //-------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    TH1* histogram_memLR_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_memLR_noSmearing_2genuineBJets_signal, "MC truth",
      histogram_memLR_jetSmearing_2genuineBJets_signal, "E_{b} smearing",
      histogram_memLR_metSmearing_2genuineBJets_signal, "#rho smearing",
      //histogram_memLR_jet_and_metSmearing_2genuineBJets_signal, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.37, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.37, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset_wRatio,
      true, yMin_memLR_wRatio, yMax_memLR, 1. - 0.29, 1. + 0.29, "dN/dP", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_signal.pdf");

    TH1* histogram_memLR_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_memLR_noSmearing_2genuineBJets_background, "MC truth",
      histogram_memLR_jetSmearing_2genuineBJets_background, "E_{b} smearing",
      histogram_memLR_metSmearing_2genuineBJets_background, "#rho smearing",
      //histogram_memLR_jet_and_metSmearing_2genuineBJets_background, "E_{b} + #rho smearing",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.37, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.37, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset_wRatio,
      true, yMin_memLR_wRatio, yMax_memLR, 1. - 0.29, 1. + 0.29, "dN/dP", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_background.pdf");

    TGraph* graph_ROC_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_noSmearing_2genuineBJets",
      histogram_memLR_noSmearing_2genuineBJets_signal, histogram_memLR_noSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_jetSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_jetSmearing_2genuineBJets",
      histogram_memLR_jetSmearing_2genuineBJets_signal, histogram_memLR_jetSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_metSmearing_2genuineBJets",
      histogram_memLR_metSmearing_2genuineBJets_signal, histogram_memLR_metSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_jet_and_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_jet_and_metSmearing_2genuineBJets",
      histogram_memLR_jet_and_metSmearing_2genuineBJets_signal, histogram_memLR_jet_and_metSmearing_2genuineBJets_background, true);

    showGraphs_wRatio(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY_wRatio,
      graph_ROC_noSmearing_2genuineBJets_logScale, "MC truth",
      graph_ROC_jetSmearing_2genuineBJets_logScale, "E_{b} smearing",
      graph_ROC_metSmearing_2genuineBJets_logScale, "#rho smearing",
      //graph_ROC_jet_and_metSmearing_2genuineBJets_logScale, "E_{b} + #rho smearing",
      nullptr, "",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showGraphs_legendOptions,
      0.065, 0.23, 0.69, 0.52, 0.26, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset_wRatio,
      true, 2.1e-4, 9.9e0, 1. - 0.29, 1. + 0.29, "Background Rate", showGraphs_yAxisOffset_wRatio, 
      "hh_bbwwMEM_dilepton_effectOfSmearing_ROC.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jetSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_signal = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kSignal_lo, histogramNames[kLR]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_memLR_missingBJet_noSmearing_2genuineBJets_signal, "MC truth",
      histogram_memLR_missingBJet_jetSmearing_2genuineBJets_signal, "E_{b} smearing",
      histogram_memLR_missingBJet_metSmearing_2genuineBJets_signal, "#rho smearing",
      //histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_signal, "E_{b} + #rho smearing",
      nullptr, "", 
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.37, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.37, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P_{m}", showHistograms_xAxisOffset_wRatio,
      true, yMin_memLR_wRatio, yMax_memLR, 1. - 0.29, 1. + 0.29, "dN/dP_{m}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_missingBJet_signal.pdf");

    TH1* histogram_memLR_missingBJet_noSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jetSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][false], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[false][true], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_background = loadHistogram(inputFile, 
      directories_part1[true][true], directories_part2_missingBJet[false], kBackground_lo, histogramNames[kLR]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_memLR_missingBJet_noSmearing_2genuineBJets_background, "MC truth",
      histogram_memLR_missingBJet_jetSmearing_2genuineBJets_background, "E_{b} smearing",
      histogram_memLR_missingBJet_metSmearing_2genuineBJets_background, "#rho smearing",
      //histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_background, "E_{b} + #rho smearing",
      nullptr, "", 
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      //0.054, 0.37, 0.63, 0.34, 0.28, showHistograms_legendOptions,
      0.065, 0.37, 0.64, 0.52, 0.26, showHistograms_legendOptions,
      labelText_background, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P_{m}", showHistograms_xAxisOffset_wRatio,
      true, yMin_memLR_wRatio, yMax_memLR, 1. - 0.29, 1. + 0.29, "dN/dP_{m}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_effectOfSmearing_memLR_missingBJet_background.pdf");

    TGraph* graph_ROC_missingBJet_noSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_noSmearing_2genuineBJets",
      histogram_memLR_missingBJet_noSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_noSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_missingBJet_jetSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_jetSmearing_2genuineBJets",
      histogram_memLR_missingBJet_jetSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_jetSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_missingBJet_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_metSmearing_2genuineBJets",
      histogram_memLR_missingBJet_metSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_metSmearing_2genuineBJets_background, true);
    TGraph* graph_ROC_missingBJet_jet_and_metSmearing_2genuineBJets_logScale = compGraphROC("graph_ROC_missingBJet_jet_and_metSmearing_2genuineBJets",
      histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_signal, histogram_memLR_missingBJet_jet_and_metSmearing_2genuineBJets_background, true);

    showGraphs_wRatio(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY_wRatio,
      graph_ROC_missingBJet_noSmearing_2genuineBJets_logScale, "MC truth",
      graph_ROC_missingBJet_jetSmearing_2genuineBJets_logScale, "E_{b} smearing",
      graph_ROC_missingBJet_metSmearing_2genuineBJets_logScale, "#rho smearing",
      //graph_ROC_missingBJet_jet_and_metSmearing_2genuineBJets_logScale, "E_{b} + #rho smearing",
      nullptr, "",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      //0.054, 0.23, 0.63, 0.34, 0.28, showGraphs_legendOptions,
      0.065, 0.23, 0.69, 0.52, 0.26, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset_wRatio,
      true, 2.1e-4, 9.9e0, 1. - 0.49, 1. + 0.49, "Background Rate", showGraphs_yAxisOffset_wRatio, 
      "hh_bbwwMEM_dilepton_effectOfSmearing_ROC_missingBJet.pdf");
    //-------------------------------------------------------------------------------------------------
  }

  if ( makePlots_effectOfHigherOrders ) {
    //-------------------------------------------------------------------------------------------------
    TH1* histogram_probS_lo_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_nlo_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_nlo, histogramNames[kProbSignal]);
    TH1* histogram_probS_lo_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbSignal]);
    TH1* histogram_probS_nlo_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_nlo, histogramNames[kProbSignal]);

    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probS_lo_signal, "LO",
      histogram_probS_nlo_signal, "NLO",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.060, 0.26, 0.72, 0.27, 0.18, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}",showHistograms_xAxisOffset_wRatio,
      true, yMin_probS_wRatio, yMax_probS, 1. - 0.29, 1. + 0.29, "dN/dlog w_{0}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_lo_vs_nlo_probS_signal.pdf");
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probS_lo_background, "LO",
      histogram_probS_nlo_background, "NLO",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.060, 0.26, 0.72, 0.27, 0.18, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probS, xMin_probS, xMax_probS, "log w_{0}",showHistograms_xAxisOffset_wRatio,
      true, yMin_probS_wRatio, yMax_probS, 1. - 0.29, 1. + 0.29, "dN/dlog w_{0}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_lo_vs_nlo_probS_background.pdf");

    TH1* histogram_probB_lo_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_nlo_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_nlo, histogramNames[kProbBackground]);
    TH1* histogram_probB_lo_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kProbBackground]);
    TH1* histogram_probB_nlo_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_nlo, histogramNames[kProbBackground]);

    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probB_lo_signal, "LO",
      histogram_probB_nlo_signal, "NLO",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.060, 0.26, 0.72, 0.27, 0.18, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}",showHistograms_xAxisOffset_wRatio,
      true, yMin_probB_wRatio, yMax_probB, 1. - 0.29, 1. + 0.29, "dN/dlog w_{1}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_lo_vs_nlo_probB_signal.pdf");
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_probB_lo_background, "LO",
      histogram_probB_nlo_background, "NLO",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.060, 0.26, 0.72, 0.27, 0.18, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_probB, xMin_probB, xMax_probB, "log w_{1}",showHistograms_xAxisOffset_wRatio,
      true, yMin_probB_wRatio, yMax_probB, 1. - 0.29, 1. + 0.29, "dN/dlog w_{1}", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_lo_vs_nlo_probB_background.pdf");

    TH1* histogram_memLR_lo_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_lo, histogramNames[kLR]);
    TH1* histogram_memLR_nlo_signal = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kSignal_nlo, histogramNames[kLR]);
    TH1* histogram_memLR_lo_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_lo, histogramNames[kLR]);
    TH1* histogram_memLR_nlo_background = loadHistogram(inputFile, 
      directories_part1[false][false], directories_part2[false][false], kBackground_nlo, histogramNames[kLR]);
  
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_memLR_lo_signal, "LO",
      histogram_memLR_nlo_signal, "NLO",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.060, 0.23, 0.72, 0.27, 0.18, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset_wRatio,
      true, yMin_memLR_wRatio, yMax_memLR, 1. - 0.29, 1. + 0.29, "dN/dP", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_lo_vs_nlo_memLR_signal.pdf");
    showHistograms_wRatio(
      showHistograms_canvasSizeX, showHistograms_canvasSizeY_wRatio,
      histogram_memLR_lo_background, "LO",
      histogram_memLR_nlo_background, "NLO",
      nullptr, "",
      nullptr, "",
      showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
      showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
      0.060, 0.23, 0.72, 0.27, 0.18, showHistograms_legendOptions,
      labelText_signal, 0.055,
      0.1800, 0.9525, 0.2900, 0.0900,
      numBinsX_memLR, xMin_memLR, xMax_memLR, "P", showHistograms_xAxisOffset_wRatio,
      true, yMin_memLR_wRatio, yMax_memLR, 1. - 0.29, 1. + 0.29, "dN/dP", showHistograms_yAxisOffset_wRatio,
      "hh_bbwwMEM_dilepton_lo_vs_nlo_memLR_background.pdf");

    TGraph* graph_ROC_lo_logScale = compGraphROC("graph_ROC_lo",
      histogram_memLR_lo_signal, histogram_memLR_lo_background, true);
    TGraph* graph_ROC_nlo_logScale = compGraphROC("graph_ROC_nlo",
      histogram_memLR_nlo_signal, histogram_memLR_nlo_background, true);

    //showGraphs_wRatio(
    //  showGraphs_canvasSizeX, showGraphs_canvasSizeY_wRatio,
    //  graph_ROC_lo_logScale, "LO",
    //  graph_ROC_nlo_logScale, "NLO",
    //  nullptr, "",
    //  nullptr, "",
    //  showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
    //  showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
    //  0.054, 0.23, 0.78, 0.21, 0.15, showGraphs_legendOptions,
    //  labelText_signal_vs_background, 0.040,
    //  0.1600, 0.9525, 0.2900, 0.0600,
    //  10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset_wRatio,
    //  true, 2.1e-4, 9.9e0, 1. - 0.49, 1. + 0.49, "Background Rate", showGraphs_yAxisOffset_wRatio, 
    //  "hh_bbwwMEM_dilepton_lo_vs_nlo_ROC.pdf");
    showGraphs(
      showGraphs_canvasSizeX, showGraphs_canvasSizeY,
      graph_ROC_lo_logScale, "LO",
      graph_ROC_nlo_logScale, "NLO",
      nullptr, "",
      nullptr, "",
      showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
      showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
      0.055, 0.23, 0.79, 0.33, 0.15, showGraphs_legendOptions,
      labelText_signal_vs_background, 0.040,
      0.1600, 0.9525, 0.2900, 0.0600,
      10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
      true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
      "hh_bbwwMEM_dilepton_lo_vs_nlo_ROC.pdf");
    //-------------------------------------------------------------------------------------------------
  }

  delete inputFile;
}
