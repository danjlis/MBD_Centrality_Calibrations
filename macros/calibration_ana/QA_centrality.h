
/* Header for the QA_centrality class */

#ifndef QA_CENTRALITY_H
#define QA_CENTRALITY_H

#include <string>
#include <string.h>
#include <iostream>
#include <filesystem>
#include <vector>

#include "TF1.h"
#include "TPad.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"

class QA_centrality{

 public:

  QA_centrality(int silent = 0, int debugger = 0);

  ~QA_centrality() {};

  struct QA_Info
  {
    int runnumber;
    int quality;
    int nevents;
    double ZDC_percentage;
    double min_bias;
    double vertex_cut;
    double vertex;
    double vertex_sigma;
    double scale_factor;
    double charge_sum_mean_ns;
    double charge_sum_mean[4];
    double glauber_mu[2];
    double glauber_k[2];
    double glauber_chi2[2];
    double glauber_npart[2];
    double glauber_trig_eff[2];
    double glauber_trig_eff_err[2];
    double glauber_chi2_forced[2];
    double glauber_npart_forced[2];
    double glauber_trig_eff_forced[2];
    double glauber_trig_eff_forced_err[2];
    double centrality_divs_chi2[4];

  };

  void QA_ReferenceRun(int runnumber);
  void Start_QA_Centrality(int runnumber);
  double NBD_getValue(int n, double mu, double k);

  void QA_MBDChannels(const int runnumber);
  void QA_ZDCCheck(const int runnumber);
  void QA_MakeChargeSum(const int runnumber);
  void QA_MakeCentralityCalibrations(const int runnumber, const bool doVertexSelection = 1, const bool use_shifted = false, const bool use_balanced = false);
  void QA_CentralityCheck(const int runnumber);

  void Print_QA_Info(bool use_table = false);

  void SetRunNumber(int run) { qa_info.runnumber = run;}
  void SetReferenceRun(int ref_run) { reference_run = ref_run ; }
  void SetChargeCut(int cut) { cthresh = cut;}
  void SetVertexCut(int cut) { z_cut = cut;}
  void SetZDCCut(int cut) { zdc_cut = cut;}
  void SetMbdNSCut(int south_cut, int north_cut) { charge_sum_south_cut = south_cut; charge_sum_north_cut = north_cut; }
  void SetTrigEffMUK(double te, double mu, double k ) { trigeff = te; set_mu = mu; set_k = k; }

 protected:

  struct QA_Info qa_info;
  char *env_p = nullptr;
  int silence = 0;
  int debug = 0;

  bool hasZDC = true;
  int reference_run = 23696;
  double cthresh = 0.25;
  double z_cut = 60.;
  double charge_sum_north_cut = 10.;
  double charge_sum_south_cut = 150.;
  double zdc_cut = 40;
  double trigeff = .91;
  double set_mu = 0.0;
  double set_k = 0.0;
  int mbd_hit_cut = 2;
};
#endif
