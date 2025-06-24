
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
#include "TH1D.h"
#include "TGraphErrors.h"
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
    double one_ZDC_percentage;
    double ZDC_background;
    double noZDC_background;
    double min_bias;
    double vertex_cut;
    double vertex;
    int fillpattern;;
    double vertex_sigma;
    double scale_factor;
    double charge_sum_mean_ns;
    double charge_sum_mean[4];
    double glauber_mu[4];
    double glauber_k[4];
    double glauber_chi2[4];
    double glauber_npart[4];
    double glauber_ecc2[4];
    double glauber_ecc3[4];
    double glauber_b[4];
    double glauber_trig_eff[4];
    double glauber_trig_eff_err[4];
    double centrality_check_chi2[4];
    double ecc2;
    double ecc3;

  };

  void QA_ReferenceRun(int runnumber);
  void QA_RunCentralityCheck(int runnumber, int scaled = 1);
  void QA_MC(std::string mc_generator);
  void Start_QA_Centrality(int runnumber);
  //  double NBD_getValue(int n, double mu, double k);
  //double NBD_getValue(double n, double mu, double k);
  //double NBDGlauberConv(double *x, double *par);

  void QA_MBDChannels(const int runnumber);
  void QA_ZDCCheck(const int runnumber);
  void QA_MakeChargeSum(const int runnumber);
  void QA_MakeDivisions(const int runnumber, const int doVertexScaled = 1);
  void QA_MakeCentralityCalibrations(const int runnumber, const bool doVertexScaled = false, const bool use_shifted = false);
  void QA_CentralityCheck(const int runnumber, const int scaled = 1);

  void Print_QA_Info(bool use_table = false);

  void SetRunNumber(int run) { qa_info.runnumber = run;}
  void SetNEvents(int events) { nevents = events;}
  void SetReferenceRun(int ref_run) { reference_run = ref_run ; }
  void SetChargeCut(float cut) { cthresh = cut;}
  void SetVertexCut(float cut) { z_cut = cut;}
  void SetCountBefore(bool count) { countbefore = count; }
  void SetNDivs(int nd) { ndivs = nd;}
  void SetDivs(int nd) { divs = nd;}
  void SetZDCCut(int cut) { zdc_cut = cut;}
  void SetMbdNSCut(int south_cut, int north_cut) { charge_sum_south_cut = south_cut; charge_sum_north_cut = north_cut; }
  void SetTrigEffMUK(double te, double mu, double k ) { trigeff = te; set_mu = mu; set_k = k; }
  void setForceZDC(bool a) { forceZDC = a; }
  void setSysName(const std::string sys) {
    sysname = sys;
    systematics = true;
  }
  void setNTupleFile(const std::string f) { ntuple_file = f; }
  void setNTupleName(const std::string f) { ntuple_name = f; }
  void setHistoFile(const std::string f) { histo_file = f; }
 protected:

  int loadTree();
  TFile *m_file = nullptr;
  TTree *m_ttree = nullptr;
  std::string sysname = "";
  bool systematics = false;
  struct QA_Info qa_info;

  bool countbefore = true;
  char *env_p = nullptr;
  char *env_tree = nullptr;
  int silence = 0;
  int debug = 0;
  bool forceZDC = false;
  bool hasZDC = true;
  int reference_run = 0;
  bool isSim = false;
  int ndivs = 20;
  int divs = 95;
  float m_maxsumcut = 2100;
  int nevents = 0;
  double cthresh = 0.5;
  double z_cut = 60.;
  double charge_sum_north_cut = 10.;
  double charge_sum_south_cut = 150.;
  double zdc_cut = 60;
  double trigeff = .94;
  double set_mu = 0.0;
  double set_k = 0.0;
  int mbd_hit_cut = 2;

  std::string histo_file = "lime_auau_hists_nominal.root";
  std::string ntuple_file = "glau_auau_ntuple_nominal.root";
  std::string ntuple_name = "nt_Au_Au";
  std::map<int, std::string> mc_map{};

  float zdc_sum[2];
  float mbd_sum[2];
  float zdc_energy[6];
  float mbd_charge[128];
  float mbd_charge_sum[2];
  float mbd_time[128];
  int minbias;
  float mbd_vertex;
  float mbd_time_zero;
  ULong64_t bunch_number;
  ULong64_t gl1_scaled;
  ULong64_t gl1_live;
};
#endif
