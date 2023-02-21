namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes

  TString spectra_file = "/home/gs627/LEEana/merge.root";
  TString flux_Xs_directory = "/home/gs627/LEEana/hist_rootfiles/XsFlux/";
  TString detector_directory = "/home/gs627/LEEana/hist_rootfiles/DetVar/";
  TString mc_directory = "/home/gs627/LEEana/mc_stat/";


  int channels_observation = 3;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse

  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 17;
  //int syst_cov_flux_Xs_end   = 19;//cov_18.root is uncorrelated reweighting and cov_19.root is correlated
 
  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  int syst_cov_mc_stat_end   = 1;
   

  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files
 
  
  ////////// display graphics flag

  bool flag_display_graphics = 1;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs          = 1;
  bool flag_syst_detector         = 1;  // ==0 No DET SYS (turn off)
  bool flag_syst_additional       = 1;
  bool flag_syst_mc_stat          = 1;
  bool flag_syst_reweight         = 0;  // Must turn to 1 for reweighting to be applied
  bool flag_syst_reweight_cor     = 0;  // Must turn to 1 for reweighting to be applied

  double Lee_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 1;
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF         = 0;

  bool flag_GoF_output2file_default_0 = 1;

  bool flag_lookelsewhere             = 1;
  
  bool flag_both_numuCC            = 0;// 1
  bool flag_CCpi0_FC_by_numuCC     = 0;// 2
  bool flag_CCpi0_PC_by_numuCC     = 0;// 3
  bool flag_NCpi0_by_numuCC        = 0;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = 0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = 0;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = 0;// 7
  bool flag_nueCC_FC_by_all        = 0;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = 0;

  ////////// MicroBooNE suggested

  bool flag_chi2_data_H0 = 0;
  bool flag_dchi2_H0toH1 = 0;
  
  ////////// Advanced tools
  
  ///// void TLee::Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed)
  ///// do the fitting on the spectra and cov_total after constraint ?
  bool flag_Lee_minimization_after_constraint = 0;// hardcoded, only for the standard 7-ch fitting
  
}
