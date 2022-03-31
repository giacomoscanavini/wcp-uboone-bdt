#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TCut.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"
#include "pfeval.h"

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

namespace LEEana{
  // this is for the real data, for fake data this should be 1 ...
  double em_charge_scale = 0.95;
  //double em_charge_scale = 1.0;

  // correct reco neutrino energy ...
  double get_reco_Enu_corr(KineInfo& kine, bool flag_data);
  
  double get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name="kine_reco_Enu");
  
  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  double get_weight(TString weight_name, EvalInfo& eval);

  int get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  
  // generic neutrino cuts
  // TCut generic_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >15";
  bool is_generic(EvalInfo& info);
  
  // preselection cuts
  // TCut preselect_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length > 0";
  bool is_preselection(EvalInfo& info); 
  
  // nueCC cuts
  // TCut nueCC_cut = "numu_cc_flag >=0 && nue_score > 7.0";
  bool is_nueCC(TaggerInfo& tagger_info);
  bool is_loosenueCC(TaggerInfo& tagger_info);
  
  bool is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  
  // numuCC cuts
  // TCut numuCC_cut = "numu_cc_flag >=0 && numu_score > 0.9";
  bool is_numuCC(TaggerInfo& tagger_info);
  bool is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data);
  bool is_numuCC_cutbased(TaggerInfo& tagger_info);
  
  // pio cuts (with and without vertex)
  // TCut pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_pi0(KineInfo& kine, bool flag_data);
  
  // must be with vertex ...
  // TCut cc_pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_cc_pi0(KineInfo& kine, bool flag_data);
  
  // NC cuts
  // TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
  bool is_NC(TaggerInfo& tagger_info);
  bool is_NCpio_bdt(TaggerInfo& tagger_info);
  bool is_NCdelta_bdt(TaggerInfo& tagger_info);

  // TCut FC_cut = "match_isFC==1";
  // TCut PC_cut = "match_isFC==0";  
  bool is_FC(EvalInfo& eval);
  bool is_filter_shower(KineInfo& kine);
  
  // TCut truth_nueCC_inside = "abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1";
  // TCut truth_numuCC_inside = "abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1";
  // TCut truth_NCpi0_inside = "truth_isCC==0 && truth_NprimPio>0 && truth_vtxInside==1";
  bool is_truth_nueCC_inside(EvalInfo& eval);
  bool is_truth_numuCC_inside(EvalInfo& eval);
  bool is_truth_NCpi0_inside(EvalInfo& eval, PFevalInfo& pfeval);

  int mcc8_pmuon_costheta_bin(float pmuon, float costh);
}


double LEEana::get_reco_Enu_corr(KineInfo& kine, bool flag_data){
  double reco_Enu_corr = 0;
  if (kine.kine_reco_Enu > 0){
    if (flag_data){
      for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
    if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
      reco_Enu_corr +=  kine.kine_energy_particle->at(j) * em_charge_scale;
    }else{
      reco_Enu_corr +=  kine.kine_energy_particle->at(j);
    }
    //  std::cout << "p: " << kine.kine_energy_particle->at(j) << " " << kine.kine_energy_info->at(j) << " " << kine.kine_particle_type->at(j) << " " << kine.kine_energy_included->at(j) << std::endl;
      }
      reco_Enu_corr += kine.kine_reco_add_energy;
      return reco_Enu_corr;
    }
  }
  return kine.kine_reco_Enu;
}

double LEEana::get_weight(TString weight_name, EvalInfo& eval){

  double addtl_weight = 1.0;

  // CV correction from numuCC cross section data
  // if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1){
  //   if (eval.truth_nuEnergy>200 && eval.truth_nuEnergy<=540) addtl_weight = 1.28043; 
  //   else if (eval.truth_nuEnergy>540 && eval.truth_nuEnergy<=705) addtl_weight = 1.21158;
  //   else if (eval.truth_nuEnergy>705 && eval.truth_nuEnergy<=805) addtl_weight = 1.19091;
  //   else if (eval.truth_nuEnergy>805 && eval.truth_nuEnergy<=920) addtl_weight = 1.17733;
  //   else if (eval.truth_nuEnergy>920 && eval.truth_nuEnergy<=1050) addtl_weight = 1.13983;
  //   else if (eval.truth_nuEnergy>1050 && eval.truth_nuEnergy<=1200) addtl_weight = 1.07864;
  //   else if (eval.truth_nuEnergy>1200 && eval.truth_nuEnergy<=1375) addtl_weight = 1.00722;
  //   else if (eval.truth_nuEnergy>1375 && eval.truth_nuEnergy<=1570) addtl_weight = 0.93857;
  //   else if (eval.truth_nuEnergy>1570 && eval.truth_nuEnergy<=2050) addtl_weight = 0.886241;
  //   else if (eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) addtl_weight = 0.858724;
  //   else if (eval.truth_nuEnergy>4000) addtl_weight = 0.858724;
  // }
  // std::cout << "energy: " << eval.truth_nuEnergy << " addtl_weight: " << addtl_weight << std::endl;
  // end of data correction

  if (weight_name == "cv_spline"){
    return addtl_weight*eval.weight_cv * eval.weight_spline;
  }else if (weight_name == "cv_spline_cv_spline"){
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "unity" || weight_name == "unity_unity"){
    return 1;
  }else if (weight_name == "lee_cv_spline"){
    return (eval.weight_lee * addtl_weight*eval.weight_cv * eval.weight_spline);
  }else if (weight_name == "lee_cv_spline_lee_cv_spline"){
    return pow(eval.weight_lee * addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "lee_cv_spline_cv_spline" || weight_name == "cv_spline_lee_cv_spline"){
    return eval.weight_lee * pow(addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "spline"){
    return eval.weight_spline;
  }else if (weight_name == "spline_spline"){
    return pow(eval.weight_spline,2);
  }else if (weight_name == "lee_spline"){
    return (eval.weight_lee * eval.weight_spline);
  }else if (weight_name == "lee_spline_lee_spline"){
    return pow(eval.weight_lee * eval.weight_spline,2);
  }else if (weight_name == "lee_spline_spline" || weight_name == "spline_lee_spline"){
    return eval.weight_lee * pow( eval.weight_spline,2);
  }else{
    std::cout <<"Unknown weights: " << weight_name << std::endl;
  }
     
  return 1;
}

double LEEana::get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data , TString var_name){
  
  
  if (var_name == "truth_nuEnergy"){      return eval.truth_nuEnergy;

  }else if (var_name == "kine_reco_Enu"){ 
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){ return get_reco_Enu_corr(kine, flag_data);
    }else{ return -1000.; }
  }else if (var_name == "kine_pio_energy_1"){   return kine.kine_pio_energy_1;
  }else if (var_name == "kine_pio_energy_2"){   return kine.kine_pio_energy_2;
  }else if (var_name == "kine_pio_theta_1"){    return kine.kine_pio_theta_1;
  }else if (var_name == "kine_pio_theta_2"){    return kine.kine_pio_theta_2;
  }else if (var_name == "kine_pio_costheta_1"){ return kine.kine_pio_theta_1;
  }else if (var_name == "kine_pio_costheta_2"){ return kine.kine_pio_theta_2;
  }else if (var_name == "kine_pio_phi_1"){      return kine.kine_pio_phi_1/180.*3.1415926;
  }else if (var_name == "kine_pio_phi_2"){      return kine.kine_pio_phi_2/180.*3.1415926;

  }else if (var_name == "kine_pio_energy_high"){ 
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      if(kine.kine_pio_energy_1 >= kine.kine_pio_energy_2) return kine.kine_pio_energy_1;
      else return kine.kine_pio_energy_2;
    }else return -1000.;
  }else if (var_name == "kine_pio_energy_low"){ 
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      if(kine.kine_pio_energy_1 <= kine.kine_pio_energy_2) return kine.kine_pio_energy_1;
      else return kine.kine_pio_energy_2;
    }else return -1000.;
  }else if (var_name == "kine_pio_costheta_high"){    
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      if(kine.kine_pio_energy_1 >= kine.kine_pio_energy_2) return kine.kine_pio_theta_1;
      else return kine.kine_pio_theta_2;
    }else return -1000.;
  }else if (var_name == "kine_pio_costheta_low"){    
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      if(kine.kine_pio_energy_1 <= kine.kine_pio_energy_2) return kine.kine_pio_theta_1;
      else return kine.kine_pio_theta_2;
    }else return -1000.;
  }else if (var_name == "kine_pio_phi_high"){      
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      if(kine.kine_pio_energy_1 >= kine.kine_pio_energy_2) return kine.kine_pio_phi_1/180.*3.1415926;
      else return kine.kine_pio_phi_2/180.*3.1415926;
    }else return -1000.;
  }else if (var_name == "kine_pio_phi_low"){      
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      if(kine.kine_pio_energy_1 <= kine.kine_pio_energy_2) return kine.kine_pio_phi_1/180.*3.1415926;
      else return kine.kine_pio_phi_2/180.*3.1415926;
    }else return -1000.;


  }else if (var_name == "pi0_mass"){
    if (kine.kine_pio_mass > 0){
      //TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      //TLorentzVector pio = p1 + p2;
      if (flag_data){
        if (kine.kine_pio_mass * em_charge_scale > 12.){return kine.kine_pio_mass * em_charge_scale;}
        else{ return -1000.; }
      }else{
        if (kine.kine_pio_mass > 12.){return kine.kine_pio_mass;}
        else{ return -1000.; }
      }
    }else{ return -1000.; }

  }else if (var_name == "pi0_costheta"){
    TLorentzVector p1(kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
    TLorentzVector p2(kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
    TLorentzVector pio = p1 + p2;
    return pio.CosTheta();
  
  }else if (var_name == "pi0_phi"){
    TLorentzVector p1(kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
    TLorentzVector p2(kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
    TLorentzVector pio = p1 + p2;
    return pio.Phi();
  
  }else if (var_name == "pi0_energy"){
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      double pi0_mass = 135;
      double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
      return pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
    }else{ return -1000.; }
  
  }else if (var_name == "pi0_total_energy"){
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      double pi0_mass = 135;
      double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
      return pi0_mass + pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
    }else{ return -1000.; }

  }else if (var_name = "pi0_momentum"){
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      double pi0_mass = 135;
      double alpha = (kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
      double pi0_energy = pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
      double pi0_total_energy = pi0_mass + pi0_energy;
      return sqrt(pi0_total_energy*pi0_total_energy - pi0_mass*pi0_mass);
    }else{ return -1000.; }

  }else if (var_name = "pi0_costheta_CM"){
    if(kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){
      double pi0_mass = 135;
      double alpha = (kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
      double pi0_energy = pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
      double pi0_total_energy = pi0_mass + pi0_energy;
      double pi0_momentum = sqrt(pi0_total_energy*pi0_total_energy - pi0_mass*pi0_mass);
      return pi0_total_energy/fabs(pi0_momentum) * alpha;
    }
    else{ return -1000.; }

  }else if (var_name == "numu_cc_3_max_muon_length"){ return tagger.numu_cc_3_max_muon_length;
  }else if (var_name == "lem_shower_num_segs"){       return tagger.lem_shower_num_segs;
  }else if (var_name == "kine_pio_energy_1"){         return kine.kine_pio_energy_1;
  }else if (var_name == "br1_2_max_length"){          return tagger.br1_2_max_length;
  }else if (var_name == "brm_acc_direct_length"){     return tagger.brm_acc_direct_length;
  }else if (var_name == "kine_pio_vtx_dis"){          return kine.kine_pio_vtx_dis;
  }else if (var_name == "anc_flag_main_outside"){     return tagger.anc_flag_main_outside;
  }else if (var_name == "stw_1_dis"){                 return tagger.stw_1_dis;
  }else if (var_name == "kine_pio_dis_1"){            return tagger.kine_pio_dis_1;
  }else if (var_name == "numu_cc_3_track_length"){    return tagger.numu_cc_3_max_length;
  }else if (var_name == "brm_connected_length"){      return tagger.brm_connected_length;
  }else if (var_name == "nue_score"){                 return (tagger.nue_score<=15.99?tagger.nue_score:15.99);
  }else if (var_name == "nc_pio_score"){              return tagger.nc_pio_score;
  }else if (var_name == "nc_delta_score"){            return tagger.nc_delta_score;
  }else if (var_name == "numu_score"){                return tagger.numu_score;
  }else if (var_name == "shower_energy"){             return tagger.mip_energy;
  }else if (var_name == "shower_angle_beam"){         return tagger.mip_angle_beam;
  }else if (var_name == "shower_angle_vertical"){     return tagger.spt_angle_vertical;
  }else if (var_name == "shwvtx_nuvtx_dis"){          return sqrt(pow(pfeval.reco_nuvtxX-pfeval.reco_showervtxX,2)+pow(pfeval.reco_nuvtxY-pfeval.reco_showervtxY,2)+pow(pfeval.reco_nuvtxZ-pfeval.reco_showervtxZ,2)); 
  }else if (var_name == "median_dQdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    return vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
  }else if (var_name == "median_dEdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    return median_dedx; // MeV/cm
  }else if (var_name == "reco_showervtxX"){       return pfeval.reco_showervtxX;
  }else if (var_name == "reco_nuvtxX"){           return pfeval.reco_nuvtxX;
  }else if (var_name == "reco_nuvtxY"){           return pfeval.reco_nuvtxY;
  }else if (var_name == "reco_nuvtxZ"){           return pfeval.reco_nuvtxZ;
  }else if (var_name == "reco_nuvtxU"){           return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) - pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "reco_nuvtxV"){           return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) + pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "mip_quality_n_tracks"){  return tagger.mip_quality_n_tracks;
  }else if (var_name == "mip_quality_n_showers"){ return tagger.mip_quality_n_showers;
  }else if (var_name == "gap_n_bad"){             return tagger.gap_n_bad;
  }else if (var_name == "muon_KE"){               return pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
  }else if (var_name == "reco_Emuon"){            return pfeval.reco_muonMomentum[3]*1000; // GeV --> MeV  
  }else if (var_name == "muon_momentum"){
      float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
      return (TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66));
  }else if (var_name == "muon_costheta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
  return TMath::Cos(muonMomentum.Theta());
      else
  return -2;
  }else if (var_name == "muon_theta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
  return muonMomentum.Theta()*180./TMath::Pi();
      else
  return -1000;
  }else if (var_name == "muon_phi"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
  return muonMomentum.Phi()/TMath::Pi()*180.;
      else
  return -1000;
  }else if (var_name == "proton_KE"){
      return pfeval.reco_protonMomentum[3]*1000.-938.27; // GeV--> MeV
  }else if (var_name == "proton_theta"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Theta()/TMath::Pi()*180.;
  }else if (var_name == "proton_phi"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Phi()/TMath::Pi()*180.;
  }else if (var_name == "Ehadron"){
    if (pfeval.reco_muonMomentum[3]>0)
      return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    else
      return -1000;
    //  }else if (var_name == "Ehadron"){
      /* Float_t Ehadron = kine.kine_reco_Enu; */
      /* for(size_t i=0; i<kine.kine_energy_particle->size(); i++) */
      /* { */
      /*     int pdgcode = kine.kine_particle_type->at(i); */
      /*     if(abs(pdgcode)==13) Ehadron = Ehadron - kine.kine_energy_particle->at(i) - 105.658; */ 
      /*     //if(abs(pdgcode)==11) Ehadron = Ehadron - kine.kine_energy_particle->at(i); */ 
      /* } */
    // return kine.kine_reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
  }else if (var_name == "Q2"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
    //  }else if (var_name == "Q2"){
    // Float_t Enu = kine.kine_reco_Enu;
    //Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    //Float_t Ehadron = Enu - Emu;
    //Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    //TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    //Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    //return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
  }else if (var_name == "x_Bjorken"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
    //  }else if (var_name == "x_Bjorken"){
    // Float_t Enu = kine.kine_reco_Enu;
    // Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    // Float_t Ehadron = Enu - Emu;
    // Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    // TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    // Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    // return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
  }else if (var_name == "N_tracks"){
      int N_tracks = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==11) continue;
          if(kine.kine_energy_particle->at(i)<10) continue;
          if(abs(pdgcode)==13 || abs(pdgcode)==211){
            N_tracks += 1;
          }
          else if(kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_tracks += 1; 
          }
      }
      return N_tracks;
  }else if (var_name == "N_showers"){
      int N_showers = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)!=11) continue;
          if(kine.kine_energy_particle->at(i)>10) N_showers += 1;
      }
      return N_showers;
  }else if (var_name == "EhadShwrFrac"){
      double EhadShwr=0, EhadTot=0;

      if (pfeval.reco_muonMomentum[3]>0) {
        EhadTot = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
        for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
          if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
            EhadShwr +=  kine.kine_energy_particle->at(j);
          }

        }
        return EhadShwr/ EhadTot;

      }
      else return -1;
  }else if (var_name == "reco_mcc8_pmuoncosth_Enu"){
    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon momentum
    float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66;
    float pmuon = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66) / 1000.0;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());
    // flattened Enu
    int indx = mcc8_pmuon_costheta_bin(pmuon, costh); // index of pmuon-costh
    double reco_Enu = get_reco_Enu_corr(kine, flag_data) / 100.0;
    if (reco_Enu<0) return -10000;
    else if (reco_Enu>25.0) return 10000; // overflow bin
    else return (indx-1)*25.0 + reco_Enu;
  }
  else{
    std::cout << "No such variable: " << var_name << std::endl;
    exit(EXIT_FAILURE);
  }
  return -1;
}

int LEEana::get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  for (auto it = map_cut_xs_bin.begin(); it != map_cut_xs_bin.end(); it++){
    TString cut_name = it->first;
    int number = it->second;

    double Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
    double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

    double truth_pi0_momentum = -1000.;
    if(pfeval.truth_pio_energy_1 > 0. && pfeval.truth_pio_energy_2 > 0.){
      double pi0_mass = 135;
      double alpha = fabs(pfeval.truth_pio_energy_1 - pfeval.truth_pio_energy_2)/(pfeval.truth_pio_energy_1 + pfeval.truth_pio_energy_2);
      double pi0_energy = pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(pfeval.truth_pio_angle/180.*3.1415926)))-1);
      double pi0_total_energy = pi0_mass + pi0_energy;
      truth_pi0_momentum = sqrt(pi0_total_energy*pi0_total_energy - pi0_mass*pi0_mass);
    }

    //float costheta_binning[10] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};   //fine binning
    //float costheta_binning[5]  = {-1,         .27,      .62,      .86,      1};   //coarse binning
    float costheta_binning[3]    = {-1,                   .62,                1};   //very coarse binning
    TLorentzVector muonMomentum(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);

    if (cut_file == 1){
      if (cut_name == "numuCC.inside.Enu.le.300"){                if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400"){          if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500"){          if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400.gt.300"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 && eval.truth_nuEnergy>300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500.gt.400"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 && eval.truth_nuEnergy>400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.600.gt.500"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=600 && eval.truth_nuEnergy>500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.700.gt.600"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=700 && eval.truth_nuEnergy>600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.800.gt.700"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=800 && eval.truth_nuEnergy>700) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.900.gt.800"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=900 && eval.truth_nuEnergy>800) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1000.gt.900"){  if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1000 && eval.truth_nuEnergy>900) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1100.gt.1000"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1100 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1100"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1100) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1000"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1500.gt.1200"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1500 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2100.gt.1500"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2100 && eval.truth_nuEnergy>1500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1400.gt.1200"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1400 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1600.gt.1400"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1600 && eval.truth_nuEnergy>1400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2000.gt.1600"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2000 && eval.truth_nuEnergy>1600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2500.gt.2000"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2500 && eval.truth_nuEnergy>2000) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2500"){         if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2500) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2100"){         if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2100) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.1500"){         if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>1500) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 2) {
      if (cut_name == "numuCC.inside.Emuon.le.100"){                if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=100 && Emuon>0) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.200.gt.100"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=200 && Emuon>100) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.300.gt.200"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=300 && Emuon>200) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.400.gt.300"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=400 && Emuon>300) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.500.gt.400"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=500 && Emuon>400) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.600.gt.500"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=600 && Emuon>500) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.700.gt.600"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=700 && Emuon>600) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.800.gt.700"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=800 && Emuon>700) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.900.gt.800"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=900 && Emuon>800) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1000.gt.900"){  if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1000 && Emuon>900) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1200.gt.1000"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1200 && Emuon>1000) return number;
      }else if (cut_name == "numuCC.inside.Emuon.gt.1200"){          if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1200) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 3) {
      if (cut_name == "numuCC.inside.Ehadron.le.100"){                if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.200.gt.100"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=200 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.300.gt.200"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=300 && Ehadron>200) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.400.gt.300"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=400 && Ehadron>300) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.500.gt.400"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=500 && Ehadron>400) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.600.gt.500"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=600 && Ehadron>500) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.700.gt.600"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=700 && Ehadron>600) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.800.gt.700"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=800 && Ehadron>700) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.900.gt.800"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=900 && Ehadron>800) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1000.gt.900"){  if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1000 && Ehadron>900) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.gt.1000"){         if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1000) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }       
    }
    else if (cut_file == 4){
      if (cut_name == "NCpi0.inside.Enu.le.540.gt.200"){          if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.705.gt.540"){    if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.805.gt.705"){    if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=805 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.920.gt.805"){    if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>805) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.1050.gt.920"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=1050 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.1200.gt.1050"){  if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1050) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.1375.gt.1200"){  if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=1375 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.1570.gt.1375"){  if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1375) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.2050.gt.1570"){  if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "NCpi0.inside.Enu.le.4000.gt.2050"){  if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 5) {
      if (cut_name == "numuCC.inside.Enu.le.4000.gt.200"){ if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 6){
      if (cut_name == "nueCC.inside.Enu.le.1200.gt.200"){         if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.4000.gt.1200"){  if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>1200) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    //very coarse angle binning
    else if (cut_file == 7){
      if (number==-1) { std::cout << "cut_name, number = " << cut_name << ", " << number << std::endl; }
      if       (cut_name == "numuCC.inside.Emuon.theta0.le.226.gt.106"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.296.gt.226"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.386.gt.296"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.505.gt.386"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.577.gt.505"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.659.gt.577"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.753.gt.659"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.861.gt.753"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.984.gt.861"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.1285.gt.984"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984  && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.2506.gt.1285"){  if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.226.gt.106"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.296.gt.226"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.386.gt.296"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.505.gt.386"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.577.gt.505"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.659.gt.577"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.753.gt.659"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.861.gt.753"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.984.gt.861"){    if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.1285.gt.984"){   if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984  && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.2506.gt.1285"){  if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }
    }
    else if (cut_file == 8) {
      if (cut_name == "NCpi0BDT.inside.Enu.le.8000.gt.275"){ if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=8000 && eval.truth_nuEnergy>275) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 9) {
      if        (cut_name == "NCpi0BDT.inside.Ppi0.le.100.gt.0"){     if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=100 && truth_pi0_momentum>0) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.150.gt.100"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=150 && truth_pi0_momentum>100) return number;  
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.200.gt.150"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=200 && truth_pi0_momentum>150) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.250.gt.200"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=250 && truth_pi0_momentum>200) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.300.gt.250"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=300 && truth_pi0_momentum>250) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.400.gt.300"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=400 && truth_pi0_momentum>300) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.500.gt.400"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=500 && truth_pi0_momentum>400) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.600.gt.500"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=600 && truth_pi0_momentum>500) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.800.gt.600"){   if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=800 && truth_pi0_momentum>600) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.1000.gt.800"){  if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=1000 && truth_pi0_momentum>800) return number;
      }else if  (cut_name == "NCpi0BDT.inside.Ppi0.le.1500.gt.1000"){ if (eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=1500 && truth_pi0_momentum>1000) return number;
      }else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
  }
  return -1;
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){

  float reco_Enu = get_reco_Enu_corr(kine, flag_data);
  double Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
  double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

  double truth_pi0_momentum = -1000.;
  if(pfeval.truth_pio_energy_1 > 0. && pfeval.truth_pio_energy_2 > 0.){
    double pi0_mass = 135;
    double alpha = fabs(pfeval.truth_pio_energy_1 - pfeval.truth_pio_energy_2)/(pfeval.truth_pio_energy_1 + pfeval.truth_pio_energy_2);
    double pi0_energy = pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(pfeval.truth_pio_angle/180.*3.1415926)))-1);
    double pi0_total_energy = pi0_mass + pi0_energy;
    truth_pi0_momentum = sqrt(pi0_total_energy*pi0_total_energy - pi0_mass*pi0_mass);
  }
  
  bool flag_truth_inside = false; // in the active volume
  if (eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3 &&  eval.truth_vtxY >-115.0 && eval.truth_vtxY<=117.0 && eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <=1036.4) flag_truth_inside = true;

  // definition of additional cuts
  std::map<std::string, bool> map_cuts_flag;
  if(is_far_sideband(kine, tagger, flag_data)) map_cuts_flag["farsideband"] = true; 
  else map_cuts_flag["farsideband"] = false; 
  
  if(is_near_sideband(kine, tagger, flag_data)) map_cuts_flag["nearsideband"] = true; 
  else map_cuts_flag["nearsideband"] = false; 
 
  if(is_nueCC(tagger)) map_cuts_flag["nueCC"] = true;
  else map_cuts_flag["nueCC"] = false;

  if(is_loosenueCC(tagger)) map_cuts_flag["loosenueCC"] = true;
  else map_cuts_flag["loosenueCC"] = false;
  
  if(is_generic(eval)) map_cuts_flag["generic"] = true;
  else map_cuts_flag["generic"] = false;

  if(eval.truth_nuEnergy <=400) map_cuts_flag["LowEintnueCC"] = true;
  else map_cuts_flag["LowEintnueCC"] = false;
  
  if (!(eval.truth_nuEnergy <=400)) map_cuts_flag["antiLowEintnueCC"] = true;
  else map_cuts_flag["antiLowEintnueCC"] = false;

  if(eval.truth_nuEnergy<=400) map_cuts_flag["LowEnu"] = true;
  else map_cuts_flag["LowEnu"] = false;
  
  if(!(eval.truth_nuEnergy<=400)) map_cuts_flag["antiLowEnu"] = true;
  else map_cuts_flag["antiLowEnu"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["badmatch"] = true;
  else map_cuts_flag["badmatch"] = false;
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["numuCCinFV"] = true;
  else map_cuts_flag["numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["RnumuCCinFV"] = true;
  else map_cuts_flag["RnumuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==-14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["AnumuCCinFV"] = true;
  else map_cuts_flag["AnumuCCinFV"] = false;

  // Xs related cuts ...
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 ) map_cuts_flag["XsnumuCCinFV"] = true;
  else map_cuts_flag["XsnumuCCinFV"] = false;
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200) map_cuts_flag["Xs_Enu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Enu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon > 105.7 && Emuon<=2506) map_cuts_flag["Xs_Emu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Emu_numuCCinFV"] = false;


  //NCpi0BDT X-section

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && eval.truth_nuEnergy<=8000 && eval.truth_nuEnergy > 275) map_cuts_flag["Xs_Enu_NCpi0BDTinFV"] = true;
  else map_cuts_flag["Xs_Enu_NCpi0BDTinFV"] = false;  

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0 && truth_pi0_momentum<=1500 && truth_pi0_momentum > 0) map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"] = true;
  else map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"] = false; 



  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron > 30 && Ehadron <=2500) map_cuts_flag["Xs_Ehad_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Ehad_numuCCinFV"] = false;

  // xs breakdown mode
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCinFV"] = true;
  else map_cuts_flag["XsecNumuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0) map_cuts_flag["XsecNC"] = true;
  else map_cuts_flag["XsecNC"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["XsecCosmic"] = true;
  else map_cuts_flag["XsecCosmic"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==1 && !(eval.truth_nuPdg==14 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200)) map_cuts_flag["XsecBkgCC"] = true;
  else map_cuts_flag["XsecBkgCC"] = false;

  // finish Xs related cuts ...
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200) map_cuts_flag["Xs_Enu_nueCCinFV"] = true;
  else map_cuts_flag["Xs_Enu_nueCCinFV"] = false;
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["nueCCinFV"] = true;
  else map_cuts_flag["nueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["RnueCCinFV"] = true;
  else map_cuts_flag["RnueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==-12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["AnueCCinFV"] = true;
  else map_cuts_flag["AnueCCinFV"] = false;
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["NCinFV"] = true;
  else map_cuts_flag["NCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==0) map_cuts_flag["outFV"] = true;
  else map_cuts_flag["outFV"] = false;
      
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0) map_cuts_flag["CCpi0inFV"] = true;
  else map_cuts_flag["CCpi0inFV"] = false;
      
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0) map_cuts_flag["NCpi0inFV"] = true;
  else map_cuts_flag["NCpi0inFV"] = false;

  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCMEC"] = true;
  else map_cuts_flag["CCMEC"] = false;
  
  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCMEC"] = true;
  else map_cuts_flag["NCMEC"] = false;
  
  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCQE"] = true;
  else map_cuts_flag["CCQE"] = false;

  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCQE"] = true;
  else map_cuts_flag["NCQE"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCRES"] = true;
  else map_cuts_flag["CCRES"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCRES"] = true;
  else map_cuts_flag["NCRES"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCDIS"] = true;
  else map_cuts_flag["CCDIS"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCDIS"] = true;
  else map_cuts_flag["NCDIS"] = false;
  
  if(pfeval.truth_nuScatType!=10 && pfeval.truth_nuScatType!=1 && pfeval.truth_nuScatType!=3 && pfeval.truth_nuScatType!=4 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["OTHER"] = true;
  else map_cuts_flag["OTHER"] = false;

  // figure out additional cuts and flag_data ...
  bool flag_add = true;
  if(add_cut == "all") flag_add = true;
  else if( (flag_data && (add_cut=="farsideband" || add_cut=="nearsideband" || add_cut=="nueCC" || add_cut=="generic" || add_cut=="loosenueCC")) || !flag_data ){ 
    std::istringstream sss(add_cut.Data());
    for(std::string line; std::getline(sss, line, '_');){
      if(map_cuts_flag.find(line)!=map_cuts_flag.end()){ 
        flag_add *= map_cuts_flag[line];
      }
      else{
        std::cout<<"ERROR: add_cut "<<line<<" not defined!\n";
        exit(EXIT_FAILURE);
      }
    } 
  }
  else{ 
    std::cout<<"ERROR: add_cut "<<add_cut<<" of channel "<< ch_name <<" is not assigned to sample "<<flag_data<<" [1: data; 0: mc]\n";
    std::cout<<"Please modify inc/WCPLEEANA/cuts.h\n";
    exit(EXIT_FAILURE);
  }
 
  if (!flag_add) return false;

  bool flag_generic = is_generic(eval);
  bool flag_numuCC = is_numuCC(tagger);
  // bool flag_numuCC = is_numuCC(tagger) and (is_far_sideband(kine, tagger, flag_data) or is_near_sideband(kine, tagger, flag_data) );
  bool flag_numuCC_tight = is_numuCC_tight(tagger, pfeval);
  bool flag_numuCC_1mu0p = is_numuCC_1mu0p(tagger, kine, pfeval);
  bool flag_numuCC_lowEhad = is_numuCC_lowEhad(tagger, kine, pfeval, flag_data);
  bool flag_numuCC_cutbased = is_numuCC_cutbased(tagger);
  bool flag_nueCC = is_nueCC(tagger);
  bool flag_nueCC_loose = is_loosenueCC(tagger);
  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);
  bool flag_FC = is_FC(eval);
  bool flag_filter_showers = is_filter_shower(kine);
  bool flag_0p = is_0p(tagger, kine, pfeval);
  // bool flag_ncpio_bdt = is_NCpio_bdt(tagger) && (!flag_0p);
  bool flag_ncpio_bdt = is_NCpio_bdt(tagger);
  bool flag_ncdelta_bdt = is_NCdelta_bdt(tagger);

  // Total efficiency evaluation 
  bool flag_truth_NCpi0_inside = is_truth_NCpi0_inside(eval, pfeval);

  //float costheta_binning[10] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};   // PeLEE binning
  //float costheta_binning[7]  = {-1,         .27,      .62, .76, .86, .94, 1};   // coarse binning
  float costheta_binning[3]    = {-1,                   .62,                1};   //very coarse binning
  TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
  
  if (ch_name == "LEE_FC_nueoverlay"  || ch_name == "nueCC_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext" || ch_name == "BG_nueCC_FC_dirt" || ch_name =="nueCC_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "LEE_PC_nueoverlay" || ch_name == "nueCC_PC_nueoverlay" ){
    // nueCC PC
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext" || ch_name == "BG_nueCC_PC_dirt" || ch_name == "nueCC_PC_bnb"){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  
  
  ////////////////////////////////////////////////////////
  // GENERIC All events
  ////////////////////////////////////////////////////////
  }else if (ch_name == "generic_bnb" || ch_name == "generic_ext" || ch_name == "generic_dirt"  || ch_name == "generic_overlay"
       || ch_name == "generic1_bnb"  || ch_name == "generic1_ext"  || ch_name == "generic1_dirt"  || ch_name == "generic1_overlay"
       || ch_name == "generic2_bnb"  || ch_name == "generic2_ext"  || ch_name == "generic2_dirt"  || ch_name == "generic2_overlay"
       ){
    if (flag_generic) return true;
    else return false;
  }else if (ch_name == "generic_bnb_FC" || ch_name == "generic_ext_FC" || ch_name == "generic_dirt_FC"  || ch_name == "generic_overlay_FC"
       || ch_name == "generic1_bnb_FC"  || ch_name == "generic1_ext_FC"  || ch_name == "generic1_dirt_FC"  || ch_name == "generic1_overlay_FC"
       || ch_name == "generic2_bnb_FC"  || ch_name == "generic2_ext_FC"  || ch_name == "generic2_dirt_FC"  || ch_name == "generic2_overlay_FC"
       ){
    if (flag_generic && flag_FC) return true;
    else return false;
  }else if (ch_name == "generic_bnb_PC" || ch_name == "generic_ext_PC" || ch_name == "generic_dirt_PC"  || ch_name == "generic_overlay_PC"
       || ch_name == "generic1_bnb_PC"  || ch_name == "generic1_ext_PC"  || ch_name == "generic1_dirt_PC"  || ch_name == "generic1_overlay_PC"
       || ch_name == "generic2_bnb_PC"  || ch_name == "generic2_ext_PC"  || ch_name == "generic2_dirt_PC"  || ch_name == "generic2_overlay_PC"
       ){
    if (flag_generic && (!flag_FC)) return true;
    else return false;
  
  ////////////////////////////////////////////////////////
  // FC+PC NCpi0 selection events
  ////////////////////////////////////////////////////////
  }else if (ch_name == "NCpi0BDT_bnb" || ch_name == "NCpi0BDT_ext" || ch_name == "NCpi0BDT_dirt"  || ch_name == "NCpi0BDT_overlay"
       || ch_name == "NCpi0BDT1_bnb" || ch_name == "NCpi0BDT1_ext"  || ch_name == "NCpi0BDT1_dirt"  || ch_name == "NCpi0BDT1_overlay"
       || ch_name == "NCpi0BDT2_bnb"  || ch_name == "NCpi0BDT2_ext"  || ch_name == "NCpi0BDT2_dirt"  || ch_name == "NCpi0BDT2_overlay"
       || ch_name == "NCpi0BDT3_bnb"  || ch_name == "NCpi0BDT3_ext"  || ch_name == "NCpi0BDT3_dirt"  || ch_name == "NCpi0BDT3_overlay"
       || ch_name == "NCpi0BDT4_bnb"  || ch_name == "NCpi0BDT4_ext"  || ch_name == "NCpi0BDT4_dirt"  || ch_name == "NCpi0BDT4_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers) return true;
    else return false;
  // FC+PC 0p NCpi0 selection events
  }else if (ch_name == "NCpi0BDT_0p_bnb" || ch_name == "NCpi0BDT_0p_ext" || ch_name == "NCpi0BDT_0p_dirt"  || ch_name == "NCpi0BDT_0p_overlay"
       || ch_name == "NCpi0BDT1_0p_bnb" || ch_name == "NCpi0BDT1_0p_ext"  || ch_name == "NCpi0BDT1_0p_dirt"  || ch_name == "NCpi0BDT1_0p_overlay"
       || ch_name == "NCpi0BDT2_0p_bnb"  || ch_name == "NCpi0BDT2_0p_ext"  || ch_name == "NCpi0BDT2_0p_dirt"  || ch_name == "NCpi0BDT2_0p_overlay"
       || ch_name == "NCpi0BDT3_0p_bnb"  || ch_name == "NCpi0BDT3_0p_ext"  || ch_name == "NCpi0BDT3_0p_dirt"  || ch_name == "NCpi0BDT3_0p_overlay"
       || ch_name == "NCpi0BDT4_0p_bnb"  || ch_name == "NCpi0BDT4_0p_ext"  || ch_name == "NCpi0BDT4_0p_dirt"  || ch_name == "NCpi0BDT4_0p_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p) return true;
    else return false;
  // FC+PC Np NCpi0 selection events
  }else if (ch_name == "NCpi0BDT_Np_bnb" || ch_name == "NCpi0BDT_Np_ext" || ch_name == "NCpi0BDT_Np_dirt"  || ch_name == "NCpi0BDT_Np_overlay"
       || ch_name == "NCpi0BDT1_Np_bnb" || ch_name == "NCpi0BDT1_Np_ext"  || ch_name == "NCpi0BDT1_Np_dirt"  || ch_name == "NCpi0BDT1_Np_overlay"
       || ch_name == "NCpi0BDT2_Np_bnb"  || ch_name == "NCpi0BDT2_Np_ext"  || ch_name == "NCpi0BDT2_Np_dirt"  || ch_name == "NCpi0BDT2_Np_overlay"
       || ch_name == "NCpi0BDT3_Np_bnb"  || ch_name == "NCpi0BDT3_Np_ext"  || ch_name == "NCpi0BDT3_Np_dirt"  || ch_name == "NCpi0BDT3_Np_overlay"
       || ch_name == "NCpi0BDT4_Np_bnb"  || ch_name == "NCpi0BDT4_Np_ext"  || ch_name == "NCpi0BDT4_Np_dirt"  || ch_name == "NCpi0BDT4_Np_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p)) return true;
    else return false;
  
  ////////////////////////////////////////////////////////
  // FC only NCpi0 selection events 
  ////////////////////////////////////////////////////////
  }else if (ch_name == "NCpi0BDT_FC_bnb" || ch_name == "NCpi0BDT_FC_ext" || ch_name == "NCpi0BDT_FC_dirt"  || ch_name == "NCpi0BDT_FC_overlay" || ch_name == "BG_NCpi0BDT_FC_ext" || ch_name == "BG_NCpi0BDT_FC_dirt"
       || ch_name == "NCpi0BDT1_FC_bnb" || ch_name == "NCpi0BDT1_FC_ext"  || ch_name == "NCpi0BDT1_FC_dirt"  || ch_name == "NCpi0BDT1_FC_overlay"
       || ch_name == "NCpi0BDT2_FC_bnb" || ch_name == "NCpi0BDT2_FC_ext"  || ch_name == "NCpi0BDT2_FC_dirt"  || ch_name == "NCpi0BDT2_FC_overlay"
       || ch_name == "NCpi0BDT3_FC_bnb" || ch_name == "NCpi0BDT3_FC_ext"  || ch_name == "NCpi0BDT3_FC_dirt"  || ch_name == "NCpi0BDT3_FC_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC) return true;
    else return false;
  // FC 0p NCpi0 selection events
  }else if (ch_name == "NCpi0BDT_FC_0p_bnb" || ch_name == "NCpi0BDT_FC_0p_ext" || ch_name == "NCpi0BDT_FC_0p_dirt"  || ch_name == "NCpi0BDT_FC_0p_overlay" || ch_name == "BG_NCpi0BDT_FC_0p_ext" || ch_name == "BG_NCpi0BDT_FC_0p_dirt"
       || ch_name == "NCpi0BDT1_FC_0p_bnb" || ch_name == "NCpi0BDT1_FC_0p_ext"  || ch_name == "NCpi0BDT1_FC_0p_dirt"  || ch_name == "NCpi0BDT1_FC_0p_overlay"
       || ch_name == "NCpi0BDT2_FC_0p_bnb" || ch_name == "NCpi0BDT2_FC_0p_ext"  || ch_name == "NCpi0BDT2_FC_0p_dirt"  || ch_name == "NCpi0BDT2_FC_0p_overlay"
       || ch_name == "NCpi0BDT3_FC_0p_bnb" || ch_name == "NCpi0BDT3_FC_0p_ext"  || ch_name == "NCpi0BDT3_FC_0p_dirt"  || ch_name == "NCpi0BDT3_FC_0p_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC && flag_0p) return true;
    else return false;
  // FC Np NCpi0 selection events
  }else if (ch_name == "NCpi0BDT_FC_Np_bnb" || ch_name == "NCpi0BDT_FC_Np_ext" || ch_name == "NCpi0BDT_FC_Np_dirt"  || ch_name == "NCpi0BDT_FC_Np_overlay" || ch_name == "BG_NCpi0BDT_FC_Np_ext" || ch_name == "BG_NCpi0BDT_FC_Np_dirt"
       || ch_name == "NCpi0BDT1_FC_Np_bnb" || ch_name == "NCpi0BDT1_FC_Np_ext"  || ch_name == "NCpi0BDT1_FC_Np_dirt"  || ch_name == "NCpi0BDT1_FC_Np_overlay"
       || ch_name == "NCpi0BDT2_FC_Np_bnb" || ch_name == "NCpi0BDT2_FC_Np_ext"  || ch_name == "NCpi0BDT2_FC_Np_dirt"  || ch_name == "NCpi0BDT2_FC_Np_overlay"
       || ch_name == "NCpi0BDT3_FC_Np_bnb" || ch_name == "NCpi0BDT3_FC_Np_ext"  || ch_name == "NCpi0BDT3_FC_Np_dirt"  || ch_name == "NCpi0BDT3_FC_Np_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC && (!flag_0p)) return true;
    else return false;

  ////////////////////////////////////////////////////////
  // PC only NCpi0 selection events 
  ////////////////////////////////////////////////////////
  }else if (ch_name == "NCpi0BDT_PC_bnb" || ch_name == "NCpi0BDT_PC_ext" || ch_name == "NCpi0BDT_PC_dirt"  || ch_name == "NCpi0BDT_PC_overlay" || ch_name == "BG_NCpi0BDT_PC_ext" || ch_name == "BG_NCpi0BDT_PC_dirt"
       || ch_name == "NCpi0BDT1_PC_bnb" || ch_name == "NCpi0BDT1_PC_ext"  || ch_name == "NCpi0BDT1_PC_dirt"  || ch_name == "NCpi0BDT1_PC_overlay"
       || ch_name == "NCpi0BDT2_PC_bnb" || ch_name == "NCpi0BDT2_PC_ext"  || ch_name == "NCpi0BDT2_PC_dirt"  || ch_name == "NCpi0BDT2_PC_overlay"
       || ch_name == "NCpi0BDT3_PC_bnb" || ch_name == "NCpi0BDT3_PC_ext"  || ch_name == "NCpi0BDT3_PC_dirt"  || ch_name == "NCpi0BDT3_PC_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC)) return true;
    else return false;
  // PC 0p NCpi0 selection events
  }else if (ch_name == "NCpi0BDT_PC_0p_bnb" || ch_name == "NCpi0BDT_PC_0p_ext" || ch_name == "NCpi0BDT_PC_0p_dirt"  || ch_name == "NCpi0BDT_PC_0p_overlay" || ch_name == "BG_NCpi0BDT_PC_0p_ext" || ch_name == "BG_NCpi0BDT_PC_0p_dirt"
       || ch_name == "NCpi0BDT1_PC_0p_bnb" || ch_name == "NCpi0BDT1_PC_0p_ext"  || ch_name == "NCpi0BDT1_PC_0p_dirt"  || ch_name == "NCpi0BDT1_PC_0p_overlay"
       || ch_name == "NCpi0BDT2_PC_0p_bnb" || ch_name == "NCpi0BDT2_PC_0p_ext"  || ch_name == "NCpi0BDT2_PC_0p_dirt"  || ch_name == "NCpi0BDT2_PC_0p_overlay"
       || ch_name == "NCpi0BDT3_PC_0p_bnb" || ch_name == "NCpi0BDT3_PC_0p_ext"  || ch_name == "NCpi0BDT3_PC_0p_dirt"  || ch_name == "NCpi0BDT3_PC_0p_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC) && flag_0p) return true;
    else return false;
  // PC Np NCpi0 selection events
  }else if (ch_name == "NCpi0BDT_PC_Np_bnb" || ch_name == "NCpi0BDT_PC_Np_ext" || ch_name == "NCpi0BDT_PC_Np_dirt"  || ch_name == "NCpi0BDT_PC_Np_overlay" || ch_name == "BG_NCpi0BDT_PC_Np_ext" || ch_name == "BG_NCpi0BDT_PC_Np_dirt"
       || ch_name == "NCpi0BDT1_PC_Np_bnb" || ch_name == "NCpi0BDT1_PC_Np_ext"  || ch_name == "NCpi0BDT1_PC_Np_dirt"  || ch_name == "NCpi0BDT1_PC_Np_overlay"
       || ch_name == "NCpi0BDT2_PC_Np_bnb" || ch_name == "NCpi0BDT2_PC_Np_ext"  || ch_name == "NCpi0BDT2_PC_Np_dirt"  || ch_name == "NCpi0BDT2_PC_Np_overlay"
       || ch_name == "NCpi0BDT3_PC_Np_bnb" || ch_name == "NCpi0BDT3_PC_Np_ext"  || ch_name == "NCpi0BDT3_PC_Np_dirt"  || ch_name == "NCpi0BDT3_PC_Np_overlay"
       ){
    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC) && (!flag_0p)) return true;
    else return false;

  ////////////////////////////////////////////////////////
  // numuCC XSEC 
  //////////////////////////////////////////////////////// 
  }else if (ch_name == "numuCC_FC_bnb" || ch_name == "BG_numuCC_FC_ext" || ch_name == "BG_numuCC_FC_dirt"
    ){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb" || ch_name == "BG_numuCC_PC_ext" || ch_name == "BG_numuCC_PC_dirt"
    ){
    if (flag_numuCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC_background_PC_overlay"
    ){
    if (ch_name == "numuCC_signal_FC_overlay"){           if (flag_numuCC && flag_FC    && map_cuts_flag["XsnumuCCinFV"])     return true;
    }else if (ch_name == "numuCC_signal_PC_overlay"){     if (flag_numuCC && (!flag_FC) && map_cuts_flag["XsnumuCCinFV"])     return true;
    }else if (ch_name == "numuCC_background_FC_overlay"){ if (flag_numuCC && flag_FC    && (!map_cuts_flag["XsnumuCCinFV"]))  return true;
    }else if (ch_name == "numuCC_background_PC_overlay"){ if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["XsnumuCCinFV"]))  return true;
    }
    return false;

  ////////////////////////////////////////////////////////
  // numuCC XSEC - Enu
  //////////////////////////////////////////////////////// 
  } else if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC_background_Enu_PC_overlay"
       || ch_name == "numuCC_signal_Enu_overlay" || ch_name == "numuCC_background_Enu_overlay"
      ){
    if (ch_name == "numuCC_signal_Enu_FC_overlay"){           if (flag_numuCC && flag_FC    &&  map_cuts_flag["Xs_Enu_numuCCinFV"])     return true;
    }else if (ch_name == "numuCC_signal_Enu_PC_overlay"){     if (flag_numuCC && (!flag_FC) &&  map_cuts_flag["Xs_Enu_numuCCinFV"])     return true;
    }else if (ch_name == "numuCC_background_Enu_FC_overlay"){ if (flag_numuCC && flag_FC    &&  (!map_cuts_flag["Xs_Enu_numuCCinFV"]))  return true;
    }else if (ch_name == "numuCC_background_Enu_PC_overlay"){ if (flag_numuCC && (!flag_FC) &&  (!map_cuts_flag["Xs_Enu_numuCCinFV"]))  return true;
    }else if (ch_name == "numuCC_signal_Enu_overlay"){        if (flag_numuCC &&                map_cuts_flag["Xs_Enu_numuCCinFV"])     return true;
    }else if (ch_name == "numuCC_background_Enu_overlay"){    if (flag_numuCC &&                (!map_cuts_flag["Xs_Enu_numuCCinFV"]))  return true;
    }
    return false;

  ////////////////////////////////////////////////////////
  // NCpi0BDT XSEC - Enu
  //////////////////////////////////////////////////////// 
  } else if (ch_name == "NCpi0BDT_signal_Enu_FC_overlay" || ch_name == "NCpi0BDT_signal_Enu_PC_overlay" || ch_name == "NCpi0BDT_background_Enu_FC_overlay" || ch_name == "NCpi0BDT_background_Enu_PC_overlay" 
       || ch_name == "NCpi0BDT_signal_Enu_overlay" || ch_name == "NCpi0BDT_background_Enu_overlay"
      ){
    if (ch_name == "NCpi0BDT_signal_Enu_FC_overlay"){           if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC     &&  map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_signal_Enu_PC_overlay"){     if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC)  &&  map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_FC_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC     &&  (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_PC_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC)  &&  (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_signal_Enu_overlay"){        if (flag_ncpio_bdt && flag_generic && flag_filter_showers &&                 map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_overlay"){    if (flag_ncpio_bdt && flag_generic && flag_filter_showers &&                 (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }
    return false;

  ////////////////////////////////////////////////////////
  // NCpi0BDT 0p XSEC - Enu 
  //////////////////////////////////////////////////////// 
  } else if (ch_name == "NCpi0BDT_signal_Enu_FC_0p_overlay" || ch_name == "NCpi0BDT_signal_Enu_PC_0p_overlay" || ch_name == "NCpi0BDT_background_Enu_FC_0p_overlay" || ch_name == "NCpi0BDT_background_Enu_PC_0p_overlay" 
       || ch_name == "NCpi0BDT_signal_Enu_0p_overlay" || ch_name == "NCpi0BDT_background_Enu_0p_overlay"
      ){
    if (ch_name == "NCpi0BDT_signal_Enu_FC_0p_overlay"){           if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p && flag_FC     && map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_signal_Enu_PC_0p_overlay"){     if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p && (!flag_FC)  && map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_FC_0p_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p && flag_FC     && (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_PC_0p_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p && (!flag_FC)  && (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_signal_Enu_0p_overlay"){        if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p &&                map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_0p_overlay"){    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_0p &&                (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }
    return false;

  ////////////////////////////////////////////////////////
  // NCpi0BDT Np XSEC - Enu 
  //////////////////////////////////////////////////////// 
  } else if (ch_name == "NCpi0BDT_signal_Enu_FC_Np_overlay" || ch_name == "NCpi0BDT_signal_Enu_PC_Np_overlay" || ch_name == "NCpi0BDT_background_Enu_FC_Np_overlay" || ch_name == "NCpi0BDT_background_Enu_PC_Np_overlay" 
       || ch_name == "NCpi0BDT_signal_Enu_Np_overlay" || ch_name == "NCpi0BDT_background_Enu_Np_overlay"
      ){
    if (ch_name == "NCpi0BDT_signal_Enu_FC_Np_overlay"){           if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p) && flag_FC     && map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_signal_Enu_PC_Np_overlay"){     if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p) && (!flag_FC)  && map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_FC_Np_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p) && flag_FC     && (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_PC_Np_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p) && (!flag_FC)  && (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_signal_Enu_Np_overlay"){        if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p) &&                map_cuts_flag["Xs_Enu_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Enu_Np_overlay"){    if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_0p) &&                (!map_cuts_flag["Xs_Enu_NCpi0BDTinFV"]))  return true;
    }
    return false;

  ////////////////////////////////////////////////////////
  // NCpi0BDT XSEC - Ppi0
  //////////////////////////////////////////////////////// 
  } else if (ch_name == "NCpi0BDT_signal_Ppi0_FC_overlay" || ch_name == "NCpi0BDT_signal_Ppi0_PC_overlay" || ch_name == "NCpi0BDT_background_Ppi0_FC_overlay" || ch_name == "NCpi0BDT_background_Ppi0_PC_overlay" 
       || ch_name == "NCpi0BDT_signal_Ppi0_overlay" || ch_name == "NCpi0BDT_background_Ppi0_overlay"
      ){
    if (ch_name == "NCpi0BDT_signal_Ppi0_FC_overlay"){           if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC     &&  map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_signal_Ppi0_PC_overlay"){     if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC)  &&  map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Ppi0_FC_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && flag_FC     &&  (!map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_background_Ppi0_PC_overlay"){ if (flag_ncpio_bdt && flag_generic && flag_filter_showers && (!flag_FC)  &&  (!map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"]))  return true;
    }else if (ch_name == "NCpi0BDT_signal_Ppi0_overlay"){        if (flag_ncpio_bdt && flag_generic && flag_filter_showers &&                 map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"])     return true;
    }else if (ch_name == "NCpi0BDT_background_Ppi0_overlay"){    if (flag_ncpio_bdt && flag_generic && flag_filter_showers &&                 (!map_cuts_flag["Xs_Ppi0_NCpi0BDTinFV"]))  return true;
    }
    return false;

  }else{ std::cout << "Not sure what cut: " << ch_name << std::endl;
  }
  return false;
}

bool LEEana::is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  bool flag_numuCC = is_numuCC(tagger);
  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  if ((reco_Enu>=800 && tagger.nue_score >=0) || (tagger.nue_score<=0 && (flag_numuCC || (flag_pi0 && flag_NC) ))) flag = true;
  return flag;
}

bool LEEana::is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  if (reco_Enu < 800 && tagger.nue_score>0 && (reco_Enu>=600 || tagger.nue_score<=7)) flag = true;
  return flag ;
}

bool LEEana::is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  if (reco_Enu < 600 && tagger.nue_score>7) flag = true;
  return flag;
}

bool LEEana::is_truth_nueCC_inside(EvalInfo& eval){
  bool flag = false;
  if (fabs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) flag = true;
  return flag;
}

bool LEEana::is_truth_numuCC_inside(EvalInfo& eval){
   bool flag = false;
  if (fabs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1) flag = true;
  return flag;
}

bool LEEana::is_truth_NCpi0_inside(EvalInfo& eval, PFevalInfo& pfeval){
   bool flag = false;
  if (eval.truth_isCC==0 && pfeval.truth_NprimPio>0 && eval.truth_vtxInside==1) flag = true;
  return flag;
}

bool LEEana::is_FC(EvalInfo& eval){
  if (eval.match_isFC){ return true;
  }else{ return false;
  }
}

bool LEEana::is_filter_shower(KineInfo& kine){
  if (kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.){ return true;
  }else{ return false;
  }
}

bool LEEana::is_cc_pi0(KineInfo& kine, bool flag_data){
  bool flag = false;
  if (flag_data){
    if (kine.kine_pio_mass>0){
      //TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      //TLorentzVector pio = p1 + p2;
      //pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;
      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300) flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300) flag = true;
  }
  return flag;
}

bool LEEana::is_pi0(KineInfo& kine, bool flag_data){
  bool flag = false;
  if (flag_data){
    if (kine.kine_pio_mass>0){
      //TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      //TLorentzVector pio = p1 + p2;
      //pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;
      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300) flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300) flag = true;
  }
  return flag;
}

bool LEEana::is_NCpio_bdt(TaggerInfo& tagger_info){
  bool flag = false;
  // if (tagger_info.nc_pio_score > 1.68 && tagger_info.numu_cc_flag >=0) flag = true;
  if (tagger_info.nc_pio_score > 1.816 && tagger_info.numu_cc_flag >=0) flag = true;
  return flag;
}
bool LEEana::is_NCdelta_bdt(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.nc_delta_score > 2.61 && tagger_info.numu_cc_flag >=0) flag = true;
  return flag;
}

bool LEEana::is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0) flag = true;
  return flag;
}

bool LEEana::is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9) flag = true;
  return flag;
}

bool LEEana::is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0) flag = true;
  return flag;
}

bool LEEana::is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.numu_cc_flag>=0){ 
    // 1 lepton <=1 proton 0 charged pion
    // 1 lepton guaranteed by numu cc flag
    // using pi0 flag to remove pi0 component in channel definition
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
      int pdgcode = kine.kine_particle_type->at(i);
      if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm? 
      if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV 
    }
    if(Nproton==0) flag = true;
  } 
  return flag;
}

bool LEEana::is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.numu_cc_flag>=0){ 
    // 1 lepton <=1 proton 0 charged pion
    // 1 lepton guaranteed by numu cc flag
    // using pi0 flag to remove pi0 component in channel definition
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++) {
        int pdgcode = kine.kine_particle_type->at(i);
        if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm? 
        if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV 
    }
    if(Npion==0) flag = true;
  } 
  return flag;
}

bool LEEana::is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){ 
    // 1 lepton <=1 proton 0 charged pion
    // 1 lepton guaranteed by numu cc flag
    // using pi0 flag to remove pi0 component in channel definition
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
        int pdgcode = kine.kine_particle_type->at(i);
        if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm? 
        if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV 
    }
    if(Nproton==0) flag = true;
  } 
  return flag;
}


bool LEEana::is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data){
  bool flag = false;
  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){
    double reco_Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Ehadron = reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
    if(Ehadron<200){ flag = true;
    }
  }
  return flag;
}

bool LEEana::is_numuCC_cutbased(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag==1 && tagger_info.cosmict_flag==0)  flag = true;
  return flag;
}

bool LEEana::is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;
  // default 7.0
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    //  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score <= 7.0 && tagger_info.nue_score > 0)
    flag = true;
  return flag;
}

bool LEEana::is_loosenueCC(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 4.0) flag = true;
  return flag;
}

bool LEEana::is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);
  flag = flag && (eval.stm_clusterlength > 15);
  return flag;
}

bool LEEana::is_preselection(EvalInfo& eval){
  bool flag = false;
  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){ tmp_match_found = eval.match_found_asInt;
  }
  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) flag = true;
  return flag;
}

int LEEana::mcc8_pmuon_costheta_bin(float pmuon, float costh){
  if (costh>=-1 and costh<-0.5) {
    if (pmuon>=0 and pmuon<0.18) return 1;
    else if (pmuon>=0.18 and pmuon<0.30) return 2;
    else if (pmuon>=0.30 and pmuon<0.45) return 3;
    else if (pmuon>=0.45 and pmuon<0.77) return 4;
    else if (pmuon>=0.77 and pmuon<2.5) return 5;
    else return -10000;
  }
  else if (costh>=-0.5 and costh<0){
    if (pmuon>=0 and pmuon<0.18) return 6;
    else if (pmuon>=0.18 and pmuon<0.30) return 7;
    else if (pmuon>=0.30 and pmuon<0.45) return 8;
    else if (pmuon>=0.45 and pmuon<0.77) return 9;
    else if (pmuon>=0.77 and pmuon<2.5) return 10;
    else return -10000;
  }
  else if (costh>0 and costh<0.27){
    if (pmuon>=0 and pmuon<0.18) return 11;
    else if (pmuon>=0.18 and pmuon<0.30) return 12;
    else if (pmuon>=0.30 and pmuon<0.45) return 13;
    else if (pmuon>=0.45 and pmuon<0.77) return 14;
    else if (pmuon>=0.77 and pmuon<2.5) return 15;
    else return -10000;    
  }
  else if (costh>=0.27 and costh<0.45){
    if (pmuon>=0 and pmuon<0.30) return 16;
    else if (pmuon>=0.30 and pmuon<0.45) return 17;
    else if (pmuon>=0.45 and pmuon<0.77) return 18;
    else if (pmuon>=0.77 and pmuon<2.5) return 19;
    else return -10000;    
  }
  else if (costh>=0.45 and costh<0.62){
    if (pmuon>=0 and pmuon<0.30) return 20;
    else if (pmuon>=0.30 and pmuon<0.45) return 21;
    else if (pmuon>=0.45 and pmuon<0.77) return 22;
    else if (pmuon>=0.77 and pmuon<2.5) return 23;
    else return -10000;   
  }
  else if (costh>=0.62 and costh<0.76){
    if (pmuon>=0 and pmuon<0.30) return 24;
    else if (pmuon>=0.30 and pmuon<0.45) return 25;
    else if (pmuon>=0.45 and pmuon<0.77) return 26;
    else if (pmuon>=0.77 and pmuon<2.5) return 27;
    else return -10000;  
  }
  else if (costh>=0.76 and costh<0.86){
    if (pmuon>=0 and pmuon<0.30) return 28;
    else if (pmuon>=0.30 and pmuon<0.45) return 29;
    else if (pmuon>=0.45 and pmuon<0.77) return 30;
    else if (pmuon>=0.77 and pmuon<1.28) return 31;
    else if (pmuon>=1.28 and pmuon<2.5) return 32;
    else return -10000;
  }
  else if (costh>=0.86 and costh<0.94){
    if (pmuon>=0 and pmuon<0.30) return 33;
    else if (pmuon>=0.30 and pmuon<0.45) return 34;
    else if (pmuon>=0.45 and pmuon<0.77) return 35;
    else if (pmuon>=0.77 and pmuon<1.28) return 36;
    else if (pmuon>=1.28 and pmuon<2.5) return 37;
    else return -10000;     
  }
  else if (costh>=0.94 and costh<1.00){
    if (pmuon>=0 and pmuon<0.30) return 38;
    else if (pmuon>=0.30 and pmuon<0.45) return 39;
    else if (pmuon>=0.45 and pmuon<0.77) return 40;
    else if (pmuon>=0.77 and pmuon<1.28) return 41;
    else if (pmuon>=1.28 and pmuon<2.5) return 42;
    else return -10000;   
  }
  else return -10000;
} 

#endif