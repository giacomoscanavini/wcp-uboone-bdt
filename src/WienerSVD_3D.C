#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TDecompSVD.h"

#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <string>

const int nbins1_2D    = 1; //1 Enu bins
const int nbins2_2D    = 4; //4 costheta bins
const int nbins3_2D    = 5; //5 momentum bins
const int nbins_tot_2D = 20;

const int nbins1_3D    = 1; //1 Enu bins
const int nbins2_3D    = 4; //4 costheta bins
const int nbins3_3D    = 5; //5 momentum bins
const int nbins_tot_3D = 20;

std::vector<int> nbins_2D = { nbins1_2D, nbins2_2D, nbins3_2D };
std::vector<int> nbins_3D = { nbins1_3D, nbins2_3D, nbins3_3D };

//maps from signal bin index to Enu slice (not used in 2D)
int map_signal_dim1_2D[nbins_tot_2D] = {0,0,0,0,0, 0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0};

//maps from signal bin index to costheta slice
int map_signal_dim2_2D[nbins_tot_2D] = { 0,0,0,0,0,
                                         1,1,1,1,1,
                                         2,2,2,2,2,
                                         3,3,3,3,3 };

//maps from signal bin index to momentum min slice
int map_signal_dim3_2D[nbins_tot_2D] = { 0,1,2,3,4,
                                         0,1,2,3,4,
                                         0,1,2,3,4,
                                         0,1,2,3,4 };



//maps from signal bin index to Enu slice
int map_signal_dim1_3D[nbins_tot_3D] = {0,0,0,0,0, 0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0};

//maps from signal bin index to costheta slice
int map_signal_dim2_3D[nbins_tot_3D] = { 0,0,0,0,0,
                                         1,1,1,1,1,
                                         2,2,2,2,2,
                                         3,3,3,3,3 };

//maps from signal bin index to momentum slice
int map_signal_dim3_3D[nbins_tot_3D] = { 0,1,2,3,4,
                                         0,1,2,3,4,
                                         0,1,2,3,4,
                                         0,1,2,3,4 };

//maps from {Enu, costheta, momentum} indices to signal bin index
std::vector<std::vector<std::vector<int>>> map_grid_signal_bin_2D = { { {  0,  1,  2,  3,  4, },
                                                                        {  5,  6,  7,  8,  9, },
                                                                        { 10, 11, 12, 13, 14, },
                                                                        { 15, 16, 17, 18, 19 } } };

//maps from {Enu, costheta, momentum} indices to signal bin index
std::vector<std::vector<std::vector<int>>> map_grid_signal_bin_3D = { { {  0,  1,  2,  3,  4, },
                                                                        {  5,  6,  7,  8,  9, },
                                                                        { 10, 11, 12, 13, 14, },
                                                                        { 15, 16, 17, 18, 19 } } };
 

//gets the index for the next (or previous) bin along a given dimension.  Returns -1 if a boundary is hit.
int get_next_bin (int ndim, int start_bin, int dim, bool direction) {
  if (start_bin==-1) { return -1; }

  int nbins_dim = -1;
  std::vector<int> bin_3D = {};
  if      (ndim==2) {
    nbins_dim = nbins_2D[dim];
    bin_3D = { map_signal_dim1_2D[start_bin], map_signal_dim2_2D[start_bin], map_signal_dim3_2D[start_bin] };
  }
  else if (ndim==3) {
    nbins_dim = nbins_3D[dim];
    bin_3D = { map_signal_dim1_3D[start_bin], map_signal_dim2_3D[start_bin], map_signal_dim3_3D[start_bin] };
  }
//std::cout << "bin_3D = {" << bin_3D[0] << ", " << bin_3D[1] << ", " << bin_3D[2] << "}" << std::endl;
  while (true) {
    //increase/decrease bin[dim] and check whether you are out of bounds
    bin_3D[dim] = bin_3D[dim] + 1*direction - 1*(!direction);
    if (bin_3D[dim]<0 || bin_3D[dim]>=nbins_dim) { return -1; }

    //check whether the bin is different, otherwise continue
    int newbin = -1;
    if      (ndim==2) { newbin = map_grid_signal_bin_2D[bin_3D[0]][bin_3D[1]][bin_3D[2]]; }
    else if (ndim==3) { newbin = map_grid_signal_bin_3D[bin_3D[0]][bin_3D[1]][bin_3D[2]]; }
    if (newbin != start_bin) { return newbin; }
  }
  return -1;
}

//returns how many bins exist along the given dimension and direction
int get_nbins_next (int ndim, int start_bin, int dim, bool direction) {
  if (start_bin==-1) { return -1; }
    std::cout << "ndim, start_bin, dim, dir = " << ndim << ",  " << start_bin << ",  " << dim << ",  " << direction << std::endl;

  int current_bin = start_bin;
  int nbins = -1;  
  while(current_bin!=-1) {
    nbins++;
std::cout << "current bin = " << current_bin << std::endl;
    current_bin = get_next_bin(ndim,current_bin,dim,direction);
  }
  return nbins;
}

/*
  void compute_2nd_derivatives (TMatrixD* C, int start_bin, int j, int k) {

      int prev_bin_j = get_next_bin(start_bin, j, false);
      int next_bin_j = get_next_bin(start_bin, j, true);
      int next_bin_k = get_next_bin(start_bin, k, true);

      Double_t epsilon = 1e-6;

    //second derivative in the same direction
    if (j==k) {
      if (prev_bin_j!=-1) {
        C_temp(start_bin,prev_bin_j) = C_temp(start_bin,prev_bin_j) + 1;
        C_temp(start_bin,start_bin ) = C_temp(start_bin,start_bin ) - 1 + epsilon;
      }
      if (next_bin_j!=-1) {
        C_temp(start_bin,next_bin_j) = C_temp(start_bin,next_bin_j) + 1;
        C_temp(start_bin,start_bin ) = C_temp(start_bin,start_bin ) - 1 + epsilon;
      }
    //second derivative in two separate directions
    } else {

      C_temp(start_bin,start_bin) = C_temp(start_bin,start_bin) + 1;
      if (next_bin_j!=-1) {
        C_temp(start_bin,next_bin_j) = C_temp(start_bin,next_bin_j) - 1 + epsilon;
        //std::cout << "start_bin, next_bin_j = " << start_bin << ",     " << next_bin_j << std::endl;
        int next_bin_j_next_bin_k = get_next_bin(next_bin_j, k, true);

        //std::vector<int> bin_3D = { map_signal_dim1_3D[next_bin_j], map_signal_dim2_3D[next_bin_j], map_signal_dim3_3D[next_bin_j] };
        //std::cout << "looking for bin after " << next_bin_j << ",  (" << bin_3D[0] << ", " << bin_3D[1] << ", " << bin_3D[2] << ")" << std::endl;
        if (next_bin_j_next_bin_k!=-1) {
          //std::cout << "found!" << std::endl;
          C_temp(start_bin,next_bin_j_next_bin_k) = C_temp(start_bin,next_bin_j_next_bin_k) + 0.5;
        }
      }
      if (next_bin_k!=-1) {
        C_temp(start_bin,next_bin_k) = C_temp(start_bin,next_bin_k) - 1 + epsilon;
        int next_bin_k_next_bin_j = get_next_bin(next_bin_k, j, true);
        if (next_bin_k_next_bin_j!=-1) {
          C_temp(start_bin,next_bin_k_next_bin_j) = C_temp(start_bin,next_bin_k_next_bin_j) + 0.5;
        }
      }
    }

  }
*/
//main function to compute C2 regularization matrix using the full 3D bin structure
TMatrixD C_3D (int derivative, int ndim) {

  Double_t epsilon = 1e-2;
  int nbins_tot = -1;
  if      (ndim==2) { nbins_tot = nbins_tot_2D; }
  else if (ndim==3) { nbins_tot = nbins_tot_3D; }

  TMatrixDSym C_sym(nbins_tot);
  for (int i=0;i<nbins_tot;i++) {
    for (int j=0;j<nbins_tot;j++) {
      C_sym(i,j) = 0;
    }
  }

  TMatrixD C_temp_dim1(nbins_tot,nbins_tot);
  TMatrixD C_temp_dim2(nbins_tot,nbins_tot);
  TMatrixD C_temp_dim3(nbins_tot,nbins_tot);
  TMatrixD C_temp_dim1_2(nbins_tot,nbins_tot);
  TMatrixD C_temp_dim2_2(nbins_tot,nbins_tot);
  TMatrixD C_temp_dim3_2(nbins_tot,nbins_tot);
  std::vector<TMatrixD> C_temp_vec;
  std::vector<TMatrixD> C_temp_2_vec;
  C_temp_vec.push_back(C_temp_dim1);
  C_temp_vec.push_back(C_temp_dim2);
  C_temp_vec.push_back(C_temp_dim3);
  C_temp_2_vec.push_back(C_temp_dim1_2);
  C_temp_2_vec.push_back(C_temp_dim2_2);
  C_temp_2_vec.push_back(C_temp_dim3_2);

  //take derivative along Enu, theta, and Pmu
  int dim_max = 3;
  int dim_min = dim_max-ndim;
  for (int dim=dim_min;dim<dim_max;dim++) {

    //TMatrixD C_temp(nbins_tot,nbins_tot);

    //iterate over 36 (2D) or 138 (3D) truth signal bins and compute 3rd derivative for each
    for (int start_bin=0;start_bin<nbins_tot;start_bin++) {
      //set matrix to 0 initially
      for (int j=0;j<nbins_tot;j++) { C_temp_vec[dim](start_bin,j) = 0; }

      int prev_bin  = get_next_bin(ndim,start_bin, dim, false);
      int prev_bin2 = get_next_bin(ndim,prev_bin,  dim, false);
      int next_bin  = get_next_bin(ndim,start_bin, dim, true);
      int next_bin2 = get_next_bin(ndim,next_bin,  dim, true);

      int nbins_prev = get_nbins_next(ndim,start_bin,dim,false);
      int nbins_next = get_nbins_next(ndim,start_bin,dim,true);

      if (derivative==2) {
        //C2 matrix
        if (nbins_prev==0) {
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) - 1 + epsilon;
          C_temp_vec[dim](start_bin,next_bin)  = C_temp_vec[dim](start_bin,next_bin)  + 1;
        } else if (nbins_next==0) {
          C_temp_vec[dim](start_bin,prev_bin)  = C_temp_vec[dim](start_bin,prev_bin)  + 1;
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) - 1 + epsilon;
        } else {
          C_temp_vec[dim](start_bin,prev_bin)  = C_temp_vec[dim](start_bin,prev_bin)  + 1;
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) - 2 + epsilon;
          C_temp_vec[dim](start_bin,next_bin)  = C_temp_vec[dim](start_bin,next_bin)  + 1;
        }
      } else if (derivative==3) {
        //C3 matrix
        if (nbins_prev>=2 && nbins_next>=2) {
          C_temp_vec[dim](start_bin,prev_bin2) = C_temp_vec[dim](start_bin,prev_bin2) - 1;
          C_temp_vec[dim](start_bin,prev_bin)  = C_temp_vec[dim](start_bin,prev_bin)  + 2;
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) + epsilon;
          C_temp_vec[dim](start_bin,next_bin)  = C_temp_vec[dim](start_bin,next_bin)  - 2;
          C_temp_vec[dim](start_bin,next_bin2) = C_temp_vec[dim](start_bin,next_bin2) + 1;
        } else if (nbins_prev==1 && nbins_next==1) {
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) + epsilon;
        } else if (nbins_next==0 || (nbins_next==1 && nbins_prev>=2)) {
          C_temp_vec[dim](start_bin,prev_bin2) = C_temp_vec[dim](start_bin,prev_bin2) - 1;
          C_temp_vec[dim](start_bin,prev_bin)  = C_temp_vec[dim](start_bin,prev_bin)  + 2;
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) - 1 + epsilon;
        } else if (nbins_prev==0 || (nbins_prev==1 && nbins_next>=2)) {
          C_temp_vec[dim](start_bin,start_bin) = C_temp_vec[dim](start_bin,start_bin) + 1 + epsilon;
          C_temp_vec[dim](start_bin,next_bin)  = C_temp_vec[dim](start_bin,next_bin)  - 2;
          C_temp_vec[dim](start_bin,next_bin2) = C_temp_vec[dim](start_bin,next_bin2) + 1;
        } else {
          std::cout << "error here" << std::endl;
        }
      }
    }

    //square C_temp for each dimension and then add to C_sym
    TMatrixD C_temp_t(nbins_tot,nbins_tot);
    C_temp_t.Transpose(C_temp_vec[dim]);
    TMatrixD C_temp_2 = C_temp_t * C_temp_vec[dim];
    for (int i=0;i<nbins_tot;i++) {
      for (int j=0;j<nbins_tot;j++) {
        C_temp_2_vec[dim](i,j) = C_temp_2(i,j);
        C_sym(i,j) = C_sym(i,j) + C_temp_2(i,j);
      }
    }

  }

  //Decompose into C = V*D*V^T
  TMatrixDSymEigen C_eigen(C_sym);
  TMatrixD V = C_eigen.GetEigenVectors();
  TMatrixD V_t(nbins_tot,nbins_tot);
  V_t.Transpose(V);

 //take square root of D
  TVectorD V_eigen = C_eigen.GetEigenValues();
  TMatrixD D(nbins_tot,nbins_tot);
  for (int i=0;i<nbins_tot;i++) {
    for (int j=0;j<nbins_tot;j++) { D(i,j) = 0; }
    D(i,i) = TMath::Sqrt(V_eigen[i]);
  }

  TMatrixD C = V * D * V_t;

/*
  TFile *file = new TFile("/uboone/data/users/lcoopert/LEE/LEEana_xs_3D_feb24/wiener_svd/debug_wienerSVD_3D_full.root","RECREATE");
  file->cd();
  C_temp_vec[0].Write("C_dim1");
  C_temp_vec[1].Write("C_dim2");
  C_temp_vec[2].Write("C_dim3");
  C_temp_2_vec[0].Write("C_dim1_2");
  C_temp_2_vec[1].Write("C_dim2_2");
  C_temp_2_vec[2].Write("C_dim3_2");
  C_sym.Write("C_2");
  C.Write("C");
  file->Write();
  file->Close();
*/
  return C;
}