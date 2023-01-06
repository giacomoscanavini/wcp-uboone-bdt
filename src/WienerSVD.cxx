// Core implementation of Wiener-SVD
// Author: Hanyu WEI,   March 17, 2017
// Version 2.0

#include "WCPLEEANA/WienerSVD.h"
#include "WienerSVD_3D.C" // 2D cross section extraction

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include <iostream>
#include "TMath.h"

using namespace std;

bool is_on_edge(std::vector<int> edges, int element) {
    int size = edges.size();
    for (int i=0;i<size;i++) { if (edges[i] == element) { return true; } }
    return false;
}

int get_slice (std::vector<int> edges, int index) {
    int nbins = edges.size()-1;
    for (int i=0;i<nbins;i++) { if (index >= edges[i] && index < edges[i+1]) { return i; } }
    return -1;
}

TMatrixD Matrix_C(Int_t n, Int_t type)
{
    //Type 0: Unit matrix
    //Type I: First derivative matrix
    //Type II: Second derivative matrix
    Double_t epsilon = 1e-6; // needed for 2nd derivative matrix inversion
    Double_t epsilon2 = 1e-2; // needed for 3rd derivative matrix inversion

    std::cout << "type = " << type << std::endl;
    if      (type==232){ std::cout << "type = " << type << std::endl; }

    std::vector<int> dim_edges { 0, n };
    //2nd derivative
    if      (type==232) { return C_3D(2,2);} // (derivative, dimension)
    else if (type==233) { return C_3D(2,3);} // bugged doens't handle dim = 1 
    //3rd derivative
    else if (type==332) { return C_3D(3,2);}
    else if (type==333) { return C_3D(3,3);} // bugged doens't handle dim = 1
    //else if (type==233) { return C2_3D(3);} // Second derivative for 3D
    else if (type==22 || type==32) { dim_edges = { 0, 5, 10, 15, 20 }; }  //2D edges between 1D dimension slices 

    else if (type==23 || type==33) { dim_edges = { 0, 5, 10, 15, 20 }; }  //3D edges between 1D dimension slices

    std::cout << "n, type, dim_edges.count = " << n << ",   " << type << ",   " << dim_edges.size() << std::endl;

    TMatrixD C(n, n);
    for(Int_t i=0; i<n; i++){
        for(Int_t j=0; j<n; j++){
            C(i, j) = 0;
            if(type == 0){
                if(i == j) C(i, j) = 1;
            }

            else if(type == 1){
                if( j-i == 1 ) C(i, j) = 1;
                if(i==j){
                    C(i, j) = -1;
                }
            }

            //else if(type == 2){
            //    if(TMath::Abs(i-j) == 1) C(i, j) = 1;
            //    if(i==j){
            //        C(i, j) = -2+epsilon;
            //        if(i==0 || i==n-1){
            //            C(i, j) = -1+epsilon;
            //        }
            //    }
            //}

            // Upgrade to type 2 breaking the matrix with multiple channels (type==20)
            else if(type == 2 || type==20){
                if(TMath::Abs(i-j) == 1) C(i, j) = 1;
                if(i==j){
                    C(i, j) = -2+epsilon;
                    if(i==0 || i==n-1){
                        C(i, j) = -1+epsilon;
                    }
                }
                if(type==20 && i==n/2){
                    if(j==i) C(i,j)=-1;
                    if(j==i-1) C(i,j)=0;
                }
                if(type==20 && i==n/2-1){
                    if(j==i) C(i,j)=-1;
                    if(j==i+1) C(i,j)=0;
                }
            }

            else if(type==23)
            {
                bool on_edge_offdiag = ( is_on_edge(dim_edges,i) && j==(i-1)) || ( is_on_edge(dim_edges,j) && i==(j-1));
                //if (TMath::Abs(i-j) == 1 && on_edge_offdiag) { std::cout << "on_edge_offdiag.  i,j = " << i << ",  " << j << std::endl; }
                if(TMath::Abs(i-j) == 1 && !on_edge_offdiag) C(i, j) = 1;
                if(i==j) 
                {
                    bool on_edge = is_on_edge(dim_edges,i) || is_on_edge(dim_edges,i+1);
                    //if (on_edge) { std::cout << "on_edge.          i, j = " << i << ",  " << j << std::endl; }
                    if (on_edge) { C(i, j) = -1+epsilon; }
                    else         { C(i, j) = -2+epsilon; }
                }
            }

            else if(type==32 || type==33) // third derivative matrix
            {
                if (get_slice(dim_edges,i) != get_slice(dim_edges,j)) { C(i,j) = 0; continue; }

                bool ihigh = i>j;
                if      (std::abs(i-j) == 2) { C(i, j) = !ihigh - ihigh; }
                else if (std::abs(i-j) == 1) {
                    C(i, j) = 2*(ihigh - !ihigh);
                    if((is_on_edge(dim_edges,i-1) && ihigh) || (is_on_edge(dim_edges,j+1) && !ihigh)) { C(i, j) = 0; }
                }
                else if(i==j) {
                    C(i, j) = 0 + epsilon2;
                    if(is_on_edge(dim_edges,i)   || is_on_edge(dim_edges,i-1)) { C(i, j) = C(i, j) + 1; }
                    if(is_on_edge(dim_edges,i+1) || is_on_edge(dim_edges,i+2)) { C(i, j) = C(i, j) - 1; }
                }
            }

            else{
                if(i-j == 2){ 
                    C(i, j) = -1;
                }
                if(i-j == 1){
                    C(i, j) = 2;
                    if(i==1){
                        C(i, j) = 0; 
                    }
                }
                if(i-j == -1){
                    C(i, j) = -2;
                    if(i==n-2){
                        C(i, j) = 0;
                    }
                }
                if(i-j == -2){
                    C(i, j) = 1;
                }
                if(i==j){
                    C(i, j) = 0 + epsilon2;
                    if(i==0 || i==1){
                        C(i, j) = 1 + epsilon2;
                    }
                    if(i==n-1 || i==n-2){
                        C(i, j) = -1 + epsilon2;
                    }
                }
            }
        }
    }

    return C;
}



TVectorD WienerSVD(TMatrixD Response, TVectorD Signal, TVectorD Measure, TMatrixD Covariance, Int_t C_type, Float_t Norm_type, TMatrixD& AddSmear, TVectorD& WF, TMatrixD& UnfoldCov, TMatrixD& covRotation, TMatrixD& covRotation_t, Float_t flag_WienerFilter)
{
    Int_t m = Response.GetNrows(); // measure, M
    Int_t n = Response.GetNcols(); // signal, S
    // Decomposition of Covariance Matrix to get orthogonal Q to rotate the current frame, 
    // then make the uncertainty for each bin equal to 1
    TDecompSVD decV(Covariance);
    TMatrixD Q0 (TMatrixD::kTransposed, decV.GetV());
    TVectorD err0 = decV.GetSig();

    TMatrixD err(m, m);
    for(Int_t i=0; i<m; i++)
    {
        for(Int_t j=0; j<m; j++)
        {
            err(i, j) = 0;
            if(i == j)
            {
                if(err0(i)) err(i, j) = 1./TMath::Sqrt( err0(i) );
                else err(i, j) = 0;
            }
        }
    }

    TMatrixD Q = err*Q0;
    // transform Measure and Response
    TVectorD M = Q*Measure;
    TMatrixD R = Q*Response;

    // For addtion of smoothness matrix, e.g. 2nd derivative C
    TMatrixD C0(n, n);
    C0 = Matrix_C(n, C_type);
    TMatrixD normsig(n, n);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            normsig(i, j) = 0;
            if(i==j) normsig(i, j) = 1./TMath::Power(Signal(i), Norm_type);
        }
    }
    C0 = C0*normsig;

    TMatrixD C = C0;
    C0.Invert();
    TMatrixD C_inv = C0;
    Signal = C*Signal;
    R = R*C_inv;
  
    // SVD decomposition of R 
    TDecompSVD udv(R);
    TMatrixD U = udv.GetU();
    TMatrixD U_t (TMatrixD::kTransposed, U);
    TMatrixD V = udv.GetV();
    TMatrixD V_t (TMatrixD::kTransposed, V);
    TVectorD D = udv.GetSig();
    // for matrix D transverse
    TMatrixD D_t(n, m);
    for(Int_t i=0; i<n; i++)
    {
        for(Int_t j=0; j<m; j++)
        {
            D_t(i,j) = 0;
            if(i==j)
            {
                D_t(i, j) = D(i);
            }
        }
    }

    TVectorD S = V_t*Signal;
    // Wiener Filter 
    TMatrixD W(n, n);
    TMatrixD W0(n, n);
    for(Int_t i=0; i<n; i++)
    {
        for(Int_t j=0; j<n; j++)
        {
            W(i, j) = 0;
            W0(i, j) = 0;
            if(i == j)
            {
                if(flag_WienerFilter!=0){
                //W(i, j) = 1./(D(i)*D(i)+2e-7); //S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                //WF(i) = D(i)*D(i)*W(i, j);//S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                W(i, j) = flag_WienerFilter*S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 ); // flag_WienerFilter -> normalization
                WF(i) = D(i)*D(i)*W(i, j);
                W0(i, j) = WF(i); 
                }
                else{
                W(i, j) = 1.0/D(i)/D(i);
                WF(i) = 1.0; // = W(i, j)*D(i)*D(i);
                W0(i, j) = WF(i);
                }
            }
        }
    }

    TVectorD unfold = C_inv*V*W*D_t*U_t*M;
    AddSmear = C_inv*V*W0*V_t*C;
    // covariance matrix of the unfolded spectrum
    //TMatrixD covRotation = C_inv*V*W*D_t*U_t*Q;
    //TMatrixD covRotation_t (TMatrixD::kTransposed, covRotation); 
    covRotation = C_inv*V*W*D_t*U_t*Q;
    covRotation_t.Transpose(covRotation);
    UnfoldCov = covRotation*Covariance*covRotation_t; 


    return unfold;
}
