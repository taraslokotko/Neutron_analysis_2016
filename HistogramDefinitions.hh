#ifndef HISTOGRAM_DEFINITIONS_HH
#define HISTOGRAM_DEFINITIONS_HH

#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <map>

// Namespace for histograms
namespace Histos {
    // Track properties
    TH2D* histPID = new TH2D("histPID", ";p/Z (MeV/Z); dEdx (ADC/mm)", 800, -500, 2500, 800, 0, 2000);
    TH2D* histPID_va = new TH2D("histPID_va", ";p/Z (MeV/Z); dEdx (ADC/mm)", 800, -500, 2500, 800, 0, 2000);
    TH2D* histPID_va_left = new TH2D("histPID_va_left", ";p/Z (MeV/Z); dEdx (ADC/mm)", 800, -500, 2500, 800, 0, 2000);
    TH2D* histPID_va_right = new TH2D("histPID_va_right", ";p/Z (MeV/Z); dEdx (ADC/mm)", 800, -500, 2500, 800, 0, 2000);

    // Angular distributions
    TH2D* histEmitAngle = new TH2D("histEmitAngle", ";#Theta (Deg.); #Phi (Deg.)", 90, 0, 90, 180, -180, 180);
    TH2D* histEmitAngle_left = new TH2D("histEmitAngle_left", ";#Theta (Deg.); #Phi (Deg.)", 90, 0, 90, 180, -180, 180);
    TH2D* histEmitAngle_right = new TH2D("histEmitAngle_right", ";#Theta (Deg.); #Phi (Deg.)", 90, 0, 90, 180, -180, 180);

    auto histEmitAngle_CM_pion_minus = new TH2D("histEmitAngle_CM_pion_minus",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);
    auto histEmitAngle_CM_pion_plus = new TH2D("histEmitAngle_CM_pion_plus",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);
    auto histEmitAngle_CM_proton = new TH2D("histEmitAngle_CM_proton",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);
    auto histEmitAngle_CM_deuteron = new TH2D("histEmitAngle_CM_deuteron",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);
    auto histEmitAngle_CM_triton = new TH2D("histEmitAngle_CM_triton",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);
    auto histEmitAngle_CM_3he = new TH2D("histEmitAngle_CM_3he",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);
    auto histEmitAngle_CM_4he = new TH2D("histEmitAngle_CM_4he",";#Theta_{CM} (Deg.); #Phi_{CM} (Deg.)",180,0,180,180,0,360);

    auto histEmitAngle_Theta_LAB_vs_CM_pion_minus = new TH2D("histEmitAngle_Theta_LAB_vs_CM_pion_minus",";#Theta_{CM} (Deg.); #Theta_{LAB} (Deg.)",180,-360,360,180,-360,360);
    auto histEmitAngle_Phi_LAB_vs_CM_pion_minus = new TH2D("histEmitAngle_Phi_LAB_vs_CM_pion_minus",";#Phi_{CM} (Deg.); #Phi_{LAB} (Deg.)",180,-360,360,180,-360,360);

    // Vertex distributions
    TH1D* histva_dist = new TH1D("histva_dist", ";Counts; Distance [mm]", 50, 0, 25);
    TH1D* histva_dist_left = new TH1D("histva_dist_left", ";Counts; Distance [mm]", 50, 0, 25);
    TH1D* histva_dist_right = new TH1D("histva_dist_right", ";Counts; Distance [mm]", 50, 0, 25);

    // Particle distributions
    TH1D* histMultiplicity_real = new TH1D("histMultiplicity_real", ";Counts; Multiplicity", 100, 0, 100);
    TH2D* histPt_vs_rapidity = new TH2D("histPt_vs_rapidity", ";Rapidity; Pt [MeV]", 100, 0, 1, 500, 0, 1000);
    TH2D* histPt_vs_rapidity_pion_minus = new TH2D("histPt_vs_rapidity_pion_minus", ";Rapidity; Pt [MeV]", 100, 0, 2, 500, 0, 500);
    TH2D* histPt_vs_rapidity_pion_plus = new TH2D("histPt_vs_rapidity_pion_plus", ";Rapidity; Pt [MeV]", 100, 0, 2, 500, 0, 500);
    TH2D* histPt_vs_rapidity_proton = new TH2D("histPt_vs_rapidity_proton", ";Rapidity; Pt [MeV]", 100, 0, 1, 500, 0, 1000);
    TH2D* histPt_vs_rapidity_deuteron = new TH2D("histPt_vs_rapidity_deuteron", ";Rapidity; Pt [MeV]", 100, 0, 1, 500, 0, 1000);
    TH2D* histPt_vs_rapidity_triton = new TH2D("histPt_vs_rapidity_triton", ";Rapidity; Pt [MeV]", 100, 0, 1, 500, 0, 1000);
    TH2D* histPt_vs_rapidity_3he = new TH2D("histPt_vs_rapidity_3he", ";Rapidity; Pt [MeV]", 100, 0, 1, 500, 0, 1000);
    TH2D* histPt_vs_rapidity_4he = new TH2D("histPt_vs_rapidity_4he", ";Rapidity; Pt [MeV]", 100, 0, 1, 500, 0, 1000);

    //Normalized spectra
    TH2D* histPt_vs_rapidity_norm_pion_minus = new TH2D("histPt_vs_rapidity_norm_pion_minus", ";Rapidity norm; Pt [MeV]", 100,-3.5,3.5,100,0,500);
    TH2D* histPt_vs_rapidity_norm_pion_plus = new TH2D("histPt_vs_rapidity_norm_pion_plus", ";Rapidity norm; Pt [MeV]", 100,-3.5,3.5,100,0,500);
    TH2D* histPt_vs_rapidity_norm_proton = new TH2D("histPt_vs_rapidity_norm_proton", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500);
    TH2D* histPt_vs_rapidity_norm_deuteron = new TH2D("histPt_vs_rapidity_norm_deuteron", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500);
    TH2D* histPt_vs_rapidity_norm_triton = new TH2D("histPt_vs_rapidity_norm_triton", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500);
    TH2D* histPt_vs_rapidity_norm_3he = new TH2D("histPt_vs_rapidity_norm_3he", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500);
    TH2D* histPt_vs_rapidity_norm_4he = new TH2D("histPt_vs_rapidity_norm_4he", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500);

    auto histPt_vs_rapidity_norm_pion_minus_left = new TH2D("histPt_vs_rapidity_norm_pion_minus_left", ";Rapidity norm; Pt [MeV]", 100,-3.5,3.5,100,0,500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_pion_plus_left = new TH2D("histPt_vs_rapidity_norm_pion_plus_left", ";Rapidity norm; Pt [MeV]", 100,-3.5,3.5,100,0,500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_proton_left = new TH2D("histPt_vs_rapidity_norm_proton_left", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_deuteron_left = new TH2D("histPt_vs_rapidity_norm_deuteron_left", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_triton_left = new TH2D("histPt_vs_rapidity_norm_triton_left", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_3he_left = new TH2D("histPt_vs_rapidity_norm_3he_left", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_4he_left = new TH2D("histPt_vs_rapidity_norm_4he_left", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra

    auto histPt_vs_rapidity_norm_pion_minus_right = new TH2D("histPt_vs_rapidity_norm_pion_minus_right", ";Rapidity norm; Pt [MeV]", 100,-3.5,3.5,100,0,500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_pion_plus_right = new TH2D("histPt_vs_rapidity_norm_pion_plus_right", ";Rapidity norm; Pt [MeV]", 100,-3.5,3.5,100,0,500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_proton_right = new TH2D("histPt_vs_rapidity_norm_proton_right", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_deuteron_right = new TH2D("histPt_vs_rapidity_norm_deuteron_right", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_triton_right = new TH2D("histPt_vs_rapidity_norm_triton_right", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_3he_right = new TH2D("histPt_vs_rapidity_norm_3he_right", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra
    auto histPt_vs_rapidity_norm_4he_right = new TH2D("histPt_vs_rapidity_norm_4he_right", ";Rapidity; Pt [MeV]", 100,-3.5,3.5,300,0,1500); //Rapidity vs Pt spectra

    // Beam histograms
    TH2D* histBeamAngle = new TH2D("histBeamAngle",";#Theta (mRad); #Phi (mRad)",80,-40,40,100,-100,0);

    // Gather all histograms into a vector for easy iteration
    std::vector<TH1*> allHistos = {
	Histos::histPID, Histos::histPID_va, Histos::histPID_va_left, Histos::histPID_va_right,
    Histos::histEmitAngle, Histos::histEmitAngle_left, Histos::histEmitAngle_right,
    Histos::histva_dist, Histos::histva_dist_left, Histos::histva_dist_right,
    Histos::histPt_vs_rapidity, Histos::histPt_vs_rapidity_pion_minus, Histos::histPt_vs_rapidity_pion_plus,
    Histos::histPt_vs_rapidity_proton, Histos::histPt_vs_rapidity_deuteron, Histos::histPt_vs_rapidity_triton,
    Histos::histPt_vs_rapidity_3he, Histos::histPt_vs_rapidity_4he, histBeamAngle, histPt_vs_rapidity_norm_pion_minus, histPt_vs_rapidity_norm_pion_plus,
	histPt_vs_rapidity_norm_proton, histPt_vs_rapidity_norm_deuteron, histPt_vs_rapidity_norm_triton, histPt_vs_rapidity_norm_3he, histPt_vs_rapidity_norm_4he,
    histEmitAngle_CM_pion_minus, histEmitAngle_CM_pion_plus, histEmitAngle_CM_proton, histEmitAngle_CM_deuteron, histEmitAngle_CM_triton, histEmitAngle_CM_3he, histEmitAngle_CM_4he,
    histPt_vs_rapidity_norm_pion_minus_left, histPt_vs_rapidity_norm_pion_plus_left, histPt_vs_rapidity_norm_pion_minus_right, histPt_vs_rapidity_norm_pion_plus_right, histPt_vs_rapidity_norm_proton_left, 
    histPt_vs_rapidity_norm_deuteron_left, histPt_vs_rapidity_norm_triton_left, histPt_vs_rapidity_norm_3he_left, histPt_vs_rapidity_norm_4he_left, histPt_vs_rapidity_norm_proton_right,
    histPt_vs_rapidity_norm_deuteron_right, histPt_vs_rapidity_norm_triton_right, histPt_vs_rapidity_norm_3he_right, histPt_vs_rapidity_norm_4he_right, histMultiplicity_real, 
    histEmitAngle_Theta_LAB_vs_CM_pion_minus, histEmitAngle_Phi_LAB_vs_CM_pion_minus
    };

    TH1* GetHistogram(const std::string& name) {
        for (auto& hist : allHistos) {
            if (hist && name == hist->GetName()) {
                return hist;
            }
        }
        std::cerr << "Histogram " << name << " not found!" << std::endl;
        return nullptr;
    }

}
// Function to write all histograms to a file
void WriteHistogramsToFile(TFile *outputFile) {
    outputFile->cd();
    for (auto& hist : Histos::allHistos) {
        if (hist) hist->Write();
    }
}
#endif // HISTOGRAM_DEFINITIONS_HH

