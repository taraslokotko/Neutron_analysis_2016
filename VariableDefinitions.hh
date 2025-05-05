// TreeDefinitions.hh
#ifndef TREE_DEFINITIONS_HH
#define TREE_DEFINITIONS_HH

#include <vector>
#include <TLorentzVector.h>
#include <TTree.h>

// Declare variables
//Beam variables
double Beam_Z = 0;
double Beam_AoQ = 0;
double BDC_X = 0;
double BDC_Y = 0;
double RotationAngleATargetPlane = 0;
double RotationAngleBTargetPlane = 0;
double XTargetPlane = 0;
double YTargetPlane = 0;
double BeamEnergyTargetPlane = 0;
double BeamVelocityTargetPlane = 0;
double Beam_energy = 0, Beam_momenta = 0, Beam_momenta_X = 0, Beam_momenta_Y = 0, Beam_momenta_Z = 0;

//Spirt tracks variables
int numRecoTrack = 0;
int numVATrack = 0;
double dedx_value = 0;
double momentum_mag = 0;
double reco_track_length = 0;
double reco_eff_curvature = 0;
int reco_ClusterNum = 0;
int reco_NumRowClusters = 0;
int reco_NumLayerClusters = 0;
double momentum_X = 0, momentum_Y = 0, momentum_Z = 0;
double MomentumTheta = 0, MomentumPhi = 0;
double dist = 0;
int Event_ID = 0;
int Track_ID = 0;
double Chi2 = 0;
double NDF = 0;
double Chi2_norm = 0;


//VA tracks variables
vector <double> *track_pos_X = new vector <double>;
vector <double> *track_pos_Y = new vector <double>;
vector <double> *track_pos_Z = new vector <double>;
vector <double> *va_Charge = new vector <double>;
vector <double> *va_gen_fit_charge = new vector <double>;
vector <double> *va_p = new vector <double>;
vector <double> *va_dedx = new vector <double>;
vector <double> *va_momentum_mag = new vector <double>;
vector <double> *va_momentum_x = new vector <double>;
vector <double> *va_momentum_y = new vector <double>;
vector <double> *va_momentum_z = new vector <double>;
vector <double> *va_rapidity = new vector <double>;
vector <double> *va_Pt = new vector <double>;
vector <double> *va_rapidity_beam_ref = new vector <double>;
vector <double> *va_Pt_beam_ref = new vector <double>;
vector <double> *va_rapidity_CM_ref = new vector <double>;
vector <double> *va_Pt_CM_ref = new vector <double>;
vector <double> *va_dist = new vector <double>;
vector <double> *va_MomentumPhi = new vector <double>;
vector <double> *va_MomentumTheta = new vector <double>;
vector <double> *va_MomentumPhi_pion_minus_CM = new vector <double>;
vector <double> *va_MomentumTheta_pion_minus_CM = new vector <double>;
vector <double> *va_MomentumPhi_pion_plus_CM = new vector <double>;
vector <double> *va_MomentumTheta_pion_plus_CM = new vector <double>;

vector <double> *va_MomentumPhi_proton_CM = new vector <double>;
vector <double> *va_MomentumTheta_proton_CM = new vector <double>;
vector <double> *va_MomentumPhi_deuteron_CM = new vector <double>;
vector <double> *va_MomentumTheta_deuteron_CM = new vector <double>;
vector <double> *va_MomentumPhi_triton_CM = new vector <double>;
vector <double> *va_MomentumTheta_triton_CM = new vector <double>;
vector <double> *va_MomentumPhi_He3_CM = new vector <double>;
vector <double> *va_MomentumTheta_He3_CM = new vector <double>;
vector <double> *va_MomentumPhi_He4_CM = new vector <double>;   
vector <double> *va_MomentumTheta_He4_CM = new vector <double>;




//kinematics variables for particles
//std::vector<TLorentzVector>* pion_minus_Vector = new std::vector<TLorentzVector>();
TLorentzVector beam_Vector, target_Vector, beam_Vector_beam_ref, beam_Vector_CM_ref;
TLorentzVector system_Vector, system_Vector_beam_ref, system_Vector_CM_ref;
TLorentzVector particle_Vector_pion_minus, particle_Vector_pion_plus, particle_Vector_proton, particle_Vector_deuteron, particle_Vector_triton, particle_Vector_3he, particle_Vector_4he;

TLorentzVector particle_Vector_pion_minus_beam_ref, particle_Vector_pion_plus_beam_ref, particle_Vector_proton_beam_ref, particle_Vector_deuteron_beam_ref, particle_Vector_triton_beam_ref, particle_Vector_3he_beam_ref, particle_Vector_4he_beam_ref;
TLorentzVector particle_Vector_pion_minus_CM_ref, particle_Vector_pion_plus_CM_ref, particle_Vector_proton_CM_ref, particle_Vector_deuteron_CM_ref, particle_Vector_triton_CM_ref, particle_Vector_3he_CM_ref, particle_Vector_4he_CM_ref;
double beam_rapidity, beam_rapidity_beam_ref, beam_rapidity_CM_ref;
double system_rapidity, system_rapidity_beam_ref, system_rapidity_CM_ref;

//Constants
const double amu = 931.5; //MeV/c2
const double mass_pion = 139.57039; //MeV/c2
const double mass_proton = 938.2723; //MeV/c2
const double mass_deuteron = 1875.612928; //MeV/c2
const double mass_triton = 2808.921136; //MeV/c2
const double mass_3he = 2808.391611; //MeV/c2
const double mass_4he = 3727.379; //MeV/c2
const double mass_124Xe = 123.8765*amu;
const double mass_136Xe = 135.8778*amu;
const double mass_124Sn = 123.878*amu;
const double mass_112Sn = 111.8776*amu;

// Declare the tree
TTree* Spirit_tree = nullptr;

// Function to initialize the tree and its branches
void InitializeTree() {
    Spirit_tree = new TTree("Spirit_tree", "Physics Analysis Tree");
    Spirit_tree->Branch("Event_ID", &Event_ID);
    Spirit_tree->Branch("Beam_Z", &Beam_Z);
    Spirit_tree->Branch("Beam_AoQ", &Beam_AoQ);
    Spirit_tree->Branch("BDC_X", &BDC_X);
    Spirit_tree->Branch("BDC_Y", &BDC_Y);
    Spirit_tree->Branch("RotationAngleATargetPlane", &RotationAngleATargetPlane);
    Spirit_tree->Branch("RotationAngleBTargetPlane", &RotationAngleBTargetPlane);
    Spirit_tree->Branch("XTargetPlane", &XTargetPlane);
    Spirit_tree->Branch("YTargetPlane", &YTargetPlane);
    Spirit_tree->Branch("BeamEnergyTargetPlane", &BeamEnergyTargetPlane);
    Spirit_tree->Branch("BeamVelocityTargetPlane", &BeamVelocityTargetPlane);

    Spirit_tree->Branch("numRecoTrack", &numRecoTrack);
    Spirit_tree->Branch("numVATrack", &numVATrack);
    Spirit_tree->Branch("dedx", &dedx_value);
    Spirit_tree->Branch("momentum_mag", &momentum_mag);
    Spirit_tree->Branch("reco_track_length", &reco_track_length);
    Spirit_tree->Branch("reco_eff_curvature", &reco_eff_curvature);
    Spirit_tree->Branch("reco_ClusterNum", &reco_ClusterNum);
    Spirit_tree->Branch("reco_NumRowClusters", &reco_NumRowClusters);
    Spirit_tree->Branch("reco_NumLayerClusters", &reco_NumLayerClusters);
    
    Spirit_tree->Branch("va_Charge", &va_Charge);
    Spirit_tree->Branch("va_dedx", &va_dedx);
    Spirit_tree->Branch("va_momentum_mag", &va_momentum_mag);
    Spirit_tree->Branch("va_gen_fit_charge", &va_gen_fit_charge);
    Spirit_tree->Branch("va_p", &va_p);
    Spirit_tree->Branch("va_momentum_mag", &va_momentum_mag);
    Spirit_tree->Branch("va_momentum_x", &va_momentum_x);
    Spirit_tree->Branch("va_momentum_y", &va_momentum_y);
    Spirit_tree->Branch("va_momentum_z", &va_momentum_z);
    Spirit_tree->Branch("va_rapidity", &va_rapidity);
    Spirit_tree->Branch("va_Pt", &va_Pt);
    Spirit_tree->Branch("va_rapidity_beam_ref", &va_rapidity_beam_ref);
    Spirit_tree->Branch("va_Pt_beam_ref", &va_Pt_beam_ref);
    Spirit_tree->Branch("va_rapidity_CM_ref", &va_rapidity_CM_ref);
    Spirit_tree->Branch("va_Pt_CM_ref", &va_Pt_CM_ref);
    Spirit_tree->Branch("va_dist", &va_dist);
    Spirit_tree->Branch("va_MomentumPhi", &va_MomentumPhi);
    Spirit_tree->Branch("va_MomentumTheta", &va_MomentumTheta);
    Spirit_tree->Branch("va_MomentumPhi_pion_minus_CM", &va_MomentumPhi_pion_minus_CM);
    Spirit_tree->Branch("va_MomentumTheta_pion_minus_CM", &va_MomentumTheta_pion_minus_CM);
    Spirit_tree->Branch("va_MomentumPhi_pion_plus_CM", &va_MomentumPhi_pion_plus_CM);
    Spirit_tree->Branch("va_MomentumTheta_pion_plus_CM", &va_MomentumTheta_pion_plus_CM);

    //Spirit_tree->Branch("pion_minus_Vector", &pion_minus_Vector);

    Spirit_tree->Branch("Beam_4Vect","TLorentzVector", &beam_Vector);
    Spirit_tree->Branch("Beam_4Vect_beam_ref","TLorentzVector", &beam_Vector_beam_ref);
    Spirit_tree->Branch("Beam_4Vect_CM_ref","TLorentzVector", &beam_Vector_CM_ref);

    Spirit_tree->Branch("Pion_minus_4Vect","TLorentzVector", &particle_Vector_pion_minus);
    Spirit_tree->Branch("Pion_plus_4Vect","TLorentzVector", &particle_Vector_pion_plus);
    Spirit_tree->Branch("Proton_4Vect","TLorentzVector", &particle_Vector_proton);
    Spirit_tree->Branch("Deuteron_4Vect","TLorentzVector", &particle_Vector_deuteron);
    Spirit_tree->Branch("Triton_4Vect","TLorentzVector", &particle_Vector_triton);
    Spirit_tree->Branch("3He_4Vect","TLorentzVector", &particle_Vector_3he);
    Spirit_tree->Branch("4He_4Vect","TLorentzVector", &particle_Vector_4he);

    Spirit_tree->Branch("Pion_minus_4Vect_beam_ref","TLorentzVector", &particle_Vector_pion_minus_beam_ref);
    Spirit_tree->Branch("Pion_plus_4Vect_beam_ref","TLorentzVector", &particle_Vector_pion_plus_beam_ref);
    Spirit_tree->Branch("Proton_4Vect_beam_ref","TLorentzVector", &particle_Vector_proton_beam_ref);
    Spirit_tree->Branch("Deuteron_4Vect_beam_ref","TLorentzVector", &particle_Vector_deuteron_beam_ref);
    Spirit_tree->Branch("Triton_4Vect_beam_ref","TLorentzVector", &particle_Vector_triton_beam_ref);
    Spirit_tree->Branch("3He_4Vect_beam_ref","TLorentzVector", &particle_Vector_3he_beam_ref);
    Spirit_tree->Branch("4He_4Vect_beam_ref","TLorentzVector", &particle_Vector_4he_beam_ref);

    Spirit_tree->Branch("Pion_minus_4Vect_CM_ref","TLorentzVector", &particle_Vector_pion_minus_CM_ref);
    Spirit_tree->Branch("Pion_plus_4Vect_CM_ref","TLorentzVector", &particle_Vector_pion_plus_CM_ref);
    Spirit_tree->Branch("Proton_4Vect_CM_ref","TLorentzVector", &particle_Vector_proton_CM_ref);
    Spirit_tree->Branch("Deuteron_4Vect_CM_ref","TLorentzVector", &particle_Vector_deuteron_CM_ref);
    Spirit_tree->Branch("Triton_4Vect_CM_ref","TLorentzVector", &particle_Vector_triton_CM_ref);
    Spirit_tree->Branch("3He_4Vect_CM_ref","TLorentzVector", &particle_Vector_3he_CM_ref);
    Spirit_tree->Branch("4He_4Vect_CM_ref","TLorentzVector", &particle_Vector_4he_CM_ref);
}


void ClearVectors() {
    va_Charge->clear();
    va_gen_fit_charge->clear();
    va_p->clear();
    va_dedx->clear();
    va_momentum_mag->clear();
    va_momentum_x->clear();
    va_momentum_y->clear();
    va_momentum_z->clear();
    va_rapidity->clear();
    va_Pt->clear();
    va_rapidity_beam_ref->clear();
    va_Pt_beam_ref->clear();
    va_rapidity_CM_ref->clear();
    va_Pt_CM_ref->clear(); 
    va_dist->clear();
    va_MomentumPhi->clear();
    va_MomentumTheta->clear();
    va_MomentumPhi_pion_minus_CM->clear();
    va_MomentumTheta_pion_minus_CM->clear();
    va_MomentumPhi_pion_plus_CM->clear();
    va_MomentumTheta_pion_plus_CM->clear();

    beam_Vector.SetPxPyPzE(NAN,NAN,NAN,NAN);
    target_Vector.SetPxPyPzE(NAN,NAN,NAN,NAN);
    system_Vector.SetPxPyPzE(NAN,NAN,NAN,NAN);
    system_Vector_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    beam_Vector_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);

    particle_Vector_pion_minus.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_pion_plus.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_proton.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_deuteron.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_triton.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_3he.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_4he.SetPxPyPzE(NAN,NAN,NAN,NAN);

    particle_Vector_pion_minus_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_pion_plus_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_proton_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_deuteron_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_triton_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_3he_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_4he_beam_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);

    particle_Vector_pion_minus_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_pion_plus_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_proton_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_deuteron_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_triton_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_3he_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
    particle_Vector_4he_CM_ref.SetPxPyPzE(NAN,NAN,NAN,NAN);
}


#endif // TREE_DEFINITIONS_HH

