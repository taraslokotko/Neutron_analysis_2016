#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TClonesArray.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <iostream>
#include <sys/stat.h> // For file size check
#include "/data/RB230064/common/spiritroot/Hokusai/SpiRITROOT_2024Spring/ana/ST_ClusterNum_DB.hh"
#include "/data/RB230064/common/spiritroot/Hokusai/SpiRITROOT_2024Spring/ana/STPIDAnalysisTask.hh"
#include "/data/RB230064/common/spiritroot/Hokusai/SpiRITROOT_2024Spring/ana/STConcReaderTask.hh"
#include "/data/RB230064/common/spiritroot/Hokusai/SpiRITROOT_2024Spring/ana/ST_ProduceDB_ClusterNum.hh"
#include "/data/RB230064/common/spiritroot/Hokusai/SpiRITROOT_2024Spring/ana/ST_TrackCut.hh"
#include "/data/RB230064/common/spiritroot/Hokusai/SpiRITROOT_2024Spring/ana/STAnaLinkDef.h"
#include "HistogramDefinitions.hh"
#include "VariableDefinitions.hh"
#include "ProcessSpectra.h"

Int_t fClustCut = 15;
Int_t fTrackCut = 40;
Double_t fPOCACut = 5;
Double_t fVertZPosCut = -32.5347;
Double_t fVertZSigCut = 2.2037;
Double_t fTargetXWidth = 30;
Double_t fThreshold = 5;

void Analysis_PID_neutron(int run_number_initial = 1163, int run_number_final = 1296, bool server_type = 0 )
{
        /*server_type = 0 corresponds to HOKUSAI server, = 1 corresponds to FISHTANK*/

        auto tree = new TChain("cbmsim");

        std::vector<int> runNumbers; // Add or remove run numbers as needed
        for (int run = run_number_initial; run <= run_number_final; ++run)
        {
                runNumbers.push_back(run);
        }
        std::vector<std::array<TString, 90>> filenames(runNumbers.size());

        if(server_type == 0)
        {
                for (size_t j = 0; j < runNumbers.size(); ++j)
                {
                        for (int i = 0; i < 90; ++i)
                        {
                                filenames[j][i] = Form("/data/RB230064/common_2024/reconstruction_master/data/run%d_s%02d.reco.2024.root", runNumbers[j], i);
                                tree->AddFile(filenames[j][i]);
                        }
                }
        }
        else
        {
                for (size_t j = 0; j < runNumbers.size(); ++j)
                {
                        for (int i = 0; i < 6; ++i)
                        {
                                if(run_number_initial>=2272 && run_number_final<=2509)
                                {
                                        filenames[j][i] = Form("/mnt/spirit/analysis/user/tsangc/SpiRITROOT/macros/data/Sn108/run%d_s%02d.reco.develop.1988.bf2b00e.conc.root", runNumbers[j], i);
                                        tree->AddFile(filenames[j][i]);
                                }
                        }
                }
        }

        TClonesArray *RecoTrackArray = nullptr;
        TClonesArray *HelixArray = nullptr;
        TClonesArray *VertexArray = nullptr;
        TClonesArray *HitArray = nullptr;
        TClonesArray *VATrackArray = nullptr;
        TClonesArray *BeamData = nullptr;
        TClonesArray *BDCData = nullptr;

        // tree -> SetBranchAddress("STHelixTrack", &HelixArray);
        tree->SetBranchAddress("STRecoTrack", &RecoTrackArray);
        tree->SetBranchAddress("VATracks", &VATrackArray);
        tree->SetBranchAddress("STVertex", &VertexArray);
        // tree -> SetBranchAddress("STHit", &HitArray); //Don't see this variable in the tree structure, probably wrtting the data to the file was scipted?
        tree->SetBranchAddress("STBeamInfo", &BeamData);
        tree->SetBranchAddress("BDCVertex", &BDCData);

        // SPIRIT cuts:
        TFile *cutImportProtons = new TFile("/data/RB230064/common_2024/PIDcut/ProtonCutG.root", "READ");
        TFile *cutImportDeuterons = new TFile("/data/RB230064/common_2024/PIDcut/DeuteronCutG.root", "READ");
        TFile *cutImportTritons = new TFile("/data/RB230064/common_2024/PIDcut/TritonCutG.root", "READ");
        TFile *cutImportHe3 = new TFile("/data/RB230064/common_2024/PIDcut/He3CutG.root", "READ");
        TFile *cutImportHe4 = new TFile("/data/RB230064/common_2024/PIDcut/AlphaCutG.root", "READ");
        TFile *cutImportHe6 = new TFile("/data/RB230064/common_2024/PIDcut/He6CutG.root", "READ");
        TFile *cutImportLi6 = new TFile("/data/RB230064/common_2024/PIDcut/Li6CutG.root", "READ");
        TFile *cutImportLi7 = new TFile("/data/RB230064/common_2024/PIDcut/Li7CutG.root", "READ");
        TFile *cutImportPion_minus = new TFile("/data/RB230064/common_2024/PIDcut/Pion_minusCutG.root", "READ");
        TFile *cutImportPion_plus = new TFile("/data/RB230064/common_2024/PIDcut/Pion_plusCutG.root", "READ");

        TCutG *cutProton, *cutDeuteron, *cutTriton, *cutHe3, *cutHe4, *cutHe6, *cutLi6, *cutLi7, *cutpion_plus, *cutpion_minus;
        cutImportProtons->GetObject("ProtonCutG", cutProton);
        cutImportDeuterons->GetObject("DeuteronCutG", cutDeuteron);
        cutImportTritons->GetObject("TritonCutG", cutTriton);
        cutImportHe3->GetObject("He3CutG", cutHe3);
        cutImportHe4->GetObject("AlphaCutG", cutHe4);
        cutImportHe6->GetObject("He6CutG", cutHe6);
        cutImportLi6->GetObject("Li6CutG", cutLi6);
        cutImportLi7->GetObject("Li7CutG", cutLi7);
        cutImportPion_plus->GetObject("Pion_plusCutG", cutpion_plus);
        cutImportPion_minus->GetObject("Pion_minusCutG", cutpion_minus);

        // TFile *rootfile = new TFile("./test.root","RECREATE");
        TFile *rootfile = new TFile(Form("./Data_124Xe_exp_runs_%d_%d.root", run_number_initial, run_number_final), "RECREATE");
        rootfile->cd();
        // TTree *Spirit_tree = new TTree("Spirit_tree","analyzed tree");

        InitializeTree();

        Int_t nEvents = tree->GetEntries();
        //nEvents = 1000; // uncomment for testing
        cout << "EvtNum: " << nEvents << endl;

        for (int iEvent = 0; iEvent < nEvents; iEvent++)
        {
                if (iEvent % 10000 == 0)
                {
                        double ratio = (static_cast<double>(iEvent) / nEvents) * 100;
                        cout << "Evt: " << iEvent << "; Percentage complete: " << ratio << " % " << endl;
                }
                tree->GetEntry(iEvent);
                Event_ID = iEvent;

                ClearVectors();

                if (VertexArray->GetEntries() == 0)
                {
                        continue;
                }
                STVertex *aVertex = (STVertex *)VertexArray->At(0);

                // auto IsCollisionVertex_correct = aVertex->IsCollisionVertex();
                //    if(aVertex->IsTargetVertex() == 0) { continue; }
                // auto IsCollisionVertex_correct = aVertex->IsTargetVertex();
                // cout<<"IsCollisionVertex_correct = "<<IsCollisionVertex_correct<<endl;

                if ((aVertex->GetPos().Z() < fVertZPosCut - fVertZSigCut * 5 || aVertex->GetPos().Z() > fVertZPosCut + fVertZSigCut * 5))
                        continue;
                if ((aVertex->GetPos().X() < -fTargetXWidth / 2. || aVertex->GetPos().X() > fTargetXWidth / 2.))
                        continue;

                numRecoTrack = RecoTrackArray->GetEntries();
                numVATrack = VATrackArray->GetEntries();
                char first;
                // cout<<"numRecoTrack = "<<numRecoTrack<<endl;
                // cout<<"numVATrack = "<<numVATrack<<endl;

                // Cut on the numner of reco tracks, so the multiplicity
                //    if(numRecoTrack<15){continue; } //So far - opened cut, so all the tracks are taken

                // Beam Data:
                STBeamInfo *beam_event = (STBeamInfo *)BeamData;
                Beam_Z = beam_event->fBeamZ;
                Beam_AoQ = beam_event->fBeamAoQ;
                Beam_energy = beam_event->fBeamEnergyTargetPlane;
                RotationAngleATargetPlane = beam_event->fRotationAngleATargetPlane;
                RotationAngleBTargetPlane = beam_event->fRotationAngleBTargetPlane;
                XTargetPlane = beam_event->fXTargetPlane;
                YTargetPlane = beam_event->fYTargetPlane;
                BeamEnergyTargetPlane = beam_event->fBeamEnergyTargetPlane;
                BeamVelocityTargetPlane = beam_event->fBeamVelocityTargetPlane;

                RotationAngleATargetPlane = -45.63 / 1000; // Test. Inputting average angle for 1197 run
                RotationAngleBTargetPlane = 0.;

                // Inputting manually Energy and Momenta for the beam at the center of the target
                // 124Xe exp
                Beam_energy = 319.5416 * 124;                                             // Kinetic beam energy; MeV
                Beam_momenta = TMath::Sqrt(Beam_energy * (Beam_energy + 2 * mass_124Xe)); // Momenta Scalar of a beam particle
                Beam_momenta_X = Beam_momenta * TMath::Cos(RotationAngleBTargetPlane) * TMath::Sin(RotationAngleATargetPlane);
                Beam_momenta_Y = Beam_momenta * TMath::Sin(RotationAngleBTargetPlane) * TMath::Sin(RotationAngleATargetPlane);
                Beam_momenta_Z = Beam_momenta * TMath::Cos(RotationAngleATargetPlane);

                //beam_Vector.SetPxPyPzE(Beam_momenta_X, Beam_momenta_Y, Beam_momenta_Z, Beam_energy + mass_124Xe); // Momenta in Lab Frame
                beam_Vector.SetPxPyPzE(Beam_momenta_X, Beam_momenta_Y, Beam_momenta_Z, TMath::Sqrt(Beam_momenta*Beam_momenta + mass_124Xe*mass_124Xe)); // Momenta in Lab Frame
                target_Vector.SetPxPyPzE(0, 0, 0, mass_112Sn);                                                    // Momenta of Target in Lab Frame

                beam_Vector_beam_ref = beam_Vector;

                TVector3 beam_3D = beam_Vector_beam_ref.Vect();

                TVector3 newZ = beam_3D.Unit();
                TVector3 oldZ(0, 0, 1);
                TVector3 rotationAxis = oldZ.Cross(newZ);
                double rotationAngle = acos(oldZ.Dot(newZ));
                beam_Vector_beam_ref.Rotate(-rotationAngle, rotationAxis);

                system_Vector = beam_Vector_beam_ref + target_Vector;

                TVector3 System_BoostVector = system_Vector.BoostVector(); // Boost Vector to CM frame.

                system_Vector_beam_ref = system_Vector;

                system_Vector_CM_ref = system_Vector_beam_ref;
                system_Vector_CM_ref.Boost(-System_BoostVector);
                beam_Vector_CM_ref = beam_Vector_beam_ref;
                beam_Vector_CM_ref.Boost(-System_BoostVector);

                // beam_rapidity = 0.5*TMath::Log( (beam_Vector_beam_ref.E() + beam_Vector_beam_ref.Pz() )/(beam_Vector_beam_ref.E() - beam_Vector_beam_ref.Pz()));

                beam_rapidity = beam_Vector.Rapidity();
                beam_rapidity_beam_ref = beam_Vector_beam_ref.Rapidity();
                beam_rapidity_CM_ref = beam_Vector_CM_ref.Rapidity();

                // cout<<"beam_rapidity = "<<beam_rapidity<<";	beam_rapidity_beam_ref = "<<beam_rapidity_beam_ref<<";		beam_rapidity_CM_ref = "<<beam_rapidity_CM_ref<<endl;

                system_rapidity_CM_ref = 0.5 * TMath::Log((system_Vector_CM_ref.E() + system_Vector_CM_ref.Pz()) / (system_Vector_CM_ref.E() - system_Vector_CM_ref.Pz()));

                // BDC Data: not working......
                // BDCVertex* BDC_event = (BDCVertex*) BDCData;// -> At(0);
                // STVertex* BDC_event = (STVertex*) BDCData;
                // TVector3 BDC_X_mod = static_cast<STVertex*>(BDCData -> At(0))->GetPos();

                // map that stores id of reco track
                // will be used to map VATracks back to STRecoTracks

                std::map<int, STRecoTrack *> RecoToVATracks;
                if (VATrackArray)
                {
                        for (int ii = 0; ii < VATrackArray->GetEntries(); ++ii)
                        {
                                auto bdc_track = static_cast<STRecoTrack *>(VATrackArray->At(ii)); // bdc_track - track in SPIRIT TPC with vertex reconstructed with BDC projection
                                RecoToVATracks[bdc_track->GetRecoID()] = bdc_track;
                                // cout<<"bdc_track->GetRecoID() = "<<bdc_track->GetRecoID()<<endl;
                        }
                }
                int real_track_number = 0;

                for (int iTrack = 0; iTrack < numRecoTrack; iTrack++)
                {
                        STRecoTrack *aTrack_tmp = (STRecoTrack *)RecoTrackArray->At(iTrack);
                        vector<Int_t> *fHitClusterIDArray_tmp = aTrack_tmp->GetClusterIDArray();
                        reco_ClusterNum = (*fHitClusterIDArray_tmp).size();
                        reco_NumRowClusters = aTrack_tmp->GetNumRowClusters();
                        reco_NumLayerClusters = aTrack_tmp->GetNumLayerClusters();
                        if (reco_ClusterNum >= fClustCut)
                        {
                                if ((reco_NumRowClusters + reco_NumLayerClusters) >= reco_ClusterNum && (aTrack_tmp->GetPOCAVertex() - aVertex->GetPos()).Mag() <= fPOCACut)
                                {
                                        real_track_number++;
                                }
                        }
                }
                Histos::histMultiplicity_real->Fill(real_track_number);

                if (real_track_number <= fTrackCut)
                {
                        continue;
                }

                if (numRecoTrack <= fTrackCut)
                {
                        // continue;
                }

                for (int iTrack = 0; iTrack < numRecoTrack; iTrack++)
                {
                        // cout<<"Analyzing track-by-track"<<endl;
                        // What do I want?
                        // 1. Draw clusters and Pads from a certain event.
                        // 2. Draw tracks from a certain event.
                        // 3. Draw tracks fitted with BDC position.
                        // What's the difference between momentum and momentum-target-plane?

                        /*
                        aVertex - element of the class STVertex
                        aHelixTrack - element of the class STHelixTrack.
                          Cannot call STHelixTrack, maybe the data is not written into the file save the space?
                        aTrack - element of the class STRecoTrack.
                        */

                        Track_ID = iTrack;
                        STRecoTrack *aTrack = (STRecoTrack *)RecoTrackArray->At(iTrack);

                        // VATracks* aTrack = (VATracks*) VATrackArray -> At(iTrack);

                        // if(aTrackCut->IsTrackInVertex(aTrack,aVertex)==0) { continue; }
                        vector<Int_t> *fHitClusterIDArray = aTrack->GetClusterIDArray();
                        reco_ClusterNum = (*fHitClusterIDArray).size();
                        reco_NumRowClusters = aTrack->GetNumRowClusters();
                        reco_NumLayerClusters = aTrack->GetNumLayerClusters();
                        reco_track_length = aTrack->GetTrackLength();
                        reco_eff_curvature = aTrack->GetEffCurvature3();

                        Chi2 = aTrack->GetChi2();
                        NDF = aTrack->GetNDF();
                        Chi2_norm = Chi2 / NDF;

                        int HelixID = aTrack->GetHelixID();

                        double Charge = aTrack->GetCharge();
                        auto p = aTrack->GetMomentum().Mag();
                        auto dedx = aTrack->GetdEdxWithCut(0, 0.7, 0.5); // What does 0.7 mean?
                        momentum_X = aTrack->GetMomentum().X();
                        momentum_Y = aTrack->GetMomentum().Y();
                        momentum_Z = aTrack->GetMomentum().Z();

                        dedx_value = dedx;
                        momentum_mag = p / Charge;

                        // Applying proton cut on SPIRIT data:
                        // if(!cutpion_minus -> IsInside(p/Charge, dedx)){
                        //    continue;
                        //  }

                        auto Helix_track = aTrack->GetHelixTrack();
                        auto Genfit_track = aTrack->GetGenfitTrack();

                        // construct BDC data with identical vector range
                        int reco_id = aTrack->GetRecoID();
                        auto it_track = RecoToVATracks.find(reco_id);
                        // fill VA branches if data is found
                        if (it_track != RecoToVATracks.end())
                        {
                                if (reco_ClusterNum < fClustCut) continue;
                                if ((aTrack->GetPOCAVertex() - aVertex->GetPos()).Mag() > fPOCACut) continue;
                                
                                //cout<<"POCA diff = "<<(aTrack->GetPOCAVertex() - aVertex->GetPos()).Mag()<<endl;
                                //cout<<"reco_ClusterNum = "<<reco_ClusterNum<<endl;
                                //cout<<"                  "<<endl;
                                va_dist->push_back((aTrack->GetPOCAVertex() - aVertex->GetPos()).Mag());
                                auto VATrack = it_track->second;
                                va_Charge->push_back(VATrack->GetCharge());
                                va_gen_fit_charge->push_back(VATrack->GetGenfitCharge());
                                va_p->push_back(VATrack->GetMomentumTargetPlane().Mag());
                                va_dedx->push_back(VATrack->GetdEdxWithCut(0, 0.7, 0.5)); // What does 0.7 mean?
                                va_momentum_mag->push_back(va_p->back() / va_gen_fit_charge->back());
                                va_momentum_x->push_back(VATrack->GetMomentumTargetPlane().X());
                                va_momentum_y->push_back(VATrack->GetMomentumTargetPlane().Y());
                                va_momentum_z->push_back(VATrack->GetMomentumTargetPlane().Z());
                                va_MomentumPhi->push_back(VATrack->GetMomentumTargetPlane().Phi() * 180.0 / Pi());
                                va_MomentumTheta->push_back(VATrack->GetMomentumTargetPlane().Theta() * 180.0 / Pi());
                                va_Pt->push_back(VATrack->GetMomentumTargetPlane().Pt());

                                auto dEdxPointArray = VATrack->GetdEdxPointArray(); //.fPosition.x();
                                // cout <<"typeid(dEdxPointArray).name() = "<< typeid(dEdxPointArray).name() << endl;
                                for (int i = 0; i < dEdxPointArray->size(); i++)
                                {
                                        auto var = dEdxPointArray->at(i);
                                        track_pos_X->push_back(var.fPosition.x());
                                        track_pos_Y->push_back(var.fPosition.y());
                                        track_pos_Z->push_back(var.fPosition.z());
                                }
                        
                                Histos::histEmitAngle->Fill(va_MomentumTheta->back(), va_MomentumPhi->back());
                                //Remember: Below angle cuts are made in LAB FRAME.
                                //Rapidity analysis is performed in CM Frame for angles.
                                if ((va_MomentumPhi->back() > -40 && va_MomentumPhi->back() < 25) ||
                                    (va_MomentumPhi->back() > 160 && va_MomentumPhi->back() < 180) ||
                                    (va_MomentumPhi->back() > -180 && va_MomentumPhi->back() < -150))
                                {
                                        Histos::histva_dist->Fill(va_dist->back());
                                        Histos::histPID_va->Fill(va_p->back() / va_gen_fit_charge->back(), va_dedx->back());                             
                                }

                                // beam left
                                if (va_MomentumPhi->back() > -30 && va_MomentumPhi->back() < 20)
                                {
                                        Histos::histPID_va_left->Fill(va_p->back() / va_gen_fit_charge->back(), va_dedx->back());
                                        Histos::histEmitAngle_left->Fill(va_MomentumTheta->back(), va_MomentumPhi->back());
                                        Histos::histva_dist_left->Fill(va_dist->back());
                                }

                                // beam right
                                if ((va_MomentumPhi->back() > 160 && va_MomentumPhi->back() < 180) ||
                                    (va_MomentumPhi->back() > -180 && va_MomentumPhi->back() < -150))
                                {
                                        Histos::histPID_va_right->Fill(va_p->back() / va_gen_fit_charge->back(), va_dedx->back());
                                        Histos::histEmitAngle_right->Fill(va_MomentumTheta->back(), va_MomentumPhi->back());
                                        Histos::histva_dist_right->Fill(va_dist->back());
                                }

                                ProcessSpectra("pion_minus", cutpion_minus, mass_pion, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_pion_minus_CM, va_MomentumPhi_pion_minus_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);
                                ProcessSpectra("pion_plus", cutpion_plus, mass_pion, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_pion_plus_CM, va_MomentumPhi_pion_plus_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);
                                ProcessSpectra("proton", cutProton, mass_proton, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_proton_CM, va_MomentumPhi_proton_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);
                                ProcessSpectra("deuteron", cutDeuteron, mass_deuteron, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_deuteron_CM, va_MomentumPhi_deuteron_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);
                                ProcessSpectra("triton", cutTriton, mass_triton, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_triton_CM, va_MomentumPhi_triton_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);
                                ProcessSpectra("3he", cutHe3, mass_3he, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_He3_CM, va_MomentumPhi_He3_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);
                                ProcessSpectra("4he", cutHe4, mass_4he, rotationAngle, rotationAxis, System_BoostVector, beam_rapidity_CM_ref, va_MomentumTheta_He4_CM, va_MomentumPhi_He4_CM, va_Pt_CM_ref, va_Pt_beam_ref, va_rapidity_CM_ref, va_rapidity_beam_ref, VATrack);        
                        }
                }

                //Spirit_tree->Fill();
        }
        Spirit_tree->Write();
        WriteHistogramsToFile(rootfile);
        cutProton->Write();
        cutDeuteron->Write();
        cutTriton->Write();
        cutHe3->Write();
        cutHe4->Write();
        cutHe6->Write();
        cutLi6->Write();
        cutLi7->Write();
        cutpion_plus->Write();
        cutpion_minus->Write();
        rootfile->Close();

        // delete Spirit_tree;
        // delete rootfile;
}
