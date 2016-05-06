//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm5/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
// $Id: TrackingAction.cc 76464 2013-11-11 10:22:56Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
//#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* DET, EventAction* EA)
:G4UserTrackingAction(),fDetector(DET), fEventAction(EA)
{ 

runLast = -1;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  // few initialisations
  //
  if (aTrack->GetTrackID() == 1) {
    fXstartAbs = fDetector->GetxstartAbs();
    fXendAbs   = fDetector->GetxendAbs();
    fPrimaryCharge = aTrack->GetDefinition()->GetPDGCharge();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4AnalysisManager* analysisManager 
                          = G4AnalysisManager::Instance();
  G4RunManager* run       = G4RunManager::GetRunManager();
  G4int currentRun        = run->GetCurrentEvent()->GetEventID();
  G4ThreeVector position = aTrack->GetPosition();
  G4ThreeVector vertex   = aTrack->GetVertexPosition();  
  G4double charge        = aTrack->GetDefinition()->GetPDGCharge();

  G4bool transmit = ((position.x() >= fXendAbs) && (vertex.x() < fXendAbs));
  G4bool charged  = (charge != 0.);
  G4bool neutral = !charged;

  G4int id = -1;
  G4ThreeVector direction = aTrack->GetMomentumDirection();
  if (aTrack->GetKineticEnergy() > 0 && transmit)
  {
  G4double theta  = std::acos(direction.x());
  G4double et    = aTrack->GetKineticEnergy();
//  std::cout << "X,Y,Z = " << direction.x() << ", " << direction.y() << ", " << direction.z() << std::endl;
  if (charged) id = 0;
  if (neutral && aTrack->GetKineticEnergy()>1000) id = 1;
   if (id>-1) {
    if (theta > 0.0) {
      G4double et    = aTrack->GetKineticEnergy();
   //   analysisManager->FillH1(id,et);
      analysisManager->FillNtupleDColumn(id,et);
      analysisManager->FillNtupleDColumn(!id,-1);

  
    } 
  }
  id = -1;   
  if (charged) id = 0;
  else if (neutral && aTrack->GetKineticEnergy() > 1000) id = 1;
  if (id>-1) {
    if (direction.x() != 0.0) {
     analysisManager->FillH1(id,theta);
     analysisManager->FillNtupleDColumn(id+2,theta);
     analysisManager->FillNtupleDColumn(!id+2,-1);
     //std::cout << "X,Y,Z = " << direction.x() << ", " << direction.y() << ", " << direction.z() << std::endl;

     //std::cout << "The value of theta is " << theta << std::endl;
    }
  }

  id = -1;
  if (charged) id = 0;
  else if (neutral && aTrack->GetKineticEnergy() > 1000) id = 1;
  if (id>-1) {
    if (direction.x() != 0.0) {
     G4double phi = std::asin(direction.z()/std::sqrt(direction.y()*direction.y()+direction.z()*direction.z()));
     analysisManager->FillNtupleDColumn(id+4,phi);
     analysisManager->FillNtupleDColumn(!id+4,-1);
     analysisManager->FillNtupleIColumn(6,currentRun);
     analysisManager->AddNtupleRow();
     //std::cout << "The value of phi is " << phi << std::endl;

    }
  } 

  }
  runLast = currentRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

