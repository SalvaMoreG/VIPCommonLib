
#include "MyGmDataTree.h"
#include <iostream>

using namespace std;


MyGmDataTree::MyGmDataTree()
	: mychain(0)
{
	mychain = new TChain("GmDataTTree");

	Initialize();
}

MyGmDataTree::~MyGmDataTree()
{
	delete mychain;
}

void
MyGmDataTree::AddRootFile( const std::string& in_filename )
{
	cout << "adding ROOT file: " << in_filename << endl;
	if (!mychain) 
		cout << "mychain does not exist" << endl;
	else
		mychain->Add( in_filename.c_str() );
}

void
MyGmDataTree::InitTree()
{
//	mychain = tree;
//	mychain->SetMakeClass(1);

// #define EXTRASTUFF

TObjArray* mylist = mychain->GetListOfBranches();

	// EVENT
	mychain->SetBranchAddress("Event_EventID", 				&m_event_EventID);

    if(mylist->FindObject("Event_InitialPosX"))
    {
        mychain->SetBranchAddress("Event_InitialPosX", 			&m_event_InitialPosX);
        mychain->SetBranchAddress("Event_InitialPosY", 			&m_event_InitialPosY);
        mychain->SetBranchAddress("Event_InitialPosZ", 			&m_event_InitialPosZ);
    }
	mychain->SetBranchAddress("Event_InitialTime",          &m_event_InitialTime);

	if (mylist->FindObject("Event_InitialKineticEnergy"))
		mychain->SetBranchAddress("Event_InitialKineticEnergy", &m_event_InitialKineticEnergy);
	if (mylist->FindObject("Event_InitialTotalEnergy"))
		mychain->SetBranchAddress("Event_InitialTotalEnergy", 	&m_event_InitialTotalEnergy);

	// TRACK
	if (mylist->FindObject("Track_EventID"))
	{
		mychain->SetBranchAddress("Track_EventID", 		&m_track_EventID);
		mychain->SetBranchAddress("Track_TrackID", 		&m_track_TrackID);
    }
    if (mylist->FindObject("Track_ParentTrackID"))
        mychain->SetBranchAddress("Track_ParentTrackID", 	&m_track_ParentTrackID);
    if (mylist->FindObject("Track_StepNumber"))
        mychain->SetBranchAddress("Track_StepNumber", 		&m_track_StepNumber);

    if (mylist->FindObject("Track_InitialPosX"))
    {
		mychain->SetBranchAddress("Track_InitialPosX", 	&m_track_InitialPosX );
		mychain->SetBranchAddress("Track_InitialPosY", 	&m_track_InitialPosY );
		mychain->SetBranchAddress("Track_InitialPosZ", 	&m_track_InitialPosZ );
    }
    
    if (mylist->FindObject("Track_FinalPosX"))
    {    
		mychain->SetBranchAddress("Track_FinalPosX", 	&m_track_FinalPosX );
		mychain->SetBranchAddress("Track_FinalPosY", 	&m_track_FinalPosY );
		mychain->SetBranchAddress("Track_FinalPosZ", 	&m_track_FinalPosZ );
    }
    
    if (mylist->FindObject("Track_InitialKineticEnergy"))
		mychain->SetBranchAddress("Track_InitialKineticEnergy", 	&m_track_InitialKineticEnergy);
    if (mylist->FindObject("Track_FinalKineticEnergy"))
		mychain->SetBranchAddress("Track_FinalKineticEnergy", 		&m_track_FinalKineticEnergy);
    if (mylist->FindObject("Track_InitialTotalEnergy"))
		mychain->SetBranchAddress("Track_InitialTotalEnergy", 		&m_track_InitialTotalEnergy);
    if (mylist->FindObject("Track_FinalTotalEnergy"))
		mychain->SetBranchAddress("Track_FinalTotalEnergy", 		&m_track_FinalTotalEnergy);	
    if (mylist->FindObject("Track_KineticEnergyChange"))
        mychain->SetBranchAddress("Track_KineticEnergyChange", 		&m_track_KineticEnergyChange);
		
    if (mylist->FindObject("Track_InitialMomX"))
    {
        mychain->SetBranchAddress("Track_InitialMomX",  &m_track_InitialMomX );
        mychain->SetBranchAddress("Track_InitialMomY",  &m_track_InitialMomY );
        mychain->SetBranchAddress("Track_InitialMomZ",  &m_track_InitialMomZ );
    }
    if (mylist->FindObject("Track_InitialDirX"))
    {
        mychain->SetBranchAddress("Track_InitialDirX",  &m_track_InitialDirX );
        mychain->SetBranchAddress("Track_InitialDirY",  &m_track_InitialDirY );
        mychain->SetBranchAddress("Track_InitialDirZ",  &m_track_InitialDirZ );
    }
	
    if (mylist->FindObject("Track_FinalMomX"))
    {
        mychain->SetBranchAddress("Track_FinalMomX",  &m_track_FinalMomX );
        mychain->SetBranchAddress("Track_FinalMomY",  &m_track_FinalMomY );
        mychain->SetBranchAddress("Track_FinalMomZ",  &m_track_FinalMomZ );
    }
    if (mylist->FindObject("Track_FinalDirX"))
    {
        mychain->SetBranchAddress("Track_FinalDirX",  &m_track_FinalDirX );
        mychain->SetBranchAddress("Track_FinalDirY",  &m_track_FinalDirY );
        mychain->SetBranchAddress("Track_FinalDirZ",  &m_track_FinalDirZ );
    }
	
    if (mylist->FindObject("Track_InitialLogicalVolume"))
        mychain->SetBranchAddress("Track_InitialLogicalVolume", 	&m_track_InitialLogicalVolume);
	if (mylist->FindObject("Track_FinalLogicalVolume"))
		mychain->SetBranchAddress("Track_FinalLogicalVolume", 		&m_track_FinalLogicalVolume);
    if (mylist->FindObject("Track_InitialPhysicalVolume"))
		mychain->SetBranchAddress("Track_InitialPhysicalVolume", 	&m_track_InitialPhysicalVolume);
    if (mylist->FindObject("Track_FinalPhysicalVolume"))
		mychain->SetBranchAddress("Track_FinalPhysicalVolume", 		&m_track_FinalPhysicalVolume);

    if (mylist->FindObject("Track_InitialPVCopyNumber"))
		mychain->SetBranchAddress("Track_InitialPVCopyNumber", 		&m_track_InitialPVCopyNumber);
    if (mylist->FindObject("Track_FinalPVCopyNumber"))
		mychain->SetBranchAddress("Track_FinalPVCopyNumber",		&m_track_FinalPVCopyNumber);

    if (mylist->FindObject("Track_CreatorProcess"))
        mychain->SetBranchAddress("Track_CreatorProcess", 			&m_track_CreatorProcess);
    if (mylist->FindObject("Track_ParticleType"))
		mychain->SetBranchAddress("Track_ParticleType", 			&m_track_ParticleType);
    if (mylist->FindObject("Track_ParticleSubType"))
		mychain->SetBranchAddress("Track_ParticleSubType", 			&m_track_ParticleSubType);
    if (mylist->FindObject("Track_ParticleCharge"))
		mychain->SetBranchAddress("Track_ParticleCharge", 			&m_track_ParticleCharge);

	// STEP
	if(mylist->FindObject("Step_EventID")) 
		mychain->SetBranchAddress("Step_EventID", 				&m_step_EventID);
	if(mylist->FindObject("Step_TrackID")) 
		mychain->SetBranchAddress("Step_TrackID", 				&m_step_TrackID);

	if(mylist->FindObject("Step_ParentTrackID")) 
		mychain->SetBranchAddress("Step_ParentTrackID", 		&m_step_ParentTrackID);
	if(mylist->FindObject("Step_StepNumber")) 
		mychain->SetBranchAddress("Step_StepNumber", 			&m_step_StepNumber);

	if(mylist->FindObject("Step_InitialPosX")) 
		mychain->SetBranchAddress("Step_InitialPosX", 	&m_step_InitialPosX );
	if(mylist->FindObject("Step_InitialPosY")) 
		mychain->SetBranchAddress("Step_InitialPosY", 	&m_step_InitialPosY );
	if(mylist->FindObject("Step_InitialPosZ")) 
		mychain->SetBranchAddress("Step_InitialPosZ", 	&m_step_InitialPosZ );

	if (mylist->FindObject("Step_FinalPosX"))
	{
		mychain->SetBranchAddress("Step_FinalPosX", 	&m_step_FinalPosX );
		mychain->SetBranchAddress("Step_FinalPosY", 	&m_step_FinalPosY );
		mychain->SetBranchAddress("Step_FinalPosZ", 	&m_step_FinalPosZ );
	}

    if (mylist->FindObject("Step_InitialTouchablePosX"))
    {
        mychain->SetBranchAddress("Step_InitialTouchablePosX",    &m_step_InitialTouchablePosX );
        mychain->SetBranchAddress("Step_InitialTouchablePosY",    &m_step_InitialTouchablePosY );
        mychain->SetBranchAddress("Step_InitialTouchablePosZ",    &m_step_InitialTouchablePosZ );
    }

	if(mylist->FindObject("Step_InitialKineticEnergy")) 
		mychain->SetBranchAddress("Step_InitialKineticEnergy", 	&m_step_InitialKineticEnergy);
	if(mylist->FindObject("Step_FinalKineticEnergy")) 
		mychain->SetBranchAddress("Step_FinalKineticEnergy", 	&m_step_FinalKineticEnergy);

	if(mylist->FindObject("Step_InitialTotalEnergy")) 
		mychain->SetBranchAddress("Step_InitialTotalEnergy", 	&m_step_InitialTotalEnergy);
	if(mylist->FindObject("Step_FinalTotalEnergy")) 
		mychain->SetBranchAddress("Step_FinalTotalEnergy", 		&m_step_FinalTotalEnergy);
	if(mylist->FindObject("Step_KineticEnergyChange")) 
		mychain->SetBranchAddress("Step_KineticEnergyChange", 	&m_step_KineticEnergyChange);

	if(mylist->FindObject("Step_AccumulatedEnergyDeposited")) 
		mychain->SetBranchAddress("Step_AccumulatedEnergyDeposited", &m_step_AccumulatedEnergyDeposited);

	if(mylist->FindObject("Step_InitialLogicalVolume")) 
		mychain->SetBranchAddress("Step_InitialLogicalVolume", 	&m_step_InitialLogicalVolume);
	if(mylist->FindObject("Step_FinalLogicalVolume")) 
		mychain->SetBranchAddress("Step_FinalLogicalVolume", 	&m_step_FinalLogicalVolume);

	if(mylist->FindObject("Step_InitialPhysicalVolume")) 
		mychain->SetBranchAddress("Step_InitialPhysicalVolume", &m_step_InitialPhysicalVolume);
	if(mylist->FindObject("Step_FinalPhysicalVolume")) 
		mychain->SetBranchAddress("Step_FinalPhysicalVolume", 	&m_step_FinalPhysicalVolume);

	if(mylist->FindObject("Step_InitialPVCopyNumber")) 
	{
		mychain->SetBranchAddress("Step_InitialPVCopyNumber", 	&m_step_InitialPVCopyNumber);
		mychain->SetBranchAddress("Step_FinalPVCopyNumber", 	&m_step_FinalPVCopyNumber);
	}

	if(mylist->FindObject("Step_ParticleType")) 
		mychain->SetBranchAddress("Step_ParticleType",			&m_step_ParticleType);
    if(mylist->FindObject("Step_ParticleSubType")) 
		mychain->SetBranchAddress("Step_ParticleSubType",		&m_step_ParticleSubType);
    if(mylist->FindObject("Step_ParticleCharge")) 
		mychain->SetBranchAddress("Step_ParticleCharge",		&m_step_ParticleCharge);

	if(mylist->FindObject("Step_CreatorProcess")) 
		mychain->SetBranchAddress("Step_CreatorProcess", 		&m_step_CreatorProcess);
    if(mylist->FindObject("Step_InitialProcess")) 
		mychain->SetBranchAddress("Step_InitialProcess", 		&m_step_InitialProcess);
    if(mylist->FindObject("Step_FinalProcess")) 
		mychain->SetBranchAddress("Step_FinalProcess", 			&m_step_FinalProcess);

}

void
MyGmDataTree::LoadTree( int in_entry )
{
	mychain->GetEntry(in_entry);
}

void
MyGmDataTree::Initialize()
{
	// Initialize bloody pointers....
	// VERY IMPORTANT!!!!!!!!!!!!!!!!!!!! POINTERS ARE EXPECTED TO BE NULL AT INIT!!!!!!!!!!!!!
  	m_event_EventID = -123456;
  	m_event_InitialTime = 0;
	m_event_InitialTotalEnergy = 0;


	m_track_EventID = 0;
	m_track_TrackID = 0;

	m_track_ParentTrackID = 0;
	m_track_StepNumber = 0;

	m_track_InitialPosX = 0;
	m_track_InitialPosY = 0;
	m_track_InitialPosZ = 0;
	m_track_FinalPosX = 0;
	m_track_FinalPosY = 0;
	m_track_FinalPosZ = 0;
	m_track_InitialKineticEnergy = 0;
	m_track_FinalKineticEnergy = 0;
	m_track_InitialTotalEnergy = 0;
	m_track_FinalTotalEnergy = 0;

	m_track_KineticEnergyChange = 0;
	
	m_track_InitialMomX = 0;
	m_track_InitialMomY = 0;
	m_track_InitialMomZ = 0;
	m_track_InitialDirX = 0;
	m_track_InitialDirY = 0;
	m_track_InitialDirZ = 0;

	m_track_FinalMomX = 0;
	m_track_FinalMomY = 0;
	m_track_FinalMomZ = 0;
	m_track_FinalDirX = 0;
	m_track_FinalDirY = 0;
	m_track_FinalDirZ = 0;

	m_track_InitialLogicalVolume = 0;

	m_track_FinalLogicalVolume = 0;
	m_track_InitialPhysicalVolume = 0;
	m_track_FinalPhysicalVolume = 0;

	m_track_InitialPVCopyNumber = 0;
	m_track_FinalPVCopyNumber = 0;

	m_track_CreatorProcess = 0;
	m_track_ParticleType = 0;
	m_track_ParticleSubType = 0;
	m_track_ParticleCharge = 0;

	m_step_EventID = 0;
	m_step_TrackID = 0;

	m_step_ParentTrackID = 0;
	m_step_StepNumber = 0;

	m_step_InitialPosX = 0;
	m_step_InitialPosY = 0;
	m_step_InitialPosZ = 0;

	m_step_FinalPosX  = 0;
	m_step_FinalPosY  = 0;
	m_step_FinalPosZ  = 0;

    m_step_InitialTouchablePosX = 0;
    m_step_InitialTouchablePosY = 0;
    m_step_InitialTouchablePosZ = 0;

	m_step_InitialKineticEnergy = 0;
	m_step_FinalKineticEnergy = 0;

	m_step_InitialTotalEnergy = 0;
	m_step_FinalTotalEnergy = 0;
	m_step_KineticEnergyChange = 0;

	m_step_AccumulatedEnergyDeposited = 0;

	m_step_InitialLogicalVolume = 0;
	m_step_FinalLogicalVolume = 0;

	m_step_InitialPhysicalVolume = 0;
	m_step_FinalPhysicalVolume = 0;

	m_step_InitialTouchable = 0;
	m_step_InitialPVCopyNumber = 0;
	m_step_FinalPVCopyNumber = 0;

	m_step_ParticleType = 0; 
	m_step_ParticleSubType = 0;
	m_step_ParticleCharge = 0;

	m_step_CreatorProcess = 0;
	m_step_InitialProcess = 0;
	m_step_FinalProcess = 0;

	m_step_InitialLocalTime = 0;
	m_step_FinalLocalTime = 0;
}



