
#ifndef _CSA_MYGMDATATREE_H___
#define _CSA_MYGMDATATREE_H___

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <string>

class MyGmDataTree
{
public:
					MyGmDataTree();
	virtual			~MyGmDataTree();
			

// METHODS
	void 			AddRootFile( const std::string& in_filename );
	void			InitTree();
	Int_t			GetEntries() const { return (Int_t) mychain->GetEntries(); }

	void			LoadTree( int in_entry );

    int				GetNumberOfTracks() const  { return m_track_EventID->size(); }
	int				GetNumberOfSteps() const  { return m_step_EventID->size(); }

	int         	GetEventID() const { return m_event_EventID; }

// DATA------------
	// EVENT DATA
	int				m_event_EventID;
	double 			m_event_InitialPosX; 
	double			m_event_InitialPosY; 
	double			m_event_InitialPosZ; 
	double 			m_event_InitialKineticEnergy;
	double 			m_event_InitialTotalEnergy;
	double			m_event_InitialTime; 

	// TRACK DATA
    std::vector<int>*   m_track_EventID;
	std::vector<int>* 	m_track_TrackID;
	std::vector<int>* 	m_track_ParentTrackID;
	std::vector<int>* 	m_track_StepNumber;

	std::vector<double>* m_track_InitialPosX;
	std::vector<double>* m_track_InitialPosY;
	std::vector<double>* m_track_InitialPosZ;
	std::vector<double>* m_track_FinalPosX;
	std::vector<double>* m_track_FinalPosY;
	std::vector<double>* m_track_FinalPosZ;
	std::vector<double>* m_track_InitialKineticEnergy;
	std::vector<double>* m_track_FinalKineticEnergy;
	std::vector<double>* m_track_InitialTotalEnergy;
	std::vector<double>* m_track_FinalTotalEnergy;
	std::vector<double>* m_track_KineticEnergyChange;

	std::vector<double>* m_track_InitialMomX;
	std::vector<double>* m_track_InitialMomY;
	std::vector<double>* m_track_InitialMomZ;
	std::vector<double>* m_track_InitialDirX;
	std::vector<double>* m_track_InitialDirY;
	std::vector<double>* m_track_InitialDirZ;

	std::vector<double>* m_track_FinalMomX;
	std::vector<double>* m_track_FinalMomY;
	std::vector<double>* m_track_FinalMomZ;
	std::vector<double>* m_track_FinalDirX;
	std::vector<double>* m_track_FinalDirY;
	std::vector<double>* m_track_FinalDirZ;

	std::vector<std::string>* m_track_InitialLogicalVolume;
	std::vector<std::string>* m_track_FinalLogicalVolume;
	std::vector<std::string>* m_track_InitialPhysicalVolume;
	std::vector<std::string>* m_track_FinalPhysicalVolume;

	std::vector<int>* 	m_track_InitialPVCopyNumber;
	std::vector<int>*	m_track_FinalPVCopyNumber;

    std::vector<std::string>* m_track_CreatorProcess;
	std::vector<std::string>* m_track_ParticleType;
	std::vector<std::string>* m_track_ParticleSubType;
	std::vector<double>*      m_track_ParticleCharge;

	// STEP DATA
	std::vector<int>* 	m_step_EventID;
	std::vector<int>* 	m_step_TrackID;

	std::vector<int>* 	m_step_ParentTrackID;
	std::vector<int>* 	m_step_StepNumber;

	std::vector<double>* m_step_InitialPosX;
	std::vector<double>* m_step_InitialPosY;
	std::vector<double>* m_step_InitialPosZ;

	std::vector<double>* m_step_FinalPosX;
	std::vector<double>* m_step_FinalPosY;
	std::vector<double>* m_step_FinalPosZ;

    std::vector<double>* m_step_InitialTouchablePosX;
    std::vector<double>* m_step_InitialTouchablePosY;
    std::vector<double>* m_step_InitialTouchablePosZ;

    std::vector<double>* m_step_FinalTouchablePosX;
    std::vector<double>* m_step_FinalTouchablePosY;
    std::vector<double>* m_step_FinalTouchablePosZ;

	std::vector<double>* m_step_InitialKineticEnergy;
	std::vector<double>* m_step_FinalKineticEnergy;
	std::vector<double>* m_step_InitialTotalEnergy;
	std::vector<double>* m_step_FinalTotalEnergy;
	std::vector<double>* m_step_KineticEnergyChange;

	std::vector<double>* m_step_AccumulatedEnergyDeposited;

	std::vector<std::string>* m_step_InitialLogicalVolume;
	std::vector<std::string>* m_step_FinalLogicalVolume;

	std::vector<std::string>* m_step_InitialPhysicalVolume;
	std::vector<std::string>* m_step_FinalPhysicalVolume;

	std::vector<std::string>* m_step_InitialTouchable;
	std::vector<int>* 	      m_step_InitialPVCopyNumber;
	std::vector<int>* 	      m_step_FinalPVCopyNumber;

	std::vector<std::string>* m_step_ParticleType;
	std::vector<std::string>* m_step_ParticleSubType;
	std::vector<double>*      m_step_ParticleCharge;

	std::vector<std::string>* m_step_CreatorProcess;
	std::vector<std::string>* m_step_InitialProcess;
	std::vector<std::string>* m_step_FinalProcess;

	std::vector<double>* m_step_InitialLocalTime;
	std::vector<double>* m_step_FinalLocalTime;

private:

	// private methods
	void			Initialize();
	void			InitializeBloodyPointers();

	// private data
	TChain* 		mychain;

};

#endif 
