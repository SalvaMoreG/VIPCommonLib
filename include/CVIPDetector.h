#pragma once

#ifndef __CDETECTOR_H__
#define __CDETECTOR_H__

#include "CVIP3Vector.h"

enum DETECTOR_TYPE
{
	DETTYPE_NONE = 0,
	DETTYPE_SCAT,
	DETTYPE_ABS
};


class CDetector
{

public:
							CDetector();
	virtual					~CDetector();

	void					SetType(DETECTOR_TYPE in_dettype);

	// copy and assign constructors
							CDetector(const CDetector& in_obj);
	CDetector& 				operator= (const CDetector& in_obj);

	//
	void					SetSize(const C3Vector& in_size) { m_size = in_size; }
	const C3Vector&			GetSize() const { return m_size; }

	void					SetPosition( const C3Vector& in_position ) { m_position = in_position; }
	const C3Vector&			GetPosition() const { return m_position; }

	bool					IsInsideDetector( const C3Vector& in_position ) const;

	C3Vector				CreateRandomRealHitPositionInDetector() const;
	C3Vector				GetVoxelizedHitPosition( const C3Vector& in_realHitPosition ) const;

	void					SetEThreshold( const double& in_eThreshold ) { m_eThreshold = in_eThreshold; }
	double					GetEThreshold() const { return m_eThreshold; }

	void					SetVoxelSize( const C3Vector& in_voxelSize) { m_detVoxelSize = in_voxelSize; }
	const C3Vector&			GetVoxelSize() const { return m_detVoxelSize; }

private:
	// private methods

	// private data
	DETECTOR_TYPE			m_dettype;

	C3Vector 				m_size;
	C3Vector 				m_position;
	
	double 					m_eThreshold;

	C3Vector				m_detVoxelSize;	// Spatial resolution
};

#endif

