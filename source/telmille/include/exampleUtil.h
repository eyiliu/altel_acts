/*
 * exampleUtil.h
 *
 *  Created on: 6 Nov 2018
 *      Author: kleinwrt
 */

/** \file
 *  Definitions for exampleUtil(ities).
 *
 *  \author Claus Kleinwort, DESY, 2018 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2018 Deutsches Elektronen-Synchroton,
 *  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 *  This library is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version. \n\n
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details. \n\n
 *  You should have received a copy of the GNU Library General Public
 *  License along with this program (see the file COPYING.LIB for more
 *  details); if not, write to the Free Software Foundation, Inc.,
 *  675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef SRC_EXAMPLEUTIL_H_
#define SRC_EXAMPLEUTIL_H_

#include "Eigen/Dense"
#include<iostream>

/// Detector layer
/**
 *
 */
class GblDetectorLayer {
public:
	GblDetectorLayer(const std::string aName, const unsigned int aLayer,
			const int aDim, const double thickness, Eigen::Vector3d& aCenter,
			Eigen::Vector2d& aResolution, Eigen::Vector2d& aPrecision,
			Eigen::Matrix3d& measTrafo, Eigen::Matrix3d& alignTrafo);
	virtual ~GblDetectorLayer();
	void print() const;
	unsigned int getLayerID() const;
	double getRadiationLength() const;
	Eigen::Vector2d getResolution() const;
	Eigen::Vector2d getPrecision() const;
	Eigen::Vector3d getCenter() const;
	Eigen::Matrix3d getMeasSystemDirs() const;
	Eigen::Matrix3d getAlignSystemDirs() const;
	Eigen::Matrix<double, 3, 6> getRigidBodyDerGlobal(Eigen::Vector3d& position,
			Eigen::Vector3d& direction) const;
	Eigen::Matrix<double, 2, 6> getRigidBodyDerLocal(Eigen::Vector3d& position,
			Eigen::Vector3d& direction) const;
	Eigen::Matrix<double, 2, 6> getRigidBodyDerLocal_mod(Eigen::Vector3d& position,
			Eigen::Vector3d& direction) const;
	Eigen::Matrix<double, 2, 6> getRigidBodyDerLocal_ai(Eigen::Vector3d& position,
			Eigen::Vector3d& direction) const;
	Eigen::Matrix<double, 6, 6> getTrafoGlobalToLocal(Eigen::Vector3d& offset,
			Eigen::Matrix3d& rotation) const;
	Eigen::Matrix<double, 6, 6> getTrafoLocalToGlobal(Eigen::Vector3d& offset,
			Eigen::Matrix3d& rotation) const;

private:
	std::string name; ///< name
	unsigned int layer; ///< layer ID
	unsigned int measDim; ///< measurement dimension (1 or 2)
	double xbyx0; ///< normalized material thickness
	Eigen::Vector3d center; ///< center
	Eigen::Vector2d resolution; ///< measurements resolution
	Eigen::Vector2d precision; ///< measurements precision
	Eigen::Vector3d udir; ///< 1. measurement direction
	Eigen::Vector3d vdir; ///< 2. measurement direction
	Eigen::Vector3d ndir; ///< normal to measurement plane
	Eigen::Matrix3d global2meas; ///< transformation into measurement system
	Eigen::Matrix3d global2align; ///< transformation into (local) alignment system
};
#endif /* SRC_EXAMPLEUTIL_H_ */
