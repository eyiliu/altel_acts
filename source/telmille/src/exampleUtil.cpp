/*
 * exampleUtil.cpp
 *
 *  Created on: 6 Nov 2018
 *      Author: kleinwrt
 */

/** \file
 *  Example utilities.
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

#include "exampleUtil.h"

using namespace Eigen;


/// Create a detector layer.
/**
 * Create detector layer with 1D or 2D measurement (u,v)
 * \param [in] aName          name
 * \param [in] aLayer         layer ID
 * \param [in] aDim           dimension (1,2)
 * \param [in] thickness      thickness / radiation_length
 * \param [in] aCenter        center of detector (origin of local systems)
 * \param [in] aResolution    resolution vector
 * \param [in] aPrecision     diagonal of precision matrix
 * \param [in] measTrafo      matrix of row vectors defining local measurement system
 * \param [in] alignTrafo     matrix of row vectors defining local alignment system
 */
GblDetectorLayer::GblDetectorLayer(const std::string aName,
		const unsigned int aLayer, const int aDim, const double thickness,
		Eigen::Vector3d& aCenter, Eigen::Vector2d& aResolution,
		Eigen::Vector2d& aPrecision, Eigen::Matrix3d& measTrafo,
		Eigen::Matrix3d& alignTrafo) :
		name(aName), layer(aLayer), measDim(aDim), xbyx0(thickness),
                center(aCenter), resolution(aResolution), precision(aPrecision),
                global2meas(measTrafo), global2align(alignTrafo) { //TODO: better to have meas2align, instead of global2meas
	udir = global2meas.row(0);
	vdir = global2meas.row(1);
	ndir = global2meas.row(2);
}

GblDetectorLayer::~GblDetectorLayer() {
}

/// Print GblDetectorLayer.
void GblDetectorLayer::print() const {
	IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::cout << " Layer " << name << " " << layer << " : " << measDim << "D, "
			<< xbyx0 << " X0, @ " << center.transpose().format(CleanFmt)
			<< ", res " << resolution.transpose().format(CleanFmt) << ", udir "
			<< udir.transpose().format(CleanFmt) << ", vdir "
			<< vdir.transpose().format(CleanFmt) << std::endl;
}

/// Get layer ID
unsigned int GblDetectorLayer::getLayerID() const {
	return layer;
}

/// Get radiation length.
double GblDetectorLayer::getRadiationLength() const {
	return xbyx0;
}

/// Get resolution.
Eigen::Vector2d GblDetectorLayer::getResolution() const {
	return resolution;
}

/// Get precision.
Eigen::Vector2d GblDetectorLayer::getPrecision() const {
	return precision;
}

/// Get center.
Eigen::Vector3d GblDetectorLayer::getCenter() const {
	return center;
}

/// Get directions of measurement system.
/**
 * Matrix from row vectors (transformation from global to measurement system)
 */
Eigen::Matrix3d GblDetectorLayer::getMeasSystemDirs() const {
	return global2meas;;
}

/// Get directions of alignment system.
/**
 * Matrix from row vectors (transformation from global to alignment system)
 */
Eigen::Matrix3d GblDetectorLayer::getAlignSystemDirs() const {
	return global2align;;
}


/// Get rigid body derivatives in global frame.
/**
 * \param[in] position   position (of prediction or measurement)
 * \param[in] direction  track direction
 */
Matrix<double, 3, 6> GblDetectorLayer::getRigidBodyDerGlobal(
		Eigen::Vector3d& position, Eigen::Vector3d& direction) const {
// lever arms (for rotations)
	Vector3d dist = position;
// dr/dg (residual vs measurement, 1-tdir*ndir^t/tdir*ndir)
	Matrix3d drdm = Matrix3d::Identity()
			- (direction * ndir.transpose()) / (direction.transpose() * ndir);
// dg/db (measurement vs 6 rigid body parameters, global system)
	Matrix<double, 3, 6> dmdg = Matrix<double, 3, 6>::Zero();
        dmdg<<
          1., 0., 0.,  0.,     -dist(2), dist(1),
          0., 1., 0.,  dist(2), 0.,     -dist(0),
          0., 0., 1., -dist(1), dist(0), 0.;

        // CONG
        // dmdg<<
        //   -1., 0., 0.,  0.,     -dist(2), dist(1),
        //   0., -1., 0.,  dist(2), 0.,     -dist(0),
        //   0., 0., -1., -dist(1), dist(0), 0.;

// drl/dg (local residuals vs rigid body parameters)
	return global2meas * drdm * dmdg; //  drlm/dg in world system   // dr/db in meas system? (global2meas does not consider center offset)
 }


/// Get rigid body derivatives in local (alignment) frame (rotated in measurement plane).
/**
 * The orthogonal alignment frame differs from measurement frame only by rotations
 * around normal to measurement plane.
 *
 * Equivalent to:
 * \code
 * getRigidBodyDerGlobal(position, direction) * getTrafoLocalToGlobal(center, global2align)
 * \endcode
 *
 * \param[in] position   position (of prediction or measurement)
 * \param[in] direction  track direction
 */
Matrix<double, 2, 6> GblDetectorLayer::getRigidBodyDerLocal(
		Eigen::Vector3d& position, Eigen::Vector3d& direction) const {
	// track direction in local system
	Vector3d tLoc = global2align * direction;
	// local slopes
	const double uSlope = tLoc[0] / tLoc[2];
	const double vSlope = tLoc[1] / tLoc[2];
	// lever arms (for rotations)
	Vector3d dist = global2align * (position - center);
	const double uPos = dist[0];
	const double vPos = dist[1];
	// wPos = 0 (in detector plane)
	// drl/dg (local residuals (du,dv) vs rigid body parameters)
	Matrix<double, 2, 6> drdb; // drla/dg   // dr/db in align system
	drdb <<
          1.0, 0.0, -uSlope, vPos * uSlope, -uPos * uSlope, vPos,
          0.0, 1.0, -vSlope, vPos * vSlope, -uPos * vSlope, -uPos;

// local (alignment) to measurement system
	Matrix3d local2meas = global2meas * global2align.transpose();
        //Note: block<2,2>, assuming in-plane rotation only
	return local2meas.block<2, 2>(0, 0) * drdb; // drlm/dg in DetectorLayer-local-align system // dr/db in align system
}

Matrix<double, 2, 6> GblDetectorLayer::getRigidBodyDerLocal_mod(
  Eigen::Vector3d& position, Eigen::Vector3d& direction) const {
	// track direction in local system
	Vector3d tDir = global2align * direction;
	// local slopes
	const double uSlope = tDir[0] / tDir[2];
	const double vSlope = tDir[1] / tDir[2];
	// lever arms (for rotations)
	Vector3d dist = global2align * (position - center);
	const double uPos = dist[0];
	const double vPos = dist[1];
	// wPos = 0 (in detector plane)
	// drl/dg (local residuals (du,dv) vs rigid body parameters)
	Matrix<double, 2, 6> drdb; // drla/dg
	drdb <<
          -1.0, 0.0, uSlope, vPos * uSlope, -uPos * uSlope, vPos,
          0.0, -1.0, vSlope, vPos * vSlope, -uPos * vSlope, -uPos;

// local (alignment) to measurement system
	Matrix3d local2meas = global2meas * global2align.transpose();
        //Note: block<2,2>, assuming meas in-plane rotation on align system
	return local2meas.block<2, 2>(0, 0) * drdb; // drlm/dg in DetectorLayer-local-align system//
}


//Assuming rotation around original local x around alpha, --> orignal local y around beta --> original local z around gamma
//The new rotation matrix is R * R_z(gamma) * R_y(beta) * R_x(gamma) (R * rotation(x) represents rotation around the new axis by x)
/*
     // The delta translation
    Eigen::Vector3d deltaCenter =
        deltaAlignmentParam.segment<3>(0);
    // The delta Euler angles
    Eigen::Vector3d deltaEulerAngles =
        deltaAlignmentParam.segment<3>(3);

    const Acts::Vector3 newCenter = oldCenter + deltaCenter;
    Eigen::Transform3d newTransform = oldTransform;
    newTransform.translation() = newCenter;
   
    newTransform *=
        Eigen::AngleAxis3(deltaEulerAngles(2), Acts::Vector3::UnitZ());
    newTransform *=
        Eigen::AngleAxis3(deltaEulerAngles(1), Acts::Vector3::UnitY());
    newTransform *=
        Eigen::AngleAxis3(deltaEulerAngles(0), Acts::Vector3::UnitX());
*/


Matrix<double, 2, 6> GblDetectorLayer::getRigidBodyDerLocal_ai(
  Eigen::Vector3d& position, Eigen::Vector3d& direction) const {

  //Get the derivatives of local x axis, local y axis , local z axis to rotation parameters
  Matrix<double, 3, 3> matrix1;
  matrix1.col(0) = Eigen::Vector3d(0, 0, 0);
  matrix1.col(1) = Eigen::Vector3d(0, 0, -1);
  matrix1.col(2) = Eigen::Vector3d(0, 1, 0);

  Matrix<double, 3, 3> matrix2;
  matrix2.col(0) = Eigen::Vector3d(0, 0, 1);
  matrix2.col(1) = Eigen::Vector3d(0, 0, 0);
  matrix2.col(2) = Eigen::Vector3d(-1, 0, 0);
  
  Matrix<double, 3, 3> matrix3;
  matrix3.col(0) = Eigen::Vector3d(0, -1, 0);
  matrix3.col(1) = Eigen::Vector3d(1, 0, 0);
  matrix3.col(2) = Eigen::Vector3d(0, 0, 0);

  Matrix<double, 3, 3> rotationMatrix;
  rotationMatrix.col(0) = udir;
  rotationMatrix.col(1) = vdir;
  rotationMatrix.col(2) = ndir;

  Matrix<double, 3, 3> localXAxisToRotation = rotationMatrix*matrix1; 
  Matrix<double, 3, 3> localYAxisToRotation = rotationMatrix*matrix2; 
  Matrix<double, 3, 3> localZAxisToRotation = rotationMatrix*matrix3; 

  /*
  If the orignal rotation matrix is [cosu cosv 0]
                                    [sinu sinv 0]
				    [0    0    1]
				    then,
       localXAxisToRotation = [ 0  0  cosv ]
                              [ 0  0  sinv ]
			      [ 0 -1   0   ] 

       localYAxisToRotation = [ 0  0 -cosU ]
                              [ 0  0 -sinU ]
			      [ 1  0   0   ]
  
       localZAxisToRotation = [-cosv cosU 0]
                              [-sinv sinU 0]
			      [ 0  0   0   ]
  */

  Eigen::Vector3d dist = position - center;
  double directionInX = udir.dot(direction);
  double directionInY = vdir.dot(direction);
  double directionInZ = ndir.dot(direction);

  Matrix<double, 1, 6> localXToAlignmentParameters = Matrix<double, 1, 6>::Zero();
  Matrix<double, 1, 6> localYToAlignmentParameters = Matrix<double, 1, 6>::Zero();
  
  localXToAlignmentParameters.block<1,3>(0,0) = -1.0 * udir.transpose() + directionInX/directionInZ*ndir.transpose();
  localXToAlignmentParameters.block<1,3>(0,3) = dist.transpose()*localXAxisToRotation - directionInX/directionInZ*dist.transpose()*localZAxisToRotation;  
  
  localYToAlignmentParameters.block<1,3>(0,0) = -1.0 * vdir.transpose() + directionInY/directionInZ*ndir.transpose();
  localYToAlignmentParameters.block<1,3>(0,3) = dist.transpose()*localYAxisToRotation - directionInY/directionInZ*dist.transpose()*localZAxisToRotation;  
 
  Matrix<double, 2, 6> localToAlignmentParameters;
  localToAlignmentParameters.row(0) = localXToAlignmentParameters;
  localToAlignmentParameters.row(1) = localYToAlignmentParameters;

  /*
  v = u + 90 
  cosv = cos(u + π/2) = -sinu
  sinv = sin(u + π/2) = cosu
   */

  /*
   -cosu, -sinu, directionInX/directionInZ, 
   directionInX/directionInZ*(dist(0)*cosv + dist(1)*sinv),  
   -dist(2) -directionInX/directionInZ* (dist(0)*cosu + dist(1)sinu),  
   dist(0)cosv + dist(1)sinv   // dist(1)cosu - dist(0)sinu
  
   -cosv, -sinv, directionInY/directionInZ, 
   dist(2) + directionInY/directionInZ*(dist(0)*cosv + dist(1)*sinv), 
   -directionInY/directionInZ* (dist(0)*cosu + dist(1)sinu),  
   -(dist(0)cosu + dist(1)sinu)  // dist(1)cosv - dist(0)sinv  
  */

  return localToAlignmentParameters;
}





/// Get transformation for rigid body derivatives from global to local (alignment) system.
/**
 * local = rotation * (global-offset)
 *
 * \param[in] offset    offset of alignment system
 * \param[in] rotation  rotation of alignment system
 */
Matrix<double, 6, 6> GblDetectorLayer::getTrafoGlobalToLocal(
		Eigen::Vector3d& offset, Eigen::Matrix3d& rotation) const {
	// transformation global to local
	Matrix<double, 6, 6> glo2loc = Matrix<double, 6, 6>::Zero();
	Matrix3d leverArms;
	leverArms<<
          0.,         offset[2], -offset[1],
          -offset[2], 0.,         offset[0],
          offset[1], -offset[0],  0.;

	glo2loc.block<3, 3>(0, 0) = rotation;
	glo2loc.block<3, 3>(0, 3) = -rotation * leverArms;
	glo2loc.block<3, 3>(3, 3) = rotation;
	return glo2loc;
}

/// Get transformation for rigid body derivatives from local (alignment) to global system.
/**
 * local = rotation * (global-offset)
 *
 * \param[in] offset    offset of alignment system
 * \param[in] rotation  rotation of alignment system
 */
Matrix<double, 6, 6> GblDetectorLayer::getTrafoLocalToGlobal(
		Eigen::Vector3d& offset, Eigen::Matrix3d& rotation) const {
	// transformation local to global
	Matrix<double, 6, 6> loc2glo = Matrix<double, 6, 6>::Zero();
	Matrix3d leverArms;
	leverArms <<
          0.,         offset[2], -offset[1],
          -offset[2], 0.,         offset[0],
          offset[1], -offset[0],  0.;

	loc2glo.block<3, 3>(0, 0) = rotation.transpose();
	loc2glo.block<3, 3>(0, 3) = leverArms * rotation.transpose();
	loc2glo.block<3, 3>(3, 3) = rotation.transpose();
	return loc2glo;
}

