#pragma once
//=============================================================================
// This plugin example illustrates how to create a new material. 
// It requires FEBio 3.0 (or up)
//
// Author : Steve Maas
// Copyright (c) 2015 - 2020
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial, which is defined in this include file.
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
// This material class implements a neo-Hookean constitutive model. 
// Since it is a (compressible, coupled) hyper-elatic material, it needs to inherit
// from FEElasticMaterial. 
class DemirayUncoupled : public FEUncoupledMaterial
{
public:
	// The constructor is called when an instance of this class is created.
	// All classes registered by the framework must take the FEModel* as the only
	// parameter in the constructor, even if the class does not need it (which most often
	// will be the case). For material classes, the FEModel parameter is passed to the 
	// base class in the initialization list.
	DemirayUncoupled(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

private:
	// The neo-Hookean material defines two material parameters.
	// They are defined here, but the class also needs to let the framework
	// know that these parameters exist. This is necessary to define the parameters
	// in the FEBio input file. This is done by declaring the DECLARE_FECORE_CLASS() below.
	
	FEParamDouble  m_a;
	FEParamDouble  m_b;

public:
	// The following three functions are the only functions that have to be defined.
	// (Actually, the Init function is optional). They are defined as virtual member 
	// functions in the base classes and will be called by FEBio to invoke the specific
	// class implementations.

	// This function calculates the spatial (i.e. Cauchy or true) stress.
	// It takes one parameter, the FEMaterialPoint and returns a mat3ds object
	// which is a symmetric second-order tensor.
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	// This function calculates the spatial elasticity tangent tensor. 
	// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
	// which is a fourth-order tensor with major and minor symmetries.
	virtual tens4ds DevTangent(FEMaterialPoint& pt);
	virtual double DevStrainEnergyDensity(FEMaterialPoint& pt);
	//void logMat3ds(mat3ds mat, string str);


	// This macro defines that the class will define a material parameter list.
	// The material parameter list itself is defined elsewhere (e.g. in the .cpp file.)
	DECLARE_FECORE_CLASS();
};
