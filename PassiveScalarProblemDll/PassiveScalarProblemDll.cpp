// PassiveScalarProblemDll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "PassiveScalarProblemDll.h"
#include <FractionalStep_x.h>
#include <PreInformation.h>

namespace SIM {

	typedef FractionalStep_x<Parameters::DataType, Parameters::Dimension, Parameters::Order> FS;
	typedef FS* FSPtr;

	static FSPtr objPtr;

	void PassiveScalarProblemDll::Initialize() {
		objPtr = new FS();
		objPtr->init();
	}

	void PassiveScalarProblemDll::Run() {
		objPtr->stepGL();
	}

	int PassiveScalarProblemDll::Number() {
		return objPtr->part->np;
	}
	NPtr PassiveScalarProblemDll::Type() {
		return NPtr(objPtr->type());
	}
	NPtr PassiveScalarProblemDll::Position() {
		return NPtr(objPtr->position());
	}
	NPtr PassiveScalarProblemDll::Scalar() {
		return NPtr(objPtr->scalar());
	}

	void PassiveScalarProblemDll::SaveData() {
		objPtr->saveData();
	}
	void PassiveScalarProblemDll::SensorOut() {
		objPtr->sensorOut();
	}

}