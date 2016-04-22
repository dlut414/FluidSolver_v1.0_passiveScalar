/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//PassiveScalarProblemDll.h
///defination of class PassiveScalarProblemDll
#pragma once
#ifdef PASSIVESCALARPROBLEMDLL_EXPORTS
#define PASSIVESCALARPROBLEMDLL_API __declspec(dllexport)
#else
#define PASSIVESCALARPROBLEMDLL_API __declspec(dllimport)
#endif

namespace SIM {

	typedef void* NPtr;

	class PassiveScalarProblemDll {
	public:
		static PASSIVESCALARPROBLEMDLL_API void Initialize();
		static PASSIVESCALARPROBLEMDLL_API void Run();
		static PASSIVESCALARPROBLEMDLL_API int Number();
		static PASSIVESCALARPROBLEMDLL_API NPtr Type();
		static PASSIVESCALARPROBLEMDLL_API NPtr Position();
		static PASSIVESCALARPROBLEMDLL_API NPtr Scalar();
		static PASSIVESCALARPROBLEMDLL_API void SaveData();
		static PASSIVESCALARPROBLEMDLL_API void SensorOut();
	};

}