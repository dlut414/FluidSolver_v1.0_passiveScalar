/*
*/
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "Header.h"
#include "Parameter.h"
#include "Particle.h"
#include "MatSolver.h"
#include "Shifter.h"
#include "Sensor.h"
#include "PassiveScalarProblem.h"

namespace SIM {

	template <typename R, unsigned D, typename Derived>
	class Simulator {
		typedef Eigen::Matrix<R,D,1> Vec;
		typedef Eigen::Matrix<R,D,D> Mat;
		typedef Eigen::Triplet<R> Tpl;
	public:
		Simulator() { timeStep = 0; }
		~Simulator() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void operator >> (const std::string& str) const {
			saveData(str);
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " No file Para. found ! " << std::endl;
			file >> para.k >> para.rho >> para.niu >> para.dtMax >> para.cfl >> para.tt >> para.eps >> para.beta >> para.alpha >> para.c;
			std::cout << " Effective radius (times of dp)   : " << para.k << std::endl;
			std::cout << " Density (kg/m3)                  : " << para.rho << std::endl;
			std::cout << " Kinematic viscosity (m2/s)       : " << para.niu << std::endl;
			std::cout << " Maximum time step (s)            : " << para.dtMax << std::endl;
			std::cout << " CFL number                       : " << para.cfl << std::endl;
			std::cout << " Total time (s)                   : " << para.tt << std::endl;
			std::cout << " EPS                              : " << para.eps << std::endl;
			std::cout << " Beta of surface detection        : " << para.beta << std::endl;
			std::cout << " Re-meshing effective radius      : " << para.alpha << std::endl;
			std::cout << " Scaling of re-meshing strength   : " << para.c << std::endl;
			std::cout << " Reading Para. done " << std::endl;
			file.close();
		}

		void init() {
			*this << "Para.txt";
			derived().init_();
			mSol = new MatSolver<R, D>(unsigned(derived().part->np), para.eps);
			std::cout << " Particle number : " << derived().part->np << std::endl;
			sen << "Sensor.in";
			R tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			timeStep = int(derived().part->ct / para.dt);
		}

		void mainLoop() {
			auto* const part = derived().part;
			while (part->ct <= para.tt) {
				std::cout << " step ----------------------------------> " << timeStep << std::endl;
				R tmp = cfl();
				para.dt = tmp < para.dtMax ? tmp : para.dtMax;
				part->updateCell();
				derived().step();
				part->ct += para.dt;	timeStep++;
				std::cout << " time --------> " << part->ct << std::endl;
				std::cout << " dt ----------> " << para.dt << std::endl;
			}
			saveData();
		}

		R stepGL() {
			auto* const part = derived().part;
			if (part->ct > para.tt) {
				saveData();
			}
			std::cout << " step ----------------------------------> " << timeStep << std::endl;
			R tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			part->updateCell();
			derived().step();
			part->ct += para.dt;	timeStep++;
			std::cout << " time --------> " << part->ct << std::endl;
			std::cout << " dt ----------> " << para.dt << std::endl;
			return part->ct;
		}

		void sensorOut() {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			sen.writeVect(derived().part);
			sen >> convert.str();
		}

		void profileOut() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			sen.writeScal(derived().part);
			sen.profile(rt, pf);
		}

		void profileOut_avgVel2() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			auto sum = 0.;
			auto count = 0;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					sum += part->vel2[p].squaredNorm();
					count++;
				}
			}
			sum = sum / count;
			std::ofstream file("./out/" + pf + ".out", std::ofstream::app);
			file << std::setprecision(6) << std::scientific << timeStep*para.cfl << " "
				<< std::setprecision(6) << std::scientific << sum
				<< std::endl;
			file.close();
			std::cout << " Writing profile. done " << std::endl;
		}

		void saveData() const {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().part) >> ("./out/" + convert.str() + ".out");
		}
		void saveData(const std::string& str) const {
			*(derived().part) >> ("./out/" + str + ".out");
		}

		void fina() {}
		
		__forceinline const Vec* position() const {
			return derived().part->pos.data();
		}
		__forceinline const R* scalar() const {
			return derived().part->phi.data();
		}
		__forceinline const int* type() const {
			return (int*)(derived().part->type.data());
		}

	public:
		Parameter<R,D> para;
		MatSolver<R,D>* mSol;
		Shifter<R,D> shi;
		Sensor<R,D> sen;
		PassiveScalarProblem<R,2> psp;

	protected:
		void step() {}
		void convect() {}
		void visTerm_e() {}
		void visTerm_i() {}
		void presTerm_e() {}
		void presTerm_i() {}
		void makeDirchlet_v() {}

		void makeNeumann_p() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				//part->neumann[p] = para.rho* (para.niu* part->lap(part->vel2, p) + para.g)* part->bdnorm.at(p);
				part->neumann[p] = 0.;
			}
		}

		void solvMat_p() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pres[p] = mSol->x[p];
				if (part->pres[p] < -1.e5) part->pres[p] = -1.e5;
				if (part->pres[p] > 1.e5) part->pres[p] = 1.e5;
			}
		}

		void solvMat_phi() {
			auto* const part = derived().part;
			mSol->ccBiCg_augment(part->type);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->phi[p] = mSol->x[p];
				if (part->phi[p] < -1.e5) part->phi[p] = -1.e5;
				if (part->phi[p] > 1.e5) part->phi[p] = 1.e5;
			}
		}

		void solvMat_v() {
			mSol->biCg_v();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				for (auto d = 0; d < D; d++) {
					part->vel2[p][d] = mSol->u[D*p + d];
				}
			}
		}

		R cfl() {
			R umax = 0.;
			const auto* const part = derived().part;
			for (unsigned p = 0; p < part->np; p++) {
				R tmp = part->vel1[p].norm();
				if (tmp > umax) umax = tmp;
			}
			para.umax = umax;
			return para.cfl * part->dp / umax;
		}

		void calBdNoSlip() {
			derived().part->bdNoSlip();
		}
		void bdSetZero() {
			derived().part->bdSetZero();
		}

		void calCell() {
			derived().part->updateCell();
		}

		void calPnd() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pnd[p] = part->cPnd(p);
				//part->pn[p] = part->cPn(p);
				//part->nbd[p] = part->cNbd(p);
			}
		}

		void calInvMat() {
			derived().part->updateInvMat();
		}

		void makeFs() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) part->fs[p] = part->_isFs(p);
			//derived().part->updateTeam();
		}

		void pthOrderPresSpatialFilter() {
			static auto counter = 0;
			if (counter++ % 1 != 0) return;
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == BD2 || part->isFs(p)) continue;
				part->pres[p] = part->func(part->pres, p);
			}
		}
		void pthOrderVelSpatialFilter() {
			static auto counter = 0;
			if (counter++ % 30 != 0) return;
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				part->vel2[p] = part->func(part->vel1, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				part->vel1[p] = part->vel2[p];
			}
		}

		void calForVis() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				part->vort[p] = part->rot(part->vel1, p);
			}
		}

		void check() const {
			const auto* const part = derived().part;
			R velMax = std::numeric_limits<R>::min();
			R phiMax = std::numeric_limits<R>::min();
			R divMax = std::numeric_limits<R>::min();
			unsigned idv = 0, idp = 0, idd = 0;
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) continue;
				const R vel = part->vel1[p].norm();
				const R phi = part->phi[p];
				const R div = part->div(part->vel1, p);
				if (vel > velMax) {
					velMax = vel;
					idv = p;
				}
				if (abs(phi) > abs(phiMax)) {
					phiMax = phi;
					idp = p;
				}
				if (part->type[p] == BD1) continue;
				if (abs(div) > abs(divMax)) {
					divMax = div;
					idd = p;
				}
			}
			std::cout << " max vel: " << velMax << " --- id: " << idv << std::endl;
			std::cout << " max phi: " << phiMax << " --- id: " << idp << std::endl;
			std::cout << " max Div: " << divMax << " --- id: " << idd << std::endl;
		}
#if 0
		void bvpSource() {
			const auto* const part = derived().part;
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->isFs(p)) {
					mSol->b[p] = bvp->func(part->pos[p]);
					continue;
				}
				mSol->b[p] = bvp->lap(part->pos[p]);
			}
		}

		void bvpAvgError() {
			const auto* const part = derived().part;
			int n = 0;
			R err = 0.;
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				R ext = bvp->func(part->pos[p]);
				err += abs(part->pres[p] - ext);
				n++;
			}
			std::cout << " bvp --- avg Error: " << err / n << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err / n << std::endl;
			file.close();
		}

		void bvpMaxError() {
			const auto* const part = derived().part;
			R err = 0.;
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				R ext = bvp->func(part->pos[p]);
				R tmp = abs(part->pres[p] - ext);
				if (tmp > err) err = tmp;
			}
			std::cout << " bvp --- max Error: " << err << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err << std::endl;
			file.close();
		}

		void gradMaxError() {
			auto* const part = derived().part;
			for (auto p = 0; p < part->np; p++) {
				part->pres[p] = bvp->func(part->pos[p]);
			}
			R err = 0.;
			for (auto p = 0; p < part->np; p++) {
				vec ext = bvp->grad(part->pos[p]);
				R tmp = (part->grad(part->pres, p) - ext).mag();
				if (tmp > err) err = tmp;
			}
			std::cout << " |grad| --- max Error: " << err << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err << std::endl;
			file.close();
		}

		void lapMaxError() {
			auto* const part = derived().part;
			for (auto p = 0; p < part->np; p++) {
				part->pres[p] = bvp->func(part->pos[p]);
			}
			R err = 0.;
			for (auto p = 0; p < part->np; p++) {
				R ext = bvp->lap(part->pos[p]);
				R tmp = abs(part->lap(part->pres, p) - ext);
				if (tmp > err) err = tmp;
			}
			std::cout << " lap --- max Error: " << err << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err << std::endl;
			file.close();
		}
#endif

		void insertRand() {
			auto* const part = derived().part;
			R coef = 0.1;
			std::default_random_engine gen;
			std::normal_distribution<R> dis(0., 0.5);
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				vec dr = coef* part->dp* vec(dis(gen), 0., dis(gen));
				part->pos[p] += dr;
			}
		}

	protected:
		int timeStep;
	};

}