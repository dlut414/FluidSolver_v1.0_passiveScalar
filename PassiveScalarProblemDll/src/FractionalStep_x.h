/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#pragma once
#include "Simulator.h"
#include "Particle_x.h"

namespace SIM {

	template <typename R, unsigned D, unsigned P>
	class FractionalStep_x : public Simulator <R,D,FractionalStep_x<R,D,P>> {
		typedef mMath::Polynomial_A<R,D,P> PN;
		typedef mMath::Derivative_A<R,D,P> DR;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<int,D,1>	iVec;
		typedef Eigen::Matrix<R,D,1>	Vec;
		typedef Eigen::Triplet<R>		Tpl;
	public:
		FractionalStep_x() {}
		~FractionalStep_x() {}

		void step() {
			updateVelocity();

			calInvMat();
			syncPos();

			convect_s1();

			calCell();
			calInvMat();
			calForVis();
			
			calCell();
			calInvMat();
			shift();

			sync();
		}

		void updateVelocity() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == FLUID) part->vel1[p] = part->vel2[p] = psp.velocity(part->pos[p]);
			}
		}

		void visTerm_e() {
			/*Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
		}

		void visTerm_i_q1r0() {
			makeLhs_v_q1();
			makeRhs_v_q1r0();
			solvMat_v();
		}

		void visTerm_i_q2r0() {
			makeLhs_v_q2();
			makeRhs_v_q2r0();
			solvMat_v();
		}


		void convect_s1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->pos[p] += para.dt * (part->vel1[p]);
			}
		}

		void convect_s2() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->pos[p] += 0.5* para.dt * (3.*part->vel1[p] - 1.* part->vel_m1[p]);
			}
		}

		void shift() {
			shi.shiftOriginWENO(part, part->phi);
		}

		void init_() {
			part = new Particle_x<R,D,P>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k, para.beta);
			part->buildCell();
			part->b2b();
			part->b2norm();
			//part->updateTeam();
			part->init_x();
			initialPassiveScalar();
		}

	public:
		Particle_x<R,D,P>* part;

	private:
		__forceinline void initialPassiveScalar() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == FLUID) part->phi[p] = psp.scalar(part->pos[p]);
			}
		}

		__forceinline void makeLhs_v_q2() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (int d = 0; d < D; d++) {
						coef.push_back(Tpl(D*p + d, D*p + d, 1.));
					}
					continue;
				}
				auto pp = 0.;
				const auto& mm = part->invMat[p];
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (part->bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						const auto w = part->w3(dr1);
						VecP npq;
						part->poly(dr, npq);
						const auto a = mm * (w* npq);
						const auto& lp = part->pn_lap_o;
						const auto pq = -para.niu* lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						for (auto d = 0; d < D; d++) {
							coef.push_back(Tpl(D*p + d, D*q + d, pq));
						}
					}
				}
				pp += 3. / (2. * para.dt);
				for (auto d = 0; d < D; d++) {
					coef.push_back(Tpl(D*p + d, D*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeLhs_v_q1() {
			coef.clear();
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						coef.push_back(Tpl(D*p + d, D*p + d, 1.));
					}
					continue;
				}
				auto pp = 0.;
				const auto& mm = part->invMat[p];
				const auto& c = part->cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
						const auto q = part->cell->linkList[key][m];
#if BD_OPT
						if (part->bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.mag();
						if (dr1 > part->r0) continue;
						if (q == p) continue;
						const auto w = part->w3(dr1);
						VecP npq;
						part->poly(dr, npq);
						const auto a = mm * (w* npq);
						const auto& lp = part->pn_lap_o;
						const auto pq = -para.niu* lp.dot(a);
						pp -= pq;
						for (auto d = 0; d < D; d++) {
							coef.push_back(Tpl(D*p + d, D*q + d, pq));
						}
					}
				}
				pp += 1. / para.dt;
				for (auto d = 0; d < D; d++) {
					coef.push_back(Tpl(D*p + d, D*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}


		__forceinline void makeRhs_v_q2r1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[d*p + d] = part->vel1[p][d];
					}
					continue;
				}
				const auto gd = part->grad(part->pres, p);
				const auto rhs = 1. / (2.* para.dt)* (4.* part->vel1[p] - part->vel_m1[p])
					- (1. / para.rho)* gd
					+ para.g;
				for (auto d = 0; d < D; d++) {
					mSol->rhs[d*p + d] = part->rhs[d];
				}
			}
		}

		__forceinline void makeRhs_v_q1r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[d*p + d] = part->vel1[p][d];
					}
					continue;
				}
				const auto rhs = (1. / para.dt)* part->vel1[p] + para.g;
				for (auto d = 0; d < D; d++) {
					mSol->rhs[d*p + d] = part->rhs[d];
				}
			}
		}

		__forceinline void makeRhs_v_q2r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[D*p + d] = part->vel1[p][d];
					}
					continue;
				}
				const auto rhs = 1. / (2.* para.dt)* (4.* part->vel1[p] - part->vel_m1[p]) + para.g;
				for (auto d = 0; d < D; d++) {
					mSol->rhs[D*p + d] = rhs[d];
				}
			}
		}


		__forceinline void sync() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->vel_m1[p] = part->vel1[p];
				part->vel1[p] = part->vel2[p];
			}
		}
		__forceinline void syncPos() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pos_m1[p] = part->pos[p];
			}
		}

	private:
		std::vector<Tpl> coef;
	};

}