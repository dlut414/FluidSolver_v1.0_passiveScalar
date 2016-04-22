/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//Interpolation.h
///defination of class Interpolation
#pragma once

namespace SIM {

	template <typename R, int D, int P>
	class Interpolation {
		typedef Particle_x<R,D,P> Ptc;
	public:
		Interpolation() {}
		~Interpolation() {}

		template <int Stencils = 8, int Dimension = D>	struct interpolateWENO_ {};
		template <int Stencils>							struct interpolateWENO_<Stencils, 1> {
			template <typename U> static const U Gen(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {}
		};
		template <int Stencils>							struct interpolateWENO_<Stencils, 2> {
			template <typename U> static const U Gen(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {
				const auto dp = p_new - part->pos[p];
				const auto up = dp.normalized();
				static const auto alpha = 2.* M_PI / Stencils;
				Vec dir[Stencils];
				for (auto i = 0; i < Stencils; i++) {
					const auto theta = i* alpha;
					const auto ct = cos(theta);
					const auto st = sin(theta);
					dir[i] << ct*up[0] + st*up[1], ct*up[1] - st*up[0];
				}
				MatPP mm[Stencils];
				VecP vv[Stencils];
				for (auto i = 0; i < Stencils; i++) {
					mm[i] = MatPP::Zero();
					vv[i] = VecP::Zero();
				}
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						const auto w = part->w3(dr1);
						for (auto stc = 0; stc < Stencils; stc++) {
							if (dr.dot(dir[stc]) < 0) continue;
							VecP npq;
							part->poly(dr, npq);
							mm[stc] += (w* npq)* npq.transpose();
							vv[stc] += (w* npq)* (phi[q] - phi[p]);
						}
					}
				}
				R oscillationIndicator[Stencils];
				R stencilWeight[Stencils];
				R stencilWeightNorm[Stencils];
				VecP polyCoef[Stencils];
				for (auto i = 0; i < Stencils; i++) {
					MatPP inv = MatPP::Zero();
					if (abs(mm[i].determinant()) < part->eps_mat) {
						auto mm_ = mm[i].block<2, 2>(0, 0);
						if (abs(mm_.determinant()) < part->eps_mat) {
							inv = MatPP::Zero();
						}
						else inv.block<2, 2>(0, 0) = mm_.inverse();
					}
					else inv = mm[i].inverse();
					polyCoef[i] = inv * vv[i];
					const auto offset = PN::value - mMath::H<D, P>::value;
					oscillationIndicator[i] = R(0.);
					for (auto term = offset; term < PN::value; term++) {
						oscillationIndicator[i] += abs(polyCoef[i][term]);
					}
				}
				static const R epsilon = 1.e-5;
				static const int magnifier = 2;
				for (auto i = 0; i < Stencils; i++) {
					stencilWeight[i] = 1. / pow(epsilon + oscillationIndicator[i], magnifier);
				}
				R stencilWeightSum = R(0.);
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightSum += stencilWeight[i];
				}
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
				}
				VecP combinedCoef = VecP::Zero();
				for (auto i = 0; i < Stencils; i++) {
					combinedCoef += stencilWeightNorm[i] * polyCoef[i];
				}
				const auto gd = part->pn_p_o * combinedCoef;
				const auto mgd = part->pn_pp_o * combinedCoef;
				Mat hes;
				int counter = 0;
				for (auto i = 0; i < D; i++) {
					for (auto j = i; j < D; j++) {
						hes(i, j) = mgd(counter++);
					}
				}
				for (auto i = 0; i < D; i++) {
					for (auto j = 0; j < i; j++) {
						hes(i, j) = hes(j, i);
					}
				}
				const auto dpt = dp.transpose();
				auto ret = phi[p];
				ret = ret + (dpt*gd) + 0.5*dpt * hes * dp;
				return ret;
			}
		};
		template <int Stencils>							struct interpolateWENO_<Stencils, 3> {
			template <typename U> static const U Gen(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {}
		};

		__forceinline const R interpolateWENO(const std::vector<R>& phi, const unsigned& p, const Vec& p_new) {
			return interpolateWENO_<>::Gen(phi, p, p_new, this);
		}
	};

}