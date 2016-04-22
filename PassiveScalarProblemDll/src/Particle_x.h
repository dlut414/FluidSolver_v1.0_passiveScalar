/*
*/
#pragma once
#include "Header.h"
#include "Particle.h"
#include "Polynomial.h"
#include "Derivative.h"

#define UPWIND_VEL 0

namespace SIM {

	template <typename R, unsigned D, unsigned P>
	class Particle_x : public Particle<R,D,Particle_x<R,D,P>> {
		typedef mMath::Polynomial_A<R, D, P> PN;
		typedef mMath::Derivative_A<R, D, P> DR;
		typedef Eigen::Matrix<int, D, 1>	iVec;
		typedef Eigen::Matrix<R, D, 1>		Vec;
		typedef Eigen::Matrix<R, D, D>		Mat;
		typedef Eigen::Matrix<R, PN::value, 1>	VecP;
		typedef Eigen::Matrix<R, PN::value, D>	MatPD;
		typedef Eigen::Matrix<R, PN::value, PN::value> MatPP;
	public:
		Particle_x() : Particle() {}
		~Particle_x() {}
		
		__forceinline void poly(const Vec& in, VecP& out) const { PN::Run(varrho, in.data(), out.data()); }

		const Vec grad(const std::vector<R>& phi, const unsigned& p) const {
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_p_o*a);
		}

		const Mat grad(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w* npq)* (u[q] - u[p]).transpose();
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_p_o*a);
		}

		const R div(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w* npq)* (u[q] - u[p]).transpose();
				}
			}
			const auto a = invMat[p] * vv;
			auto ret = static_cast<R>(0);
			for (auto d = 0; d < D; d++) {
				ret += pn_p_o.block<1, PN::value>(d, 0) * a.block<PN::value, 1>(0, d);
			}
			return ret;
		}

		const R lap(const std::vector<R>& phi, const unsigned& p) const {
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP	npq;
					poly(dr, npq);
					vv += w * (phi[q] - phi[p])* npq;
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_lap_o*a);
		}

		const Vec lap(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w * npq) * (u[q] - [p]);
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_lap_o*a).transpose();
		}

		const R rot(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w * npq) * (u[q] - u[p]).transpose();
				}
			}
			const auto a = invMat[p] * vv;
			const Mat der = pn_p_o*a;
			switch (D) {
			case 1:
				return R(0);
				break;
			case 2:
				return R(der(0, 1) - der(1, 0));
				break;
			case 3:
				return R(0);
				break;
			default:
				return R(0);
			}
		}

		template <typename U>
		const U func(const std::vector<U>& phi, const unsigned& p) const {
			return phi[p];
		}

		template <typename U>
		const U func(const std::vector<U>& phi, const Vec& p) const {
			auto rid = 0;
			auto isNear = 0;
			auto rr = std::numeric_limits<R>::max();
			auto c = cell->iCoord(p);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
					const auto dr = pos[q] - p;
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					else {
						isNear = 1;
						if (dr1 < rr) {
							rr = dr1;
							rid = q;
						}
					}
				}
			}
			if (!isNear) return U::Zero();
			const auto dp = p - pos[rid];
			const auto dPhi = (dp.transpose()* grad(phi, rid)).transpose();
			return phi[rid] + dPhi;
		}

		const Vec func_mafl(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			Vec pLocal[6];

			for (auto i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			Vec ret[6];
			ret[3] = phi[p];
			for (auto fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = Vec::Zero();
				auto ww = static_cast<R>(0);
				auto c = cell->iCoord(pLocal[fp]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
						if (type[q] == BD2) continue;
#if BD_OPT
						if (bdOpt(q)) continue;
#endif
						const auto dr1 = (pos[q] - pLocal[fp]).norm();
						const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).norm();
						const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).norm();
						if (dr1 > re) continue;
						if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
						const auto w = w1(dr1);
						ww += w;
						ret[fp] += w * phi[q];
					}
				}
				if (abs(ww) < eps) ww = 1.;
				ret[fp] = ret[fp] / ww;
			}
			return ret[3] - (dmove.norm() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
		}

		const Vec func_mafl_mmt(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			Vec pLocal[6];

			for (auto i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			Vec ret[6];
			ret[3] = phi[p];
			for (auto fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = Vec::Zero();
				auto ww = static_cast<R>(0);
				auto c = cell->iCoord(pLocal[fp]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
						if (type[q] == BD2) continue;
#if BD_OPT
						if (bdOpt(q)) continue;
#endif
						const auto dr1 = (pos[q] - pLocal[fp]).norm();
						const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).norm();
						const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).norm();
						if (dr1 > re) continue;
						if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
						const auto w = w1(dr1);
						ww += w;
						ret[fp] += w * phi[q];
					}
				}
				if (abs(ww) < eps) continue;
				ret[fp] = ret[fp] / ww;
			}
			auto ret_mmt = ret[3] - (dmove.norm() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
			auto ret_min = ret[1];
			auto ret_max = ret[1];
			for (int fp = 1; fp <= 4; fp++) {
				if (ret[fp].squaredNorm() < ret_min.squaredNorm()) ret_min = ret[fp];
				if (ret[fp].squaredNorm() > ret_max.squaredNorm()) ret_max = ret[fp];
			}
			if (ret_mmt.squaredNorm() < ret_min.squaredNorm()) ret_mmt = ret_min* ret_mmt.norm();
			if (ret_mmt.squaredNorm() > ret_max.squaredNorm()) ret_mmt = ret_max* ret_mmt.norm();
			return ret_mmt;
		}

		const R func_lsA(const std::vector<R>& phi, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* npq)* (phi[q] - phi[p]);
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = MatPP::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
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

		const Vec func_lsA(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			auto mm = MatPP::Zero();
			auto vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* npq)* (phi[q] - phi[p]).transpose();
				}
			}
			const auto inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = MatPP::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
			const Mat hes[D];
			for (auto d = 0; d < D; d++) {
				for (auto i = 0; i < D; i++) {
					for (auto j = i; j < D; j++) {
						hes[d](i, j) = mgd.block<mMath::H<D,2>,1>(0, d);
					}
				}
				for (auto i = 0; i < D; i++) {
					for (auto j = 0; j < D; j++) {
						hes[d](i, j) = hes[d](j, i);
					}
				}
			}
			const auto dp = p_new - pos[p];
			const auto dpt = dp.transpose();
			auto ret = phi[p];
			for (auto d = 0; d < D; d++) ret[d] += (dpt*gd).transpose() + 0.5*dpt*hes[d]*dp;
			return ret;
		}

		const Vec func_lsA_upwind(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(phi[p].norm());
#else
			const auto up = dp.normalized();
#endif
			MatPP mm = MatPP::Zero();
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					if (dr.dot(up) < 0) continue;
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* npq)* (phi[q] - phi[p]).transpose();
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = MatPP::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
			Mat hes[D];
			for (auto d = 0; d < D; d++) {
				int counter = 0;
				for (auto i = 0; i < D; i++) {
					for (auto j = i; j < D; j++) {
						hes[d](i, j) = mgd(counter++, d);
					}
				}
				for (auto i = 0; i < D; i++) {
					for (auto j = 0; j < i; j++) {
						hes[d](i, j) = hes[d](j, i);
					}
				}
			}
			const auto dpt = dp.transpose();
			auto ret = phi[p];
			for (auto d = 0; d < D; d++) ret[d] += (dpt*gd)[d] + 0.5*dpt * hes[d] * dp;
			return ret;
		}

		const R func_lsA_upwind(const std::vector<R>& phi, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(phi[p].norm());
#else
			const auto up = dp.normalized();
#endif
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					if (dr.dot(up) < 0) continue;
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* npq)* (phi[q] - phi[p]);
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = MatPP::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
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
			ret = ret + (dpt*gd) +0.5*dpt * hes * dp;
			return ret;
		}

		//const Vec func_lsB(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {}

		//const Vec func_lsB_upwind(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {}
		
		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY, int Dimension = D>	struct interpolateWENO_A_ {
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_A_<StencilsX, StencilsY, Stencils, 1> {
			template <typename U> static const U Run(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {}
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_A_<StencilsX, StencilsY, Stencils, 2> {
			template <typename U> static const U Run(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R,D,P>* part) {
				const auto dp = p_new - part->pos[p];
				if (dp.norm() < part->eps) return phi[p];
				const auto up = dp.normalized();
				const auto alpha = 2.* M_PI / StencilsX;
				Vec dir[StencilsX];
				Vec ctr[Stencils];
				for (auto i = 0; i < StencilsX; i++) {
					const auto theta = i* alpha;
					const auto ct = cos(theta);
					const auto st = sin(theta);
					dir[i] << ct*up[0] + st*up[1], ct*up[1] - st*up[0];
				}
				for (auto j = 0; j < StencilsY; j++) {
					//const auto dis = part->r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
					const auto dis = part->r0* ( R(1.) - R(1.)*(j + 1) / (StencilsY) );
					for (auto i = 0; i < StencilsX; i++) {
						const auto stcId = i* StencilsY + j;
						ctr[stcId] = part->pos[p] + dis*dir[i];
					}
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
						for (auto stcId = 0; stcId < Stencils; stcId++) {
							const auto dis = (part->pos[q] - ctr[stcId]).norm();
							const auto w = part->w3(dis);
							VecP npq;
							part->poly(dr, npq);
							mm[stcId] += (w* npq)* npq.transpose();
							vv[stcId] += (w* npq)* (phi[q] - phi[p]);
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
						return phi[p];
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
				const R epsilon = 1.e-6;
				const int magnifier = 5;
				for (auto i = 0; i < Stencils; i++) {
					stencilWeight[i] = 1. / pow(epsilon + oscillationIndicator[i], magnifier);
				}
				R stencilWeightSum = R(0.);
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightSum += stencilWeight[i];
				}
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
					//if (p == 8625)std::cout << stencilWeightNorm[i] << ", ";
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
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_A_<StencilsX, StencilsY, Stencils, 3> {
			template <typename U> static const U Run(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {}
		};
		
		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY, int Dimension = D>	struct interpolateWENO_B_ {
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_B_<StencilsX, StencilsY, Stencils, 1> {
			template <typename U> static const U Run(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {}
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_B_<StencilsX, StencilsY, Stencils, 2> {
			template <typename U> static const U Run(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {
				const auto dp = p_new - part->pos[p];
				if (dp.norm() < part->eps) return phi[p];
				const auto up = dp.normalized();
				const auto alpha = 2.* M_PI / StencilsX;
				Vec dir[StencilsX];
				Vec ctr[Stencils];
				for (auto i = 0; i < StencilsX; i++) {
					const auto theta = i* alpha;
					const auto ct = cos(theta);
					const auto st = sin(theta);
					dir[i] << ct*up[0] + st*up[1], ct*up[1] - st*up[0];
				}
				for (auto j = 0; j < StencilsY; j++) {
					//const auto dis = part->r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
					const auto dis = part->r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
					for (auto i = 0; i < StencilsX; i++) {
						const auto stcId = i* StencilsY + j;
						ctr[stcId] = part->pos[p] + dis*dir[i];
					}
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
						for (auto stcX = 0; stcX < StencilsX; stcX++) {
							for (auto stcY = 0; stcY < StencilsY; stcY++) {
								const auto stcId = stcX* StencilsY + stcY;
								const auto dis = (part->pos[q] - ctr[stcId]).norm();
								const auto w = part->w3(dis);
								VecP npq;
								part->poly(dr, npq);
								mm[stcId] += (w* npq)* npq.transpose();
								vv[stcId] += (w* npq)* (phi[q] - phi[p]);
							}
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
						return phi[p];
						auto mm_ = mm[i].block<2, 2>(0, 0);
						if (abs(mm_.determinant()) < part->eps_mat) {
							inv = MatPP::Zero();
						}
						else inv.block<2, 2>(0, 0) = mm_.inverse();
					}
					else inv = mm[i].inverse();
					polyCoef[i] = inv * vv[i];
					oscillationIndicator[i] = R(0.);
				}

				for (auto i = 0; i < Stencils; i++) {
					const auto dp = part->dp;
					const auto A = polyCoef[i][0] * polyCoef[i][0];
					const auto B = polyCoef[i][1] * polyCoef[i][1];
					const auto C = polyCoef[i][2] * polyCoef[i][2];
					const auto D = polyCoef[i][3] * polyCoef[i][3];
					const auto E = polyCoef[i][4] * polyCoef[i][4];
					const auto beta1 = (dp*dp*dp)*(A + B) + (dp*dp*dp*dp*dp / 6.)*(2.*C + D + 2.*E);
					const auto beta2 = (dp*dp*dp*dp*dp)*(4.*C + D + 4.*E);
					//oscillationIndicator[i] = 4.*beta1 - (1. / 3.)*beta2;
					oscillationIndicator[i] = beta1 + beta2;
				}
				const R epsilon = 1.e-6;
				const int magnifier = 5;
				for (auto i = 0; i < Stencils; i++) {
					stencilWeight[i] = 1. / pow(epsilon + oscillationIndicator[i], magnifier);
				}
				R stencilWeightSum = R(0.);
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightSum += stencilWeight[i];
				}
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
					if (p == 8625)std::cout << stencilWeightNorm[i] << ", ";
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
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_B_<StencilsX, StencilsY, Stencils, 3> {
			template <typename U> static const U Run(const std::vector<U>& phi, const unsigned& p, const Vec& p_new, Particle_x<R, D, P>* part) {}
		};

		__forceinline const R interpolateWENO(const std::vector<R>& phi, const unsigned& p, const Vec& p_new) {
			return interpolateWENO_B_<>::Run(phi, p, p_new, this);
		}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				MatPP mm = MatPP::Zero();
				const auto c = cell->iCoord(pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (bdOpt(p, q)) continue;
#endif
						const auto dr = pos[q] - pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > r0) continue;
						const auto w = w3(dr1);
						VecP npq;
						poly(dr, npq);
						mm += (w* npq) * npq.transpose();
					}
				}
				invMat[p] = MatPP::Zero();
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
					auto mm_ = mm.block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) {
						invMat[p] = MatPP::Zero();
						continue;
					}
					invMat[p].block<2, 2>(0, 0) = mm_.inverse();
					continue;
				}
				invMat[p] = mm.inverse();
			}
		}

		template <unsigned D_ = D>
		void init_x() {}

		template <>
		void init_x<1>() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(MatPP());
			}

			varrho = 1./(1.*dp);
			Vec zero = Vec::Zero();
			DR::Run<1>(varrho, zero.data(), pn_p_o.data());
			DR::Run<2>(varrho, zero.data(), pn_pp_o.data());
			DR::Run<2>(varrho, zero.data(), pn_lap_o.data());
		}

		template <>
		void init_x<2>() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(MatPP());
			}

			varrho = 1./(1.*dp);
			Vec zero = Vec::Zero();
			DR::Run<1, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(0, 0).data());
			DR::Run<0, 1>(varrho, zero.data(), pn_p_o.block<1, PN::value>(1, 0).data());
			DR::Run<2, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(0, 0).data());
			DR::Run<1, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(1, 0).data());
			DR::Run<0, 2>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(2, 0).data());
			pn_lap_o = pn_pp_o.block<1, PN::value>(0, 0) + pn_pp_o.block<1, PN::value>(2, 0);
		}

		template <>
		void init_x<3>() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(MatPP());
			}

			varrho = 1./(1.*dp);
			Vec zero = Vec::Zero();
			DR::Run<1, 0, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(0, 0).data());
			DR::Run<0, 1, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(1, 0).data());
			DR::Run<0, 0, 1>(varrho, zero.data(), pn_p_o.block<1, PN::value>(2, 0).data());
			DR::Run<2, 0, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(0, 0).data());
			DR::Run<1, 1, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(1, 0).data());
			DR::Run<1, 0, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(2, 0).data());
			DR::Run<0, 2, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(3, 0).data());
			DR::Run<0, 1, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(4, 0).data());
			DR::Run<0, 0, 2>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(5, 0).data());
			pn_lap_o = pn_pp_o.block<1, PN::value>(0, 0) + pn_pp_o.block<1, PN::value>(3, 0) + pn_pp_o.block<1, PN::value>(5, 0);
		}

		void init_x() { init_x<>(); }

	public:
		std::vector<MatPP> invMat;

		R varrho;
		Eigen::Matrix<R,D,PN::value,Eigen::RowMajor>					pn_p_o;
		Eigen::Matrix<R,mMath::H<D,2>::value,PN::value,Eigen::RowMajor>	pn_pp_o;
		Eigen::Matrix<R,1,PN::value,Eigen::RowMajor>					pn_lap_o;
	};

}