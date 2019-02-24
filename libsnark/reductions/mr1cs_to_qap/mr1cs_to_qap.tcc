/** @file
 *****************************************************************************
 Test program that exercises the ppzkSNARK (first generator, then
 prover, then verifier) on a synthetic mr1cs instance.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <cassert>
#include <cstdio>
#include<iostream>
#include<algorithm>
#include<set>

#include<libff/algebra/fields/bigint.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>
#include <libfqfft/polynomial_arithmetic/xgcd.hpp>


#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>

#include<stdlib.h>

using namespace std;

namespace libsnark{



	/**
	 * Instance map for the mr1cs-to-QAP reduction.
	 *
	 * Namely, given a mr1cs constraint system cs, construct a QAP instance for which:
	 *   A := (A_0(z),A_1(z),...,A_m(z))
	 *   B := (B_0(z),B_1(z),...,B_m(z))
	 *   C := (C_0(z),C_1(z),...,C_m(z))
	 * where
	 *   m = number of variables of the QAP
	 * and
	 *   each A_i,B_i,C_i is expressed in the Lagrange basis.
	 */
	template<typename FieldT>
		 qap_instance<FieldT> mr1cs_to_mqap_instance_map(const mr1cs_constraint_system<FieldT> &cs)
	{
		libff::enter_block("Call to mr1cs_to_mqap_instance_map");

		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);

		std::vector<std::map<size_t, FieldT> > A_in_Lagrange_basis(cs.num_variables()+1);
		std::vector<std::map<size_t, FieldT> > B_in_Lagrange_basis(cs.num_variables()+1);
		std::vector<std::map<size_t, FieldT> > C_in_Lagrange_basis(cs.num_variables()+1);

		libff::enter_block("Compute polynomials A, B, C in Lagrange basis");
		/**
		 * add and process the constraints
		 *     input_i * 0 = 0
		 * to ensure soundness of input consistency
		 */
		for (size_t i = 0; i <= cs.num_inputs(); ++i)
		{
			A_in_Lagrange_basis[i][cs.num_constraints() + i] = FieldT::one();
		}
		/* process all other constraints */
		for (size_t i = 0; i < cs.num_constraints(); ++i)
		{
			for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j)
			{
				A_in_Lagrange_basis[cs.constraints[i].a.terms[j].index][i] +=
					cs.constraints[i].a.terms[j].coeff;
			}

			for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j)
			{
				B_in_Lagrange_basis[cs.constraints[i].b.terms[j].index][i] +=
					cs.constraints[i].b.terms[j].coeff;
			}

			for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j)
			{
				C_in_Lagrange_basis[cs.constraints[i].c.terms[j].index][i] +=
					cs.constraints[i].c.terms[j].coeff;
			}
		}
		libff::leave_block("Compute polynomials A, B, C in Lagrange basis");

		libff::leave_block("Call to mr1cs_to_mqap_instance_map");

		return qap_instance<FieldT>(domain,
				cs.num_variables(),
				domain->m,
				cs.num_inputs(),
				std::move(A_in_Lagrange_basis),
				std::move(B_in_Lagrange_basis),
				std::move(C_in_Lagrange_basis));
	}

	/**
	 * Instance map for the mr1cs-to-QAP reduction followed by evaluation of the resulting QAP instance.
	 *
	 * Namely, given a mr1cs constraint system cs and a field element t, construct
	 * a QAP instance (evaluated at t) for which:
	 *   At := (A_0(t),A_1(t),...,A_m(t))
	 *   Bt := (B_0(t),B_1(t),...,B_m(t))
	 *   Ct := (C_0(t),C_1(t),...,C_m(t))
	 *   Ht := (1,t,t^2,...,t^n)
	 *   Zt := Z(t) = "vanishing polynomial of a certain set S, evaluated at t"
	 * where
	 *   m = number of variables of the QAP
	 *   n = degree of the QAP
	 */
	template<typename FieldT>
		qap_instance_evaluation<FieldT> mr1cs_to_mqap_instance_map_with_evaluation(const mr1cs_constraint_system<FieldT> &cs,
				const FieldT &t)
		{
			libff::enter_block("Call to mr1cs_to_mqap_instance_map_with_evaluation");

			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs()-1)*(cs.num_convol()));

			std::vector<FieldT> At, Bt, Ct, Ht;

			At.resize(cs.num_variables()+1, FieldT::zero());
			Bt.resize(cs.num_variables()+1, FieldT::zero());
			Ct.resize(cs.num_variables()+1, FieldT::zero());
			Ht.reserve(domain->m+1);

			const FieldT Zt = domain->compute_vanishing_polynomial(t);

			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2 = libfqfft::get_evaluation_domain<FieldT>((domain->m)*4+1);

			/*
			FieldT mvT[] = new FieldT[cs.num_inputs()-1];
			FieldT temp = (t^domain2->m);
			mvT[0] = FieldT::one();
			for(size_t i = 1; i<cs.num_inputs()-1;i++)
			{
				mvT[i] *= temp;
			}
			*/
			libff::enter_block("Compute evaluations of A, B, C, H at t");
			const std::vector<FieldT> u = domain->evaluate_all_lagrange_polynomials(t);
			/**
			 * add and process the constraints
			 *     input_i * 0 = 0
			 * to ensure soundness of input consistency
			 */
			for (size_t i = 0; i <= cs.num_inputs(); ++i)
			{
				At[i] = u[cs.num_constraints() + i];
			}

			size_t conv_index = 0;
			/* process all other constraints */
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j)
				{
					At[cs.constraints[i].a.terms[j].index] +=
						u[i]*cs.constraints[i].a.terms[j].coeff;	
				}
				/*
				for (size_t j = 0; j < cs.constraints[i].a.terms2.size(); ++j)
				{
					At[cs.constraints[i].a2.terms[j].index] +=
						u[i]*(mvT[i]);
				}
				*/


				for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j)
				{
					Bt[cs.constraints[i].b.terms[j].index] +=
						u[i]*cs.constraints[i].b.terms[j].coeff;
				}

				/*
				for (size_t j = 0; j < cs.constraints[i].b.terms2.size(); ++j)
				{
					Bt[cs.constraints[i].b2.terms[j].index] +=
						u[i]*mvT[i];
				}
				*/

				for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j)
				{
					Ct[cs.constraints[i].c.terms[j].index] +=
						u[i]*cs.constraints[i].c.terms[j].coeff;
				}

				/*
				for (size_t j = 0; j < cs.constraints[i].c.terms2.size(); ++j)
				{
					Ct[cs.constraints[i].c2.terms[j].index] +=
						u[i]*mvT[i];
				}
				*/
				for(size_t j=0;j<cs.constraints[i].a2.terms.size();j++){
					///TODO change (j+1) to safe value like prime
					if(j==0) conv_index ++;
					for(size_t k=0;k<cs.num_convol_outputs();k++){
						At[cs.constraints[i].a2.terms[j].index] +=
							u[cs.num_constraints() + j*(conv_index) + k ] * ((j+1)^cs.constraints[i].a2.terms[j].coeff);
					}
				}
				for(size_t j=0;j<cs.constraints[i].b2.terms.size();j++){
					///TODO change (j+1) to safe value like prime
					if(j==0) conv_index ++;
					for(size_t k=0;k<cs.num_convol_outputs();k++){
						Bt[cs.constraints[i].b2.terms[j].index] +=
							u[cs.num_constraints() + j*(conv_index) + k ] * ((j+1)^cs.constraints[i].b2.terms[j].coeff);
					}
				}

				for(size_t j=0;j<cs.constraints[i].c2.terms.size();j++){
					///TODO change (j+1) to safe value like prime
					if(j==0) conv_index ++;
					for(size_t k=0;k<cs.num_convol_outputs();k++){
						Ct[cs.constraints[i].c2.terms[j].index] +=
							u[cs.num_constraints() + j*(conv_index) + k ] * ((j+1)^cs.constraints[i].c2.terms[j].coeff);
					}
				}

			}
			

			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>((domain2->m)*2);

			FieldT ti = FieldT::one();
			for (size_t i = 0; i < domain->m;i++)//domain4->m; ++i)
			{
				Ht.emplace_back(ti);
				ti *= t;
			}
			libff::leave_block("Compute evaluations of A, B, C, H at t");

			libff::leave_block("Call to mr1cs_to_mqap_instance_map_with_evaluation");

			return qap_instance_evaluation<FieldT>(domain,
					cs.num_variables(),
					domain->m,
					//domain2->m,
					//domain4->m,
					cs.num_inputs(),
					//cs.num_kernels(),
					t,
					std::move(At),
					std::move(Bt),
					std::move(Ct),
					std::move(Ht),
					Zt);
		}

	/**
	 * Witness map for the mr1cs-to-QAP reduction.
	 *
	 * The witness map takes zero knowledge into account when d1,d2,d3 are random.
	 *
	 * More precisely, compute the coefficients
	 *     h_0,h_1,...,h_n
	 * of the polynomial
	 *     H(z) := (A(z)*B(z)-C(z))/Z(z)
	 * where
	 *   A(z) := A_0(z) + \sum_{k=1}^{m} w_k A_k(z) + d1 * Z(z)
	 *   B(z) := B_0(z) + \sum_{k=1}^{m} w_k B_k(z) + d2 * Z(z)
	 *   C(z) := C_0(z) + \sum_{k=1}^{m} w_k C_k(z) + d3 * Z(z)
	 *   Z(z) := "vanishing polynomial of set S"
	 * and
	 *   m = number of variables of the QAP
	 *   n = degree of the QAP
	 *
	 * This is done as follows:
	 *  (1) compute evaluations of A,B,C on S = {sigma_1,...,sigma_n}
	 *  (2) compute coefficients of A,B,C
	 *  (3) compute evaluations of A,B,C on T = "coset of S"
	 *  (4) compute evaluation of H on T
	 *  (5) compute coefficients of H
	 *  (6) patch H to account for d1,d2,d3 (i.e., add coefficients of the polynomial (A d2 + B d1 - d3) + d1*d2*Z )
	 *
	 * The code below is not as simple as the above high-level description due to
	 * some reshuffling to save space.
	 */
	template<typename FieldT>
		qap_witness<FieldT> mr1cs_to_qap_witness_map(const mr1cs_constraint_system<FieldT> &cs,
				const mr1cs_primary_input<FieldT> &primary_input,
				const mr1cs_auxiliary_input<FieldT> &auxiliary_input,
				const FieldT &d1,
				const FieldT &d2,
				const FieldT &d3)
		{
			libff::enter_block("Call to mr1cs_to_qap_witness_map");

			/* sanity check */
			assert(cs.is_satisfied(primary_input, auxiliary_input));

			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs()-1)*(cs.num_convol()));
			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2 = libfqfft::get_evaluation_domain<FieldT>((domain->m*4)+1);
			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>((domain2->m*2));

			mr1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
			full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());

			libff::enter_block("Compute evaluation of polynomials A, B on set S");
			std::vector<FieldT> aA(domain->m, FieldT::zero()), aB(domain->m, FieldT::zero());

			/* account for the additional constraints input_i * 0 = 0 */
			for (size_t i = 0; i <= cs.num_inputs(); ++i)
			{
				aA[i+cs.num_constraints()] = (i > 0 ? full_variable_assignment[i-1] : FieldT::one());
			}
			/* account for all other constraints */
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				aA[i] += cs.constraints[i].a.evaluate(full_variable_assignment);
				//aA[i] += cs.constraints[i].a2.evaluate(full_variable_assignment);

				aB[i] += cs.constraints[i].b.evaluate(full_variable_assignment);
				//aB[i] += cs.constraints[i].b2.evaluate(full_variable_assignment);
			}
			size_t conv_index = 0;
			size_t conv_index2 = 0;
			for( size_t i=0;i<cs.num_constraints();i++){
				FieldT acc = FieldT::zero();
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					for(auto &lt : cs.constraints[i].a2){
						//if(lt.idex ==0) acc += (FieldT::one()*cs.conv_num_output()) * lt.coeff; //there is no first one in conv constraints
						acc += full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff);
					}
					//if(acc != FieldT::zero()){
						aA[conv_index*cs.num_convol_outputs() + j + cs.num_constraints() + cs.num_inputs()] += acc;
						conv_index++;
					//}
				}
				acc = FieldT::zero();
				for(size_t j=1;j<=cs.conv_num_output();j++){
					for(auto &lt : cs.constraints[i].b2){
						acc += full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff);
					}
					aB[conv_index*cs.num_convol_outputs() + j + cs.num_constraints() + cs.num_inputs()] += acc;
					conv_index2++;
					//aB[i*cs.num_convol() + j + cs.num_constraints() + cs.num_inputs()] += acc;
				}
			}
			libff::leave_block("Compute evaluation of polynomials A, B on set S");

			libff::enter_block("Compute coefficients of polynomial A");
			domain->iFFT(aA);
			libff::leave_block("Compute coefficients of polynomial A");

			libff::enter_block("Compute coefficients of polynomial B");
			domain->iFFT(aB);
			libff::leave_block("Compute coefficients of polynomial B");

			libff::enter_block("Compute ZK-patch");
			std::vector<FieldT> coefficients_for_H(domain->m+1, FieldT::zero());
#ifdef MULTICORE
#pragma omp parallel for
#endif
			/* add coefficients of the polynomial (d2*A + d1*B - d3) + d1*d2*Z */
			for (size_t i = 0; i < domain->m; ++i)
			{
				coefficients_for_H[i] = d2*aA[i] + d1*aB[i];
			}
			coefficients_for_H[0] -= d3;
			domain->add_poly_Z(d1*d2, coefficients_for_H);
			libff::leave_block("Compute ZK-patch");

			libff::enter_block("Compute evaluation of polynomial A on set T");
			domain->cosetFFT(aA, FieldT::multiplicative_generator);
			libff::leave_block("Compute evaluation of polynomial A on set T");

			libff::enter_block("Compute evaluation of polynomial B on set T");
			domain->cosetFFT(aB, FieldT::multiplicative_generator);
			libff::leave_block("Compute evaluation of polynomial B on set T");

			libff::enter_block("Compute evaluation of polynomial H on set T");
			std::vector<FieldT> &H_tmp = aA; // can overwrite aA because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
			for (size_t i = 0; i < domain->m; ++i)
			{
				H_tmp[i] = aA[i]*aB[i];
			}
			std::vector<FieldT>().swap(aB); // destroy aB

			libff::enter_block("Compute evaluation of polynomial C on set S");
			std::vector<FieldT> aC(domain->m, FieldT::zero());
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				aC[i] += cs.constraints[i].c.evaluate(full_variable_assignment);
			}

			conv_index = 0;
			for( size_t i=0;i<cs.num_constraints();i++){
				FieldT acc = FieldT::zero();
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					for(auto &lt : cs.constraints[i].c2){
						acc += full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff);
					}
					//if(acc != FieldT::zero()){
						aC[conv_index*cs.num_convol_outputs() + j + cs.num_constraints() + cs.num_inputs()] += acc;
						conv_index++;
					//}
				}
			}
			libff::leave_block("Compute evaluation of polynomial C on set S");
			libff::leave_block("Compute evaluation of polynomial C on set S");

			libff::enter_block("Compute coefficients of polynomial C");
			domain->iFFT(aC);
			libff::leave_block("Compute coefficients of polynomial C");

			libff::enter_block("Compute evaluation of polynomial C on set T");
			domain->cosetFFT(aC, FieldT::multiplicative_generator);
			libff::leave_block("Compute evaluation of polynomial C on set T");

#ifdef MULTICORE
#pragma omp parallel for
#endif
			for (size_t i = 0; i < domain->m; ++i)
			{
				H_tmp[i] = (H_tmp[i]-aC[i]);
			}

			libff::enter_block("Divide by Z on set T");
			domain->divide_by_Z_on_coset(H_tmp);
			libff::leave_block("Divide by Z on set T");

			libff::leave_block("Compute evaluation of polynomial H on set T");

			libff::enter_block("Compute coefficients of polynomial H");
			domain->icosetFFT(H_tmp, FieldT::multiplicative_generator);
			libff::leave_block("Compute coefficients of polynomial H");

			libff::enter_block("Compute sum of H and ZK-patch");
#ifdef MULTICORE
#pragma omp parallel for
#endif
			for (size_t i = 0; i < domain->m; ++i)
			{
				coefficients_for_H[i] += H_tmp[i];
			}
			libff::leave_block("Compute sum of H and ZK-patch");

			libff::leave_block("Call to mr1cs_to_qap_witness_map");

			return qap_witness<FieldT>(cs.num_variables(),
					domain->m,
					cs.num_inputs(),
					d1,
					d2,
					d3,
					full_variable_assignment,
					std::move(coefficients_for_H));
		}


	template<typename FieldT>
		void print_Poly(const vector<FieldT> c){
			for (size_t i = 0; i < c.size(); i++)
			{
				unsigned long coefficient = c[i].as_ulong();

				if (i == 0) std::cout << coefficient << " + ";
				else if (i < (c.size()-1)) std::cout << coefficient << "x^" << i << " + ";
				else std::cout << coefficient << "x^" << i << std::endl;
			}

		}

	template<typename FieldT>
		void print_point(const vector<FieldT> f){
			for (size_t i = 0; i < f.size(); i++)
			{
				printf("%ld: %ld\n", i, f[i].as_ulong());
			}
		}

	/*
	template<typename FieldT>
		vector<FieldT> mvPolycalculate(vector<vector<FieldT>> mvPoly, 
				size_t num_constraints, size_t num_z_order,
				shared_ptr<libfqfft::evaluation_domain<FieldT>> domain
				)
		{
			vector<vector<FieldT>> tempPoly(num_constraints, vector<FieldT>(num_constraints, FieldT::zero()));

			bool isEmpty = true;
			for(size_t i=0; i<num_constraints;i++){
				size_t j;
				for(j=0; j<num_z_order+1;j++){
					if(zxPoly[i][j] != FieldT::zero()){
						isEmpty = false;
						break;
					}
				}
				if(j < num_z_order+1) tempPoly[i][i] = FieldT::one();
			}

			if(isEmpty) return vector<FieldT>(1, FieldT::zero());

			for(size_t i=0; i<num_constraints;i++){
				domain->iFFT(tempPoly[i]);
			}

			//v0(x) * z(x^(4d+1))
			vector<FieldT> zExtendPoly((num_z_order*((num_constraints-1)*4+1))+1, FieldT::zero());
			for(size_t i=0;i<num_z_order+1;i++){
				zExtendPoly[i*((num_constraints-1)*4 + 1)] = zxPoly[0][i];
			}

			vector<FieldT> mulResult(1, FieldT::zero());
			libfqfft::_polynomial_multiplication(mulResult, tempPoly[0], zExtendPoly);

			//point add
			for(size_t i=0; i<num_constraints;i++){
				FieldT temp = FieldT::zero();
				for(size_t j=1; j<num_constraints;j++){
					temp+=tempPoly[j][i];
				}
				mulResult[i] += temp;
			}

			return mulResult;
		}
	template<typename FieldT>
		qap_instance_evaluation<FieldT> mr1cs_to_qap_instance_map_with_evaluation(
				const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain,
				const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2,
				const vector<vector<FieldT>> &V,
				const vector<vector<FieldT>> &W,
				const vector<vector<FieldT>> &Y,
				const size_t num_variables,
				const size_t num_inputs,
				const size_t num_constraints,
				const size_t num_z_order,
				const FieldT &t
				)
		{
			libff::enter_block("Call to mr1cs_to_qap_instance_map_with_evaluation");


			const FieldT Zt = domain->compute_vanishing_polynomial(t);


			/*
			// check mr1cs
			FieldT ckVk = FieldT::zero();
			FieldT ckWk = FieldT::zero();
			FieldT ckYk = FieldT::zero();
			FieldT omega = libff::get_root_of_unity<FieldT>(4);
			//omega = omega^16;
			for(size_t i=0; i<num_variables; i++){

			FieldT temp = FieldT::zero();
			for(size_t j=0;j<V[i].size();j++){
			temp += V[i][j] * (omega^j);
			}
			cout<<temp.as_ulong()<< " * " << wires[i].as_ulong() << " = "<< (temp*wires[i]).as_ulong() <<"\t";
			ckVk += temp;

			temp = FieldT::zero();
			for(size_t j=0;j<W[i].size();j++){
			temp += W[i][j] * (omega^j);
			}
			cout<<temp.as_ulong()<< " * " << wires[i].as_ulong() << " = "<< (temp*wires[i]).as_ulong() <<"\t";
			ckWk +=temp;

			temp = FieldT::zero();
			for(size_t j=0;j<Y[i].size();j++){
			temp += Y[i][j] * (omega^j);
			}
			cout<<temp.as_ulong()<< " * " << wires[i].as_ulong() << " = "<< (temp*wires[i]).as_ulong() <<endl;
			ckYk +=temp;
			}
			if(ckVk * ckWk == ckYk) cout<<"ok"<<endl;
			else cout<<"no"<<endl;
			 

			vector<FieldT> At(domain2->m, FieldT::zero());
			vector<FieldT> Bt(domain2->m, FieldT::zero());
			vector<FieldT> Ct(domain2->m, FieldT::zero());
			vector<FieldT> Ht;
			Ht.reserve(domain2->m);
			FieldT ti = FieldT::one();
			for(size_t i=0; i<domain2->m;i++){
				Ht.emplace_back(ti);
				ti*= t;
			}

			//At, Bt, Ct = poly(t);
			for(size_t i=0; i<num_variables; i++){
				for(size_t j=0;j<V[i].size();j++){
					At[i] += V[i][j] * Ht[j];
				}

				for(size_t j=0;j<W[i].size();j++){
					Bt[i] += W[i][j] * Ht[j];
				}

				for(size_t j=0;j<Y[i].size();j++){
					Ct[i] += Y[i][j] * Ht[j];
				}
			}
			libff::leave_block("Call to mr1cs_to_qap_instance_map_with_evaluation");

			return qap_instance_evaluation<FieldT>(domain2,
					num_variables,
					domain2->m,
					num_inputs,
					t,
					std::move(At),
					std::move(Bt),
					std::move(Ct),
					std::move(Ht),
					Zt);



		}

	template<typename FieldT>
		qap_witness<FieldT> mr1cs_to_qap_witness_map(
				const shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain2,
				vector<vector<FieldT>> &V,
				vector<vector<FieldT>> &W,
				vector<vector<FieldT>> &Y,
				vector<FieldT> &Ht,
				const FieldT Zt,
				const size_t num_variables,
				const size_t num_inputs,
				const size_t num_constraints,
				const size_t num_z_order,
				const vector<FieldT> &wires
				)
		{
			const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(num_constraints);
			const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>(domain2->m*2);

			///qap witness start

			vector<FieldT> aA(domain2->m, FieldT::zero());
			vector<FieldT> aB(domain2->m, FieldT::zero());
			vector<FieldT> aC(domain2->m, FieldT::zero());

			//V, W, Y poly -> point
			for(size_t i=0; i<num_variables; i++){
				V[i].resize(domain2->m, FieldT::zero());
				domain2->FFT(V[i]);
				W[i].resize(domain2->m, FieldT::zero());
				domain2->FFT(W[i]);
				Y[i].resize(domain2->m, FieldT::zero());
				domain2->FFT(Y[i]);
			}

			// mr1cs point * wire
			for(size_t i=0; i<domain2->m;i++){
				for(size_t j=0;j< num_variables;j++){
					aA[i] += V[j][i] * wires[j];
					aB[i] += W[j][i] * wires[j];
					aC[i] += Y[j][i] * wires[j];
				}
			}

			FieldT resA3 = FieldT::zero();
			FieldT resB3 = FieldT::zero();
			FieldT resC3 = FieldT::zero();
			// res = poly(t) using lagrange
			/*
			   for(size_t i=0;i<domain2->m;i++){
			   resA3 += aA[i] * u[i];
			   resB3 += aB[i] * u[i];
			   resC3 += aC[i] * u[i];
			   }
			 


			//A, B, C  point -> poly
			domain2->iFFT(aA);
			domain2->iFFT(aB);
			domain2->iFFT(aC);

			aA.resize(domain4->m, FieldT::zero());
			aB.resize(domain4->m, FieldT::zero());
			aC.resize(domain4->m, FieldT::zero());

			//A(x) => A(gx) same B, C
			domain4->cosetFFT(aA, FieldT::multiplicative_generator);
			domain4->cosetFFT(aB, FieldT::multiplicative_generator);
			domain4->cosetFFT(aC, FieldT::multiplicative_generator);

			vector<FieldT> tmpH(domain4->m, FieldT::zero());
			for(size_t i=0;i<domain4->m;i++){
				tmpH[i] = aA[i] * aB[i] - aC[i];
			}

			domain4->icosetFFT(aA, FieldT::multiplicative_generator);
			domain4->icosetFFT(aB, FieldT::multiplicative_generator);
			domain4->icosetFFT(aC, FieldT::multiplicative_generator);


			vector<FieldT> Zpoly(1, FieldT::one());
			vector<FieldT> xw(2, FieldT::one());
			for(size_t i=0;i<domain->m;i++){
				xw[0] = -domain->get_domain_element(i);
				libfqfft::_polynomial_multiplication(Zpoly, Zpoly, xw);
			}
			//Z(gx) points
			Zpoly.resize(domain4->m, FieldT::zero());
			domain4->cosetFFT(Zpoly, FieldT::multiplicative_generator);

			//(A*B-C(gx) point) / (Z(gx) point)
			for(size_t i=0;i<domain4->m;i++){
				tmpH[i] = tmpH[i] * Zpoly[i].inverse();
			}
			domain4->icosetFFT(tmpH, FieldT::multiplicative_generator);

			domain4->icosetFFT(Zpoly, FieldT::multiplicative_generator);

			//h(gx) -> h(x)
			FieldT resH = FieldT::zero();
			for(size_t i=0;i<domain4->m;i++){
				resH += tmpH[i] * Ht[i];
			}
			cout<<"cal resH = "<<resH.as_ulong()<<endl;

			FieldT resA2 = FieldT::zero();
			FieldT resB2 = FieldT::zero();
			FieldT resC2 = FieldT::zero();
			// res = poly(t)
			for(size_t i=0;i<domain2->m;i++){
				resA2 += aA[i] * Ht[i];
				resB2 += aB[i] * Ht[i];
				resC2 += aC[i] * Ht[i];
			}
			aA.resize(domain2->m);
			aB.resize(domain2->m);
			aC.resize(domain2->m);
			domain2->cosetFFT(aA, FieldT::multiplicative_generator);
			domain2->cosetFFT(aB, FieldT::multiplicative_generator);
			domain2->cosetFFT(aC, FieldT::multiplicative_generator);


			/*
			   FieldT resA = FieldT::zero();
			   FieldT resB = FieldT::zero();
			   FieldT resC = FieldT::zero();
			   for(size_t i=0;i< num_variables;i++){
			   resA += At[i] * wires[i];
			   resB += Bt[i] * wires[i];
			   resC += Ct[i] * wires[i];
			   }
			 


			if(((resA2 * resB2) - resC2) == (Zt * resH)) cout<< "ok2"<<endl;
			else cout <<"no2"<<endl;

			/*
			   cout<<"A = "<<resA.as_ulong() <<", A2 = "<<resA2.as_ulong()<< ", A3 = "<< resA3.as_ulong()<<endl;
			   cout << (resA == resA2) << endl;
			   cout<<"B = "<<resB.as_ulong() <<", B2 = "<<resB2.as_ulong()<< ", B3 = "<<resB3.as_ulong()<<endl;
			   cout << (resA == resA2) << endl;
			   cout<<"C = "<<resC.as_ulong() <<", C2 = "<<resC2.as_ulong()<< ", C3 = "<<resC3.as_ulong()<<endl;
			   cout << (resA == resA2) << endl;

			   cout<<"resA * resB - resC / Zt = "<< ((resA*resB - resC)*Zt.inverse()).as_ulong()<<endl;
			   cout<<"H = "<<resH.as_ulong()<<endl;
			   cout<<"Z = "<<Zt.as_ulong()<<endl;
			   cout<<"Z^-1 = "<<Zt.inverse().as_ulong()<<endl;
			 
			return qap_witness<FieldT>(num_variables,
					domain2->m,
					num_inputs,
					FieldT::zero(),
					FieldT::zero(),
					FieldT::zero(),
					wires,
					move(tmpH)
					);
		}

	template<typename FieldT>
		void evaluate_mqap(){

			size_t num_variables = 15;
			size_t num_inputs = 8;
			size_t num_constraints = 3;
			size_t num_z_order = 6;



			/*
			   1	a0	1	1	0		0	0	0		0	0	0
			   2	a1	s	0	0		0	0	0		0	0	0
			   3	a2	s^2	0	1		0	0	0		0	0	0
			   1	x0	0	0	0		1	1	0		0	0	0
			   2	x1	0	0	0		s	0	0		0	0	0
			   3	x2	0	0	0		s^2	0	0		0	0	0
			   4	x3	0	0	0		s^3	0	0		0	0	0
			   5	x4	0	0	0		s^4	0	1		0	0	0
			   1	y0	0	0	0		0	0	0		1	1	0
			   4	y1	0	0	0		0	0	0		s	0	0
			   10	y2	0	0	0		0	0	0		s^2	0	0
			   16	y3	0	0	0		0	0	0		s^3	0	0
			   22	y4	0	0	0		0	0	0		s^4	0	0
			   22	y5	0	0	0		0	0	0		s^5	0	0
			   15	y6	0	0	0		0	0	0		s^6	0	1

			 

			const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(num_constraints);
			const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2 = libfqfft::get_evaluation_domain<FieldT>(num_z_order*((num_constraints-1)*4+1) + (num_constraints-1) +1);
			const shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>(domain2->m*2);

			vector<vector<vector<FieldT>>> mr1csV (num_variables, vector<vector<FieldT>>(num_constraints, vector<FieldT>(num_z_order+1, FieldT::zero())));
			vector<vector<vector<FieldT>>> mr1csW (num_variables, vector<vector<FieldT>>(num_constraints, vector<FieldT>(num_z_order+1, FieldT::zero())));
			vector<vector<vector<FieldT>>> mr1csY (num_variables, vector<vector<FieldT>>(num_constraints, vector<FieldT>(num_z_order+1, FieldT::zero())));

			vector<FieldT> wires={
				FieldT(1), FieldT(2), FieldT(3), //a
				FieldT(1), FieldT(2), FieldT(3), FieldT(4), FieldT(5), //x
				FieldT(1), FieldT(4), FieldT(10), FieldT(16), FieldT(22), FieldT(22), FieldT(15) //y
			};

			//V points
			mr1csV[0][0][0] = FieldT::one();
			mr1csV[0][1][0] = FieldT::one();

			mr1csV[1][0][1] = FieldT::one();

			mr1csV[2][2][0] = FieldT::one();
			mr1csV[2][0][2] = FieldT::one();

			//W points
			mr1csW[3][0][0] = FieldT::one();
			mr1csW[3][1][0] = FieldT::one();

			mr1csW[4][0][1] = FieldT::one();
			mr1csW[5][0][2] = FieldT::one();
			mr1csW[6][0][3] = FieldT::one();

			mr1csW[7][2][0] = FieldT::one();
			mr1csW[7][0][4] = FieldT::one();

			//Y points
			mr1csY[8][0][0] = FieldT::one();
			mr1csY[8][1][0] = FieldT::one();

			mr1csY[9][0][1] = FieldT::one();
			mr1csY[10][0][2] = FieldT::one();
			mr1csY[11][0][3] = FieldT::one();
			mr1csY[12][0][4] = FieldT::one();
			mr1csY[13][0][5] = FieldT::one();

			mr1csY[14][0][6] = FieldT::one();
			mr1csY[14][2][0] = FieldT::one();

			libff::enter_block("Compute V, W, Y");
			cout<<"===========Calculate V==========="<<endl;
			vector<vector<FieldT>> V(num_variables, vector<FieldT>(1, FieldT::zero()));
			for(size_t i=0; i<num_variables; i++){
				V[i] = zxPolycalculate(mr1csV[i], num_constraints, num_z_order, domain);	
			}


			cout<<"===========Calculate W==========="<<endl;
			vector<vector<FieldT>> W(num_variables, vector<FieldT>(1, FieldT::zero()));
			for(size_t i=0; i<num_variables; i++){
				W[i] = zxPolycalculate(mr1csW[i], num_constraints, num_z_order, domain);	
			}

			cout<<"===========Calculate Y==========="<<endl;
			vector<vector<FieldT>> Y(num_variables, vector<FieldT>(1, FieldT::zero()));
			for(size_t i=0; i<num_variables; i++){
				Y[i] = zxPolycalculate(mr1csY[i], num_constraints, num_z_order, domain);	
			}
			libff::leave_block("Compute V, W, Y");

			FieldT t = FieldT::random_element();
			//qap_instance_evaluation<FieldT> qap_inst =  mr1cs_to_qap_instance_map_with_evaluation(domain, domain2, V, W, Y, num_variables, num_inputs, num_constraints, num_z_order, t);
			//qap_witness<FieldT> qap_wit =  mr1cs_to_qap_witness_map(qap_inst.domain, V, W, Y, qap_inst.Ht, qap_inst.Zt, num_variables, num_inputs, num_constraints, num_z_order, wires);
		}

	template<typename ppT>
		void testPoly(){
			evaluate_mqap<libff::Fr<ppT>>();
		}

	*/

}

