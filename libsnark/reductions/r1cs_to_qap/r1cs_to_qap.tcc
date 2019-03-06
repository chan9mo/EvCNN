/** @file
 *****************************************************************************

 Implementation of interfaces for a R1CS-to-QAP reduction.

 See r1cs_to_qap.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef R1CS_TO_QAP_TCC_
#define R1CS_TO_QAP_TCC_

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>

namespace libsnark {

	/**
	 * Instance map for the R1CS-to-QAP reduction.
	 *
	 * Namely, given a R1CS constraint system cs, construct a QAP instance for which:
	 *   A := (A_0(z),A_1(z),...,A_m(z))
	 *   B := (B_0(z),B_1(z),...,B_m(z))
	 *   C := (C_0(z),C_1(z),...,C_m(z))
	 * where
	 *   m = number of variables of the QAP
	 * and
	 *   each A_i,B_i,C_i is expressed in the Lagrange basis.
	 */
	template<typename FieldT>
		 qap_instance<FieldT> r1cs_to_qap_instance_map(const r1cs_constraint_system<FieldT> &cs)
	{
		libff::enter_block("Call to r1cs_to_qap_instance_map");

		//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()));

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

		size_t conv_index_a = 0;
		size_t conv_index_b = 0;
		size_t conv_index_c = 0;

		
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

			for(size_t j=0;j<cs.constraints[i].a2.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_a ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					A_in_Lagrange_basis[cs.constraints[i].a2.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].a2.terms[j].coeff.as_ulong());
				}
			}
			for(size_t j=0;j<cs.constraints[i].b2.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_b ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					B_in_Lagrange_basis[cs.constraints[i].b2.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].b2.terms[j].coeff.as_ulong());
				}
			}

			for(size_t j=0;j<cs.constraints[i].c2.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_c ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					C_in_Lagrange_basis[cs.constraints[i].c2.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
				}
			}
		}
		libff::leave_block("Compute polynomials A, B, C in Lagrange basis");

		libff::leave_block("Call to r1cs_to_qap_instance_map");

		return qap_instance<FieldT>(domain,
				cs.num_variables(),
				domain->m,
				cs.num_inputs(),
				std::move(A_in_Lagrange_basis),
				std::move(B_in_Lagrange_basis),
				std::move(C_in_Lagrange_basis));
	}

	/**
	 * Instance map for the R1CS-to-QAP reduction followed by evaluation of the resulting QAP instance.
	 *
	 * Namely, given a R1CS constraint system cs and a field element t, construct
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
		qap_instance_evaluation<FieldT> r1cs_to_qap_instance_map_with_evaluation(const r1cs_constraint_system<FieldT> &cs,
				const FieldT &t)
		{



			libff::enter_block("Call to r1cs_to_qap_instance_map_with_evaluation");

			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()+1));

			std::cout<<"#constraint : "<< cs.num_constraints()<<"\t#input : "<< cs.num_inputs()<<"\t#convolOut : "<<(cs.num_convol_outputs())<<"\t#convol : "<<(cs.num_convol())<<std::endl;

			std::vector<FieldT> At, Bt, Ct, Ht;

			At.resize(cs.num_variables()+1, FieldT::zero());
			Bt.resize(cs.num_variables()+1, FieldT::zero());
			Ct.resize(cs.num_variables()+1, FieldT::zero());
			Ht.reserve(domain->m+1);

			const FieldT Zt = domain->compute_vanishing_polynomial(t);

			libff::enter_block("Compute evaluations of A, B, C, H at t");
			const std::vector<FieldT> u = domain->evaluate_all_lagrange_polynomials(t);
			/**
			 * add and process the constraints
			 *     input_i * 0 = 0
			 * to ensure soundness of input consistency
			 */
			/*
			for(size_t j=0;j<cs.num_variables();j++){
				
				for(size_t i=0;i<cs.num_constraints();i++){
					if(i==0) std::cout<<j<<"\t";
					std::cout<<cs.constraints[i].a.terms[j].coeff.as_ulong()<<"\t";
				}
				std::cout<<std::endl;
				
			}
			*/
			//std::cout<<"At ";
			for (size_t i = 0; i <= cs.num_inputs(); ++i)
			{
				At[i] = u[cs.num_constraints() + i];
				//std::cout<<i<<", "<<cs.num_constraints()+i<<"\t";
			}
			//std::cout<<std::endl;
			/* process all other constraints */

			size_t conv_index_a = 0;
			size_t conv_index_b = 0;
			size_t conv_index_c = 0;

			
			
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				std::cout<<"At\n";
				for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j)
				{
					At[cs.constraints[i].a.terms[j].index] +=
						u[i]*cs.constraints[i].a.terms[j].coeff;
					std::cout<<"( "<<i<<", "<<cs.constraints[i].a.terms[j].index<<")="<<cs.constraints[i].a.terms[j].coeff.as_ulong()<<"\n";
				}
				//std::cout<<std::endl;
				std::cout<<"bt\n";
				for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j)
				{
					Bt[cs.constraints[i].b.terms[j].index] +=
						u[i]*cs.constraints[i].b.terms[j].coeff;
					std::cout<<"( "<<i<<", "<<cs.constraints[i].b.terms[j].index<<")="<<cs.constraints[i].b.terms[j].coeff.as_ulong()<<"\n";

				}
				std::cout<<"Ct\n";
				for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j)
				{
					Ct[cs.constraints[i].c.terms[j].index] +=
						u[i]*cs.constraints[i].c.terms[j].coeff;					
					std::cout<<"( "<<i<<", "<<cs.constraints[i].c.terms[j].index<<")="<<cs.constraints[i].c.terms[j].coeff.as_ulong()<<"\n";

				}

				bool flag = true;
				for(size_t k=0;k<cs.num_convol_outputs();k++){

					///TODO change (j+1) to safe value like prime
					FieldT temp = FieldT::one();
					for(size_t j=0;j<cs.constraints[i].a2.terms.size() != 0;j++){
						if(cs.constraints[i].a2.terms[j].index){
						
						if(flag){
							conv_index_a ++;
							flag = false;
						}
						std::cout<<"( "<<i<<", "<<cs.constraints[i].a2.terms[j].index<<")="<<cs.constraints[i].a2.terms[j].coeff.as_ulong()<<"\n";

					//for(size_t k=0;k<cs.num_convol_outputs();k++){
						//temp = (FieldT(k+1)^(cs.constraints[i].a2.terms[j].coeff.as_ulong()-1));
						At[cs.constraints[i].a2.terms[j].index] +=
							u[(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_a-1) + k ] * temp;//(FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						//	u[(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_a-1)+k] * (FieldT(k+1)^cs.constraints[i].a2.terms[j].coeff.as_ulong());
						
						std::cout<<"#in : "<<cs.num_inputs()<<" #cons : "<<cs.num_constraints()<<" #convOut : "<<cs.num_convol_outputs()<<" convIndex :"<<conv_index_a<<" k :"<<k<<std::endl;
						std::cout<<"At["<<(cs.constraints[i].a2.terms[j].index)<<
							"] = u["<<(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_a-1) + k <<
							"]*("<<k+1<<")^"<<cs.constraints[i].a2.terms[j].coeff.as_ulong()<<
							"= "<<temp.as_ulong()<<"\n";
						
						temp *= FieldT(k+1);
						}

						
					}
				}
				flag = true;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					///TODO change (j+1) to safe value like prime
					FieldT temp = FieldT::one();
					
					for(size_t j=0;j<cs.constraints[i].b2.terms.size();j++){
						if(cs.constraints[i].b2.terms[j].index){
						if(flag){
							conv_index_b ++;
							flag = false;
						}
						std::cout<<"( "<<i<<", "<<cs.constraints[i].b2.terms[j].index<<")="<<cs.constraints[i].b2.terms[j].coeff.as_ulong()<<"\n";

						//temp = (FieldT(k+1)^(cs.constraints[i].b2.terms[j].coeff.as_ulong()-1));
						Bt[cs.constraints[i].b2.terms[j].index] +=
							u[(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_b-1) + k ] * temp;//(FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						 
						//	u[(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_b-1)+k] * (FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						
						std::cout<<"Bt["<<(cs.constraints[i].b2.terms[j].index)<<
							"] = u["<<(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_b-1)+k<<
								"]*("<<k+1<<")^"<<cs.constraints[i].b2.terms[j].coeff.as_ulong()<<
								"= "<<temp.as_ulong()<<"\n";

						temp *= FieldT(k+1);
						}
					}
				}
				flag = true;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					///TODO change (j+1) to safe value like prime
					FieldT temp = FieldT::one();

					for(size_t j=0;j<cs.constraints[i].c2.terms.size();j++){
						if(cs.constraints[i].c2.terms[j].index){

						if(flag){
							conv_index_c ++;
							flag = false;
						} 
						std::cout<<"( "<<i<<", "<<cs.constraints[i].c2.terms[j].index<<")="<<cs.constraints[i].c2.terms[j].coeff.as_ulong()<<"\n";
						//temp = (FieldT(k+1)^(cs.constraints[i].c2.terms[j].coeff.as_ulong()-1));
						Ct[cs.constraints[i].c2.terms[j].index] +=
							u[(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_c-1) + k ] * temp;//(FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						
						//	u[(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_c-1)+k] * (FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());

						
						std::cout<<"Ct["<<(cs.constraints[i].c2.terms[j].index)<<
							"] = u["<<(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_c-1)+k<<
							"]*("<<j+1<<")^"<<cs.constraints[i].c2.terms[j].coeff.as_ulong()<<
							"= "<<temp.as_ulong()<<"\n";

						temp *= FieldT(k+1);
						}
					}
				}
				
			}

			FieldT ti = FieldT::one();
			for (size_t i = 0; i < domain->m+1; ++i)
			{
				Ht.emplace_back(ti);
				ti *= t;
			}
			libff::leave_block("Compute evaluations of A, B, C, H at t");

			libff::leave_block("Call to r1cs_to_qap_instance_map_with_evaluation");

			return qap_instance_evaluation<FieldT>(domain,
					cs.num_variables(),
					domain->m,
					cs.num_inputs(),
					t,
					std::move(At),
					std::move(Bt),
					std::move(Ct),
					std::move(Ht),
					Zt);
		}

	/**
	 * Witness map for the R1CS-to-QAP reduction.
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
		qap_witness<FieldT> r1cs_to_qap_witness_map(const r1cs_constraint_system<FieldT> &cs,
				const r1cs_primary_input<FieldT> &primary_input,
				const r1cs_auxiliary_input<FieldT> &auxiliary_input,
				const FieldT &d1,
				const FieldT &d2,
				const FieldT &d3)
		{
			libff::enter_block("Call to r1cs_to_qap_witness_map");

			/* sanity check */
			assert(cs.is_satisfied(primary_input, auxiliary_input));

			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()+1));
			std::cout<<"#constraint : "<< cs.num_constraints()<<"\t#input : "<< cs.num_inputs()<<"\t#convolOut : "<<(cs.num_convol_outputs())<<"\t#convol : "<<(cs.num_convol())<<std::endl;
			r1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
			full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());

			std::cout<<"variables : ";
			for(size_t i=0;i<full_variable_assignment.size();i++){
				std::cout<<full_variable_assignment[i].as_ulong()<<"\t";
			}
			std::cout<<"\n";
			std::cout<<"degree : "<<domain->m<<" sum = "<<cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()+1)<<"\n";

			libff::enter_block("Compute evaluation of polynomials A, B on set S");
			std::vector<FieldT> aA(domain->m, FieldT::zero()), aB(domain->m, FieldT::zero());

			/* account for the additional constraints input_i * 0 = 0 */
			for (size_t i = 0; i <= cs.num_inputs(); ++i)
			{
				aA[i+cs.num_constraints()] = (i > 0 ? full_variable_assignment[i-1] : FieldT::one());
				//std::cout<<"aA["<<i+cs.num_constraints()<<"] = "<<aA[i+cs.num_constraints()].as_ulong()<<std::endl;

			}
			/* account for all other constraints */
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				aA[i] += cs.constraints[i].a.evaluate(full_variable_assignment);
				std::cout<<"aA["<<i<<"] = "<<aA[i].as_ulong()<<std::endl;
				aB[i] += cs.constraints[i].b.evaluate(full_variable_assignment);
				std::cout<<"aB["<<i<<"] = "<<aB[i].as_ulong()<<std::endl;

			}


			size_t conv_index = 0;
			size_t conv_index2 = 0;
			for( size_t i=0;i<cs.num_constraints();i++){
				FieldT acc = FieldT::zero();
				bool flag = false;
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					
					FieldT temp = FieldT::one();
					for(auto &lt : cs.constraints[i].a2){
						if(lt.index){
						if(!flag){ 
							conv_index++;
							flag = true;
						}
						//temp = (FieldT(j)^(lt.coeff.as_ulong()-1));
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						//std::cout<<j<<"^"<<lt.coeff.as_ulong()<<" = "<<temp.as_ulong()<<std::endl;

						temp *= FieldT(j);
						//acc = full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff.as_ulong());

						std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						std::cout<<"j : "<<j<<" convIdx :"<<conv_index<<" #convOut : "<<cs.num_convol_outputs()<<" #cons : "<<cs.num_constraints()<<+" #in :"<<cs.num_inputs()<<std::endl;
						aA[(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						std::cout<<"aA["<<(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()<<"] = acc\n";
						}
					}
				}
				/*
				for(auto &lt : cs.constraints[i].a2){
					if(!flag){ 
						conv_index++;
						flag = true;
					}
					FieldT temp = FieldT::one();
					for(size_t j=1;j<=cs.num_convol_outputs();j++){
						//if(lt.idex ==0) acc += (FieldT::one()*cs.conv_num_output()) * lt.coeff; //there is no first one in conv constraints
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						temp *= FieldT(j);
						//Fixing temp value
						std::cout<<j<<"^"<<lt.coeff.as_ulong()<<" = "<<temp.as_ulong()<<std::endl;
						//acc = full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff.as_ulong());

						std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						aA[(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						std::cout<<"aA["<<(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()<<"] = acc\n";
					}
					//if(acc != FieldT::zero()){
					//}
					*/
				
				/*
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					std::cout<<"aA["<<(j)+cs.num_constraints()+cs.num_inputs()<<"] ="<<
						aA[(j) + cs.num_constraints() + cs.num_inputs()].as_ulong() <<"\n";
				}
				*/
			
				acc = FieldT::zero();
				flag = false;
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					
					FieldT temp = FieldT::one();
					for(auto &lt : cs.constraints[i].b2){
						if(lt.index){
						if(!flag){ 
							conv_index2++;
							flag = true;
						}
						//temp = (FieldT(j)^(lt.coeff.as_ulong()-1));
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						std::cout<<j<<"^"<<lt.coeff.as_ulong()<<" = "<<temp.as_ulong()<<std::endl;

						temp *= FieldT(j);
						//acc = full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff.as_ulong());

						std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						aB[(j)+(conv_index2-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						//std::cout<<"j : "<<j<<" convind2-1 : "<<conv_index2-1<<" convout : "<<cs.num_convol_outputs()<<" cons+in : "<<cs.num_constraints() + cs.num_inputs()<<std::endl;
						std::cout<<"aB["<<(j)+(conv_index2-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()<<"] = acc\n";
						}
					}
				}
				/*
				for(auto &lt : cs.constraints[i].b2){
					if(!flag){ 
						conv_index2++;
						flag = true;
					}
					FieldT temp = FieldT::one();
					for(size_t j=1;j<=cs.num_convol_outputs();j++){
						//acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						//temp *= FieldT(j);
						acc = full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff.as_ulong());

						std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						aB[(j)+(conv_index2-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						std::cout<<"aB["<<(j)+(conv_index2-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()<<"] = acc\n";
					}
					//aB[i*cs.num_convol() + j + cs.num_constraints() + cs.num_inputs()] += acc;
				}
				/*
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					std::cout<<"aB["<< +(j)+cs.num_constraints()+cs.num_inputs()<<"] ="<<
						aB[(j) + cs.num_constraints() + cs.num_inputs()].as_ulong() <<"\n";
				}
				*/
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
				std::cout<<"aC["<<i<<"] = "<<aC[i].as_ulong()<<std::endl;

			}

			conv_index = 0;
			for( size_t i=0;i<cs.num_constraints();i++){
				FieldT acc = FieldT::zero();
				bool flag = false;
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					
					FieldT temp = FieldT::one();
					for(auto &lt : cs.constraints[i].c2){
						if(lt.index){
						if(!flag){ 
							conv_index++;
							flag = true;
						}
						//temp = (FieldT(j)^(lt.coeff.as_ulong()-1));
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						///TODO lt.index break when output is not used to sequentially   
						//acc = full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff.as_ulong());
						//std::cout<<j<<"^"<<lt.coeff.as_ulong()<<" = "<<temp.as_ulong()<<std::endl;

						temp *= FieldT(j);
						aC[(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						std::cout<<"aC["<<(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()<<"] = acc\n";
						}
					}
				}
				/*
				for(auto &lt : cs.constraints[i].c2){
					if(!flag){
						conv_index++;
						flag = true;
					}
					FieldT temp = FieldT::one();
					for(size_t j=1;j<=cs.num_convol_outputs();j++){
						//acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						//temp *= FieldT(j);
						acc = full_variable_assignment[lt.index-1] * (FieldT(j)^lt.coeff.as_ulong());

						aC[(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						std::cout<<"aC["<<(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()<<"] = acc\n";
					}
				}
				/*
				//for debug
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					std::cout<<"aC["<<(j)+cs.num_constraints()+cs.num_inputs()<<"] ="<<
						aC[(j) + cs.num_constraints() + cs.num_inputs()].as_ulong() <<"\n";
				}
				*/

			}
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

			libff::leave_block("Call to r1cs_to_qap_witness_map");

			return qap_witness<FieldT>(cs.num_variables(),
					domain->m,
					cs.num_inputs(),
					d1,
					d2,
					d3,
					full_variable_assignment,
					std::move(coefficients_for_H));
		}



	/////NEW SYSTEM *******/////////

	/**
	 * Instance map for the R1CS-to-QAP reduction.
	 *
	 * Namely, given a R1CS constraint system cs, construct a QAP instance for which:
	 *   A := (A_0(z),A_1(z),...,A_m(z))
	 *   B := (B_0(z),B_1(z),...,B_m(z))
	 *   C := (C_0(z),C_1(z),...,C_m(z))
	 * where
	 *   m = number of variables of the QAP
	 * and
	 *   each A_i,B_i,C_i is expressed in the Lagrange basis.
	 */
	/*
	template<typename FieldT>
		 qap_instance<FieldT> r1cs_to_qap_instance_map(const r1cs_constraint_convol_system<FieldT> &cs)
	{
		libff::enter_block("Call to r1cs_convol_to_qap_instance_map");

		//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()));

		std::vector<std::map<size_t, FieldT> > A_in_Lagrange_basis(cs.num_variables()+1);
		std::vector<std::map<size_t, FieldT> > B_in_Lagrange_basis(cs.num_variables()+1);
		std::vector<std::map<size_t, FieldT> > C_in_Lagrange_basis(cs.num_variables()+1);

		libff::enter_block("Compute polynomials A, B, C in Lagrange basis");
		/**
		 * add and process the constraints
		 *     input_i * 0 = 0
		 // to ensure soundness of input consistency
		 
		for (size_t i = 0; i <= cs.num_inputs(); ++i)
		{
			A_in_Lagrange_basis[i][cs.num_constraints() + i] = FieldT::one();
		}


		size_t conv_index_a = 0;
		size_t conv_index_b = 0;
		size_t conv_index_c = 0;

		// process all other constraints 
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


			for(size_t j=0;j<cs.constraints[i].a2.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_a ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					A_in_Lagrange_basis[cs.constraints[i].a2.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].a2.terms[j].coeff.as_ulong());
				}
			}
			for(size_t j=0;j<cs.constraints[i].b2.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_b ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					B_in_Lagrange_basis[cs.constraints[i].b2.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].b2.terms[j].coeff.as_ulong());
				}
			}

			for(size_t j=0;j<cs.constraints[i].c2.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_c ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					C_in_Lagrange_basis[cs.constraints[i].c2.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
				}
			}
		}
		libff::leave_block("Compute polynomials A, B, C in Lagrange basis");

		libff::leave_block("Call to r1cs_convol_to_qap_instance_map");

		return qap_instance<FieldT>(domain,
				cs.num_variables(),
				domain->m,
				cs.num_inputs(),
				std::move(A_in_Lagrange_basis),
				std::move(B_in_Lagrange_basis),
				std::move(C_in_Lagrange_basis));
	}
	*/
	/**
	 * Instance map for the R1CS-to-QAP reduction followed by evaluation of the resulting QAP instance.
	 *
	 * Namely, given a R1CS constraint system cs and a field element t, construct
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
	/*
	template<typename FieldT>
		qap_instance_evaluation<FieldT> r1cs_to_qap_instance_map_with_evaluation(const r1cs_constraint_convol_system<FieldT> &cs,
				const FieldT &t)
		{



			libff::enter_block("Call to r1cs_convol_to_qap_instance_map_with_evaluation");

			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()));
			std::cout<<"degree m = "<<domain->m<<"\n";
			std::vector<FieldT> At, Bt, Ct, Ht;

			At.resize(cs.num_variables()+1, FieldT::zero());
			Bt.resize(cs.num_variables()+1, FieldT::zero());
			Ct.resize(cs.num_variables()+1, FieldT::zero());
			Ht.reserve(domain->m+1);

			const FieldT Zt = domain->compute_vanishing_polynomial(t);

			libff::enter_block("Compute evaluations of A, B, C, H at t");
			const std::vector<FieldT> u = domain->evaluate_all_lagrange_polynomials(t);
			/**
			 * add and process the constraints
			 *     input_i * 0 = 0
			 * to ensure soundness of input consistency
			 */
			/*
			   for(size_t j=0;j<cs.num_variables();j++){
			   for(size_t i=0;i<cs.num_constraints();i++){
			   if(i==0) std::cout<<j<<"\t";
			   std::cout<<cs.constraints[i].a.terms[j].coeff.as_ulong()<<"\t";
			   }
			   std::cout<<std::endl;
			   }
			 
			for (size_t i = 0; i <= cs.num_inputs(); ++i)
			{
				At[i] = u[cs.num_constraints() + i];
				//std::cout<<i<<", "<<cs.num_constraints()+i<<"\t";
			}
			//std::cout<<std::endl;

			size_t conv_index_a = 0;
			size_t conv_index_b = 0;
			size_t conv_index_c = 0;


			/* process all other constraints 
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				//std::cout<<"At ";
				for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j)
				{

					At[cs.constraints[i].a.terms[j].index] +=
						u[i]*cs.constraints[i].a.terms[j].coeff;
					// std::cout<<"( "<<i<<", "<<j<<"), const a term index"<<cs.constraints[i].a.terms[j].index<<"\t";
				}
				//std::cout<<std::endl;

				for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j)
				{
					Bt[cs.constraints[i].b.terms[j].index] +=
						u[i]*cs.constraints[i].b.terms[j].coeff;
				}

				for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j)
				{
					Ct[cs.constraints[i].c.terms[j].index] +=
						u[i]*cs.constraints[i].c.terms[j].coeff;
				}

				for(size_t j=0;j<cs.constraints[i].a2.terms.size();j++){
					///TODO change (j+1) to safe value like prime
					FieldT temp = FieldT::one();
					if(j==0) conv_index_a ++;
					for(size_t k=0;k<cs.num_convol_outputs();k++){

						At[cs.constraints[i].a2.terms[j].index] +=
							u[(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_a-1) + k ] * temp;//(FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						temp *= FieldT(k+1);
							//u[cs.num_inputs()+1 + (cs.num_constraints()) + (cs.num_convol_outputs())*(conv_index_a-1) + k ] * (FieldT(k+1)^cs.constraints[i].a2.terms[j].coeff.as_ulong());
						/*
						std::cout<<"At["<<(cs.constraints[i].a2.terms[j].index)<<
							"] = u["<<(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_a-1)+k<<
							"]*("<<j+1<<")^"<<cs.constraints[i].a2.terms[j].coeff.as_ulong()<<"\n";
						
					}
				}
				for(size_t j=0;j<cs.constraints[i].b2.terms.size();j++){
					///TODO change (j+1) to safe value like prime
					FieldT temp = FieldT::one();
					if(j==0) conv_index_b ++;
					for(size_t k=0;k<cs.num_convol_outputs();k++){
						Bt[cs.constraints[i].b2.terms[j].index] +=
							u[(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_b-1) + k ] * temp;//(FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						temp *= FieldT(k+1);
						/*
						std::cout<<"Bt["<<(cs.constraints[i].b2.terms[j].index)<<
							"] = u["<<(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_b-1)+k<<
							"]*("<<j+1<<")^"<<cs.constraints[i].b2.terms[j].coeff.as_ulong()<<"\n";
						
					}
				}

				for(size_t j=0;j<cs.constraints[i].c2.terms.size();j++){
					///TODO change (j+1) to safe value like prime
					FieldT temp = FieldT::one();
					if(j==0) conv_index_c ++;
					for(size_t k=0;k<cs.num_convol_outputs();k++){
						Ct[cs.constraints[i].c2.terms[j].index] +=
							u[(cs.num_inputs()+1+cs.num_constraints()) + (cs.num_convol_outputs()) *(conv_index_c-1) + k ] * temp;//(FieldT(k+1)^cs.constraints[i].c2.terms[j].coeff.as_ulong());
						temp *= FieldT(k+1);
						/*
						std::cout<<"Ct["<<(cs.constraints[i].c2.terms[j].index)<<
							"] = u["<<(cs.num_inputs()+1+cs.num_constraints())+(cs.num_convol_outputs())*(conv_index_c-1)+k<<
							"]*("<<j+1<<")^"<<cs.constraints[i].c2.terms[j].coeff.as_ulong()<<"\n";
						
					}
				}
			}

			FieldT ti = FieldT::one();
			for (size_t i = 0; i < domain->m+1; ++i)
			{
				Ht.emplace_back(ti);
				ti *= t;
			}
			libff::leave_block("Compute evaluations of A, B, C, H at t");

			libff::leave_block("Call to r1cs_convol_to_qap_instance_map_with_evaluation");

			return qap_instance_evaluation<FieldT>(domain,
					cs.num_variables(),
					domain->m,
					cs.num_inputs(),
					t,
					std::move(At),
					std::move(Bt),
					std::move(Ct),
					std::move(Ht),
					Zt);
		}
		*/

	/**
	 * Witness map for the R1CS-to-QAP reduction.
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
	/*
	template<typename FieldT>
		qap_witness<FieldT> r1cs_to_qap_witness_map(const r1cs_constraint_convol_system<FieldT> &cs,
				const r1cs_primary_input<FieldT> &primary_input,
				const r1cs_auxiliary_input<FieldT> &auxiliary_input,
				const FieldT &d1,
				const FieldT &d2,
				const FieldT &d3)
		{
			libff::enter_block("Call to r1cs_convol_to_qap_witness_map");

			/* sanity check 
			assert(cs.is_satisfied(primary_input, auxiliary_input));

			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()));

			r1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
			full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());


			libff::enter_block("Compute evaluation of polynomials A, B on set S");
			std::vector<FieldT> aA(domain->m, FieldT::zero()), aB(domain->m, FieldT::zero());

			std::cout<<"num inputs : "<<cs.num_inputs()<<" num constraints : "<< cs.num_constraints()<<"\n";
			/* account for the additional constraints input_i * 0 = 0 
			for (size_t i = 0; i <= cs.num_inputs(); ++i)
			{
				aA[i+cs.num_constraints()] = (i > 0 ? full_variable_assignment[i-1] : FieldT::one());
			}
			/* account for all other constraints 
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				aA[i] += cs.constraints[i].a.evaluate(full_variable_assignment);
				aB[i] += cs.constraints[i].b.evaluate(full_variable_assignment);
			}

			size_t conv_index = 0;
			size_t conv_index2 = 0;
			for( size_t i=0;i<cs.num_constraints();i++){
				FieldT acc = FieldT::zero();
				bool flag = false;
				for(auto &lt : cs.constraints[i].a2){
					if(!flag){ 
						conv_index++;
						flag = true;
					}
					FieldT temp = FieldT::one();
					for(size_t j=1;j<=cs.num_convol_outputs();j++){
						//if(lt.idex ==0) acc += (FieldT::one()*cs.conv_num_output()) * lt.coeff; //there is no first one in conv constraints
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						temp *= FieldT(j);
						//std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						aA[(j)+(conv_index-1)*cs.num_convol_outputs() + cs.num_constraints() + cs.num_inputs()] += acc;
						//std::cout<<"aA["<<(j-1)+(conv_index-1)*cs.num_convol_outputs()+cs.num_constraints()+cs.num_inputs()<<"] = acc\n";
					}
					//if(acc != FieldT::zero()){
					//}
				}
				/*
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					std::cout<<"aA["<<(j)+cs.num_constraints()+cs.num_inputs()<<"] ="<<
						aA[(j) + cs.num_constraints() + cs.num_inputs()].as_ulong() <<"\n";
				}
				
				acc = FieldT::zero();
				flag = false;
				for(auto &lt : cs.constraints[i].b2){
					if(!flag){ 
						conv_index2++;
						flag = true;
					}
					FieldT temp = FieldT::one();
					for(size_t j=1;j<=cs.num_convol_outputs();j++){
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						temp *= FieldT(j);
						//std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						aB[(conv_index2-1)*cs.num_convol_outputs() +(j) + cs.num_constraints() + cs.num_inputs()] += acc;
						//std::cout<<"aB["<<(conv_index2-1)*cs.num_convol_outputs()+(j-1)+cs.num_constraints()+cs.num_inputs()<<"] = acc\n";
					}
					//aB[i*cs.num_convol() + j + cs.num_constraints() + cs.num_inputs()] += acc;
				}
				/*
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					std::cout<<"aB["<< +(j)+cs.num_constraints()+cs.num_inputs()<<"] ="<<
						aB[(j) + cs.num_constraints() + cs.num_inputs()].as_ulong() <<"\n";
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
			// add coefficients of the polynomial (d2*A + d1*B - d3) + d1*d2*Z 
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
				bool flag = false;
				for(auto &lt : cs.constraints[i].c2){
					if(!flag){
						conv_index++;
						flag = true;
					}
					FieldT temp = FieldT::one();
					for(size_t j=1;j<=cs.num_convol_outputs();j++){
						acc = full_variable_assignment[lt.index-1] * temp;//(FieldT(j)^lt.coeff.as_ulong());
						temp *= FieldT(j);
						aC[(conv_index-1)*cs.num_convol_outputs()+ (j) + cs.num_constraints() + cs.num_inputs()] += acc;
						//std::cout<<"acc : "<<acc.as_ulong()<<" var["<<lt.index-1<<"]*("<<j<<"^"<<lt.coeff.as_ulong()<<")\n";
						//std::cout<<"aC["<<(conv_index-1)*cs.num_convol_outputs() +(j-1)+cs.num_constraints()+cs.num_inputs()<<"] = acc\n";
					}
				}
				/*
				//for debug
				for(size_t j=1;j<=cs.num_convol_outputs();j++){
					std::cout<<"aC["<<(j)+cs.num_constraints()+cs.num_inputs()<<"] ="<<
						aC[(j) + cs.num_constraints() + cs.num_inputs()].as_ulong() <<"\n";
				}
				/

			}
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

			libff::leave_block("Call to r1cs_convol_to_qap_witness_map");

			return qap_witness<FieldT>(cs.num_variables(),
					domain->m,
					cs.num_inputs(),
					d1,
					d2,
					d3,
					full_variable_assignment,
					std::move(coefficients_for_H));
		}
*/
} // libsnark

#endif // R1CS_TO_QP_TCC_