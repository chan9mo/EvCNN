/** @file


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
#include <libff/algebra/scalar_multiplication/multiexp.hpp>

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

		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
		//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs(0)));//*(cs.num_convol()));

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

			for(size_t j=0;j<cs.constraints[i].a1.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_a ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					A_in_Lagrange_basis[cs.constraints[i].a1.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].a1.terms[j].coeff.as_ulong());
				}
			}
			for(size_t j=0;j<cs.constraints[i].b1.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_b ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					B_in_Lagrange_basis[cs.constraints[i].b1.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].b1.terms[j].coeff.as_ulong());
				}
			}

			for(size_t j=0;j<cs.constraints[i].c1.terms.size();j++){
				///TODO change (j+1) to safe value like prime
				if(j==0) conv_index_c ++;
				for(size_t k=0;k<cs.num_convol_outputs();k++){
					C_in_Lagrange_basis[cs.constraints[i].c1.terms[j].index][i] +=
						(FieldT(j+1)^cs.constraints[i].c1.terms[j].coeff.as_ulong());
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

			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
			// size_t degree = cs.num_constraints() + cs.num_inputs() + 1;
			// for(size_t i=0;i<cs.num_convol;i++){
			// 	degree += (cs.num_convol_outputs(i))*(cs.num_convol_outputs2(i));
			// }
			// const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(degree);


			/*
			else{
				domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()+1));
			}
			*/
			//std::cout<<"#constraint : "<< cs.num_constraints()<<"\t#input : "<< cs.num_inputs()<<"\t#convolOut : "<<(cs.num_convol_outputs(0))<<"\t#convol : "<<(cs.num_convol)<<std::endl;

			std::vector<FieldT> At, Bt, Ct, Ht;


			std::cout<<"test 4^-1: "<<FieldT(4).inverse().as_ulong()<<", -1 : "<<FieldT(-1).as_ulong()<<std::endl;

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
			 * to ensure soundness of input consistency			 */
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

			// size_t conv_index_a = 0;
			// size_t conv_index_b = 0;
			// size_t conv_index_c = 0;

			// size_t conv_index_a2 = 0;
			// size_t conv_index_b2 = 0;
			// size_t conv_index_c2 = 0;

			// bool conv_flag = false;

			/*
			if(false && cs.num_convol){

				libff::enter_block("Compute temps");
				std::vector<std::vector<std::vector<FieldT>>> temp11,temp22;
				temp11.reserve(cs.num_convol);
				temp22.reserve(cs.num_convol);
				for (size_t i = 0; i < cs.num_convol; ++i)
				{
					std::vector<std::vector<FieldT>> tv1(cs.num_convol_outputs(i), std::vector<FieldT>(cs.num_convol_outputs(i), FieldT::zero()));
					std::vector<std::vector<FieldT>> tv2(cs.num_convol_outputs2(i), std::vector<FieldT>(cs.num_convol_outputs2(i), FieldT::zero()));
					temp11.push_back(tv1);
					temp22.push_back(tv2);
					for(size_t k=0;k<cs.num_convol_outputs(i);k++){
						temp11[i][k][0] = FieldT::one();
						for(size_t kk=1;kk<cs.num_convol_outputs(i);kk++){
							temp11[i][k][kk] = temp11[i][k][kk-1]*FieldT(k+1);
						}


						// for(size_t kk=0;kk<cs.num_convol_outputs(i);kk++){
						// 	cout<<temp11[i][k][kk].as_ulong()<<"\t";
						// }
						// cout<<endl;
					}
					for(size_t j=0;j<cs.num_convol_outputs2(i);j++){
						temp22[i][j][0] = FieldT::one();
						for(size_t jj=1;jj<cs.num_convol_outputs2(i);jj++){
							temp22[i][j][jj] = temp22[i][j][jj-1]*FieldT(j+1);
						}
					}
				}

				//libff::enter_block("Cal test");
				//cout<<"yH : "<<cs.num_convol_outputs(0)<<", yW :"<<cs.num_convol_outputs(0)<<endl;
				//cout<<"iH : "<<cs.num_convol_input_height(0)<<", iW :"<<cs.num_convol_input_width(0)<<endl;
				//cout<<"aterm : "<<cs.constraints[0].a1.terms.size()<<endl;

				std::vector<FieldT> test(cs.num_convol_outputs(conv_index_a)*cs.num_convol_outputs2(conv_index_a), FieldT::zero());
				for (size_t i = 0; i < cs.num_constraints(); ++i){
					if(cs.constraints[i].a1.terms.size()>1){
					for(int x=0;x<cs.num_convol_outputs(conv_index_a);x++){
						for(int j=0;j<cs.num_convol_outputs2(conv_index_a);j++){
							size_t u_index = (cs.num_inputs()+1+cs.num_constraints()) + (x*(cs.num_convol_outputs2(conv_index_a))) + j;
							for(int k=0;k<cs.num_convol_outputs(conv_index_a);k++){
								for(int l=0;l<cs.num_convol_outputs2(conv_index_a);l++){
									test[k*cs.num_convol_outputs2(conv_index_a) +l] = temp11[conv_index_a][x][k]*temp22[conv_index_a][j][l];
								}
							}
							for(int k=0;k<cs.num_convol_input_height(conv_index_a);k++){
								for(int l=0;l<cs.num_convol_input_width(conv_index_a);l++){
									size_t a_index = k*cs.num_convol_input_width(conv_index_a) + l+1;
									//cout<<x+1<<"^"<<k<<"="<<temp11[0][x][k].as_ulong()<<", "<<j+1<<"^"<<l<<"="<<temp22[0][j][l].as_ulong()<<endl;
									if(cs.constraints[i].a1.terms[a_index].index){
										At[cs.constraints[i].a1.terms[a_index].index] +=
											u[u_index] * test[k*cs.num_convol_outputs2(conv_index_a)+l];//temp11[0][i][k]*temp22[0][j][l];
										//cout<<"aidx : "<<cs.constraints[i].a1.terms[a_index].index<<", uidx :"<<u_index<<" temp : "<<test[k*cs.num_convol_outputs2(0)+l].as_ulong()<<endl;
									}
								}
							}
							for(int k=0;k<cs.num_convol_kernel_height(conv_index_a);k++){
								for(int l=0;l<cs.num_convol_kernel_width(conv_index_a);l++){
									size_t b_index = k*cs.num_convol_kernel_width(conv_index_a) + l+1;
									Bt[cs.constraints[i].b1.terms[b_index].index] +=
										u[u_index] * test[k*cs.num_convol_outputs2(conv_index_a)+l];//temp11[0][i][k]*temp22[0][j][l];
									//cout<<"k :"<<k<<",l :"<<l<<endl;
									//cout<<"bidx : "<<cs.constraints[i].b1.terms[b_index].index<<", uidx :"<<u_index<<" temp : "<<test[k*cs.num_convol_outputs2(0)+l].as_ulong()<<endl;

								}
							}
							for(int k=0;k<cs.num_convol_outputs(conv_index_a);k++){
								for(int l=0;l<cs.num_convol_outputs2(conv_index_a);l++){
									size_t c_index = k*cs.num_convol_outputs2(conv_index_a) + l+1;
									Ct[cs.constraints[i].c1.terms[c_index].index] +=
										u[u_index] * test[k*cs.num_convol_outputs2(conv_index_a)+l];//temp11[0][i][k]*temp22[0][j][l];
									//cout<<"cidx : "<<cs.constraints[i].c1.terms[c_index].index<<", uidx :"<<u_index<<" temp : "<<test[k*cs.num_convol_outputs2(0)+l].as_ulong()<<endl;
								}
							}
						}
					}
					conv_index_a++;
					}
				}
				//libff::leave_block("Cal test");
				libff::leave_block("Compute temps");
			}
			*/

			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				//std::cout<<"At\n";
				for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j)
				{
					At[cs.constraints[i].a.terms[j].index] +=
						u[i]*cs.constraints[i].a.terms[j].coeff;
					//std::cout<<"( "<<i<<", "<<cs.constraints[i].a.terms[j].index<<")="<<cs.constraints[i].a.terms[j].coeff.as_ulong()<<"\n";
				}
				//std::cout<<std::endl;
				//std::cout<<"bt\n";
				for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j)
				{
					Bt[cs.constraints[i].b.terms[j].index] +=
						u[i]*cs.constraints[i].b.terms[j].coeff;
					//std::cout<<"( "<<i<<", "<<cs.constraints[i].b.terms[j].index<<")="<<cs.constraints[i].b.terms[j].coeff.as_ulong()<<"\n";

				}
				//std::cout<<"Ct\n";
				for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j)
				{
					Ct[cs.constraints[i].c.terms[j].index] +=
						u[i]*cs.constraints[i].c.terms[j].coeff;
					//std::cout<<"( "<<i<<", "<<cs.constraints[i].c.terms[j].index<<")="<<cs.constraints[i].c.terms[j].coeff.as_ulong()<<"\n";

				}
				//if(cs.num_convol_dimensions(0) == 2){
				//std::cout<<"a1 size : "<<cs.constraints[i].a1.terms.size()<<std::endl;
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

			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1);
			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs())*(cs.num_convol()+1));
			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs(0))*(cs.num_convol_outputs2(0)));
			/*
			size_t degree = cs.num_constraints() + cs.num_inputs() + 1;
			for(size_t i=0;i<cs.num_convol;i++){
				degree += (cs.num_convol_outputs(i))*(cs.num_convol_outputs2(i));
			}
			const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(degree);
			*/

			//std::cout<<"#constraint : "<< cs.num_constraints()<<"\t#input : "<< cs.num_inputs()<<"\t#convolOut : "<<(cs.num_convol_outputs(0))<<"\t#convol : "<<(cs.num_convol)<<std::endl;
			r1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
			full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());
			/*
			std::cout<<"variables : ";
			for(size_t i=0;i<full_variable_assignment.size();i++){
				std::cout<<full_variable_assignment[i].as_ulong()<<"\t";
			}
			*/
			std::cout<<"\n";
			//std::cout<<"degree : "<<domain->m<<" sum = "<<cs.num_constraints() + cs.num_inputs() + 1 + (cs.num_convol_outputs(0))*(cs.num_convol+1)<<"\n";

			libff::enter_block("Compute evaluation of polynomials A, B on set S");
			std::vector<FieldT> aA(domain->m, FieldT::zero()), aB(domain->m, FieldT::zero());
			std::vector<FieldT> aC(domain->m, FieldT::zero());

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
				//std::cout<<"aA["<<i<<"] = "<<aA[i].as_ulong()<<std::endl;
				aB[i] += cs.constraints[i].b.evaluate(full_variable_assignment);
				//std::cout<<"aB["<<i<<"] = "<<aB[i].as_ulong()<<std::endl;

			}
			/*
			size_t conv_index_a = 0;
			size_t conv_index = 0;
			size_t conv_index2 = 0;

			if(cs.num_convol){
				libff::enter_block("Compute temps");
				std::vector<std::vector<std::vector<FieldT>>> temp11,temp22;
				temp11.reserve(cs.num_convol);
				temp22.reserve(cs.num_convol);
				for (size_t i = 0; i < cs.num_convol; ++i)
				{
					std::vector<std::vector<FieldT>> tv1(cs.num_convol_outputs(i), std::vector<FieldT>(cs.num_convol_outputs(i), FieldT::zero()));
					std::vector<std::vector<FieldT>> tv2(cs.num_convol_outputs2(i), std::vector<FieldT>(cs.num_convol_outputs2(i), FieldT::zero()));
					temp11.push_back(tv1);
					temp22.push_back(tv2);
					for(size_t k=0;k<cs.num_convol_outputs(i);k++){
						temp11[i][k][0] = FieldT::one();
						for(size_t kk=1;kk<cs.num_convol_outputs(i);kk++){
							temp11[i][k][kk] = temp11[i][k][kk-1]*FieldT(k+1);

						}

						// for(size_t kk=0;kk<cs.num_convol_outputs(i);kk++){
						// 	cout<<temp11[i][k][kk].as_ulong()<<"\t";
						// }
						// cout<<endl;
					}
					for(size_t j=0;j<cs.num_convol_outputs2(i);j++){
						temp22[i][j][0] = FieldT::one();
						for(size_t jj=1;jj<cs.num_convol_outputs2(i);jj++){
							temp22[i][j][jj] = temp22[i][j][jj-1]*FieldT(j+1);
						}
						// for(size_t jj=0;jj<cs.num_convol_outputs2(i);jj++){
						// 	cout<<temp22[i][j][jj].as_ulong()<<"\t";
						// }
						// cout<<endl;
					}
				}
				libff::leave_block("Compute temps");

				libff::enter_block("Cal test");
				std::vector<FieldT> test(cs.num_convol_outputs(0)*cs.num_convol_outputs2(0), FieldT::zero());
				for (size_t i = 0; i < cs.num_constraints(); ++i){
					if(cs.constraints[i].a1.terms.size()>1){

					FieldT acc = FieldT::zero();
					for(int x=0;x<cs.num_convol_outputs(0);x++){
						for(int j=0;j<cs.num_convol_outputs2(0);j++){
							size_t u_index = (cs.num_inputs()+1+cs.num_constraints()) + (x*(cs.num_convol_outputs2(0))) + j;
							for(int k=0;k<cs.num_convol_outputs(0);k++){
								for(int l=0;l<cs.num_convol_outputs2(0);l++){
									test[k*cs.num_convol_outputs2(0) +l] = temp11[0][x][k]*temp22[0][j][l];
								}
							}
							for(int k=0;k<cs.num_convol_input_height(0);k++){
								for(int l=0;l<cs.num_convol_input_width(0);l++){
									size_t a_index = k*cs.num_convol_input_width(0) + l+1;

									if(cs.constraints[i].a1.terms[a_index].index){
										acc = full_variable_assignment[cs.constraints[i].a1.terms[a_index].index-1] * test[k*cs.num_convol_outputs2(0)+l];
										aA[u_index] += acc;
										//cout<<"aidx : "<<cs.constraints[i].a1.terms[a_index].index<<", uidx :"<<u_index<<" temp : "<<test[k*cs.num_convol_outputs2(0)+l].as_ulong()<<endl;
									}
								}
							}
							for(int k=0;k<cs.num_convol_kernel_height(0);k++){
								for(int l=0;l<cs.num_convol_kernel_width(0);l++){
									size_t b_index = k*cs.num_convol_kernel_width(0) + l+1;
									acc = full_variable_assignment[cs.constraints[i].b1.terms[b_index].index-1] * test[k*cs.num_convol_outputs2(0)+l];
									aB[u_index] += acc;
								}
							}
							for(int k=0;k<cs.num_convol_outputs(0);k++){
								for(int l=0;l<cs.num_convol_outputs2(0);l++){
									size_t c_index = k*cs.num_convol_outputs2(0) + l+1;
									acc = full_variable_assignment[cs.constraints[i].c1.terms[c_index].index-1] * test[k*cs.num_convol_outputs2(0)+l];
									aC[u_index] += acc;
									//cout<<"cidx : "<<cs.constraints[i].c1.terms[c_index].index<<", uidx :"<<u_index<<" temp : "<<test[k*cs.num_convol_outputs2(0)+l].as_ulong()<<endl;
								}
							}
						}
					}
					}
				}
				libff::leave_block("Cal test");
			}
			*/


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
			//std::vector<FieldT> aC(domain->m, FieldT::zero());
			for (size_t i = 0; i < cs.num_constraints(); ++i)
			{
				aC[i] += cs.constraints[i].c.evaluate(full_variable_assignment);
				//std::cout<<"aC["<<i<<"] = "<<aC[i].as_ulong()<<std::endl;

			}
			//libff::leave_block("Cal Conv C");
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

// 		 	template<typename FieldT>
// std::vector<FieldT> zxPolycalculate(
// 	const r1cs_constraint_system<FieldT> &cs,
// 	r1cs_variable_assignment<FieldT> full_variable_assignment,
// 	std::vector<FieldT> aA,
// 	linear_combination<FieldT> a,
// 	linear_combination<FieldT> a1,
// 	size_t num_variables,

// 	std::vector<std::vector<FieldT>> zxPoly,
// 		size_t num_constraints, size_t num_z_order,
// 		std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain
// 		)
// {
// 	std::vector<FieldT> tempPoly(domain->m, FieldT::zero());
// 	for(size_t i=0; i<cs.num_constraints();i++){
// 		aA[i] += cs.constraints[i].a.evaluate(full_variable_assignment);


// 	}
// 	for(size_t i=0; i<cs.num_constraints();i++){
// 		if(cs.constraints[i].a1.size()){
// 			std:::vector<FieldT> tempPoly(domain->m, FieldT::zero());
// 			tempPoly[i] = FieldT::one();
// 			domain->iFFT(tempPoly[i])
// 		}
// 	}
// 	//std::vector<std::vector<FieldT>> tempPoly(num_constraints, std::vector<FieldT>(num_constraints, FieldT::zero()));
//  	bool isEmpty = true;
// 	for(size_t i=0; i<num_constraints;i++){
// 		size_t j;
// 		for(j=0; j<num_z_order+1;j++){
// 			if(zxPoly[i][j] != FieldT::zero()){
// 				isEmpty = false;
// 				break;
// 			}
// 		}
// 		if(j < num_z_order+1) tempPoly[i][i] = FieldT::one();
// 	}
//  	if(isEmpty) return std::vector<FieldT>(1, FieldT::zero());
//  	for(size_t i=0; i<num_constraints;i++){
// 		domain->iFFT(tempPoly[i]);
// 	}
//  	//v0(x) * z(x^(4d+1))
// 	std::vector<FieldT> zExtendPoly((num_z_order*((num_constraints-1)*4+1))+1, FieldT::zero());
// 	for(size_t i=0;i<num_z_order+1;i++){
// 		zExtendPoly[i*((num_constraints-1)*4 + 1)] = zxPoly[0][i];
// 	}
//  	std::vector<FieldT> mulResult(1, FieldT::zero());
// 	libfqfft::_polynomial_multiplication(mulResult, tempPoly[0], zExtendPoly);
//  	//point add
// 	for(size_t i=0; i<num_constraints;i++){
// 		FieldT temp = FieldT::zero();
// 		for(size_t j=1; j<num_constraints;j++){
// 			temp+=tempPoly[j][i];
// 		}
// 		mulResult[i] += temp;
// 	}
//  	return mulResult;
// }
	template<typename FieldT>
qap_instance_evaluation<FieldT> mr1cs_to_qap_instance_map_with_evaluation(
	const r1cs_constraint_system<FieldT> &cs,
		const FieldT &t
		)
{

	libff::enter_block("Call to mr1cs_to_qap_instance_map_with_evaluation");


	std::cout<<"cs num : "<<cs.num_constraints()<<std::endl;
	const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints() + 1); //+ cs.num_inputs() + 1);
	size_t degree = 0;
	if(cs.num_convol) degree = (2*domain->m+1)*cs.num_convol_outputs(0);
	else degree = (2*domain->m+1);
	const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2 = libfqfft::get_evaluation_domain<FieldT>(degree);
	const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>(domain2->m*2);

	std::cout<<"dx : "<<domain->m<<"dz : "<<domain2->m<<std::endl;
	std::vector<FieldT> At, Bt, Ct, Ht;

 	//vector<FieldT> At(domain2->m, FieldT::zero());
	//vector<FieldT> Bt(domain2->m, FieldT::zero());
	//vector<FieldT> Ct(domain2->m, FieldT::zero());
	//vector<FieldT> Ht;


	const FieldT Zt = domain->compute_vanishing_polynomial(t);
	std::cout<<"Zt : "<<Zt.as_ulong()<<std::endl;
	libff::enter_block("Compute evaluations of A, B, C, H at t");

	const std::vector<FieldT> u = domain->evaluate_all_lagrange_polynomials(t);
	std::cout<<"num var : "<<cs.num_variables()<<std::endl;

	At.resize(cs.num_variables()+1, FieldT::zero());
	Bt.resize(cs.num_variables()+1, FieldT::zero());
	Ct.resize(cs.num_variables()+1, FieldT::zero());
 	//At, Bt, Ct = poly(t);
	std::cout<<"size(#var +1) :"<<cs.num_variables()+1<<std::endl;

	Ht.reserve(domain4->m);
	FieldT ti = FieldT::one();
	for(size_t i=0; i<domain4->m;i++){
		Ht.emplace_back(ti);
		ti*= t;
	}

	/*
	for (size_t i = 0; i <= cs.num_inputs(); ++i)
	{
		At[i] = u[cs.num_constraints() + i];
		//std::cout<<i<<", "<<cs.num_constraints()+i<<"\t";
	}
	*/
	/*
	for (size_t i = 0; i < cs.num_constraints(); ++i)
	{
		//std::cout<<"At\n";
		for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j)
		{
			At[cs.constraints[i].a.terms[j].index] +=
				u[i]*cs.constraints[i].a.terms[j].coeff;
		}
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
	}
	*/
	for (size_t i = 0; i < cs.num_constraints(); ++i){
		if(cs.constraints[i].a1.terms.size()>1){
			//  std::cout << "cs : "<<i<<std::endl;
			//  std::cout <<"degree"<<std::endl;
			//  std::cout<< "A\n";
			for (size_t j = 1; j < cs.constraints[i].a1.terms.size(); ++j)
			{
				At[cs.constraints[i].a1.terms[j].index] +=
					u[i]*Ht[(cs.constraints[i].a1.terms[j].coeff.as_ulong()-1)*(2*domain->m+1)];//cs.constraints[i].a.terms[j].coeff;
				  //std::cout<<cs.constraints[i].a1.terms[j].index<<", "<<(cs.constraints[i].a1.terms[j].coeff.as_ulong()-1)*(2*domain->m+1)<<std::endl;
			}

			// std::cout<< "B\n";
			for (size_t j = 1; j < cs.constraints[i].b1.terms.size(); ++j)
			{
				Bt[cs.constraints[i].b1.terms[j].index] +=
					u[i]*Ht[(cs.constraints[i].b1.terms[j].coeff.as_ulong()-1)*(2*domain->m+1)];
				// std::cout<<cs.constraints[i].b1.terms[j].index<<", "<<(cs.constraints[i].b1.terms[j].coeff.as_ulong()-1)*(2*domain->m+1)<<std::endl;

			}
			// std::cout<< "C\n";
			for (size_t j = 1; j < cs.constraints[i].c1.terms.size(); ++j)
			{
				Ct[cs.constraints[i].c1.terms[j].index] +=
					u[i]*Ht[(cs.constraints[i].c1.terms[j].coeff.as_ulong()-1)*(2*domain->m+1)];
				// std::cout<<cs.constraints[i].c1.terms[j].index<<", "<<(cs.constraints[i].c1.terms[j].coeff.as_ulong()-1)*(2*domain->m+1)<<std::endl;
			}
		}
	}




	libff::leave_block("Compute evaluations of A, B, C, H at t");

	libff::leave_block("Call to mr1cs_to_qap_instance_map_with_evaluation");
 	return qap_instance_evaluation<FieldT>(domain2,
			cs.num_variables(),
			domain2->m,
			cs.num_inputs(),
			t,
			move(At),
			move(Bt),
			move(Ct),
			move(Ht),
			Zt);
 }
		template<typename FieldT>
	qap_witness<FieldT> mr1cs_to_qap_witness_map(
			//const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain2,
			const r1cs_constraint_system<FieldT> &cs,
			const r1cs_primary_input<FieldT> &primary_input,
			const r1cs_auxiliary_input<FieldT> &auxiliary_input,
			const FieldT &d1,
			const FieldT &d2,
			const FieldT &d3,
			const FieldT &t
			)
	{
		std::cout<<"t : "<<t.as_ulong()<<std::endl;
		//t = FieldT(2);
		//d1 = d2 = d3 = FieldT::zero();
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(cs.num_constraints()+1);//+ cs.num_inputs() + 1);
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2 =
						libfqfft::get_evaluation_domain<FieldT>((2*domain->m+1)*cs.num_convol_outputs(0));

		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>(domain2->m*2);
			std::cout<<"d : "<<domain2->m<<std::endl;


		r1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
		full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());



		libff::enter_block("Compute evaluation of polynomials A, B on set S");
		std::vector<FieldT> aA(domain->m, FieldT::zero()), aB(domain->m, FieldT::zero());
		std::vector<FieldT> aC(domain->m, FieldT::zero());
		/* account for the additional constraints input_i * 0 = 0 */
		// for (size_t i = 0; i <= cs.num_inputs(); ++i)
		// {
		// 	aA[i+cs.num_constraints()] = (i > 0 ? full_variable_assignment[i-1] : FieldT::one());
		// 	//std::cout<<"aA["<<i+cs.num_constraints()<<"] = "<<aA[i+cs.num_constraints()].as_ulong()<<std::endl;

		// }
		/* account for all other constraints */
		/*
		for (size_t i = 0; i < cs.num_constraints(); ++i)
		{
			aA[i] += cs.constraints[i].a.evaluate(full_variable_assignment);
			//std::cout<<"aA["<<i<<"] = "<<aA[i].as_ulong()<<std::endl;
			aB[i] += cs.constraints[i].b.evaluate(full_variable_assignment);
			//std::cout<<"aB["<<i<<"] = "<<aB[i].as_ulong()<<std::endl;
			aC[i] += cs.constraints[i].c.evaluate(full_variable_assignment);

		}
		*/
	/*
		libff::enter_block("Compute coefficients of polynomial A");
		domain->iFFT(aA);
		libff::leave_block("Compute coefficients of polynomial A");

		libff::enter_block("Compute coefficients of polynomial B");
		domain->iFFT(aB);
		libff::leave_block("Compute coefficients of polynomial B");

		libff::enter_block("Compute coefficients of polynomial C");
		domain->iFFT(aC);
		libff::leave_block("Compute coefficients of polynomial C");
	*/
		// FieldT aA1 = FieldT::zero();
		// for (auto x : aA)
		// {
		// 	aA1 += x;
		// 	std::cout << x.as_ulong() << std::endl;
		// }

		// std::cout << "aB\n";
		// FieldT aB1 = FieldT::zero();
		// for (auto x : aB)
		// {
		// 	aB1 += x;
		// 	std::cout << x.as_ulong() << std::endl;
		// }
		// std::cout << "aC\n";
		// FieldT aC1 = FieldT::zero();
		// for (auto x : aC)
		// {
		// 	aC1 += x;
		// 	std::cout << x.as_ulong() << std::endl;
		// }
		// std::cout << "aA : " << (aA1).as_ulong() << ", aB: " << (aB1).as_ulong() << ", aC : " << (aC1).as_ulong() << std::endl;

		// std::cout << "aA*aB-aC : " << (aA1 * aB1 - aC1).as_ulong() << std::endl;

		// ///qap witness start
		const std::vector<FieldT> u = domain->evaluate_all_lagrange_polynomials(t);

		std::vector<FieldT> tempMul(1, FieldT::zero());
		size_t conv_index = 0;
		libff::enter_block("calcuate A, B, C for z");
//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
		for(size_t i=0;i<cs.num_constraints();i++){
			FieldT f1 = FieldT::zero();
			FieldT f2 = FieldT::zero();
			FieldT f3 = FieldT::zero();

			if(cs.constraints[i].a1.terms.size()>1){
				//std::cout << "cs : "<<i<<std::endl;
				std::vector<FieldT> tempPoly(domain->m, FieldT::zero());
				//tempPoly has (w^0, 1) == L_0(x)
				tempPoly[i] = FieldT::one();
				//Point to Polynomial
				domain->iFFT(tempPoly);
				/* //debug check lagrange poly
				std::cout<<"L(x) Points for degree 4"<<std::endl;
				std::vector<FieldT> tempD4(domain4->m, FieldT::zero());
				for(size_t j=0; j<domain->m;j++){
					std::cout<<"j : "<<j<<std::endl;
					tempD4[j] = tempPoly[j];
					
				}

				domain4->FFT(tempD4);
				for(size_t j=0;j<domain4->m;j++){
					std::cout<<j<<", "<<tempD4[j].as_ulong()<<std::endl;
				}
				*/

				/* //debug
				//std::vector<FieldT> tempP2(tempPoly);
				//tempP2 = (L_0(x)-1)
				//tempP2[0] = tempP2[0]-1;
				//tempMul = L_0(x)(L_0(x)-1)
				//libfqfft::_polynomial_multiplication(tempMul, tempPoly, tempP2);
				//tempMul.resize(domain4->m, FieldT::zero());
				//tempMul Poly to Point
				//domain4->FFT(tempMul);
				std::cout<<"L(L-1) Points"<<std::endl;
				for(size_t j=0;j<domain4->m;j++){
					std::cout<<j<<", "<<tempMul[j].as_ulong()<<std::endl;
				}
				*/

				
				std::vector<FieldT> zExtendPoly(domain2->m+1, FieldT::zero());//(cs.num_convol_input_height(0)*((domain->m)*4+1)+1), FieldT::zero());
				for(size_t j=1;j<cs.constraints[i].a1.terms.size();j++){
					//ttt = extended degree = degree * (2d_x + 1)
					size_t ttt = (cs.constraints[i].a1.terms[j].coeff.as_ulong()-1)*(2*(domain->m)+1);
					//zExtPoly has (var * x^ttt)
					zExtendPoly[ttt] = full_variable_assignment[cs.constraints[i].a1.terms[j].index -1];
					//std::cout<<ttt<<", var["<<cs.constraints[i].a1.terms[j].index -1<<"] = "<<full_variable_assignment[cs.constraints[i].a1.terms[j].index -1].as_ulong()<<std::endl;
				}

				/* //debug for check extPoly
				 std::cout<<"zExtendPoly : ";
				 for(size_t j=0; j<zExtendPoly.size();j++){
				 	std::cout<<zExtendPoly[j].as_ulong()<< "*x^"<<j<<" ";
				 }
				 for(auto x : zExtendPoly){
				 	f1 += x;
				 }
				 std::cout<<std::endl;
				 std::cout<<"temp : ";
				 FieldT tem = FieldT::one();
				 FieldT sumt = FieldT::zero();
				 for(size_t j=0; j<tempPoly.size();j++){
				 	std::cout<<tempPoly[j].as_ulong()<< "*x^"<<j<<" ";
				 	sumt += tempPoly[j]*tem;
				 	tem *= t;
				 }
				 std::cout<<std::endl;
				 std::cout<<"L(t) : "<<u[i].as_ulong()<<", sumt : "<<sumt.as_ulong()<<std::endl;
				 */
				std::vector<FieldT> mulResult(1, FieldT::zero());
				// tempPoly = L_0(x), zExtPoly = X(z) => V(x,z) = L_0(x) * X(z=x^2d_x+1)
				libfqfft::_polynomial_multiplication(mulResult, tempPoly, zExtendPoly);

				/* // debug check result of poly mul
				 std::cout<<"mul : ";
				 for(size_t j=0; j<mulResult.size();j++){
				 	std::cout<<mulResult[j].as_ulong()<< "*x^"<<j<<" ";
				 }
				 std::cout<<std::endl;

				 std::cout<<"f1 : "<<f1.as_ulong()<<std::endl;
				 */
				for(size_t j=0;j<mulResult.size();j++){
					if(j < aA.size()){
						aA[j] += mulResult[j];
					}
					else{
						aA.push_back(mulResult[j]);
					}
				}
				/*
				  //debug check evalutation of mul poly 
				std::vector<FieldT> testA(aA);
				testA.resize(domain4->m, FieldT::zero());
				domain4->FFT(testA);
				for(size_t j=0;j<testA.size();j++){
					std::cout<<j<<", "<<testA[j].as_ulong()<<std::endl;
				}
				*/
				
			}
			if(cs.constraints[i].b1.terms.size()>1){
				std::vector<FieldT> tempPoly(domain->m, FieldT::zero());
				//tempPoly = (w^0, 1)
				tempPoly[i] = FieldT::one();
				//tempPoly point to poly
				domain->iFFT(tempPoly);

				std::vector<FieldT> zExtendPoly(domain2->m+1, FieldT::zero());//(cs.num_convol_kernel_height(0)*((domain->m)*4+1)+1), FieldT::zero());
				for(size_t j=1;j<cs.constraints[i].b1.terms.size();j++){
					size_t ttt = (cs.constraints[i].b1.terms[j].coeff.as_ulong()-1)*(2*(domain->m)+1);

					//std::cout<<(cs.constraints[i].b1.terms[j].coeff.as_ulong()-1)<<std::endl;
					zExtendPoly[(cs.constraints[i].b1.terms[j].coeff.as_ulong()-1)*(2*(domain->m)+1)] = full_variable_assignment[cs.constraints[i].b1.terms[j].index - 1];
					//std::cout<<ttt<<", var["<<cs.constraints[i].b1.terms[j].index -1<<"] = "<<full_variable_assignment[cs.constraints[i].b1.terms[j].index -1].as_ulong()<<std::endl;
				}

				/* //debug check ext poly
				 std::cout<<"zExtendPoly : ";
				 for(size_t j=0; j<zExtendPoly.size();j++){
				 	std::cout<<zExtendPoly[j].as_ulong()<< "*x^"<<j<<" ";
				 }
				
				 for(auto x : zExtendPoly){
				 	f2 += x;
				 }
				 std::cout<<"temp : ";
				 FieldT tem = FieldT::one();
				 FieldT sumt = FieldT::zero();
				 for(size_t j=0; j<tempPoly.size();j++){
				 	std::cout<<tempPoly[j].as_ulong()<< "*x^"<<j<<" ";
				 	sumt += tempPoly[j]*tem;
				 	tem *= t;
				 }
				 std::cout<<std::endl;
				 //L(t) = tempPoly(t), sumt = tempPoly(t); check above relation
				 std::cout<<"L(t) : "<<u[i].as_ulong()<<", sumt : "<<sumt.as_ulong()<<std::endl;
				 */
				
				std::vector<FieldT> mulResult(1, FieldT::zero());
				libfqfft::_polynomial_multiplication(mulResult, tempPoly, zExtendPoly);
				
				//std::cout<<"f2 : "<<f2.as_ulong()<<std::endl;
				for(size_t j=0;j<mulResult.size();j++){
					if(j < aB.size()){
						aB[j] += mulResult[j];
					}
					else{
						aB.push_back(mulResult[j]);
					}
				}
				/* //debug check evaluation
				std::vector<FieldT> testB(aB);
				testB.resize(domain4->m, FieldT::zero());
				domain4->FFT(testB);
				for(size_t j=0;j<testB.size();j++){
					std::cout<<j<<", "<<testB[j].as_ulong()<<std::endl;
				}
				*/
			}
			if(cs.constraints[i].c1.terms.size()>1){
				std::vector<FieldT> tempPoly(domain->m, FieldT::zero());
				tempPoly[i] = FieldT::one();
				domain->iFFT(tempPoly);
				//std::vector<FieldT> tempP2(tempPoly);
				//tempP2[0] = tempP2[0]-1;
				//libfqfft::_polynomial_multiplication(tempMul, tempPoly, tempP2);

				std::vector<FieldT> zExtendPoly(domain2->m+1, FieldT::zero());//(cs.num_convol_outputs(0)*((domain->m)*4+1)+1), FieldT::zero());
				//std::cout<<"d2m = "<<domain2->m<<std::endl;
				for(size_t j=1;j<cs.constraints[i].c1.terms.size();j++){
					 size_t ttt = (cs.constraints[i].c1.terms[j].coeff.as_ulong()-1)*(2*(domain->m)+1);
					zExtendPoly[(cs.constraints[i].c1.terms[j].coeff.as_ulong()-1)*(2*(domain->m)+1)] = full_variable_assignment[cs.constraints[i].c1.terms[j].index -1];
					// std::cout<<ttt<<", var["<<cs.constraints[i].c1.terms[j].index -1<<"] = "<<full_variable_assignment[cs.constraints[i].c1.terms[j].index -1].as_ulong()<<std::endl;
				}
				/*
				 std::cout<<"zExtendPoly : ";
				 for(size_t j=0; j<zExtendPoly.size();j++){
				 	std::cout<<zExtendPoly[j].as_ulong()<< "*x^"<<j<<" ";
				 }
				
				 for(auto x : zExtendPoly){
				 	f3 += x;
				 }
				 std::cout<<"f3 :"<<f3.as_ulong()<<std::endl;
				 std::cout<<"temp : ";
				 FieldT tem = FieldT::one();
				 FieldT sumt = FieldT::zero();
				 for(size_t j=0; j<tempPoly.size();j++){
				 	std::cout<<tempPoly[j].as_ulong()<< "*x^"<<j<<" ";
				 	sumt += tempPoly[j]*tem;
				 	tem *= t;
				 }
				 std::cout<<std::endl;
				 std::cout<<"L(t) : "<<u[i].as_ulong()<<", sumt : "<<sumt.as_ulong()<<std::endl;
				 */
				std::vector<FieldT> mulResult(1, FieldT::zero());
				libfqfft::_polynomial_multiplication(mulResult, tempPoly, zExtendPoly);

				//FieldT f3 = FieldT::zero();
				for(size_t j=0;j<mulResult.size();j++){
					if(j < aC.size()){
						aC[j] += mulResult[j];
					}
					else{
						aC.push_back(mulResult[j]);
					}
				}
				/*
				std::vector<FieldT> testC(aC);
				testC.resize(domain4->m, FieldT::zero());
				domain4->FFT(testC);
				for(size_t j=0;j<testC.size();j++){
					std::cout<<j<<", "<<testC[j].as_ulong()<<std::endl;
				}
				*/
			}
			 //std::cout<<"f1f2-f3 = "<<(f1*f2-f3).as_ulong()<<std::endl;
		}
		// r1cs poly * wire
		/*
		for(size_t i=0; i<num_variables;i++){
			for(size_t j=0;j< V[i].size();j++){
				aA[j] += V[i][j] * wires[i];
			}
			for(size_t j=0;j< W[i].size();j++){
				aB[j] += W[i][j] * wires[i];
			}
			for(size_t j=0;j< Y[i].size();j++){
				aC[j] += Y[i][j] * wires[i];
			}
		}
		*/

		/*
		//DEBUG Calculate V(t), W(t), Y(t)
		FieldT resA3 = FieldT::zero();
		FieldT resB3 = FieldT::zero();
		FieldT resC3 = FieldT::zero();
		//res = poly(t) using lagrange 
		for (size_t i = 0; i < domain2->m; i++)
		{
			resA3 += aA[i] * u[i];
			resB3 += aB[i] * u[i];
			resC3 += aC[i] * u[i];
		}

		//A, B, C  point -> poly
		aA.resize(domain2->m, FieldT::zero());
		aB.resize(domain2->m, FieldT::zero());
		aC.resize(domain2->m, FieldT::zero());
		domain2->iFFT(aA);
		domain2->iFFT(aB);
		domain2->iFFT(aC);
		*/
		 std::vector<FieldT> Ht;
		//FieldT t = FieldT(2);
		 Ht.reserve(domain4->m);
		 FieldT ti = FieldT::one();
		 //std::cout<<"t : "<<t.as_ulong()<<std::endl;
		 for(size_t i=0; i<domain4->m;i++){
		 	Ht.emplace_back(ti);
		 	ti*= t;
		 }
		 /*
		 FieldT sum = FieldT::zero();
		 for(size_t i=0;i<aA.size();i++){
		 	sum += aA[i]*Ht[i];
		 }

		 FieldT sum2 = FieldT::zero();
		 for(size_t i=0;i<aB.size();i++){
		 	sum2 += aB[i]*Ht[i];
		 }


		 FieldT sum3 = FieldT::zero();
		 for(size_t i=0;i<aC.size();i++){
		 	sum3 += aC[i]*Ht[i];
		 }
		 std::cout<<"aA :"<<(sum).as_ulong()<<", aB :"<<(sum2).as_ulong()<<", aC :"<<(sum3).as_ulong()<<std::endl;
		 FieldT ABC =sum*sum2 - sum3;
		 */
		//  domain4->iFFT(tempMul);
	
		//  for(size_t j=0;j<domain4->m;j++){
		// 	ABC = ABC + (tempMul[j]*Ht[j]);
		//  }
		// std::cout<<"aA*aB-aC : "<<(ABC).as_ulong()<<std::endl;
		/*
		 //ABC poly -> point
		 domain2->FFT(aA);
		 domain2->FFT(aB);
		 domain2->FFT(aC);
		 //END DEBUG
		*/
		libff::leave_block("calcuate A, B, C for z");

		libff::enter_block("Compute ZK-patch");
 		std::vector<FieldT> coefficients_for_H(domain4->m+1, FieldT::zero());
#ifdef MULTICORE
#pragma omp parallel for
#endif
		/* add coefficients of the polynomial (d2*A + d1*B - d3) + d1*d2*Z */
		for (size_t i = 0; i < domain->m; ++i)
		{
			coefficients_for_H[i] = d2*aA[i] + d1*aB[i];
		}
		coefficients_for_H[0] -= d3;
		domain4->add_poly_Z(d1*d2, coefficients_for_H);
		libff::leave_block("Compute ZK-patch");


		libff::enter_block("A, B, C coset FFT");
		aA.resize(domain4->m, FieldT::zero());
		aB.resize(domain4->m, FieldT::zero());
		aC.resize(domain4->m, FieldT::zero());
		
		/*
		 std::vector<FieldT> tempaA(aA);
		 std::vector<FieldT> tempaB(aB);
		 std::vector<FieldT> tempaC(aC);
		 domain4->FFT(tempaA);
		 domain4->FFT(tempaB);
		 domain4->FFT(tempaC);
		 std::cout<<"A*B-C point"<<std::endl;
		 for(size_t j=0;j<domain4->m;j++){
		 	std::cout<<j<<", "<<tempaA[j].as_ulong()<<", "<<tempaB[j].as_ulong()<<", "<<tempaC[j].as_ulong()<<", "<<(tempaA[j]*tempaB[j]-tempaC[j]).as_ulong()<<std::endl;
		 }
		 */
		 

		//A(x) => A(gx) same B, C
		libff::enter_block("Compute evaluation of polynomial A on set T");
		domain4->cosetFFT(aA, FieldT::multiplicative_generator);
		libff::leave_block("Compute evaluation of polynomial A on set T");

		libff::enter_block("Compute evaluation of polynomial B on set T");
		domain4->cosetFFT(aB, FieldT::multiplicative_generator);
		libff::leave_block("Compute evaluation of polynomial B on set T");

		libff::enter_block("Compute evaluation of polynomial C on set T");
		domain4->cosetFFT(aC, FieldT::multiplicative_generator);
		libff::leave_block("Compute evaluation of polynomial C on set T");
		libff::leave_block("A, B, C coset FFT");


		libff::enter_block("Compute evaluation of polynomial H on set T");
		std::vector<FieldT> &H_tmp = aA; // can overwrite aA because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
		for (size_t i = 0; i < domain4->m; ++i)
		{
			H_tmp[i] = aA[i]*aB[i] - aC[i];
		}
		std::vector<FieldT>().swap(aB); // destroy aB

		libff::enter_block("Divide by Z on set T");
		///TODO change Z : domain Z is change to domain4 Z
		std::vector<FieldT> Zpoly(domain->m+1, FieldT::zero());
		domain->add_poly_Z(FieldT::one(), Zpoly);
		/*
		 FieldT zsum = FieldT::zero();
		 for(auto x : Zpoly){
		 	zsum += x;
		 }
		 std::cout<<"Z(1) = "<<zsum.as_ulong()<<std::endl;

		 zsum = FieldT::zero();

		 for(size_t i=0;i<Zpoly.size();i++){
		 	zsum += Zpoly[i]*Ht[i];
		 }
		 std::cout<<"Z("<<t.as_ulong()<<") = "<<zsum.as_ulong()<<std::endl;
		 const FieldT Zt = domain->compute_vanishing_polynomial(t);
		 std::cout<<"Zt : "<<Zt.as_ulong()<<std::endl;
		*/

		Zpoly.resize(domain4->m, FieldT::zero());
		
		/*
		 std::vector<FieldT> tempZ(Zpoly);
		 domain4->FFT(tempZ);
		 std::cout<<"tempZ point"<<std::endl;
		 for(size_t j=0; j<domain4->m;j++){
		 	std::cout<<j<<", "<<tempZ[j].as_ulong()<<std::endl;
		 }
		 */

		domain4->cosetFFT(Zpoly, FieldT::multiplicative_generator);
		// libff::leave_block("claculate Zt point");
		//(A*B-C(gx) point) / (Z(gx) point)
		for(size_t i=0;i<domain4->m;i++){
			H_tmp[i] = H_tmp[i] * Zpoly[i].inverse();
		}
		domain4->icosetFFT(H_tmp, FieldT::multiplicative_generator);
		libff::leave_block("Divide by Z on set T");

		libff::leave_block("Compute evaluation of polynomial H on set T");

		libff::enter_block("Compute sum of H and ZK-patch");

#ifdef MULTICORE
#pragma omp parallel for
#endif
		for (size_t i = 0; i < domain4->m; ++i)
		{
			coefficients_for_H[i] += H_tmp[i];
		}
		libff::leave_block("Compute sum of H and ZK-patch");

		/*
		 FieldT Hsum = FieldT::zero();
		 for(size_t i=0;i<coefficients_for_H.size();i++){
		 	Hsum += coefficients_for_H[i]*Ht[i];
		 }
		 */
		 //std::cout<<"H("<<t.as_ulong()<<") = "<<Hsum.as_ulong()<<std::endl;
		 //std::cout<<"ABC - H*Z = "<<(ABC - zsum*Hsum).as_ulong()<<std::endl;






	// 	/////////DEBUG////////

	// std::vector<FieldT> At, Bt, Ct;



	// At.resize(cs.num_variables()+1, FieldT::zero());
	// Bt.resize(cs.num_variables()+1, FieldT::zero());
	// Ct.resize(cs.num_variables()+1, FieldT::zero());
 	// //At, Bt, Ct = poly(t);


	// for (size_t i = 0; i < cs.num_constraints(); ++i){
	// 	if(cs.constraints[i].a1.terms.size()>1){
	// 		for (size_t j = 1; j < cs.constraints[i].a1.terms.size(); ++j)
	// 		{
	// 			At[cs.constraints[i].a1.terms[j].index] +=
	// 				u[i]*Ht[(cs.constraints[i].a1.terms[j].coeff.as_ulong())*(4*domain->m+1)];//cs.constraints[i].a.terms[j].coeff;
	// 		}

	// 		for (size_t j = 1; j < cs.constraints[i].b1.terms.size(); ++j)
	// 		{
	// 			Bt[cs.constraints[i].b1.terms[j].index] +=
	// 				u[i]*Ht[(cs.constraints[i].b1.terms[j].coeff.as_ulong())*(4*domain->m+1)];

	// 		}
	// 		for (size_t j = 1; j < cs.constraints[i].c1.terms.size(); ++j)
	// 		{
	// 			Ct[cs.constraints[i].c1.terms[j].index] +=
	// 				u[i]*Ht[(cs.constraints[i].c1.terms[j].coeff.as_ulong())*(4*domain->m+1)];
	// 		}
	// 	}
	// }
	// FieldT sumAA = FieldT::zero();
	// FieldT sumBB = FieldT::zero();
	// FieldT sumCC = FieldT::zero();
	// sumAA = At[0];
	// sumBB = Bt[0];
	// sumCC = Ct[0];
	// for(size_t i=1;i<cs.num_variables()+1;i++){
	// 	sumAA += full_variable_assignment[i-1] * At[i];
	// 	sumBB += full_variable_assignment[i-1] * Bt[i];
	// 	sumCC += full_variable_assignment[i-1] * Ct[i];
	// }
	// std::cout<<"A1 : "<<sum.as_ulong()<<", A2 : "<<sumAA.as_ulong()<<std::endl;
	// std::cout<<"B1 : "<<sum2.as_ulong()<<", B2 : "<<sumBB.as_ulong()<<std::endl;
	// std::cout<<"C1 : "<<sum3.as_ulong()<<", C2 : "<<sumCC.as_ulong()<<std::endl;
	// std::cout<<"AB-C : "<<(sumAA*sumBB-sumCC).as_ulong()<<", CC :"<<(ABC).as_ulong()<<std::endl;

	// 	//////////DEBUGEND//////




		return qap_witness<FieldT>(cs.num_variables(),
				domain4->m,
				cs.num_inputs(),
				d1,
				d2,
				d3,
				full_variable_assignment,
				std::move(coefficients_for_H));

	}
	/*
		template<typename FieldT>
	vector<FieldT> imagePolyCalculate(vector<vector<vector<FieldT>>> stPoly,
			size_t num_constraints, size_t num_s_order, size_t num_t_order,
			shared_ptr<libfqfft::evaluation_domain<FieldT>> domain
			)
	{
		const shared_ptr<libfqfft::evaluation_domain<FieldT>>domain = libfqfft::get_evaluation_domain<FieldT>(num_constraints);
		size_t order = ((((((num_constraints-1)*4+1)*num_s_order)*4+1)*num_t_order), FieldT::zero());
	//	vector<vecotr<<FieldT>> extPoly(num_constraints, vector<FieldT>(order,FieldT::zero()));
		vector<FieldT> tempST((order), FieldT::zero());
		for(size_t i =0; i<num_constraints;i++){
			for(size_t j=0; j<num_s_order;j++){
				for(size_t k=0; k<num_t_order;k++){
					if(stPoly[i][j][k] != FieldT::zero()){
						if(j!=0 && k !=0){
							vector<FieldT> tempPoly(num_constraints, FieldT::zero());
							size_t s_order = j*(4*(num_constraints-1)+1); // j*(4d_x + 1)
							size_t k_order = k*(4*(num_s_order-1)+1); //k*(4d_s+1)
							tempPoly[i] = FieldT::one();
							domain->iFFT(tempPoly[i]);
							size_t new_order = s_order+k_order;
							for(size_t x = 0; x<num_constraints;x++){
								tempST[x*new_order] += tempPoly[x];
							}
						}
					}
				}
			}
		}
		return tempST;
	}
	*/
} // libsnark

#endif // R1CS_TO_QP_TCC_
