/** @file
 *****************************************************************************
 Test program that exercises the ppzkSNARK (first generator, then
 prover, then verifier) on a synthetic R1CS instance.

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
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/examples/run_r1cs_ppzksnark.hpp>
//#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/examples/R1CS.hpp>

#include<stdlib.h>


using namespace libsnark;

using namespace std;

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

	template<typename FieldT>
vector<FieldT> zxPolycalculate(vector<vector<FieldT>> zxPoly, 
		size_t num_constraints, size_t num_z_order,
		shared_ptr<libfqfft::evaluation_domain<FieldT>> domain
		)
{
	//print points
	/*
	   for(size_t i=0;i<num_constraints;i++){
	   cout<<i<<" : ";
	   for(size_t j=0; j<num_z_order;j++){
	   cout<<zxPoly[i][j].as_ulong()<< " s^"<<j<< " + ";
	   }
	   cout<<endl;
	   }
	 */

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
	//cout<<"zEx(x) : ";
	//print_Poly(zExtendPoly);

	vector<FieldT> mulResult(1, FieldT::zero());
	libfqfft::_polynomial_multiplication(mulResult, tempPoly[0], zExtendPoly);
	//cout<<"v * z : "<<endl;
	//print_Poly(mulResult);

	//point add
	for(size_t i=0; i<num_constraints;i++){
		FieldT temp = FieldT::zero();
		for(size_t j=1; j<num_constraints;j++){
			temp+=tempPoly[j][i];
		}
		mulResult[i] += temp;
	}

	//cout<<"mulResult+v Poly"<<endl;
	//print_Poly(mulResult);

	/*
	   mulResult.resize(domain2->m, FieldT::zero());
	   domain2->FFT(mulResult);
	   cout<<"mulResult +v point"<<endl;
	   print_point(mulResult);
	 */
	return mulResult;
}

template<typename FieldT>
void evaluate_mqap(){

	size_t num_variabels = 15;
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

	 */

	shared_ptr<libfqfft::evaluation_domain<FieldT> > domain = libfqfft::get_evaluation_domain<FieldT>(num_constraints);
	shared_ptr<libfqfft::evaluation_domain<FieldT> > domain2 = libfqfft::get_evaluation_domain<FieldT>(num_z_order*((num_constraints-1)*4+1) + (num_constraints-1) +1);

	vector<vector<vector<FieldT>>> mr1csV (num_variabels, vector<vector<FieldT>>(num_constraints, vector<FieldT>(num_z_order+1, FieldT::zero())));
	vector<vector<vector<FieldT>>> mr1csW (num_variabels, vector<vector<FieldT>>(num_constraints, vector<FieldT>(num_z_order+1, FieldT::zero())));
	vector<vector<vector<FieldT>>> mr1csY (num_variabels, vector<vector<FieldT>>(num_constraints, vector<FieldT>(num_z_order+1, FieldT::zero())));

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
	vector<vector<FieldT>> V(num_variabels, vector<FieldT>(1, FieldT::zero()));
	for(size_t i=0; i<num_variabels; i++){
		V[i] = zxPolycalculate(mr1csV[i], num_constraints, num_z_order, domain);	
	}


	cout<<"===========Calculate W==========="<<endl;
	vector<vector<FieldT>> W(num_variabels, vector<FieldT>(1, FieldT::zero()));
	for(size_t i=0; i<num_variabels; i++){
		W[i] = zxPolycalculate(mr1csW[i], num_constraints, num_z_order, domain);	
	}

	cout<<"===========Calculate Y==========="<<endl;
	vector<vector<FieldT>> Y(num_variabels, vector<FieldT>(1, FieldT::zero()));
	for(size_t i=0; i<num_variabels; i++){
		Y[i] = zxPolycalculate(mr1csY[i], num_constraints, num_z_order, domain);	
	}
	libff::leave_block("Compute V, W, Y");

	const FieldT t = FieldT::zero();//random_element();

	const FieldT Zt = domain->compute_vanishing_polynomial(t);

	const vector<FieldT> u = domain2->evaluate_all_lagrange_polynomials(t);

	FieldT ckVk = FieldT::zero();
	FieldT ckWk = FieldT::zero();
	FieldT ckYk = FieldT::zero();

	/*
	   FieldT omega = libff::get_root_of_unity<FieldT>(4);
	   for(size_t i=0;i<num_variabels;i++){
	   FieldT sum = FieldT::zero();
	   FieldT sum1 = FieldT::zero();
	   FieldT sum2 = FieldT::zero();
	   for(size_t j=0; j<V[i].size();j++){
	   sum += V[i][j] ;
	   }
	   omega = libff::get_root_of_unity<FieldT>(2);
	   for(size_t j=0; j<V[i].size();j++){
	   sum1 += V[i][j] *(omega^j);
	   }
	   omega = libff::get_root_of_unity<FieldT>(4);
	   for(size_t j=0; j<V[i].size();j++){
	   sum2 += V[i][j] *(omega^j);
	   }
	   cout << "V(0) = "<<sum.as_ulong()<< "\tVt(1) = "<<sum1.as_ulong()<<"\tVt(2) = "<< sum2.as_ulong()<<endl;
	   }
	 */
	FieldT omega = libff::get_root_of_unity<FieldT>(4);
	//omega = omega^16;
	for(size_t i=0; i<num_variabels; i++){

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
	Ht.reserve(domain2->m+1);
	FieldT ti = FieldT::one();
	for(size_t i=0; i<domain2->m+1;i++){
		Ht.emplace_back(ti);
		ti*= t;
	}

	//At, Bt, Ct = poly(t);
	for(size_t i=0; i<num_variabels; i++){
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

	/////kegen qap instance end


	///qap witness start

	vector<FieldT> aA(domain2->m, FieldT::zero());
	vector<FieldT> aB(domain2->m, FieldT::zero());
	vector<FieldT> aC(domain2->m, FieldT::zero());

	//V, W, Y poly -> point
	for(size_t i=0; i<num_variabels; i++){
		V[i].resize(domain2->m, FieldT::zero());
		domain2->FFT(V[i]);
		W[i].resize(domain2->m, FieldT::zero());
		domain2->FFT(W[i]);
		Y[i].resize(domain2->m, FieldT::zero());
		domain2->FFT(Y[i]);
	}

	// r1cs point * wire
	for(size_t i=0; i<domain2->m;i++){
		for(size_t j=0;j< num_variabels;j++){
			aA[i] += V[j][i] * wires[j];
			aB[i] += W[j][i] * wires[j];
			aC[i] += Y[j][i] * wires[j];
		}
	}

	FieldT resA3 = FieldT::zero();
	FieldT resB3 = FieldT::zero();
	FieldT resC3 = FieldT::zero();
	// res = poly(t) using lagrange
	for(size_t i=0;i<domain2->m;i++){
		resA3 += aA[i] * u[i];
		resB3 += aB[i] * u[i];
		resC3 += aC[i] * u[i];
	}
	//cout<<"before iFFT = points"<<endl;
	//print_point(aA);


	// aA point -> poly
	domain2->iFFT(aA);
	//cout<<"after iFFT = poly"<<endl;
	//print_Poly(aA);

	FieldT A0 = FieldT::zero();
	for(size_t i=0;i<domain2->m;i++){
		A0 += aA[i];	
	}
	cout<<"A0 = "<<A0.as_ulong()<<endl;

	domain2->iFFT(aB);
	//print_Poly(aB);

	domain2->iFFT(aC);
	//print_Poly(aC);

	vector<FieldT> tempAB(1, FieldT::zero());
	libfqfft::_polynomial_multiplication(tempAB, aA, aB);
	vector<FieldT> tempABC(1, FieldT::zero());
	libfqfft::_polynomial_subtraction(tempABC, tempAB, aC);
	shared_ptr<libfqfft::evaluation_domain<FieldT> > domain4 = libfqfft::get_evaluation_domain<FieldT>((domain2->m*2));
	tempABC.resize(domain4->m, FieldT::zero());
	domain4->FFT(tempABC);
	cout<<"A*B-C points"<<endl;
	print_point(tempABC);


	FieldT resA2 = FieldT::zero();
	FieldT resB2 = FieldT::zero();
	FieldT resC2 = FieldT::zero();
	// res = poly(t)
	for(size_t i=0;i<domain2->m;i++){
		resA2 += aA[i] * Ht[i];
		resB2 += aB[i] * Ht[i];
		resC2 += aC[i] * Ht[i];
	}
	domain2->cosetFFT(aA, FieldT::multiplicative_generator);
	domain2->cosetFFT(aB, FieldT::multiplicative_generator);
	domain2->cosetFFT(aC, FieldT::multiplicative_generator);

	cout<<"mult gen = "<<FieldT::multiplicative_generator.as_ulong()<<endl;
	//cout<<"after cosetFFT = points"<<endl;
	//print_point(aA);
	/*

	cout<<"after iFFT = poly"<<endl;
	domain2->iFFT(aA);
	print_Poly(aA);
	domain2->iFFT(aB);
	print_Poly(aB);
	domain2->iFFT(aC);
	print_Poly(aC);


	domain2->FFT(aA);
	domain2->FFT(aB);
	domain2->FFT(aC);
	*/
	vector<FieldT> &H_tmp = aA;

	for(size_t i=0;i<domain2->m;i++){
		H_tmp[i] = (aA[i] * aB[i]) - aC[i];
	}
	domain->divide_by_Z_on_coset(H_tmp);

	domain2->icosetFFT(H_tmp, FieldT::multiplicative_generator);
	cout << "H poly : ";
	print_Poly(H_tmp);
	//cout<<"H0 + 2 = "<<(H_tmp[0]+FieldT(2)).as_ulong()<<endl;

	cout<<"H points"<<endl;
	domain2->FFT(H_tmp);
	print_point(H_tmp);

	domain2->iFFT(H_tmp);

	vector<FieldT> Zpoly(1, FieldT::one());
	vector<FieldT> xw(2, FieldT::one());
	omega = libff::get_root_of_unity<FieldT>(4);
	cout<<"omega = "<<omega.as_ulong()<<endl;
	FieldT mulTemp = FieldT::one();
	for(size_t i=0;i<domain->m;i++){
		xw[0] = -mulTemp;
		libfqfft::_polynomial_multiplication(Zpoly, Zpoly, xw);
		mulTemp *= omega;
	}
	cout<<"Zpoly : ";
	print_Poly(Zpoly);


	vector<FieldT> gcd(1, FieldT::zero());
	vector<FieldT> uu(1, FieldT::zero());
	vector<FieldT> vv(1, FieldT::zero());
	libfqfft::_polynomial_xgcd(tempABC, Zpoly, gcd, uu, vv);
	cout<<"gcd : ";
	print_Poly(gcd);
	cout<<"u : ";
	print_Poly(uu);
	cout<<"v : ";
	print_Poly(vv);

	FieldT Ztt = FieldT::zero();
	for(size_t i=0;i<domain->m;i++){
		Ztt += Zpoly[i]*Ht[i];
	}
	cout<<"Zt = "<<Ztt.as_ulong()<<endl;
	Ztt = FieldT::zero();
	omega = libff::get_root_of_unity<FieldT>(4);
	for(size_t i=0;i<domain->m;i++){
		Ztt += Zpoly[i]*(omega^i);
	}
	cout<<"Z(w) = "<<Ztt.as_ulong()<<endl;

	vector<FieldT> q, r;
	libfqfft::_polynomial_division(q, r, tempABC, Zpoly);
	cout << "A*B-C / Z = ";
	print_Poly(q);

	cout<<"remain = ";
	print_Poly(r);

	FieldT H3 = FieldT::zero();
	for(size_t i=0;i<q.size();i++){
		H3 += q[i] * Ht[i];
	}

	cout<<"H3 = "<<H3.as_ulong()<<endl;


	cout<<"Z points "<<endl;
	shared_ptr<libfqfft::evaluation_domain<FieldT> > domain3 = libfqfft::get_evaluation_domain<FieldT>(domain->m+1);
	domain3->FFT(Zpoly);
	print_point(Zpoly);

	domain3->iFFT(Zpoly);
	Zpoly.resize(domain2->m, FieldT::zero());


	FieldT resA = FieldT::zero();
	FieldT resB = FieldT::zero();
	FieldT resC = FieldT::zero();
	for(size_t i=0;i< num_variabels;i++){
		resA += At[i] * wires[i];
		resB += Bt[i] * wires[i];
		resC += Ct[i] * wires[i];
	}
	FieldT resH = FieldT::zero();
	for(size_t i=0;i<domain2->m;i++){
		resH += H_tmp[i] * Ht[i];
	}

	if(((resA * resB) - resC) == (Zt * resH)) cout<< "ok2"<<endl;
	else cout <<"no2"<<endl;

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
	cout<<"-Z = "<<(-Zt).as_ulong()<<endl;

}

template<typename ppT>
void testPoly(){
	evaluate_mqap<libff::Fr<ppT>>();
}

	template<typename ppT>
void test_r1cs_ppzksnark(size_t num_constraints,
		size_t input_size)
{
	libff::print_header("(enter) Test R1CS ppzkSNARK");

	const bool test_serialization = true;
	r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_binary_input<libff::Fr<ppT> >(num_constraints, input_size);
	const bool bit = run_r1cs_ppzksnark<ppT>(example, test_serialization);
	assert(bit);

	libff::print_header("(leave) Test R1CS ppzkSNARK");
}

int main()
{
	default_r1cs_ppzksnark_pp::init_public_params();
	libff::start_profiling();

	testPoly<default_r1cs_ppzksnark_pp>();
	//test_r1cs_ppzksnark<default_r1cs_ppzksnark_pp>(10000, 100);
}

