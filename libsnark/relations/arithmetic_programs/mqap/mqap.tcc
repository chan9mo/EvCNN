/** @file
*****************************************************************************

Implementation of interfaces for a QAP ("Quadratic Arithmetic Program").

See qap.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef MQAP_TCC_
#define MQAP_TCC_

#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libfqfft/evaluation_domain/evaluation_domain.hpp>

namespace libsnark {

template<typename FieldT>
mqap_instance<FieldT>::mqap_instance(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
                                   const size_t num_variables,
                                   const size_t degree,
                                   const size_t degree2,
                                   const size_t num_inputs,
                                   const size_t num_kernels,
                                   const std::vector<std::map<size_t, FieldT> > &A_in_Lagrange_basis,
                                   const std::vector<std::map<size_t, FieldT> > &B_in_Lagrange_basis,
                                   const std::vector<std::map<size_t, FieldT> > &C_in_Lagrange_basis) :
    num_variables_(num_variables),
    degree_(degree),
    degree_2_(degree2),
    num_inputs_(num_inputs),
    num_kernels_(num_kernels),
    domain(domain),
    A_in_Lagrange_basis(A_in_Lagrange_basis),
    B_in_Lagrange_basis(B_in_Lagrange_basis),
    C_in_Lagrange_basis(C_in_Lagrange_basis)
{
}

template<typename FieldT>
mqap_instance<FieldT>::mqap_instance(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
                                   const size_t num_variables,
                                   const size_t degree,
                                   const size_t degree2,
                                   const size_t num_inputs,
                                   const size_t num_kernels,
                                   std::vector<std::map<size_t, FieldT> > &&A_in_Lagrange_basis,
                                   std::vector<std::map<size_t, FieldT> > &&B_in_Lagrange_basis,
                                   std::vector<std::map<size_t, FieldT> > &&C_in_Lagrange_basis) :
    num_variables_(num_variables),
    degree_(degree),
    degree_2_(degree2),
    num_inputs_(num_inputs),
    num_kernels_(num_kernels),
    domain(domain),
    A_in_Lagrange_basis(std::move(A_in_Lagrange_basis)),
    B_in_Lagrange_basis(std::move(B_in_Lagrange_basis)),
    C_in_Lagrange_basis(std::move(C_in_Lagrange_basis))
{
}

template<typename FieldT>
size_t mqap_instance<FieldT>::num_variables() const
{
    return num_variables_;
}

template<typename FieldT>
size_t mqap_instance<FieldT>::degree() const
{
    return degree_;
}

template<typename FieldT>
size_t mqap_instance<FieldT>::num_inputs() const
{
    return num_inputs_;
}

template<typename FieldT>
bool mqap_instance<FieldT>::is_satisfied(const mqap_witness<FieldT> &witness) const
{
    const FieldT t = FieldT::random_element();

    std::vector<FieldT> At(this->num_variables()+1, FieldT::zero());
    std::vector<FieldT> Bt(this->num_variables()+1, FieldT::zero());
    std::vector<FieldT> Ct(this->num_variables()+1, FieldT::zero());
    std::vector<FieldT> Ht(this->degree()+1);

    const FieldT Zt = this->domain->compute_vanishing_polynomial(t);

    const std::vector<FieldT> u = this->domain->evaluate_all_lagrange_polynomials(t);

    for (size_t i = 0; i < this->num_variables()+1; ++i)
    {
        for (auto &el : A_in_Lagrange_basis[i])
        {
            At[i] += u[el.first] * el.second;
        }

        for (auto &el : B_in_Lagrange_basis[i])
        {
            Bt[i] += u[el.first] * el.second;
        }

        for (auto &el : C_in_Lagrange_basis[i])
        {
            Ct[i] += u[el.first] * el.second;
        }
    }

    FieldT ti = FieldT::one();
    for (size_t i = 0; i < this->degree()+1; ++i)
    {
        Ht[i] = ti;
        ti *= t;
    }

    const mqap_instance_evaluation<FieldT> eval_mqap_inst(this->domain,
                                                        this->num_variables(),
                                                        this->degree(),
                                                        this->num_inputs(),
                                                        t,
                                                        std::move(At),
                                                        std::move(Bt),
                                                        std::move(Ct),
                                                        std::move(Ht),
                                                        Zt);
    return eval_mqap_inst.is_satisfied(witness);
}

template<typename FieldT>
mqap_instance_evaluation<FieldT>::mqap_instance_evaluation(
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain2,
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain4,
                                                         const size_t num_variables,
                                                         const size_t degree,
                                                         const size_t degree2,
                                                         const size_t degree4,
                                                         const size_t num_inputs,
                                                         const size_t num_kernels,
                                                         const FieldT &t,
                                                         const std::vector<FieldT> &At,
                                                         const std::vector<FieldT> &Bt,
                                                         const std::vector<FieldT> &Ct,
                                                         const std::vector<FieldT> &Ht,
                                                         const FieldT &Zt) :
    num_variables_(num_variables),
    degree_(degree),
    degree2_(degree2),
    degree4_(degree4),
    num_inputs_(num_inputs),
    num_kernels_(num_kernels),
    domain(domain),
    domain2(domain2),
    domain4(domain4),
    t(t),
    At(At),
    Bt(Bt),
    Ct(Ct),
    Ht(Ht),
    Zt(Zt)
{
}

template<typename FieldT>
mqap_instance_evaluation<FieldT>::mqap_instance_evaluation(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain2,
		const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain4,
                                                         const size_t num_variables,
                                                         const size_t degree,
                                                         const size_t degree2,
														 const size_t degree4,
                                                         const size_t num_inputs,
                                                         const size_t num_kernels,
                                                         const FieldT &t,
                                                         std::vector<FieldT> &&At,
                                                         std::vector<FieldT> &&Bt,
                                                         std::vector<FieldT> &&Ct,
                                                         std::vector<FieldT> &&Ht,
                                                         const FieldT &Zt) :
    num_variables_(num_variables),
    degree_(degree),
    degree2_(degree2),
    degree4_(degree4),
    num_inputs_(num_inputs),
    num_kernels_(num_kernels),
    domain(domain),
    t(t),
    At(std::move(At)),
    Bt(std::move(Bt)),
    Ct(std::move(Ct)),
    Ht(std::move(Ht)),
    Zt(Zt)
{
}

template<typename FieldT>
size_t mqap_instance_evaluation<FieldT>::num_variables() const
{
    return num_variables_;
}

template<typename FieldT>
size_t mqap_instance_evaluation<FieldT>::degree() const
{
    return degree_;
}

template<typename FieldT>
size_t mqap_instance_evaluation<FieldT>::degree2() const
{
    return degree2_;
}

template<typename FieldT>
size_t mqap_instance_evaluation<FieldT>::degree4() const
{
    return degree4_;
}

template<typename FieldT>
size_t mqap_instance_evaluation<FieldT>::num_inputs() const
{
    return num_inputs_;
}

template<typename FieldT>
size_t mqap_instance_evaluation<FieldT>::num_kernels() const
{
    return num_kernels_;
}

template<typename FieldT>
bool mqap_instance_evaluation<FieldT>::is_satisfied(const mqap_witness<FieldT> &witness) const
{

    if (this->num_variables() != witness.num_variables())
    {
		std::cout<<"wit num error"<<std::endl;
        return false;
    }

    if (this->degree() != witness.degree())
    {
		std::cout<<"degree error"<<std::endl;
        return false;
    }

    if (this->degree2() != witness.degree2())
    {
		std::cout<<"degree2 error"<<std::endl;
        return false;
    }

    if (this->num_inputs() != witness.num_inputs())
    {
		std::cout<<"input num error"<<std::endl;
        return false;
    }

    if (this->num_kernels() != witness.num_kernels())
    {
		std::cout<<"input kernel error"<<std::endl;
        return false;
    }

    size_t new_degree = this->degree2()*((this->degree()-1)*4+1)+(this->degree()-1)+1;
    if ((new_degree*2) != witness.coefficients_for_ABCs.size())
    {
		std::cout<<"wit abc size error"<<std::endl;
        return false;
    }

    if ((new_degree*2) != witness.coefficients_for_H.size())
    {
		std::cout<<"H size error"<<std::endl;
        //return false;
    }

    ///TODO think At, Ht size
    if (this->At.size() != this->num_variables()+1 || this->Bt.size() != this->num_variables()+1 || this->Ct.size() != this->num_variables()+1)
    {
		std::cout<<"a or b or c size error"<<std::endl;
        //return false;
    }

    if (this->Ht.size() != this->degree()+1)
    {
		std::cout<<"h degree size error"<<std::endl;
        //return false;
    }

    if (this->Zt != this->domain->compute_vanishing_polynomial(this->t))
    {
		std::cout<<"zt size error"<<std::endl;
        //return false;
    }

    FieldT ans_A = this->At[0]; //+ witness.d1*this->Zt;
    FieldT ans_B = this->Bt[0]; //+ witness.d2*this->Zt;
    FieldT ans_C = this->Ct[0]; //+ witness.d3*this->Zt;
    FieldT ans_H = FieldT::zero();

	/*origin
    ans_A = ans_A + libff::inner_product<FieldT>(this->At.begin()+1,
                                                 this->At.begin()+1+this->num_variables(),
                                                 witness.coefficients_for_ABCs.begin(),
                                                 witness.coefficients_for_ABCs.begin()+this->num_variables());
    ans_B = ans_B + libff::inner_product<FieldT>(this->Bt.begin()+1,
                                                 this->Bt.begin()+1+this->num_variables(),
                                                 witness.coefficients_for_ABCs.begin(),
                                                 witness.coefficients_for_ABCs.begin()+this->num_variables());
    ans_C = ans_C + libff::inner_product<FieldT>(this->Ct.begin()+1,
                                                 this->Ct.begin()+1+this->num_variables(),
                                                 witness.coefficients_for_ABCs.begin(),
                                                 witness.coefficients_for_ABCs.begin()+this->num_variables());
    ans_H = ans_H + libff::inner_product<FieldT>(this->Ht.begin(),
                                                 this->Ht.begin()+this->degree()+1,
                                                 witness.coefficients_for_H.begin(),
                                                 witness.coefficients_for_H.begin()+this->degree()+1);
	*/
    ans_A = libff::inner_product<FieldT>(this->At.begin(),
                                                 this->At.begin()+this->num_variables(),
                                                 witness.coefficients_for_ABCs.begin(),
                                                 witness.coefficients_for_ABCs.begin()+this->num_variables());
    ans_B = libff::inner_product<FieldT>(this->Bt.begin(),
                                                 this->Bt.begin()+this->num_variables(),
                                                 witness.coefficients_for_ABCs.begin(),
                                                 witness.coefficients_for_ABCs.begin()+this->num_variables());
    ans_C = libff::inner_product<FieldT>(this->Ct.begin(),
                                                 this->Ct.begin()+this->num_variables(),
                                                 witness.coefficients_for_ABCs.begin(),
                                                 witness.coefficients_for_ABCs.begin()+this->num_variables());
    ans_H = libff::inner_product<FieldT>(this->Ht.begin(),
                                                 this->Ht.begin()+this->degree(),
                                                 witness.coefficients_for_H.begin(),
                                                 witness.coefficients_for_H.begin()+this->degree());
	/*
	for(size_t i=0; i<this->degree();i++){
		if(this->Ht[i] != FieldT::zero()){
			std::cout<<"H["<<i<<"]="<<witness.coefficients_for_H[i].as_ulong()<<"*Ht["<<i<<"]="<<this->Ht[i].as_ulong()<<(witness.coefficients_for_H[i]*this->Ht[i]).as_ulong()<<std::endl;
		}
	}
	std::cout<<"degree = "<<this->degree()<<std::endl;
	std::cout <<"wit coeff H[0] = "<<witness.coefficients_for_H[0].as_ulong()<<std::endl;
	std::cout<<"AB-C = "<<(ans_A*ans_B-ans_C).as_ulong()<<std::endl;
	std::cout<<"H = "<<ans_H.as_ulong()<<std::endl;
	std::cout<<"Z = "<<this->Zt.as_ulong()<<std::endl;
	std::cout<<"HZ = "<<(this->Zt*ans_H).as_ulong()<<std::endl;
	*/
    if (ans_A * ans_B - ans_C != ans_H * this->Zt)
    {
		std::cout<<"divise error"<<std::endl;
        return false;
    }

    return true;
}

template<typename FieldT>
mqap_witness<FieldT>::mqap_witness(const size_t num_variables,
                                 const size_t degree,
                                 const size_t degree2,
                                 const size_t num_inputs,
                                 const size_t num_kernels,
                                 const FieldT &d1,
                                 const FieldT &d2,
                                 const FieldT &d3,
                                 const std::vector<FieldT> &coefficients_for_ABCs,
                                 const std::vector<FieldT> &coefficients_for_H) :
    num_variables_(num_variables),
    degree_(degree),
    degree2_(degree2),
    num_inputs_(num_inputs),
    num_kernels_(num_kernels),
    d1(d1),
    d2(d2),
    d3(d3),
    coefficients_for_ABCs(coefficients_for_ABCs),
    coefficients_for_H(coefficients_for_H)
{
}

template<typename FieldT>
mqap_witness<FieldT>::mqap_witness(const size_t num_variables,
                                 const size_t degree,
                                 const size_t num_inputs,
                                 const FieldT &d1,
                                 const FieldT &d2,
                                 const FieldT &d3,
                                 const std::vector<FieldT> &coefficients_for_ABCs,
                                 std::vector<FieldT> &&coefficients_for_H) :
    num_variables_(num_variables),
    degree_(degree),
    degree2_(degree2),
    num_inputs_(num_inputs),
    num_kernels_(num_kernels),
    d1(d1),
    d2(d2),
    d3(d3),
    coefficients_for_ABCs(coefficients_for_ABCs),
    coefficients_for_H(std::move(coefficients_for_H))
{
}


template<typename FieldT>
size_t mqap_witness<FieldT>::num_variables() const
{
    return num_variables_;
}

template<typename FieldT>
size_t mqap_witness<FieldT>::degree() const
{
    return degree_;
}

template<typename FieldT>
size_t mqap_witness<FieldT>::degree2() const
{
    return degree2_;
}

template<typename FieldT>
size_t mqap_witness<FieldT>::num_inputs() const
{
    return num_inputs_;
}

template<typename FieldT>
size_t mqap_witness<FieldT>::num_kernels() const
{
    return num_kernels_;
}

	template<typename FieldT>
mqap_instance_evaluation<FieldT> mr1cs_to_qap_instance_map_with_evaluation(
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
	// check r1cs
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
	 */
 	vector<FieldT> At(domain2->m, FieldT::zero());
	vector<FieldT> Bt(domain2->m, FieldT::zero());
	vector<FieldT> Ct(domain2->m, FieldT::zero());
	vector<FieldT> Ht;
	Ht.reserve(domain4->m);
	FieldT ti = FieldT::one();
	for(size_t i=0; i<domain4->m;i++){
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
 	return mqap_instance_evaluation<FieldT>(domain2,
			num_variables,
			domain2->m,
			num_inputs,
			t,
			move(At),
			move(Bt),
			move(Ct),
			move(Ht),
			Zt);
 }

} // libsnark

//testtesttest
//testtest

#endif // mqap_TCC_
