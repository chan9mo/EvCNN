#ifndef R1CS_legosnark_TCC_
#define R1CS_legosnark_TCC_

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>

#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#ifdef MULTICORE
#include <omp.h>
#endif

#include <libsnark/knowledge_commitment/kc_multiexp.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>


namespace libsnark
{

template <typename ppT>
bool r1cs_legosnark_proving_key<ppT>::operator==(const r1cs_legosnark_proving_key<ppT> &other) const
{
    return (this->F_g1 == other.F_g1 &&
            this->T_g1 == other.T_g1 &&
            this->R_g1 == other.R_g1);
}
template <typename ppT>
std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_proving_key<ppT> &pk)
{
    out << pk.F_g1 << OUTPUT_NEWLINE;
    out << pk.T_g1 << OUTPUT_NEWLINE;
    out << pk.R_g1 << OUTPUT_NEWLINE;

    return out;
}
template <typename ppT>
std::istream &operator>>(std::istream &in, const r1cs_legosnark_proving_key<ppT> &pk)
{
    in >> pk.F_g1;
    libff::consume_OUTPUT_NEWLINE(in);
    in >> pk.T_g1;
    libff::consume_OUTPUT_NEWLINE(in);
    in >> pk.R_g1;
    libff::consume_OUTPUT_NEWLINE(in);

    return in;
}

template <typename ppT>
bool r1cs_legosnark_verification_key<ppT>::operator==(const r1cs_legosnark_verification_key<ppT> &other) const
{
    return (this->u_g2 == other.u_g2 &&
            this->v_g2 == other.v_g2 &&
            this->w_g2 == other.w_g2
            );
}
template <typename ppT>
std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_verification_key<ppT> &vk)
{
    out << vk.u_g2 << OUTPUT_NEWLINE;
    out << vk.v_g2 << OUTPUT_NEWLINE;
    out << vk.w_g2 << OUTPUT_NEWLINE;

    return out;
}
template <typename ppT>
std::istream &operator>>(std::istream &in, const r1cs_legosnark_verification_key<ppT> &vk)
{
    in >> vk.u_g2;
    libff::consume_OUTPUT_NEWLINE(in);
    in >> vk.v_g2;
    libff::consume_OUTPUT_NEWLINE(in);
    in >> vk.w_g2;
    libff::consume_OUTPUT_NEWLINE(in);

    return in;
}

template <typename ppT>
bool r1cs_legosnark_proof<ppT>::operator==(const r1cs_legosnark_proof<ppT> &other) const
{
    return (this->T_sum_G1 == other.T_sum_G1 &&
            this->R_sum_G1 == other.R_sum_G1);
}

template <typename ppT>
std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_proof<ppT> &proof)
{
    out << proof.T_sum_G1 << OUTPUT_NEWLINE;
    out << proof.R_sum_G1 << OUTPUT_NEWLINE;

    return out;
}

template <typename ppT>
std::istream &operator>>(std::istream &in, r1cs_legosnark_proof<ppT> &proof)
{
    in >> proof.T_sum_G1;
    libff::consume_OUTPUT_NEWLINE(in);
    in >> proof.R_sum_G1;
    libff::consume_OUTPUT_NEWLINE(in);

    return in;
}

// template <typename ppT>
// r1cs_legosnark_hash<ppT> r1cs_legosnark_hashing(const r1cs_legosnark_pp<ppT> &pp, const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input){

//     libff::G1<ppT> sigma_g1 = libff::G1<ppT>::zero();
//     sigma_g1 = pp.hfal_H[0];

//     for(size_t i = 0; i < primary_input.size(); i++){
//         sigma_g1 = primary_input[i] * pp.hfal_H[i+1] + sigma_g1;
//     }
//     r1cs_legosnark_hash<ppT> hash = r1cs_legosnark_hash<ppT>(std::move(sigma_g1));

// 	hash.print_size();

//     return hash;
// }

template <typename ppT>
r1cs_legosnark_keypair<ppT> r1cs_legosnark_generator(accumulation_vector<libff::G1<ppT> > &gamma_H, accumulation_vector<libff::G1<ppT> > &gamma_F, size_t len, size_t len2){
    libff::Fr<ppT> k1 = libff::Fr<ppT>::random_element();
    libff::Fr<ppT> k2 = libff::Fr<ppT>::random_element();
    libff::Fr<ppT> a = libff::Fr<ppT>::random_element();
    libff::G1_vector<ppT> P_vector;
    P_vector.reserve(len+1);

    P_vector.emplace_back((k1*gamma_H.first) + (k2* gamma_F.first));
    for(size_t i = 0; i < len; i++){
		//H_vector.emplace_back(libff::G1<ppT>::one());
		P_vector.emplace_back((k1 * gamma_H.rest.values[i]) + (k2 * gamma_F.rest.values[i]));
    }
    r1cs_legosnark_proving_key<ppT> ek = r1cs_legosnark_proving_key<ppT>(std::move(P_vector));
    libff::G2<ppT> g2a = a*libff::G2<ppT>::one();
    
    r1cs_legosnark_verification_key<ppT> vk 
        = r1cs_legosnark_verification_key<ppT>(
                                                        std::move(k1 * g2a),
                                                        std::move(k2 * g2a),
                                                        std::move(g2a)
                                                       );
    
    ek.print_size();
	vk.print_size();

    return r1cs_legosnark_keypair<ppT>(std::move(ek), std::move(vk));


}

template <typename ppT>
r1cs_legosnark_proof<ppT> r1cs_legosnark_prover(const r1cs_legosnark_proving_key<ppT> &pk, const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input, size_t len){
    libff::G1<ppT> Proof_sum = libff::G1<ppT>::zero();
    Proof_sum = pk.P_g1[0];
    for(size_t i = 0; i <  len; i++){ 
        // std::cout<<i<<", "<<(primary_input[i]).as_ulong()<<std::endl;
        Proof_sum = (primary_input[i] * pk.P_g1[i+1]) + Proof_sum;
    }

    r1cs_legosnark_proof<ppT> proof
        = r1cs_legosnark_proof<ppT>(std::move(Proof_sum));

	proof.print_size();

    return proof;
}

template <typename ppT>
bool r1cs_legosnark_verifier(const r1cs_legosnark_verification_key<ppT> &vk, 
                                    const libff::G1<ppT> &cm1, 
                                    const libff::G1<ppT> &cm2,
                                    const r1cs_legosnark_proof<ppT> &proof){


    libff::GT<ppT> left = ppT::reduced_pairing(proof.proof_G1, vk.a_g2);
    libff::GT<ppT> right = ppT::reduced_pairing(cm1, vk.C1_g2);
    libff::GT<ppT> right2 = ppT::reduced_pairing(cm2, vk.C2_g2);

    return (left == (right *  right2));
    //     printf("pass!!\n");
    
    
    // const libff::G1_precomp<ppT> proof_g_Proof_precomp = ppT::precompute_G1(proof.proof_G1);
    // const libff::G1_precomp<ppT> acc_precomp = ppT::precompute_G1(cm1);
    // const libff::G1_precomp<ppT> acc_precomp2 = ppT::precompute_G1(cm2);

    // const libff::G2_precomp<ppT> vk_C1_g2_precomp = ppT::precompute_G2(vk.C1_g2);
    // const libff::G2_precomp<ppT> vk_C2_g2_precomp = ppT::precompute_G2(vk.C2_g2);
    // const libff::G2_precomp<ppT> vk_a_g2_precomp = ppT::precompute_G2(vk.a_g2);

    // libff::Fqk<ppT> vf_right_one = ppT::miller_loop(proof_g_Proof_precomp, vk_a_g2_precomp);

    // libff::Fqk<ppT> vf_left_one = ppT::double_miller_loop(
    //     acc_precomp, vk_C1_g2_precomp,
    //     acc_precomp2, vk_C2_g2_precomp
    // );

    // libff::Fqk<ppT> vf_left_one2 = ppT::double_miller_loop(
    //     acc_precomp2, vk_C1_g2_precomp,
    //     acc_precomp, vk_C2_g2_precomp
    // ); 

    // vk.C1_g2.print();
    // vk.C2_g2.print();
    // vk.a_g2.print();

    // printf("vf_right_one<<\n");
    // vf_right_one.print();
    // printf("vf_left_one<<\n");
    // vf_left_one.print();
    // printf("vf_left_one2<<\n");
    // vf_left_one2.print();
    // return (vf_right_one == vf_left_one);
}
}// libsnark

#endif
