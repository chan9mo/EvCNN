/** @file
*****************************************************************************

Implementation of Kate Polynomial Commitment for R1CS.

See kzg10.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef KZG_TCC_
#define KZG_TCC_

#include <iostream>
#include <sstream>
#include <type_traits>

#include <libff/common/profiling.hpp>
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>
#include <libsnark/common/data_structures/accumulation_vector.hpp>

#include <libsnark/common/default_types/kzg_pp.hpp>


namespace libsnark {

template<typename ppT>
commitkey<ppT> kzg_setup(int t)
{
     libff::enter_block("Call to kzg_setup");

    /* Generate generator g, randomness a */
    libff::G1<ppT> generator = libff::G1<ppT>::random_element();
    // G2 is for the vefiryeval (reduced_pairing) computation
    libff::G2<ppT> generator2 = libff::G2<ppT>::random_element(); 
    libff::Fr<ppT> a = libff::Fr<ppT>::random_element();

    /* Generate t-SDH tuple : G1 */
    libff::enter_block("Generating t-SDH tuple: G1");

    libff::enter_block("Generating G1 MSM window table");
    
    const size_t g1_scalar_count = t;
    const size_t g1_scalar_size = libff::Fr<ppT>::size_in_bits();
    const size_t g1_window_size = libff::get_exp_window_size<libff::G1<ppT>>(g1_scalar_count);

    libff::print_indent(); printf("* G1 window: %zu\n", g1_window_size);
    libff::window_table<libff::G1<ppT>> g1table = libff::get_window_table(g1_scalar_size, 
                                                                          g1_window_size, 
                                                                          generator);

    libff::leave_block("Generating G1 MSM window table");

    libff::enter_block("Randomness A: Windowed Exponentiation");
    libff::Fr<ppT> exp_a3 = libff::Fr<ppT>::one();
    libff::G1_vector<ppT> newg1tuple;
    newg1tuple.reserve(t + 1);

    newg1tuple[0] = generator;

    for(size_t j = 1; j < t; j++) {
    exp_a3 = exp_a3 * a;
    newg1tuple[j] = windowed_exp(g1_scalar_size, g1_window_size, g1table, exp_a3);
    }
    libff::leave_block("Randomness A: Windowed Exponentiation");

    libff::leave_block("Generating t-SDH tuple: G1");

    /* Generate t-SDH tuple : G2 */
    libff::enter_block("Generating t-SDH tuple: G2");

    libff::Fr<ppT> exp_a2 = libff::Fr<ppT>::one();

    libff::G2_vector<ppT> g2tuple;
    g2tuple.reserve(2);

    g2tuple[0] = generator2;
    exp_a2 = exp_a2 * a;
    g2tuple[1] = exp_a2 * generator2;

    libff::leave_block("Generating t-SDH tuple: G2");

    /* Output as a commitment key */
    libff::leave_block("Call to kzg_setup");

    commitkey<ppT> tuple = commitkey<ppT>(std::move(newg1tuple), 
                                          std::move(g2tuple), 
                                          std::move(g1table));
    return tuple;

}

template<typename ppT>
libff::G1<ppT> kzg_commit(commitkey<ppT> &ck, libff::Fr_vector<ppT> &poly, int t)
{
    libff::enter_block("Call to kzg_commit");

    libff::enter_block("Commit at G1");

    libff::enter_block("Rearrange Polynomial for Multi_exp");
    libff::Fr_vector<ppT> poly_fast; poly_fast.reserve(t);

    for(size_t i = 0; i < t; i++) {
        poly_fast[i] = poly[t - i - 1];
    }

    libff::leave_block("Rearrange Polynomial for Multi_exp");
    
    const size_t chunks = 1;

    libff::G1<ppT> commit1 = libff::multi_exp<libff::G1<ppT>, 
                                              libff::Fr<ppT>, 
                                              libff::multi_exp_method_BDLO12>
        (ck.g1.begin(), ck.g1.begin() + t, 
         poly_fast.begin(), poly_fast.begin() + t, 
         chunks);

    libff::leave_block("Commit at G1");

    libff::leave_block("Call to kzg_commit");

    return commit1;
}

template<typename ppT>
libff::Fr<ppT> kzg_hash(libff::G1<ppT> &commit_a)
{
    libff::enter_block("Call to kzg_hash"); //print, print_coordinates are diff.

    libff::enter_block("Extract Commit A, B, C's Coord[2]");

    std::string a_x = commit_a.coord[2].toString(10);

    libff::leave_block("Extract Commit A, B, C's Coord[2]");

    libff::enter_block("Hash: SHA256");

    libff::enter_block("String - Type Conversion -> Bit_vector");
    std::vector<size_t> randomness;
    size_t gamma = 0;

    for (int i = 0; i < a_x.size(); i++) {
    gamma = (size_t)a_x.at(i) - 48;
    randomness.push_back((size_t)gamma);
    }
    const size_t repack = 1; //what is the use of this repack variable?

    libff::Fr_vector<default_kzg_pp> abc_x = libff::pack_int_vector_into_field_element_vector<libff::Fr<default_kzg_pp>>(randomness, repack);
    libff::bit_vector input = libff::convert_field_element_vector_to_bit_vector<libff::Fr<default_kzg_pp>>(abc_x);
    libff::leave_block("String - Type Conversion -> Bit_vector");
   
    libff::bit_vector hash_result = sha256_two_to_one_hash_gadget<libff::Fr<default_kzg_pp>>::get_hash(input);
    libff::leave_block("Hash: SHA256");

    libff::leave_block("Call to kzg_hash");

    libff::Fr<ppT> point = libff::convert_bit_vector_to_field_element<libff::Fr<default_kzg_pp>>(hash_result);
    return point;
}

template<typename ppT>
witness<ppT> kzg_witness(commitkey<ppT> &ck, libff::Fr_vector<ppT> &poly, libff::Fr<ppT> &point, int t)
{
    libff::enter_block("Call to kzg_witness");
    int i;

    /* Evaluate Polynomial */
    libff::enter_block("Evaluate Polynomial");

    libff::Fr<ppT> eval = libff::Fr<ppT>::zero();
    libff::Fr<ppT> temp = libff::Fr<ppT>::one();

    for(i = 1; i < t + 1; i++) {
    eval += poly[t - i] * temp;
    temp *= point;
    }

    libff::leave_block("Evaluate Polynomial");

    libff::enter_block("Constant Update : poly(x) - poly(i)");

    poly[t - 1] = poly[t - 1] - eval;

    libff::leave_block("Constant Update : poly(x) - poly(i)");


    /* Divisor: (x - point) */
    libff::enter_block("Compute Divisor[2]: stands for polynomial (x - point)");

    libff::Fr_vector<ppT> divisor(2);
    divisor[0] = libff::Fr<ppT>::one();

    libff::Fr<ppT> minus = libff::Fr<ppT>::zero();
    minus = minus - point;
    divisor[1] = minus;

    libff::leave_block("compute divisor[2]: stands for polynomial (x - point)");

    //division Algorithm.
    libff::enter_block("Divide Algorithm: poly(x) - poly(i) / (x - i)");
    libff::Fr_vector<ppT> psi(t);

     for(i = 0; i < t - 1; i++) {
        psi[i] = poly[i];
        poly[i] = poly[i] - (psi[i] * divisor[0]);
        poly[i + 1] = poly[i + 1] - psi[i] * divisor[1];
    }

    libff::leave_block("Divide Algorithm: poly(x) - poly(i) / (x - i)");

    /* compute w = g ^ psi(a) */
    libff::enter_block("Compute w = g ^ psi(a): G1");

    libff::enter_block("Rearrange Polynomial for Multi_exp");

    libff::Fr_vector<ppT> psi_fast(t);

    for(size_t i = 0; i < t - 1; i++) {
        psi_fast[i] = psi[t - i - 2];
    }
    
    libff::leave_block("Rearrange Polynomial for Multi_exp");

    const size_t chunks = 1;

    libff::G1<ppT> w1 = libff::multi_exp<libff::G1<ppT>, 
                                         libff::Fr<ppT>, 
                                         libff::multi_exp_method_BDLO12>
        (ck.g1.begin(), ck.g1.begin() + t, 
        psi_fast.begin(), psi_fast.begin() + t, 
        chunks);

    libff::leave_block("Compute w = g ^ psi(a): G1");

    /* For Non-Interactive: Put evaluation as Group-1 Element */
    libff::G1<ppT> eval_g1 = eval * ck.g1[0];
    
    libff::leave_block("Call to kzg_witness");

    /* Output as a witness */
    witness<ppT> wit = witness<ppT>(std::move(point), 
                                    std::move(eval_g1), 
                                    std::move(w1));
    return wit;
}

template<typename ppT>
libff::Fr<ppT> kzg_evaluate(libff::Fr_vector<ppT> &poly, libff::Fr<ppT> &point, int t)
{
    libff::enter_block("Call to kzg_evaluate");

    libff::Fr<ppT> eval = libff::Fr<ppT>::zero();
    libff::Fr<ppT> temp = libff::Fr<ppT>::one();

    for(size_t i = 1; i < t + 1; i++) {
    eval += poly[t - i] * temp;
    temp *= point;
    }
    libff::leave_block("Call to kzg_evaluate");
    return eval;
}

template<typename ppT>
bool kzg_vfyeval(commitkey<ppT> &ck, libff::G1<ppT> &commit, witness<ppT> &witness)
{
    libff::enter_block("Call to kzg_vfyeval");

    /* LEFT SIDE: e(C, g) */
    libff::enter_block("Compute LEFT : e(C, g)");

    libff::GT<ppT> left1 = ppT::reduced_pairing(commit, ck.g2[0]);

    libff::leave_block("Compute LEFT : e(C, g)");

    /* RIGHT SIDE: e(w, g ^ (a-i)) * e(g ^ eval, g) */
    libff::enter_block("Compute RIGHT : e(w, g ^ (a-i)) * e(g ^ eval, g)");

    /*right1, 3: e(w, g ^ (a-i))*/
    libff::enter_block("Compute e(w, g ^ (a-i))");
    libff::Fr<ppT> zero = libff::Fr<ppT>::zero();

    /*g ^ (-i)*/
    libff::enter_block("Compute g ^ (-i)");

    libff::Fr<ppT> num = zero - witness.point;
    libff::G1<ppT> num1 = num * ck.g1[0];
    libff::G2<ppT> num2 = num * ck.g2[0];

    libff::leave_block("Compute g ^ (-i)");

    /* e(w, g ^ (- i) * g ^ a = g ^ (a - i)) */
    libff::GT<ppT> right1 = ppT::reduced_pairing(witness.w1, num2 + ck.g2[1]);

    libff::leave_block("Compute e(w, g ^ (a-i))");

    /* right2: e(g ^ eval, g) */ //eval which side? doesnt matter.
    libff::enter_block("Compute e(g ^ eval, g");

    libff::GT<ppT> right2 = ppT::reduced_pairing(witness.eval, ck.g2[0]);

    libff::leave_block("Compute e(g ^ eval, g");

    /* RIGHT: e(w, g ^ (a-i)) * e(g ^ eval, g), */

    libff::GT<ppT> right = right1 * right2;
    libff::leave_block("Compute RIGHT : e(w, g ^ (a-i)) * e(g ^ eval, g)");

    libff::enter_block("Verification: LEFT =? RIGHT");
    bool verifyresult;

    if (left1 == right) {
        verifyresult = true;
    } else {
        verifyresult = false;
    }

    libff::leave_block("Verification: LEFT =? RIGHT");

    libff::leave_block("Call to kzg_vfyeval");

    return verifyresult;
}

} //libsnark
#endif // KZG_TCC_