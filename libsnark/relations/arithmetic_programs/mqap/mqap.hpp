/** @file
 *****************************************************************************

 Declaration of interfaces for a QAP ("Quadratic Arithmetic Program").

 QAPs are defined in \[GGPR13].

 References:

 \[GGPR13]:
 "Quadratic span programs and succinct NIZKs without PCPs",
 Rosario Gennaro, Craig Gentry, Bryan Parno, Mariana Raykova,
 EUROCRYPT 2013,
 <http://eprint.iacr.org/2012/215>

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MQAP_HPP_
#define MQAP_HPP_

#include <map>
#include <memory>

#include <libfqfft/evaluation_domain/evaluation_domain.hpp>

namespace libsnark {

/* forward declaration */
template<typename FieldT>
class mqap_witness;

/**
 * A QAP instance.
 *
 * Specifically, the datastructure stores:
 * - a choice of domain (corresponding to a certain subset of the field);
 * - the number of variables, the degree, and the number of inputs; and
 * - coefficients of the A,B,C polynomials in the Lagrange basis.
 *
 * There is no need to store the Z polynomial because it is uniquely
 * determined by the domain (as Z is its vanishing polynomial).
 */
template<typename FieldT>
class mqap_instance {
private:
    size_t num_variables_;
    size_t degree_;
	size_t degree2_;
	size_t degree4_;
    size_t num_inputs_;
    size_t num_kernels_;

public:
    std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain;

    std::vector<std::map<size_t, FieldT> > A_in_Lagrange_basis;
    std::vector<std::map<size_t, FieldT> > B_in_Lagrange_basis;
    std::vector<std::map<size_t, FieldT> > C_in_Lagrange_basis;

    mqap_instance(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
                 const size_t num_variables,
                 const size_t degree,
				 const size_t degree2, // other varible degree
                 const size_t num_inputs,
                 const size_t num_kernels,
                 const std::vector<std::map<size_t, FieldT> > &A_in_Lagrange_basis,
                 const std::vector<std::map<size_t, FieldT> > &B_in_Lagrange_basis,
                 const std::vector<std::map<size_t, FieldT> > &C_in_Lagrange_basis);

    mqap_instance(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
                 const size_t num_variables,
                 const size_t degree,
                 const size_t degree2,
                 const size_t num_inputs,
                 const size_t num_kernels,
                 std::vector<std::map<size_t, FieldT> > &&A_in_Lagrange_basis,
                 std::vector<std::map<size_t, FieldT> > &&B_in_Lagrange_basis,
                 std::vector<std::map<size_t, FieldT> > &&C_in_Lagrange_basis);

    mqap_instance(const mqap_instance<FieldT> &other) = default;
    mqap_instance(mqap_instance<FieldT> &&other) = default;
    mqap_instance& operator=(const mqap_instance<FieldT> &other) = default;
    mqap_instance& operator=(mqap_instance<FieldT> &&other) = default;

    size_t num_variables() const;
    size_t degree() const;
    size_t degree2() const;
    size_t degree4() const;
    size_t num_inputs() const;
    size_t num_kernels() const;

    bool is_satisfied(const mqap_witness<FieldT> &witness) const;
};

/**
 * A QAP instance evaluation is a QAP instance that is evaluated at a field element t.
 *
 * Specifically, the datastructure stores:
 * - a choice of domain (corresponding to a certain subset of the field);
 * - the number of variables, the degree, and the number of inputs;
 * - a field element t;
 * - evaluations of the A,B,C (and Z) polynomials at t;
 * - evaluations of all monomials of t;
 * - counts about how many of the above evaluations are in fact non-zero.
 */
template<typename FieldT>
class mqap_instance_evaluation {
private:
    size_t num_variables_;
    size_t degree_;
    size_t degree2_;
    size_t degree4_;
    size_t num_inputs_;
public:
    std::shared_ptr<libfqfft::evaluation_domain<FieldT> > domain;

    FieldT t;

    std::vector<FieldT> At, Bt, Ct, Ht;

    FieldT Zt;

    mqap_instance_evaluation(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
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
                            const FieldT &Zt);
    mqap_instance_evaluation(const std::shared_ptr<libfqfft::evaluation_domain<FieldT> > &domain,
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
                            const FieldT &Zt);

    mqap_instance_evaluation(const mqap_instance_evaluation<FieldT> &other) = default;
    mqap_instance_evaluation(mqap_instance_evaluation<FieldT> &&other) = default;
    mqap_instance_evaluation& operator=(const mqap_instance_evaluation<FieldT> &other) = default;
    mqap_instance_evaluation& operator=(mqap_instance_evaluation<FieldT> &&other) = default;

    size_t num_variables() const;
    size_t degree() const;
    size_t degree2() const;
    size_t degree4() const;
    size_t num_inputs() const;

    bool is_satisfied(const mqap_witness<FieldT> &witness) const;
};

/**
 * A QAP witness.
 */
template<typename FieldT>
class mqap_witness {
private:
    size_t num_variables_;
    const size_t degree,
    const size_t degree2,
    const size_t degree4,
    const size_t num_inputs,
    const size_t num_kernels,

public:
    FieldT d1, d2, d3;

    std::vector<FieldT> coefficients_for_ABCs;
    std::vector<FieldT> coefficients_for_H;

    mqap_witness(const size_t num_variables,
                const size_t degree,
                const size_t degree2,
                const size_t num_inputs,
                const size_t num_kernels,
                const FieldT &d1,
                const FieldT &d2,
                const FieldT &d3,
                const std::vector<FieldT> &coefficients_for_ABCs,
                const std::vector<FieldT> &coefficients_for_H);

    mqap_witness(const size_t num_variables,
                const size_t degree,
                const size_t degree2,
                const size_t num_inputs,
                const size_t num_kernels,
                const FieldT &d1,
                const FieldT &d2,
                const FieldT &d3,
                const std::vector<FieldT> &coefficients_for_ABCs,
                std::vector<FieldT> &&coefficients_for_H);

    mqap_witness(const mqap_witness<FieldT> &other) = default;
    mqap_witness(mqap_witness<FieldT> &&other) = default;
    mqap_witness& operator=(const mqap_witness<FieldT> &other) = default;
    mqap_witness& operator=(mqap_witness<FieldT> &&other) = default;

    size_t num_variables() const;
    size_t degree() const;
    size_t degree2() const;
    size_t num_inputs() const;
    size_t num_kernels() const;
};

} // libsnark

#include <libsnark/relations/arithmetic_programs/mqap/mqap.tcc>

#endif // QAP_HPP_
