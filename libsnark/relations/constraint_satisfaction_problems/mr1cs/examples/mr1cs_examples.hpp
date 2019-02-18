/** @file
 *****************************************************************************

 Declaration of interfaces for a R1CS example, as well as functions to sample
 R1CS examples with prescribed parameters (according to some distribution).

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef R1CS_EXAMPLES_HPP_
#define R1CS_EXAMPLES_HPP_

#include <libsnark/relations/constraint_satisfaction_problems/mr1cs/mr1cs.hpp>

namespace libsnark {

/**
 * A R1CS example comprises a R1CS constraint system, R1CS input, and R1CS witness.
 */
template<typename FieldT>
struct mr1cs_example {
    mr1cs_constraint_system<FieldT> constraint_system;
    mr1cs_primary_input<FieldT> primary_input;
    mr1cs_auxiliary_input<FieldT> auxiliary_input;

    mr1cs_example<FieldT>() = default;
    mr1cs_example<FieldT>(const mr1cs_example<FieldT> &other) = default;
    mr1cs_example<FieldT>(const mr1cs_constraint_system<FieldT> &constraint_system,
                         const mr1cs_primary_input<FieldT> &primary_input,
                         const mr1cs_auxiliary_input<FieldT> &auxiliary_input) :
        constraint_system(constraint_system),
        primary_input(primary_input),
        auxiliary_input(auxiliary_input)
    {};
    mr1cs_example<FieldT>(mr1cs_constraint_system<FieldT> &&constraint_system,
                         mr1cs_primary_input<FieldT> &&primary_input,
                         mr1cs_auxiliary_input<FieldT> &&auxiliary_input) :
        constraint_system(std::move(constraint_system)),
        primary_input(std::move(primary_input)),
        auxiliary_input(std::move(auxiliary_input))
    {};
};

template<typename FieldT>
mr1cs_example<FieldT> generate_mr1cs_example_with_image_convol(const size_t num_kernel,
                                                            const size_t num_w, const size_t num_h);


template<typename FieldT>
mr1cs_example<FieldT> generate_mr1cs_example_with_convol(const size_t num_a,
                                                            const size_t num_x);

/**
 * Generate a R1CS example such that:
 * - the number of constraints of the R1CS constraint system is num_constraints;
 * - the number of variables of the R1CS constraint system is (approximately) num_constraints;
 * - the number of inputs of the R1CS constraint system is num_inputs;
 * - the R1CS input consists of ``full'' field elements (typically require the whole log|Field| bits to represent).
 */
template<typename FieldT>
mr1cs_example<FieldT> generate_mr1cs_example_with_field_input(const size_t num_constraints,
                                                            const size_t num_inputs);

/**
 * Generate a R1CS example such that:
 * - the number of constraints of the R1CS constraint system is num_constraints;
 * - the number of variables of the R1CS constraint system is (approximately) num_constraints;
 * - the number of inputs of the R1CS constraint system is num_inputs;
 * - the R1CS input consists of binary values (as opposed to ``full'' field elements).
 */
template<typename FieldT>
mr1cs_example<FieldT> generate_mr1cs_example_with_binary_input(const size_t num_constraints,
                                                             const size_t num_inputs);

} // libsnark

#include <libsnark/relations/constraint_satisfaction_problems/mr1cs/examples/mr1cs_examples.tcc>

#endif // R1CS_EXAMPLES_HPP_
