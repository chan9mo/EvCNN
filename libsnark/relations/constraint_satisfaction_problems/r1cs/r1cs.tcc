/** @file
 *****************************************************************************

 Declaration of interfaces for:
 - a R1CS constraint,
 - a R1CS variable assignment, and
 - a R1CS constraint system.

 See r1cs.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef R1CS_TCC_
#define R1CS_TCC_

#include <algorithm>
#include <cassert>
#include <set>

#include <libff/algebra/fields/bigint.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

namespace libsnark {

template<typename FieldT>
r1cs_constraint<FieldT>::r1cs_constraint(const linear_combination<FieldT> &a,
                                         const linear_combination<FieldT> &b,
                                         const linear_combination<FieldT> &c) :
    a(a), b(b), c(c)
{
}

template<typename FieldT>
r1cs_constraint<FieldT>::r1cs_constraint(const std::initializer_list<linear_combination<FieldT> > &A,
                                         const std::initializer_list<linear_combination<FieldT> > &B,
                                         const std::initializer_list<linear_combination<FieldT> > &C)
{
    for (auto lc_A : A)
    {
        a.terms.insert(a.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B)
    {
        b.terms.insert(b.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C)
    {
        c.terms.insert(c.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }
}

template<typename FieldT>
r1cs_constraint<FieldT>::r1cs_constraint(const linear_combination<FieldT> &a,
                                         const linear_combination<FieldT> &b,
                                         const linear_combination<FieldT> &c,
                                         const linear_combination<FieldT> &a1,
                                         const linear_combination<FieldT> &b1,
                                         const linear_combination<FieldT> &c1) :
    a(a), b(b), c(c), a1(a1), b1(b1), c1(c1)
{
}

template<typename FieldT>
r1cs_constraint<FieldT>::r1cs_constraint(const std::initializer_list<linear_combination<FieldT> > &A,
                                         const std::initializer_list<linear_combination<FieldT> > &B,
                                         const std::initializer_list<linear_combination<FieldT> > &C,
                                         const std::initializer_list<linear_combination<FieldT> > &A1,
                                         const std::initializer_list<linear_combination<FieldT> > &B1,
                                         const std::initializer_list<linear_combination<FieldT> > &C1)
{
    for (auto lc_A : A)
    {
        a.terms.insert(a.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B)
    {
        b.terms.insert(b.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C)
    {
        c.terms.insert(c.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }

    for (auto lc_A : A1)
    {
        a1.terms.insert(a1.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B1)
    {
        b1.terms.insert(b1.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C1)
    {
        c1.terms.insert(c1.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }
}

template<typename FieldT>
r1cs_constraint<FieldT>::r1cs_constraint(const linear_combination<FieldT> &a,
                                         const linear_combination<FieldT> &b,
                                         const linear_combination<FieldT> &c,
                                         const linear_combination<FieldT> &a1,
                                         const linear_combination<FieldT> &b1,
                                         const linear_combination<FieldT> &c1,
                                         const linear_combination<FieldT> &a2,
                                         const linear_combination<FieldT> &b2,
                                         const linear_combination<FieldT> &c2) :
    a(a), b(b), c(c), a1(a1), b1(b1), c1(c1), a2(a2), b2(b2), c2(c2)
{
}

template<typename FieldT>
r1cs_constraint<FieldT>::r1cs_constraint(const std::initializer_list<linear_combination<FieldT> > &A,
                                         const std::initializer_list<linear_combination<FieldT> > &B,
                                         const std::initializer_list<linear_combination<FieldT> > &C,
                                         const std::initializer_list<linear_combination<FieldT> > &A1,
                                         const std::initializer_list<linear_combination<FieldT> > &B1,
                                         const std::initializer_list<linear_combination<FieldT> > &C1,
                                         const std::initializer_list<linear_combination<FieldT> > &A2,
                                         const std::initializer_list<linear_combination<FieldT> > &B2,
                                         const std::initializer_list<linear_combination<FieldT> > &C2)
{
    for (auto lc_A : A)
    {
        a.terms.insert(a.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B)
    {
        b.terms.insert(b.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C)
    {
        c.terms.insert(c.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }

    for (auto lc_A : A1)
    {
        a1.terms.insert(a1.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B1)
    {
        b1.terms.insert(b1.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C1)
    {
        c1.terms.insert(c1.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }

    for (auto lc_A : A2)
    {
        a2.terms.insert(a2.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B2)
    {
        b2.terms.insert(b2.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C2)
    {
        c2.terms.insert(c2.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }
}

template<typename FieldT>
bool r1cs_constraint<FieldT>::operator==(const r1cs_constraint<FieldT> &other) const
{
    return (this->a == other.a &&
            this->b == other.b &&
            this->c == other.c &&
            this->a2 == other.a1 &&
            this->b2 == other.b1 &&
            this->c2 == other.c1 &&
            this->a2 == other.a2 &&
            this->b2 == other.b2 &&
            this->c2 == other.c2);
}

template<typename FieldT>
std::ostream& operator<<(std::ostream &out, const r1cs_constraint<FieldT> &c)
{
    out << c.a;
    out << c.b;
    out << c.c;
    out << c.a1;
    out << c.b1;
    out << c.c1;
    out << c.a2;
    out << c.b2;
    out << c.c2;

    return out;
}

template<typename FieldT>
std::istream& operator>>(std::istream &in, r1cs_constraint<FieldT> &c)
{
    in >> c.a;
    in >> c.b;
    in >> c.c;
    in >> c.a1;
    in >> c.b1;
    in >> c.c1;
    in >> c.a2;
    in >> c.b2;
    in >> c.c2;

    return in;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_inputs() const
{
    return primary_input_size;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_variables() const
{
    return primary_input_size + auxiliary_input_size;
}


template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_constraints() const
{
    return constraints.size();
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol() const
{
    return convol_size;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_outputs() const
{
    return convol_outputs_size;
}
template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_outputs2() const
{
    return convol_outputs_size2;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_input_height() const
{
    return convol_input_height;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_input_width() const
{
    return convol_input_width;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_kernel_height() const
{
    return convol_kernel_height;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_kernel_width() const
{
    return convol_kernel_width;
}

template<typename FieldT>
size_t r1cs_constraint_system<FieldT>::num_convol_dimensions() const
{
    return convol_dimensions;
}

template<typename FieldT>
bool r1cs_constraint_system<FieldT>::is_valid() const
{
    if (this->num_inputs() > this->num_variables()) return false;

    for (size_t c = 0; c < constraints.size(); ++c)
    {
        if (!(constraints[c].a.is_valid(this->num_variables()) &&
              constraints[c].b.is_valid(this->num_variables()) &&
              constraints[c].c.is_valid(this->num_variables())))
        {
            return false;
        }
    }

    return true;
}

template<typename FieldT>
void dump_r1cs_constraint(const r1cs_constraint<FieldT> &constraint,
                          const r1cs_variable_assignment<FieldT> &full_variable_assignment,
                          const std::map<size_t, std::string> &variable_annotations)
{
    printf("terms for a:\n"); constraint.a.print_with_assignment(full_variable_assignment, variable_annotations);
    printf("terms for b:\n"); constraint.b.print_with_assignment(full_variable_assignment, variable_annotations);
    printf("terms for c:\n"); constraint.c.print_with_assignment(full_variable_assignment, variable_annotations);
}

template<typename FieldT>
bool r1cs_constraint_system<FieldT>::is_satisfied(const r1cs_primary_input<FieldT> &primary_input,
                                                  const r1cs_auxiliary_input<FieldT> &auxiliary_input) const
{
    assert(primary_input.size() == num_inputs());
    assert(primary_input.size() + auxiliary_input.size() == num_variables());

    r1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
    full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());

    for (size_t c = 0; c < constraints.size(); ++c)
    {
        const FieldT ares = constraints[c].a.evaluate(full_variable_assignment);
        const FieldT bres = constraints[c].b.evaluate(full_variable_assignment);
        const FieldT cres = constraints[c].c.evaluate(full_variable_assignment);

        if (!(ares*bres == cres))
        {
#ifdef DEBUG
            auto it = constraint_annotations.find(c);
            printf("constraint %zu (%s) unsatisfied\n", c, (it == constraint_annotations.end() ? "no annotation" : it->second.c_str()));
            printf("<a,(1,x)> = "); ares.print();
            printf("<b,(1,x)> = "); bres.print();
            printf("<c,(1,x)> = "); cres.print();
            printf("constraint was:\n");
            dump_r1cs_constraint(constraints[c], full_variable_assignment, variable_annotations);
#endif // DEBUG
            return false;
        }
    }

    return true;
}

template<typename FieldT>
void r1cs_constraint_system<FieldT>::add_constraint(const r1cs_constraint<FieldT> &c)
{
    constraints.emplace_back(c);
}

template<typename FieldT>
void r1cs_constraint_system<FieldT>::add_constraint(const r1cs_constraint<FieldT> &c, const std::string &annotation)
{
#ifdef DEBUG
    constraint_annotations[constraints.size()] = annotation;
#endif
    constraints.emplace_back(c);
}

template<typename FieldT>
void r1cs_constraint_system<FieldT>::swap_AB_if_beneficial()
{
    libff::enter_block("Call to r1cs_constraint_system::swap_AB_if_beneficial");

    libff::enter_block("Estimate densities");
    libff::bit_vector touched_by_A(this->num_variables() + 1, false), touched_by_B(this->num_variables() + 1, false);

    for (size_t i = 0; i < this->constraints.size(); ++i)
    {
        for (size_t j = 0; j < this->constraints[i].a.terms.size(); ++j)
        {
            touched_by_A[this->constraints[i].a.terms[j].index] = true;
        }

        for (size_t j = 0; j < this->constraints[i].b.terms.size(); ++j)
        {
            touched_by_B[this->constraints[i].b.terms[j].index] = true;
        }
    }

    size_t non_zero_A_count = 0, non_zero_B_count = 0;
    for (size_t i = 0; i < this->num_variables() + 1; ++i)
    {
        non_zero_A_count += touched_by_A[i] ? 1 : 0;
        non_zero_B_count += touched_by_B[i] ? 1 : 0;
    }

    if (!libff::inhibit_profiling_info)
    {
        libff::print_indent(); printf("* Non-zero A-count (estimate): %zu\n", non_zero_A_count);
        libff::print_indent(); printf("* Non-zero B-count (estimate): %zu\n", non_zero_B_count);
    }
    libff::leave_block("Estimate densities");

    if (non_zero_B_count > non_zero_A_count)
    {
        libff::enter_block("Perform the swap");
        for (size_t i = 0; i < this->constraints.size(); ++i)
        {
            std::swap(this->constraints[i].a, this->constraints[i].b);
        }
        libff::leave_block("Perform the swap");
    }
    else
    {
        libff::print_indent(); printf("Swap is not beneficial, not performing\n");
    }

    libff::leave_block("Call to r1cs_constraint_system::swap_AB_if_beneficial");
}

template<typename FieldT>
bool r1cs_constraint_system<FieldT>::operator==(const r1cs_constraint_system<FieldT> &other) const
{
    return (this->constraints == other.constraints &&
            this->primary_input_size == other.primary_input_size &&
            this->auxiliary_input_size == other.auxiliary_input_size);
}

template<typename FieldT>
std::ostream& operator<<(std::ostream &out, const r1cs_constraint_system<FieldT> &cs)
{
    out << cs.primary_input_size << "\n";
    out << cs.auxiliary_input_size << "\n";

    out << cs.num_constraints() << "\n";
    for (const r1cs_constraint<FieldT>& c : cs.constraints)
    {
        out << c;
    }

    return out;
}

template<typename FieldT>
std::istream& operator>>(std::istream &in, r1cs_constraint_system<FieldT> &cs)
{
    in >> cs.primary_input_size;
    in >> cs.auxiliary_input_size;

    cs.constraints.clear();

    size_t s;
    in >> s;

    char b;
    in.read(&b, 1);

    cs.constraints.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        r1cs_constraint<FieldT> c;
        in >> c;
        cs.constraints.emplace_back(c);
    }

    return in;
}

template<typename FieldT>
void r1cs_constraint_system<FieldT>::report_linear_constraint_statistics() const
{
#ifdef DEBUG
    for (size_t i = 0; i < constraints.size(); ++i)
    {
        auto &constr = constraints[i];
        bool a_is_const = true;
        for (auto &t : constr.a.terms)
        {
            a_is_const = a_is_const && (t.index == 0);
        }

        bool b_is_const = true;
        for (auto &t : constr.b.terms)
        {
            b_is_const = b_is_const && (t.index == 0);
        }

        if (a_is_const || b_is_const)
        {
            auto it = constraint_annotations.find(i);
            printf("%s\n", (it == constraint_annotations.end() ? FMT("", "constraint_%zu", i) : it->second).c_str());
        }
    }
#endif
}

template<typename FieldT>
void r1cs_constraint_system<FieldT>::add_convol_constraint(const size_t num_inputs, const size_t num_kernels)
{
    linear_combination<FieldT> A, B, C, A1, B1, C1;
    for(size_t i=0;i<num_kernels;i++){
        A1.add_term(i+1, i);
    }
    for(size_t i=0;i<num_inputs;i++){
        B1.add_term(num_kernels+i+1, i);
    }
    for(size_t i=0; i<num_kernels+num_inputs-1;i++){
        C1.add_term(num_kernels+num_inputs+i+1, i);
	}

    convol_outputs_size = num_kernels+num_inputs-1;//num_outputs;
    convol_size++;

    constraints.emplace_back(r1cs_constraint<FieldT>(A, B, C, A1, B1, C1));
}
/*
///////NEW SYSTEM////
template<typename FieldT>
r1cs_constraint_convol<FieldT>::r1cs_constraint_convol(const linear_combination<FieldT> &a,
                                         const linear_combination<FieldT> &b,
                                         const linear_combination<FieldT> &c) :
    a(a), b(b), c(c)
{
}

template<typename FieldT>
r1cs_constraint_convol<FieldT>::r1cs_constraint_convol(const std::initializer_list<linear_combination<FieldT> > &A,
                                         const std::initializer_list<linear_combination<FieldT> > &B,
                                         const std::initializer_list<linear_combination<FieldT> > &C)
{
    for (auto lc_A : A)
    {
        a.terms.insert(a.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B)
    {
        b.terms.insert(b.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C)
    {
        c.terms.insert(c.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }
}

template<typename FieldT>
r1cs_constraint_convol<FieldT>::r1cs_constraint_convol(const linear_combination<FieldT> &a,
                                         const linear_combination<FieldT> &b,
                                         const linear_combination<FieldT> &c,
                                         const linear_combination<FieldT> &a2,
                                         const linear_combination<FieldT> &b2,
                                         const linear_combination<FieldT> &c2) :
    a(a), b(b), c(c), a2(a2), b2(b2), c2(c2)
{
}

template<typename FieldT>
r1cs_constraint_convol<FieldT>::r1cs_constraint_convol(const std::initializer_list<linear_combination<FieldT> > &A,
                                         const std::initializer_list<linear_combination<FieldT> > &B,
                                         const std::initializer_list<linear_combination<FieldT> > &C,
                                         const std::initializer_list<linear_combination<FieldT> > &A2,
                                         const std::initializer_list<linear_combination<FieldT> > &B2,
                                         const std::initializer_list<linear_combination<FieldT> > &C2)
{
    for (auto lc_A : A)
    {
        a.terms.insert(a.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B)
    {
        b.terms.insert(b.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C)
    {
        c.terms.insert(c.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }

    for (auto lc_A : A2)
    {
        a2.terms.insert(a2.terms.end(), lc_A.terms.begin(), lc_A.terms.end());
    }
    for (auto lc_B : B2)
    {
        b2.terms.insert(b2.terms.end(), lc_B.terms.begin(), lc_B.terms.end());
    }
    for (auto lc_C : C2)
    {
        c2.terms.insert(c2.terms.end(), lc_C.terms.begin(), lc_C.terms.end());
    }
}


template<typename FieldT>
bool r1cs_constraint_convol<FieldT>::operator==(const r1cs_constraint_convol<FieldT> &other) const
{
    return (this->a == other.a &&
            this->b == other.b &&
            this->c == other.c &&
            this->a2 == other.a2 &&
            this->b2 == other.b2 &&
            this->c2 == other.c2);
}

template<typename FieldT>
std::ostream& operator<<(std::ostream &out, const r1cs_constraint_convol<FieldT> &c)
{
    out << c.a;
    out << c.b;
    out << c.c;
    out << c.a2;
    out << c.b2;
    out << c.c2;

    return out;
}

template<typename FieldT>
std::istream& operator>>(std::istream &in, r1cs_constraint_convol<FieldT> &c)
{
    in >> c.a;
    in >> c.b;
    in >> c.c;
    in >> c.a2;
    in >> c.b2;
    in >> c.c2;

    return in;
}

template<typename FieldT>
size_t r1cs_constraint_convol_system<FieldT>::num_inputs() const
{
    return primary_input_size;
}

template<typename FieldT>
size_t r1cs_constraint_convol_system<FieldT>::num_variables() const
{
    return primary_input_size + auxiliary_input_size;
}


template<typename FieldT>
size_t r1cs_constraint_convol_system<FieldT>::num_constraints() const
{
    return constraints.size();
}

template<typename FieldT>
size_t r1cs_constraint_convol_system<FieldT>::num_convol() const
{
    return convol_size;
}

template<typename FieldT>
size_t r1cs_constraint_convol_system<FieldT>::num_convol_outputs() const
{
    return convol_outputs_size;
}


template<typename FieldT>
bool r1cs_constraint_convol_system<FieldT>::is_valid() const
{
    if (this->num_inputs() > this->num_variables()) return false;

    for (size_t c = 0; c < constraints.size(); ++c)
    {
        if (!(constraints[c].a.is_valid(this->num_variables()) &&
              constraints[c].b.is_valid(this->num_variables()) &&
              constraints[c].c.is_valid(this->num_variables())))
        {
            return false;
        }
    }

    return true;
}

template<typename FieldT>
void dump_r1cs_constraint_convol(const r1cs_constraint_convol<FieldT> &constraint,
                          const r1cs_variable_assignment<FieldT> &full_variable_assignment,
                          const std::map<size_t, std::string> &variable_annotations)
{
    printf("terms for a:\n"); constraint.a.print_with_assignment(full_variable_assignment, variable_annotations);
    printf("terms for b:\n"); constraint.b.print_with_assignment(full_variable_assignment, variable_annotations);
    printf("terms for c:\n"); constraint.c.print_with_assignment(full_variable_assignment, variable_annotations);
}

template<typename FieldT>
bool r1cs_constraint_convol_system<FieldT>::is_satisfied(const r1cs_primary_input<FieldT> &primary_input,
                                                  const r1cs_auxiliary_input<FieldT> &auxiliary_input) const
{
    assert(primary_input.size() == num_inputs());
    assert(primary_input.size() + auxiliary_input.size() == num_variables());

    r1cs_variable_assignment<FieldT> full_variable_assignment = primary_input;
    full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());

    for (size_t c = 0; c < constraints.size(); ++c)
    {
        const FieldT ares = constraints[c].a.evaluate(full_variable_assignment);
        const FieldT bres = constraints[c].b.evaluate(full_variable_assignment);
        const FieldT cres = constraints[c].c.evaluate(full_variable_assignment);

        if (!(ares*bres == cres))
        {
#ifdef DEBUG
            auto it = constraint_annotations.find(c);
            printf("constraint %zu (%s) unsatisfied\n", c, (it == constraint_annotations.end() ? "no annotation" : it->second.c_str()));
            printf("<a,(1,x)> = "); ares.print();
            printf("<b,(1,x)> = "); bres.print();
            printf("<c,(1,x)> = "); cres.print();
            printf("constraint was:\n");
            dump_r1cs_constraint(constraints[c], full_variable_assignment, variable_annotations);
#endif // DEBUG
            return false;
        }
    }

    return true;
}

template<typename FieldT>
void r1cs_constraint_convol_system<FieldT>::add_constraint(const r1cs_constraint_convol<FieldT> &c)
{
    constraints.emplace_back(c);
}

template<typename FieldT>
void r1cs_constraint_convol_system<FieldT>::add_constraint(const r1cs_constraint_convol<FieldT> &c, const std::string &annotation)
{
#ifdef DEBUG
    constraint_annotations[constraints.size()] = annotation;
#endif
    constraints.emplace_back(c);
}

template<typename FieldT>
void r1cs_constraint_convol_system<FieldT>::add_convol_constraint(const size_t num_inputs, const size_t num_kernels)
{
    linear_combination<FieldT> A, B, C, A2, B2, C2;
    for(size_t i=0;i<num_kernels;i++){
        A2.add_term(i+1, i);
    }
    for(size_t i=0;i<num_inputs;i++){
        B2.add_term(num_kernels+i+1, i);
    }
    for(size_t i=0; i<num_kernels+num_inputs-1;i++){
        C2.add_term(num_kernels+num_inputs+i+1, i);
	}

    constraints.emplace_back(r1cs_constraint_convol<FieldT>(A, B, C, A2, B2, C2));
}

template<typename FieldT>
void r1cs_constraint_convol_system<FieldT>::swap_AB_if_beneficial()
{
    libff::enter_block("Call to r1cs_constraint_convol_system::swap_AB_if_beneficial");

    libff::enter_block("Estimate densities");
    libff::bit_vector touched_by_A(this->num_variables() + 1, false), touched_by_B(this->num_variables() + 1, false);

    for (size_t i = 0; i < this->constraints.size(); ++i)
    {
        for (size_t j = 0; j < this->constraints[i].a.terms.size(); ++j)
        {
            touched_by_A[this->constraints[i].a.terms[j].index] = true;
        }

        for (size_t j = 0; j < this->constraints[i].b.terms.size(); ++j)
        {
            touched_by_B[this->constraints[i].b.terms[j].index] = true;
        }
    }

    size_t non_zero_A_count = 0, non_zero_B_count = 0;
    for (size_t i = 0; i < this->num_variables() + 1; ++i)
    {
        non_zero_A_count += touched_by_A[i] ? 1 : 0;
        non_zero_B_count += touched_by_B[i] ? 1 : 0;
    }

    if (!libff::inhibit_profiling_info)
    {
        libff::print_indent(); printf("* Non-zero A-count (estimate): %zu\n", non_zero_A_count);
        libff::print_indent(); printf("* Non-zero B-count (estimate): %zu\n", non_zero_B_count);
    }
    libff::leave_block("Estimate densities");

    if (non_zero_B_count > non_zero_A_count)
    {
        libff::enter_block("Perform the swap");
        for (size_t i = 0; i < this->constraints.size(); ++i)
        {
            std::swap(this->constraints[i].a, this->constraints[i].b);
        }
        libff::leave_block("Perform the swap");
    }
    else
    {
        libff::print_indent(); printf("Swap is not beneficial, not performing\n");
    }

    libff::leave_block("Call to r1cs_constraint_convol_system::swap_AB_if_beneficial");
}

template<typename FieldT>
bool r1cs_constraint_convol_system<FieldT>::operator==(const r1cs_constraint_convol_system<FieldT> &other) const
{
    return (this->constraints == other.constraints &&
            this->primary_input_size == other.primary_input_size &&
            this->auxiliary_input_size == other.auxiliary_input_size);
}

template<typename FieldT>
std::ostream& operator<<(std::ostream &out, const r1cs_constraint_convol_system<FieldT> &cs)
{
    out << cs.primary_input_size << "\n";
    out << cs.auxiliary_input_size << "\n";

    out << cs.num_constraints() << "\n";
    for (const r1cs_constraint<FieldT>& c : cs.constraints)
    {
        out << c;
    }

    return out;
}

template<typename FieldT>
std::istream& operator>>(std::istream &in, r1cs_constraint_convol_system<FieldT> &cs)
{
    in >> cs.primary_input_size;
    in >> cs.auxiliary_input_size;
    in >> cs.convol_outputs_size;
    in >> cs.convol_size;

    cs.constraints.clear();

    size_t s;
    in >> s;

    char b;
    in.read(&b, 1);

    cs.constraints.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        r1cs_constraint_convol<FieldT> c;
        in >> c;
        cs.constraints.emplace_back(c);
    }

    return in;
}

template<typename FieldT>
void r1cs_constraint_convol_system<FieldT>::report_linear_constraint_statistics() const
{
#ifdef DEBUG
    for (size_t i = 0; i < constraints.size(); ++i)
    {
        auto &constr = constraints[i];
        bool a_is_const = true;
        for (auto &t : constr.a.terms)
        {
            a_is_const = a_is_const && (t.index == 0);
        }

        bool b_is_const = true;
        for (auto &t : constr.b.terms)
        {
            b_is_const = b_is_const && (t.index == 0);
        }

        if (a_is_const || b_is_const)
        {
            auto it = constraint_annotations.find(i);
            printf("%s\n", (it == constraint_annotations.end() ? FMT("", "constraint_%zu", i) : it->second).c_str());
        }
    }
#endif
}
*/

} // libsnark
#endif // R1CS_TCC_
