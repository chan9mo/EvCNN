#ifndef R1CS_LEGOSNARK_HPP_
#define R1CS_LEGOSNARK_HPP_

#include <memory>

#include <libff/algebra/curves/public_params.hpp>

#include <libsnark/common/data_structures/accumulation_vector.hpp>
#include <libsnark/knowledge_commitment/knowledge_commitment.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/r1cs_gg_ppzksnark_params.hpp>

namespace libsnark{

	template <typename ppT>
	class r1cs_legosnark_proving_key;

	template <typename ppT>
	std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_proving_key<ppT> &pk);

	template <typename ppT>
	std::istream &operator>>(std::istream &in, r1cs_legosnark_proving_key<ppT> &pk);

	template <typename ppT>
	class r1cs_legosnark_proving_key{
	public:
		libff::G1_vector<ppT> P_g1;

		r1cs_legosnark_proving_key(){};
		r1cs_legosnark_proving_key<ppT> &operator=(const r1cs_legosnark_proving_key<ppT> &other) = default;
		r1cs_legosnark_proving_key(const r1cs_legosnark_proving_key<ppT> &other) = default;
		r1cs_legosnark_proving_key(r1cs_legosnark_proving_key<ppT> &&other) = default;
		r1cs_legosnark_proving_key(const libff::G1_vector<ppT> &P_g1) : 
										P_g1(P_g1)
										{};

		size_t P_g1_size_in_bits() const
		{
			return P_g1.size()*libff::G1<ppT>::size_in_bits();
		}
		size_t size_in_bits() const
		{
			return P_g1_size_in_bits();
		}
		void print_size() const
		{
			libff::print_indent(); printf("* LEGO PK size in bits: %zu\n",this->size_in_bits());
		}

    	bool operator==(const r1cs_legosnark_proving_key<ppT> &other) const;
    	friend std::ostream& operator<< <ppT>(std::ostream &out, const r1cs_legosnark_proving_key<ppT> &pk);
    	friend std::istream& operator>> <ppT>(std::istream &in, r1cs_legosnark_proving_key<ppT> &pk);
	};

	template <typename ppT>
	class r1cs_legosnark_verification_key;

	template <typename ppT>
	std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_verification_key<ppT> &vk);

	template <typename ppT>
	std::istream &operator>>(std::istream &in, r1cs_legosnark_verification_key<ppT> &vk);

	template <typename ppT>
	class r1cs_legosnark_verification_key{
	public:
		libff::G2<ppT> C1_g2; // g2^(k1 * a)
		libff::G2<ppT> C2_g2; // g2^(k2 * a)
		libff::G2<ppT> a_g2; // g2^a

		r1cs_legosnark_verification_key() = default;
		r1cs_legosnark_verification_key(const libff::G2<ppT> &&C1_g2,
												 const libff::G2<ppT> &&C2_g2,
												 const libff::G2<ppT> &&a_g2) :
			C1_g2(std::move(C1_g2)),
			C2_g2(std::move(C2_g2)),
			a_g2(std::move(a_g2))
		{};

		size_t C1_g2_size() const
		{
			return 2;
		}
		size_t C2_g2_size() const
		{
			return 2;
		}
		size_t a_g2_size() const
		{
			return 2;
		}
		size_t size_in_bits() const
		{
			return C1_g2_size() * libff::G2<ppT>::size_in_bits() 
				+ C2_g2_size() * libff::G2<ppT>::size_in_bits() 
				+ a_g2_size() * libff::G2<ppT>::size_in_bits();
		}
		void print_size() const
		{
			libff::print_indent(); printf("* LEGO VK size in bits: %zu\n",this->size_in_bits());
		}
    	bool operator==(const r1cs_legosnark_verification_key<ppT> &other) const;
    	friend std::ostream& operator<< <ppT>(std::ostream &out, const r1cs_legosnark_verification_key<ppT> &pk);
    	friend std::istream& operator>> <ppT>(std::istream &in, r1cs_legosnark_verification_key<ppT> &pk);
	};

	template <typename ppT>
	class r1cs_legosnark_keypair;

	template <typename ppT>
	std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_keypair<ppT> &KEYPAIR);

	template <typename ppT>
	std::istream &operator>>(std::istream &in, r1cs_legosnark_keypair<ppT> &KEYPAIR);
	
	template <typename ppT>
	class r1cs_legosnark_keypair
	{
	  public:
		r1cs_legosnark_proving_key<ppT> pk;
		r1cs_legosnark_verification_key<ppT> vk;

		r1cs_legosnark_keypair(const r1cs_legosnark_keypair<ppT> &other):
			pk(std::move(other.pk)),
			vk(std::move(other.vk))
		{};
		r1cs_legosnark_keypair(r1cs_legosnark_proving_key<ppT> &&pk,
									r1cs_legosnark_verification_key<ppT> &&vk) : 
				pk(std::move(pk)),
				vk(std::move(vk))
		{};

		bool operator==(const r1cs_legosnark_keypair<ppT> &other) const;
		friend std::ostream &operator<<<ppT>(std::ostream &out, const r1cs_legosnark_keypair<ppT> &KEYPAIR);
		friend std::istream &operator>><ppT>(std::istream &in, r1cs_legosnark_keypair<ppT> &KEYPAIR);
	};

    /*
	template <typename ppT>
	class r1cs_HFAL_XP1_ppzksnark_hash;

	template <typename ppT>
	std::ostream &operator<<(std::ostream &out, const r1cs_HFAL_XP1_ppzksnark_hash<ppT> &proof);

	template <typename ppT>
	std::istream &operator>>(std::istream &in, r1cs_HFAL_XP1_ppzksnark_hash<ppT> &proof);
	
	template <typename ppT>
	class r1cs_HFAL_XP1_ppzksnark_hash{
	  public:
		libff::G1<ppT> HASH;

		r1cs_HFAL_XP1_ppzksnark_hash();
		r1cs_HFAL_XP1_ppzksnark_hash(libff::G1<ppT> &&HASH) :
			HASH(std::move(HASH))
		{};

		size_t HASH_size() const
		{
			return 1;
		}
		size_t size_in_bits() const
		{
			return HASH_size() * libff::G1<ppT>::size_in_bits();
		}
		void print_size() const
		{
			libff::print_indent(); printf("* XP1 hash size in bits: %zu\n", this->size_in_bits());
		}

		bool operator==(const r1cs_HFAL_XP1_ppzksnark_hash<ppT> &other) const;
		friend std::ostream &operator<<<ppT>(std::ostream &out, const r1cs_HFAL_XP1_ppzksnark_hash<ppT> &HASH);
		friend std::istream &operator>><ppT>(std::istream &in, r1cs_HFAL_XP1_ppzksnark_hash<ppT> &HASH);
	};
    */

	template <typename ppT>
	class r1cs_legosnark_proof;

	template <typename ppT>
	std::ostream &operator<<(std::ostream &out, const r1cs_legosnark_proof<ppT> &proof);

	template <typename ppT>
	std::istream &operator>>(std::istream &in, r1cs_legosnark_proof<ppT> &proof);

	template <typename ppT>
	class r1cs_legosnark_proof
	{
	  public:
	  	libff::G1<ppT> proof_G1;

		r1cs_legosnark_proof();
		r1cs_legosnark_proof(libff::G1<ppT> &&proof_G1) :
			proof_G1(std::move(proof_G1))
		{};

		size_t proof_G1_size() const
		{
			return 1;
		}
		size_t size_in_bits() const
		{
			return proof_G1_size() * libff::G1<ppT>::size_in_bits();
		}
		void print_size() const
		{
			libff::print_indent(); printf("* LEGO Proof size in bits: %zu\n", this->size_in_bits());
		}


		bool operator==(const r1cs_legosnark_proof<ppT> &other) const;
		friend std::ostream &operator<<<ppT>(std::ostream &out, const r1cs_legosnark_proof<ppT> &proof);
		friend std::istream &operator>><ppT>(std::istream &in, r1cs_legosnark_proof<ppT> &proof);
	};

	// template<typename ppT>
	// r1cs_HFAL_XP1_ppzksnark_hash<ppT> r1cs_HFAL_XP1_ppzksnark_hashing(const r1cs_HFAL_XP1_ppzksnark_pp<ppT> &pp,
	// 														   const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input);

	template<typename ppT>
	r1cs_legosnark_keypair<ppT> r1cs_legosnark_generator(accumulation_vector<libff::G1<ppT> > &gamma_H, accumulation_vector<libff::G1<ppT> > &gamma_F, size_t len);
	
	template<typename ppT>
	r1cs_legosnark_proof<ppT> r1cs_legosnark_prover(const r1cs_legosnark_proving_key<ppT> &pk,
																  const r1cs_gg_ppzksnark_primary_input<ppT> &r1cs_gg_ppzksnark_primary_input, size_t len);

	template<typename ppT>
	bool r1cs_legosnark_verifier(const r1cs_legosnark_verification_key<ppT> &vk, 
                                    const libff::G1<ppT> &cm1, 
                                    const libff::G1<ppT> &cm2,
                                    const r1cs_legosnark_proof<ppT> &proof);


}
#include <libsnark/zk_proof_systems/ppzksnark/convol_snark/r1cs_legosnark.tcc>

#endif


