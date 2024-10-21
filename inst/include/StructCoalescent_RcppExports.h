// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_StructCoalescent_RCPPEXPORTS_H_GEN_
#define RCPP_StructCoalescent_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace StructCoalescent {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("StructCoalescent", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("StructCoalescent", "_StructCoalescent_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in StructCoalescent");
            }
        }
    }

    inline List DTALikelihoodC(NumericMatrix ED, NumericMatrix fit_mig_mat, NumericVector node_indices) {
        typedef SEXP(*Ptr_DTALikelihoodC)(SEXP,SEXP,SEXP);
        static Ptr_DTALikelihoodC p_DTALikelihoodC = NULL;
        if (p_DTALikelihoodC == NULL) {
            validateSignature("List(*DTALikelihoodC)(NumericMatrix,NumericMatrix,NumericVector)");
            p_DTALikelihoodC = (Ptr_DTALikelihoodC)R_GetCCallable("StructCoalescent", "_StructCoalescent_DTALikelihoodC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_DTALikelihoodC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(fit_mig_mat)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List DemeDecompC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
        typedef SEXP(*Ptr_DemeDecompC)(SEXP,SEXP,SEXP);
        static Ptr_DemeDecompC p_DemeDecompC = NULL;
        if (p_DemeDecompC == NULL) {
            validateSignature("List(*DemeDecompC)(NumericMatrix,int,NumericVector)");
            p_DemeDecompC = (Ptr_DemeDecompC)R_GetCCallable("StructCoalescent", "_StructCoalescent_DemeDecompC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_DemeDecompC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(n_deme)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List MTT_transition_kernelC(NumericMatrix ED, NumericMatrix bit_rates, NumericVector node_indices, NumericVector eigen_vals, NumericMatrix eigen_vecs, NumericMatrix inverse_vecs) {
        typedef SEXP(*Ptr_MTT_transition_kernelC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_MTT_transition_kernelC p_MTT_transition_kernelC = NULL;
        if (p_MTT_transition_kernelC == NULL) {
            validateSignature("List(*MTT_transition_kernelC)(NumericMatrix,NumericMatrix,NumericVector,NumericVector,NumericMatrix,NumericMatrix)");
            p_MTT_transition_kernelC = (Ptr_MTT_transition_kernelC)R_GetCCallable("StructCoalescent", "_StructCoalescent_MTT_transition_kernelC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_MTT_transition_kernelC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(bit_rates)), Shield<SEXP>(Rcpp::wrap(node_indices)), Shield<SEXP>(Rcpp::wrap(eigen_vals)), Shield<SEXP>(Rcpp::wrap(eigen_vecs)), Shield<SEXP>(Rcpp::wrap(inverse_vecs)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline NumericMatrix FitMigMatC(NumericMatrix bit_mm, NumericVector coal_rate) {
        typedef SEXP(*Ptr_FitMigMatC)(SEXP,SEXP);
        static Ptr_FitMigMatC p_FitMigMatC = NULL;
        if (p_FitMigMatC == NULL) {
            validateSignature("NumericMatrix(*FitMigMatC)(NumericMatrix,NumericVector)");
            p_FitMigMatC = (Ptr_FitMigMatC)R_GetCCallable("StructCoalescent", "_StructCoalescent_FitMigMatC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_FitMigMatC(Shield<SEXP>(Rcpp::wrap(bit_mm)), Shield<SEXP>(Rcpp::wrap(coal_rate)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

    inline NumericMatrix BitMigMatC(NumericMatrix fit_mm, NumericVector coal_rate) {
        typedef SEXP(*Ptr_BitMigMatC)(SEXP,SEXP);
        static Ptr_BitMigMatC p_BitMigMatC = NULL;
        if (p_BitMigMatC == NULL) {
            validateSignature("NumericMatrix(*BitMigMatC)(NumericMatrix,NumericVector)");
            p_BitMigMatC = (Ptr_BitMigMatC)R_GetCCallable("StructCoalescent", "_StructCoalescent_BitMigMatC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_BitMigMatC(Shield<SEXP>(Rcpp::wrap(fit_mm)), Shield<SEXP>(Rcpp::wrap(coal_rate)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

    inline List NodeCountC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
        typedef SEXP(*Ptr_NodeCountC)(SEXP,SEXP,SEXP);
        static Ptr_NodeCountC p_NodeCountC = NULL;
        if (p_NodeCountC == NULL) {
            validateSignature("List(*NodeCountC)(NumericMatrix,int,NumericVector)");
            p_NodeCountC = (Ptr_NodeCountC)R_GetCCallable("StructCoalescent", "_StructCoalescent_NodeCountC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NodeCountC(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(n_deme)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline NumericVector NodeIndicesC(NumericMatrix ED) {
        typedef SEXP(*Ptr_NodeIndicesC)(SEXP);
        static Ptr_NodeIndicesC p_NodeIndicesC = NULL;
        if (p_NodeIndicesC == NULL) {
            validateSignature("NumericVector(*NodeIndicesC)(NumericMatrix)");
            p_NodeIndicesC = (Ptr_NodeIndicesC)R_GetCCallable("StructCoalescent", "_StructCoalescent_NodeIndicesC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NodeIndicesC(Shield<SEXP>(Rcpp::wrap(ED)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline double SC_like_C(NumericMatrix ED, NumericVector coal_rate, NumericMatrix bit_mig_mat, NumericVector node_indices) {
        typedef SEXP(*Ptr_SC_like_C)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_SC_like_C p_SC_like_C = NULL;
        if (p_SC_like_C == NULL) {
            validateSignature("double(*SC_like_C)(NumericMatrix,NumericVector,NumericMatrix,NumericVector)");
            p_SC_like_C = (Ptr_SC_like_C)R_GetCCallable("StructCoalescent", "_StructCoalescent_SC_like_C");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_SC_like_C(Shield<SEXP>(Rcpp::wrap(ED)), Shield<SEXP>(Rcpp::wrap(coal_rate)), Shield<SEXP>(Rcpp::wrap(bit_mig_mat)), Shield<SEXP>(Rcpp::wrap(node_indices)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline arma::mat sample_path(const int a, const int b, const double t0, const double t1, const arma::mat& Q, const arma::mat& P) {
        typedef SEXP(*Ptr_sample_path)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_path p_sample_path = NULL;
        if (p_sample_path == NULL) {
            validateSignature("arma::mat(*sample_path)(const int,const int,const double,const double,const arma::mat&,const arma::mat&)");
            p_sample_path = (Ptr_sample_path)R_GetCCallable("StructCoalescent", "_StructCoalescent_sample_path");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_path(Shield<SEXP>(Rcpp::wrap(a)), Shield<SEXP>(Rcpp::wrap(b)), Shield<SEXP>(Rcpp::wrap(t0)), Shield<SEXP>(Rcpp::wrap(t1)), Shield<SEXP>(Rcpp::wrap(Q)), Shield<SEXP>(Rcpp::wrap(P)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

}

#endif // RCPP_StructCoalescent_RCPPEXPORTS_H_GEN_
