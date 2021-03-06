// this file mainly refers to libint2's example

#include "jlcxx/jlcxx.hpp"
#include <libint2.hpp>
#include <iostream>
#include <thread>
#include <mutex>
#include <Eigen/Dense>

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

static std::mutex mtx;

size_t nbasis(const std::vector<libint2::Shell> &shells)
{
    size_t n = 0;
    for (const auto &shell : shells)
        n += shell.size();
    return n;
}

size_t max_nprim(const std::vector<libint2::Shell> &shells)
{
    size_t n = 0;
    for (auto shell : shells)
        n = std::max(shell.nprim(), n);
    return n;
}

int max_l(const std::vector<libint2::Shell> &shells)
{
    int l = 0;
    for (auto shell : shells)
        for (auto c : shell.contr)
            l = std::max(c.l, l);
    return l;
}

std::vector<libint2::Shell> construct_shells(int Nshell, int *l, double *alpha, double *coeff, double *cart)
{
    std::vector<libint2::Shell> shells;
    for (int i = 0; i < Nshell; ++i)
    {
        double *a = alpha + 3 * i;
        double *d = coeff + 3 * i;
        double *r = cart + 3 * i;
        shells.push_back({{a[0], a[1], a[2]},
                          {{l[i], false, {d[0], d[1], d[2]}}},
                          {{r[0], r[1], r[2]}}});
    }
    return shells;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}



void Ve(double *ve, int Nshell, int *l, double *alpha, double *d, double *R) {
  const auto shells = construct_shells(Nshell, l, alpha, d, R);
  const auto n = nbasis(shells);
  auto index = [&n](int i, int j, int k, int l) {
        return n * (n * (n * l + k) + j) + i;
  };
  std::memset(ve, 0, n*n*n*n*sizeof(double));
  libint2::Engine engine(libint2::Operator::coulomb, max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  const auto& buf = engine.results();

  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1];
    auto n1 = shells[s1].size();

    for(auto s2=0; s2!=shells.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for(auto s3=0; s3!=shells.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        for(auto s4=0; s4!=shells.size(); ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue;
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  ve[index(bf1, bf4, bf3, bf2)] = buf_1234[f1234];
                }
              }
            }
          }
        }
      }
    }
  }
}



// copy from libint2's example
void compute_2body_fock(double *veff, double* rho, int Nshell, int *l, double *alpha, double *d, double *R) {
  const auto shells = construct_shells(Nshell, l, alpha, d, R);
  const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> D(rho, n, n);

  // construct the 2-electron repulsion integrals engine
  libint2::Engine engine(libint2::Operator::coulomb, max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  const auto& buf = engine.results();

  // The problem with the simple Fock builder is that permutational symmetries of the Fock,
  // density, and two-electron integrals are not taken into account to reduce the cost.
  // To make the simple Fock builder efficient we must rearrange our computation.
  // The most expensive step in Fock matrix construction is the evaluation of 2-e integrals;
  // hence we must minimize the number of computed integrals by taking advantage of their permutational
  // symmetry. Due to the multiplicative and Hermitian nature of the Coulomb kernel (and realness
  // of the Gaussians) the permutational symmetry of the 2-e ints is given by the following relations:
  //
  // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
  //
  // (here we use chemists' notation for the integrals, i.e in (ab|cd) a and b correspond to
  // electron 1, and c and d -- to electron 2).
  //
  // It is easy to verify that the following set of nested loops produces a permutationally-unique
  // set of integrals:
  // foreach a = 0 .. n-1
  //   foreach b = 0 .. a
  //     foreach c = 0 .. a
  //       foreach d = 0 .. (a == c ? b : c)
  //         compute (ab|cd)
  //
  // The only complication is that we must compute integrals over shells. But it's not that complicated ...
  //
  // The real trick is figuring out to which matrix elements of the Fock matrix each permutationally-unique
  // (ab|cd) contributes. STOP READING and try to figure it out yourself. (to check your answer see below)

  // loop over permutationally-unique set of shells
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();   // number of basis functions in this shell

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for(auto s3=0; s3<=s1; ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for(auto s4=0; s4<=s4_max; ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // ANSWER
          // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
          //    F(a,b) += (ab|cd) * D(c,d)
          //    F(c,d) += (ab|cd) * D(a,b)
          //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
          //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
          //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
          //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
          // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
          //    i.e. the number of the integrals/sets equivalent to it
          // 3) the end result must be symmetrized
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1234[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                  G(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                  G(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                  G(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                  G(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                  G(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                }
              }
            }
          }

        }
      }
    }
  }

  // symmetrize the result
  for(int i=0; i<n; ++i)
  {
    for(int j=0; j<n; ++j)
    {
      veff[i*n + j] = (G(i,j) + G(j,i))/4.;
    }
  }
}

JLCXX_MODULE define_julia_module(jlcxx::Module &mod)
{
    mod.method("initialize", [](){libint2::initialize();});
    mod.method("Ve", &Ve);
    mod.method("Veff", &compute_2body_fock);
}