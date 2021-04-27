#pragma once

//#include<dune/geometry/quadraturerules.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

/** Lokalni operator za zadaću linearne elastičnosti
 *
 *
 * \tparam BCType tip rubnog uvjeta
 * \tparam FEM skalarni prostor konačnih elemenata
 */

template<typename BCType, typename FEM>
class ElasticityLocalOperator : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianVolume       <ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianApplyBoundary<ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianBoundary     <ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // Zastavice koje signaliziraju da na svakom elementu treba zvati:
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };

  ElasticityLocalOperator(const BCType& bctype_, // boundary cond.type
                          double mu_, double lambda_, double g_vert_, double rho_,
                          unsigned int intorder_=4) :
    bctype( bctype_ ), mu(mu_), lambda(lambda_), g_vert(g_vert_), rho(rho_), intorder( intorder_ )
  {}

  // volume integral depending on test and ansatz functions
  // eg = element
  // lfsu = lokalni prostor funkcija za rješenje
  // lfsv =  lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja
  // r    = lokalni rezidual
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    const int dim = EG::Geometry::coorddimension;

    using Gradient =  Dune::FieldVector<double,dim>;
    using size_type = typename LFSU::Traits::SizeType;

    // Kvadraturna formula na elementu.
    const auto & rule = Dune::PDELab::quadratureRule(eg.geometry(), intorder);

    // skalarni prostori
    using namespace Dune::Indices;
    auto const lfsu0 = lfsu.child(_0);
    auto const lfsu1 = lfsu.child(_1);
    auto const lfsu2 = lfsu.child(_2);

       for(auto const & qp : rule){
           // Uzimamo samo skalarne bazne funkcije jer su bazne funkcije za svaku komponentu iste.
           auto& phi0 = cache.evaluateFunction(qp.position(), lfsu0.finiteElement().localBasis());
           // Izračunajmo sve komponenete pomaka u integracijskoj točki.
           double u_0=0.0;
           double u_1=0.0;
           double u_2=0.0;
           for (size_type i=0; i<lfsu0.size(); ++i) u_0 += x(lfsu0,i)*phi0[i];
           for (size_type i=0; i<lfsu1.size(); ++i) u_1 += x(lfsu1,i)*phi0[i];
           for (size_type i=0; i<lfsu2.size(); ++i) u_2 += x(lfsu2,i)*phi0[i];

           // Gradijenti skalarnih baznih funkcija na ref elementu.
           auto const & js0 = cache.evaluateJacobian(qp.position(), lfsu0.finiteElement().localBasis());
           // Gradijenti skalarnih baznih funkcija na fizičkom elementu
           const auto &jac = eg.geometry().jacobianInverseTransposed(qp.position());
           std::vector<Gradient> gradphi0(lfsu0.size());
           for (size_type i=0; i<lfsu0.size(); i++)
               jac.mv(js0[i][0],gradphi0[i]);  // gradphi0[i] = jac * js0[i][0]

           // Gradijent komponenti rješenja (pomaka).
           Gradient gradu_0(0.0);
           Gradient gradu_1(0.0);
           Gradient gradu_2(0.0);
           for (size_type i=0; i<lfsu0.size(); ++i)
               gradu_0.axpy(x(lfsu0,i), gradphi0[i]); // gradu_0 += x(lfsu0,i)*gradphi0[i]
           for (size_type i=0; i<lfsu1.size(); ++i)
               gradu_1.axpy(x(lfsu1,i), gradphi0[i]); // gradu_1 += x(lfsu1,i)*gradphi0[i]
           for (size_type i=0; i<lfsu2.size(); ++i)
               gradu_2.axpy(x(lfsu2,i), gradphi0[i]); // gradu_2 += x(lfsu2,i)*gradphi0[i]

           Gradient f;
           f[0] = 0.0;
           f[1] = 0.0;
           f[2] = -9.81 * rho; // kgm/s^2


           double divu = gradu_0[0] + gradu_1[1] + gradu_2[2];

           double eu[3][3];
           eu[0][0] = 0.5*( gradu_0[0] + gradu_0[0] );
           eu[0][1] = 0.5*( gradu_0[1] + gradu_1[0] );
           eu[0][2] = 0.5*( gradu_0[2] + gradu_2[0] );
           eu[1][0] = 0.5*( gradu_1[0] + gradu_0[1] );
           eu[1][1] = 0.5*( gradu_1[1] + gradu_1[1] );
           eu[1][2] = 0.5*( gradu_1[2] + gradu_2[1] );
           eu[2][0] = 0.5*( gradu_2[0] + gradu_0[2] );
           eu[2][1] = 0.5*( gradu_2[1] + gradu_1[2] );
           eu[2][2] = 0.5*( gradu_2[2] + gradu_2[2] );

           double factor = qp.weight()*eg.geometry().integrationElement(qp.position());

           for (size_type i=0; i<lfsu0.size(); ++i)
             r.accumulate(lfsu0, i, (  2*mu * ( eu[0][0]*gradphi0[i][0] + eu[0][1]*gradphi0[i][1] + eu[0][2]*gradphi0[i][2] )
                                      + lambda *( divu * gradphi0[i][0] )
                                      - f[0]*phi0[i]
                                    ) * factor);

           for (size_type i=0; i<lfsu1.size(); ++i)
             r.accumulate(lfsu1, i, (  2*mu * ( eu[1][0]*gradphi0[i][0] + eu[1][1]*gradphi0[i][1] + eu[1][2]*gradphi0[i][2] )
                                      + lambda * ( divu * gradphi0[i][1] )
                                      - f[1]*phi0[i]
                                    ) * factor);

           for (size_type i=0; i<lfsu2.size(); ++i)
             r.accumulate(lfsu2, i, (  2*mu * ( eu[2][0]*gradphi0[i][0] + eu[2][1]*gradphi0[i][1] + eu[2][2]*gradphi0[i][2] )
                                      + lambda * ( divu * gradphi0[i][2] )
                                      - f[2]*phi0[i]
                                    ) * factor);

         }
  }

  // boundary integral
  // ig = intersection (= stranica elementa)
  // lfsu_s = lokalni prostor funkcija na stranici za rješenje
  // lfsv_s = lokalni prostor funkcija na stranici za test funkciju
  // x_s    = vektor koeficijenata rješenja (na stranici)
  // r_s    = rezidual (na stranici)
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    const int dim = IG::coorddimension;
    // skalarni prostori
    auto const & lfsu0 = lfsu_s.template child<0>();
    auto const & lfsu1 = lfsu_s.template child<1>();
    auto const & lfsu2 = lfsu_s.template child<2>();

    // kvadraturna formula na stranici
    auto const & rule = Dune::PDELab::quadratureRule(ig.geometry(), intorder);

    for (auto const & qp : rule)
      {
        if ( bctype.isDirichlet( ig, qp.position() ) )
          continue;

        auto globalpos = ig.geometry().global(qp.position());

        // Opterećenje imamo samo na gornjoj granici y=2.
        if(globalpos[2] < 2.0 - 1E-8)
              continue;

        auto local = ig.geometryInInside().global(qp.position());

        // Imamo onoliko vektora baznih funkcija koliko ima komponenata, ali sve su "iste"
        auto& phi0 = cache.evaluateFunction(local,lfsu0.finiteElement().localBasis());

        auto factor = qp.weight()*ig.geometry().integrationElement(qp.position());

        Dune::FieldVector<double,dim> flux;
        flux[0] = 0.0;
        flux[1] = 0.0;
        flux[2] = -g_vert;
/*
        for (size_t i=0; i<lfsu0.size(); ++i)
          r_s.accumulate(lfsu0,i, flux[0] * phi0[i] * factor);

        for (size_t i=0; i<lfsu1.size(); ++i)
          r_s.accumulate(lfsu1,i, flux[1] * phi0[i] * factor);
*/
        for (size_t i=0; i<lfsu2.size(); ++i)
          r_s.accumulate(lfsu2,i, flux[2] * phi0[i] * factor);
      }
  }

private:
  const BCType & bctype;
  double mu;
  double lambda;
  double g_vert;
  double rho;
  unsigned int intorder;
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};



