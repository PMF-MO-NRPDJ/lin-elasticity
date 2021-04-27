#pragma once

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>


/*  Klasa koja određuje vrijednost Dirichletovog rubnog uvjeta i 
    njegovo proširenje na čitavu domenu. 
    Template parametri:
       GV = GridView
       RF = Range Field Type (tip kojim su predstavljeni elementi slike funkcije)
       
       Ovo je rubni uvjet ua skalarnu komponentu.
    */
template<typename GV>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, BCExtension<GV> > {
  const GV& gv;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
  BCExtension (const GV& gv_) : gv(gv_) {}

  // uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno. 
  // e     = element 
  // xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
  // y      = izračunata Dirichletova vrijednost
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    //auto x = e.geometry().global(xlocal);
    y = 0.0;
    return;
  }

  inline const GV& getGridView () {return gv;}
};


// Rubni uvjet. Ovo je klasa za skalarnu komponentu
// što znači da pretpostavljamo isti tip uvjeta za sve komponente.
template <typename GV>
class BCTypeParam : public Dune::PDELab::DirichletConstraintsParameters
{
  const GV& gv;
public:
  BCTypeParam(const GV& gv_) : gv(gv_) {}

public:
  //  intersection = stranica elementa (u 3D) ili brid elementa (u 2D)
  //  coord        = lokalne koordinate točke na "intersectionu" koja se ispituje
  //  povratna vrijednost: true ako je točka na Dirichletovoj granici
  //                       false ako nije. 
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
    auto xg = intersection.geometry().global( coord );
    if( xg[0]  < 1E-6 )
        return true; //  lijeva strana, x=0, je učvršćena!

    return false;  
  }

  inline const GV& getGridView(){ return gv; }
};
