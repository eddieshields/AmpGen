#include "AmpGen/Mixing.h"


using namespace AmpGen;

Mixing::Mixing(const EventType& type, const AmpGen::MinuitParameterSet& mps)
{
  // Define PDF.
  INFO( type );
  pdf_dir = CoherentSum( type, mps );

  // Define conj. PDF.
  EventType eventType = type.conj();
  INFO( typeCnj );
  MinuitParameterSet mpsCnj( mps )
  AddCPConjugate(mpsCnj);
  pdf_cnj = CoherentSum( typeCnj, mps );
}

complex_t Mixing::gp(const real_t& t) const
{
  return std::cosh( _z * ( t/2 ) );
}

complex_t Mixing::gm(const real_t& t) const
{
  return std::sinh( _z * ( t/2 ) );
}

real_t Mixing::prob_unnormalised( const Event& evt ) const
{
  return std::norm( _pdf_dir.getVal( evt ) * gp( evt[ evt.size() -1 ] ) + _pdf_cnj.getVal( evt ) * gm( evt[ evt.size() -1 ] ) );
}

void Mixing::prepare()
{
  pdf_dir.prepare();
  pdf_cnj.prepare();
}

void Mixing::reset( bool resetEvents )
{
  pdf_dir.reset();
  pdf_cnj.reset();
}

void Mixing::setEvents( EventList& list )
{
  pdf_dir.setEvents( list );
  pdf_cnj.setEvents( list );
}


