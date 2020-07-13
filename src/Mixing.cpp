#include "AmpGen/Mixing.h"


using namespace AmpGen;

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
  _pdf_dir.prepare();
  _pdf_cnj.prepare();
}

void Mixing::reset( bool resetEvents )
{
  _pdf_dir.reset();
  _pdf_cnj.reset();
}

void Mixing::setEvents( EventList& list )
{
  _pdf_dir.setEvents( list );
  _pdf_cnj.setEvents( list );
}


