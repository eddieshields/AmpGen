#ifndef AMPGEN_MIXING_H
#define AMPGEN_MIXING_H

#include <stddef.h>
#include <complex>

#include "AmpGen/CoherentSum.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Types.h"
#include "AmpGen/MinuitParameterSet.h"


namespace AmpGen
{
  class Mixing
  {
  public:
    Mixing();
    Mixing( const EventType& type, const AmpGen::MinuitParameterSet& mps );

    ~Mixing() {};

    complex_t gp( const real_t& t ) const;
    complex_t gm( const real_t& t ) const;

    real_t prob( const Event& evt )              const;
    real_t prob_unnormalised( const Event& evt ) const;

    void setX( const real_t& x )    { _x = x; }
    void setY( const real_t& y )    { _y = y; }
    void setZ( const complex_t& z ) { _z = z; }

    // PDF methods.
    void prepare();
    void reset( bool resetEvents = false );
    void setEvents( EventList& list );


  private:
  CoherentSum pdf_dir;
  CoherentSum pdf_cnj;

  real_t _x = {0.004};
  real_t _y = {0.006};
  complex_t _z = complex_t(_x, _y);
  };
} // namespace AmpGen
#endif