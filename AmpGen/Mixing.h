#ifndef AMPGEN_MIXING_H
#define AMPGEN_MIXING_H

#include <stddef.h>
#include <complex>


#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  template<class PDF_DIR, class PDF_CNJ> 
  class Mixing
  {
  public:
    Mixing();
    Mixing( PDF_DIR& pdf_dir, PDF_CNJ& pdf_cnj )
    {
      _pdf_dir = pdf_dir;
      _pdf_cnj = pdf_cnj;
    };
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
  PDF_DIR _pdf_dir;
  PDF_CNJ _pdf_cnj;

  real_t _x = {0.004};
  real_t _y = {0.006};
  complex_t _z = complex_t(_x, _y);
  };
} // namespace AmpGen
#endif