#ifndef AMPGEN_MIXING_H
#define AMPGEN_MIXING_H

#include <stddef.h>
#include <complex>

#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Types.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/CoherentSum.h"

namespace AmpGen
{

  class Mixing : public CoherentSum
  {
  public:
    Mixing();
    Mixing( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix = "" );
    ~Mixing() {};

    complex_t gp( const real_t& t ) const;
    complex_t gm( const real_t& t ) const;

    real_t prob_unnormalised( const Event& evt ) const;

    complex_t getVal   ( const Event& evt ) const;
    complex_t getValDir( const Event& evt ) const;
    complex_t getValCnj( const Event& evt ) const;

    void setX( const real_t& x )    { _x = x; }
    void setY( const real_t& y )    { _y = y; }
    void setZ( const complex_t& z ) { _z = z; }

  private:
    size_t           m_sizeDir;
    size_t           m_sizeCnj;
    EventType        m_evtTypeCnj;
    AmplitudeRules   m_rulesCnj;

    real_t    _x = {0.004};
    real_t    _y = {0.006};
    complex_t _z = complex_t(_x, _y);
  };
} // namespace AmpGen
#endif