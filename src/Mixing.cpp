#include "AmpGen/Mixing.h"

#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/ThreadPool.h"

using namespace AmpGen;

Mixing::Mixing( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix  )
  : CoherentSum  ( true, type, mps, prefix )
  , m_evtTypeCnj (type.conj(1,0)) // Conjugate of D0.
{
  auto amplitudes      = m_rules.getMatchingRules( m_evtType, prefix);
  if( amplitudes.size() == 0 ){
    WARNING("The defined amplitudes don't seem to be able to be able to generate eventType: " << type);
  }
  auto amplitudesCnj  = m_rules.getMatchingRules( m_evtTypeCnj, prefix );
  if ( amplitudes.size() == 0 ){
    WARNING("The define amplitudes don't seem to be able to be able to generate eventType: " << type.conj(1,0));
  }
  m_sizeDir = amplitudes.size();
  for (auto& amp : amplitudesCnj ) amplitudes.push_back( amp );
  for( auto& amp : amplitudes )    INFO( amp.first.decayDescriptor() );
  // for( auto& amp : amplitudesCnj ) INFO( amp.first.decayDescriptor() );
  m_matrixElements.resize( amplitudes.size() );
  m_normalisations.resize( m_matrixElements.size(), m_matrixElements.size() );
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , std::thread::hardware_concurrency(), "Number of threads to use" );
  ThreadPool tp(nThreads);
  for(size_t i = 0; i < m_matrixElements.size(); ++i){
    if ( i < amplitudes.size() ) {
      tp.enqueue( [i,this,&mps,&amplitudes]{ 
      m_matrixElements[i] = TransitionMatrix<complex_t>( amplitudes[i].first,    amplitudes[i].second,    mps, this->m_evtType.getEventFormat(),    this->m_dbThis);
      CompilerWrapper().compile( m_matrixElements[i].amp, this->m_objCache); 
      } );
    }
    else {
      tp.enqueue( [i,this,&mps,&amplitudesCnj]{ 
      m_matrixElements[i] = TransitionMatrix<complex_t>( amplitudesCnj[i].first, amplitudesCnj[i].second, mps, this->m_evtTypeCnj.getEventFormat(), this->m_dbThis);
      CompilerWrapper().compile( m_matrixElements[i].amp, this->m_objCache); 
      } );
    }
  }
  m_isConstant = false;
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
  // | A |^2 = | g+( t )*A + q/p*g-( t )*Abar |^2
  return std::norm( gp( evt[evt.size()-1] )*getValDir( evt ) + gm( evt[evt.size()-1] )*getValCnj( evt ) );
}

complex_t Mixing::getVal( const Event& evt ) const
{
  complex_t value( 0., 0. ), dir( 0., 0. ), cnj( 0., 0. );
  for ( size_t i = 0; i < m_sizeDir; ++i )
    dir += m_matrixElements[i].coefficient * evt.getCache( m_matrixElements[i].addressData );
  for ( size_t i = m_sizeDir; i < m_matrixElements.size(); i++ )
    cnj += m_matrixElements[i].coefficient * evt.getCache( m_matrixElements[i].addressData );
  
  value = dir + cnj;
  return value;
}

complex_t Mixing::getValDir( const Event& evt ) const
{
  complex_t value( 0., 0. );
  for ( size_t i = 0; i < m_sizeDir; ++i )
    value += m_matrixElements[i].coefficient * evt.getCache( m_matrixElements[i].addressData );
  return value;
}

complex_t Mixing::getValCnj( const Event& evt ) const
{
  complex_t value( 0., 0. );
  for ( size_t i = m_sizeDir; i < m_matrixElements.size(); ++i )
    value += m_matrixElements[i].coefficient * evt.getCache( m_matrixElements[i].addressData );
  return value;
}




