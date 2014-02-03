#ifndef _hgcsssimhit_hh_
#define _hgcsssimhit_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

#include "G4SiHit.hh"

static const float SIZE_X=200;//mm
static const float SIZE_Y=SIZE_X;
static const float CELL_SIZE_X=2.5;//mm
static const float CELL_SIZE_Y=CELL_SIZE_X;
static const unsigned N_CELLS_XY_MAX=SIZE_X/CELL_SIZE_X*SIZE_Y/CELL_SIZE_Y;
static const unsigned N_LAYERS=30;
static const unsigned GRANULARITY[N_LAYERS]={
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1,
  1,1,1,1,1
};

class HGCSSSimHit{

public:
  HGCSSSimHit():
    energy_(0),
    time_(0),
    layer_(0),
    cellid_(0),
    nGammas_(0),
    nElectrons_(0),
    nMuons_(0),
    nHadrons_(0)
  {
  };
  HGCSSSimHit(const G4SiHit & aSiHit);

  ~HGCSSSimHit(){};

  inline double energy() const {
    return energy_;
  };

  inline double time() const {
    return time_;
  };

  inline unsigned layer() const {
    return layer_;
  };

  inline void layer(const unsigned & layer){
    layer_ = layer;
  };

  inline unsigned cellid() const {
    return cellid_;
  };

  inline unsigned nGammas() const {
    return nGammas_;
  };

  inline unsigned nElectrons() const {
    return nElectrons_;
  };

  inline unsigned nMuons() const {
    return nMuons_;
  };

  inline unsigned nHadrons() const {
    return nHadrons_;
  };
  inline unsigned numberOfParticles() const {
    return nGammas_+nElectrons_+nMuons_+nHadrons_;
  };

  inline double gFrac() const {
    return nGammas_/numberOfParticles();
  };

  inline double eFrac() const {
    return nElectrons_/numberOfParticles();
  };

  inline  double muFrac() const {
    return nMuons_/numberOfParticles();
  };

  inline double hadFrac() const {
    return nHadrons_/numberOfParticles();
  };

  void Add(const G4SiHit & aSiHit);

  void encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell);

  inline bool get_x_side() const{
    return cellid_ & 0x0001;
  };

  inline bool get_y_side() const {
    return (cellid_ & 0x0100) >> 8;
  };

  inline unsigned get_x_cell() const {
    return (cellid_ & 0x00FE) >> 1;
  };

  inline unsigned get_y_cell() const {
    return (cellid_ & 0xFE00) >> 9;
  };

  inline double get_x() const {
    float sign = get_x_side() ? 1. : -1. ;
    if (sign > 0)
      return get_x_cell()*sign*CELL_SIZE_X*getGranularity()+CELL_SIZE_X*getGranularity()/2;
    else return get_x_cell()*sign*CELL_SIZE_X*getGranularity()-CELL_SIZE_X*getGranularity()/2;
  };

  inline double get_y() const {
    float sign = get_y_side() ? 1. : -1. ;
    if (sign > 0)
      return get_y_cell()*sign*CELL_SIZE_Y*getGranularity()+CELL_SIZE_Y*getGranularity()/2;
    else return get_y_cell()*sign*CELL_SIZE_Y*getGranularity()-CELL_SIZE_Y*getGranularity()/2;
  };

  inline unsigned getGranularity() const{
    return GRANULARITY[layer_];
  };

  void PrintGeometry(std::ostream & aOs) const ;

  void Print(std::ostream & aOs) const ;

private:

  double energy_;
  double time_;
  unsigned layer_;
  unsigned cellid_;
  unsigned nGammas_;
  unsigned nElectrons_;
  unsigned nMuons_;
  unsigned nHadrons_;

  ClassDef(HGCSSSimHit,1);



};


typedef std::vector<HGCSSSimHit> HGCSSSimHitVec;



#endif
