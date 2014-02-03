#ifndef _hgcssrecohit_hh_
#define _hgcssrecohit_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>

#include "HGCSSSimHit.hh"

class HGCSSRecoHit{

public:
  HGCSSRecoHit():
    energy_(0),
    adcCounts_(0),
    layer_(0),
    cellid_(0),
    noiseFrac_(0)
  {};

  HGCSSRecoHit(const HGCSSSimHit & aSimHit, const unsigned granularity=1);

  ~HGCSSRecoHit(){};

  inline double energy() const {
    return energy_;
  };

  inline void energy(const double & energy) {
    energy_ = energy;
  };

  inline unsigned adcCounts() const {
    return adcCounts_;
  };

  inline void adcCounts(const unsigned & adcCounts){
    adcCounts_ = adcCounts;
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

  inline void cellid(const unsigned & id){
    cellid_ = id;
  };

  inline unsigned fullcellid() const {
    return cellid_ | (layer_<<24);
  };

  inline double noiseFraction() const {
    return noiseFrac_;
  };

  inline void noiseFraction(const double & aFrac){
    noiseFrac_ = aFrac;
  };

  void Add(const HGCSSSimHit & aSimHit);
  
  void encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell, const unsigned granularity);

  inline bool get_x_side() const{
    return cellid_ & 0x0001;
  };

  inline bool get_y_side() const{
    return (cellid_ & 0x0100) >> 8;
  };

  inline unsigned get_x_cell() const{
    return (cellid_ & 0x00FE) >> 1;
  };

  inline unsigned get_y_cell() const{
    return (cellid_ & 0xFE00) >> 9;
  };

  inline double get_x() const{
    float sign = get_x_side() ? 1. : -1. ;
    if (sign > 0)
      return get_x_cell()*sign*CELL_SIZE_X*getGranularity()+CELL_SIZE_X*getGranularity()/2;
    else return get_x_cell()*sign*CELL_SIZE_X*getGranularity()-CELL_SIZE_X*getGranularity()/2;
  };

  inline double get_y() const{
    float sign = get_y_side() ? 1. : -1. ;
    if (sign > 0)
      return get_y_cell()*sign*CELL_SIZE_Y*getGranularity()+CELL_SIZE_Y*getGranularity()/2;
    else return get_y_cell()*sign*CELL_SIZE_Y*getGranularity()-CELL_SIZE_Y*getGranularity()/2;
  };

  inline unsigned getGranularity() const{
    return (cellid_ & 0x00FF0000) >> 16;
  };

  void Print(std::ostream & aOs) const;

private:

  double energy_;
  unsigned adcCounts_;
  unsigned layer_;
  unsigned cellid_;
  double noiseFrac_;

  ClassDef(HGCSSRecoHit,1);

};


typedef std::vector<HGCSSRecoHit> HGCSSRecoHitVec;



#endif
