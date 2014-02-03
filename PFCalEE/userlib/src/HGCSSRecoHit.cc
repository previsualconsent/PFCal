#include "HGCSSRecoHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSRecoHit::HGCSSRecoHit(const HGCSSSimHit & aSimHit, const unsigned granularity){
  energy_ = aSimHit.energy();
  adcCounts_ = 0;
  layer_ = aSimHit.layer();
  noiseFrac_ = 0;

  double x = aSimHit.get_x();
  double y = aSimHit.get_y();
  //cellid encoding:
  bool x_side = x>0 ? true : false;
  bool y_side = y>0 ? true : false;
  unsigned x_cell = static_cast<unsigned>(fabs(x)/(CELL_SIZE_X*granularity));
  unsigned y_cell = static_cast<unsigned>(fabs(y)/(CELL_SIZE_Y*granularity));

  encodeCellId(x_side,y_side,x_cell,y_cell,granularity);

}

void HGCSSRecoHit::encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell, const unsigned granularity){
  cellid_ = 
    x_side | (x_cell<<1) |
    (y_side<<8) | (y_cell<<9) |
    (granularity<<16) ;

  // std::cout << " Cross-check of encoding: cellid=" << cellid_ << std::endl
  // 	    << " x_side " << x_side << " " << get_x_side() << std::endl
  // 	    << " y_side " << y_side << " " << get_y_side() << std::endl
  // 	    << " x_cell " << x_cell << " " << get_x_cell() << std::endl
  // 	    << " y_cell " << y_cell << " " << get_y_cell() << std::endl
  //        << " granularity " << granularity << " " << getGranularity() << std::endl
  //   ;
}

void HGCSSRecoHit::Add(const HGCSSSimHit & aSimHit){
  energy_ += aSimHit.energy();
}

void HGCSSRecoHit::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Layer " << layer_ << " cellid " << cellid_ << std::endl
      << " = Energy " << energy_ << " noiseFrac " << noiseFrac_ << std::endl
      << " = Digi E " << adcCounts_ << " adcCounts." << std::endl
      << "====================================" << std::endl;

}
