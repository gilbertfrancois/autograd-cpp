/*

File         : doublet.h

Author       : G.F. Duivesteijn
Organization : Delft University of Technology, 
               Faculty of Aerospace Technology
email        : info@gilbertfrancois.com

Copyright    : Copyright 2003 G.F. Duivesteijn, All rights reserved.

Description: 
  
  The doublet object is designed for evaluating algoritmic derivatives
  (AD) in forward mode or Tangent Linear Mode (TLM). In the
  literature, this method is also called automatic or computational
  differentiation.  The tangent value for the dependent value will be
  automatically evaluated by replacing the "double" type with the
  "doublet" object.

Disclamer:

  The product is not intended for use, and you may not use or allow
  others to use the product, in connection with any application
  requiring fail-safe performance such as the operation of nuclear power
  facilities, air traffic control or navigation systems, weapons control
  systems, life support systems, or any other system whose failure could
  lead to injury, death, environmental damage or mass destruction. You
  agree that the author will have no liability of any nature, and you are
  solely responsible, for any expense, loss, injury or damage incurred
  as a result of such use of this software.

History:
  
  2003-03-23: Created.

*/
#ifndef _doublet_
#define _doublet_

#include <assert.h>

class doublet {
 public:
  double val;
  double dot;

  doublet();
  ~doublet();


  friend  const doublet sin (const doublet& a);
  friend  const doublet cos (const doublet& a);
  friend  const doublet exp (const doublet& a);
  friend  const doublet log (const doublet& a);
  friend  const doublet pow (const doublet& l, double r);
  friend  const doublet sqrt (const doublet& a);

  friend inline const doublet operator+ (const double c,   const doublet& r) {
  doublet a;
  a.val = c + r.val;
  a.dot = r.dot;
  return a;
}  
  friend inline const doublet operator+ (const doublet& l, const double c) {
  doublet a;
  a.val = l.val + c;
  a.dot = l.dot;
  return a;
}
  friend inline const doublet operator+ (const doublet& l, const doublet& r) {
  doublet a;
  a.val = l.val + r.val;
  a.dot = l.dot + r.dot;
  return a;
}
  friend inline const doublet operator- (const doublet& l) {
  doublet a;
  a.val = - l.val;
  a.dot = - l.dot;
  return a;
}
  friend inline const doublet operator- (const double c,   const doublet& r) {
  doublet a;
  a.val = c - r.val;
  a.dot = r.dot;
  return a;
}
  friend inline const doublet operator- (const doublet& l, const double c) {
  doublet a;
  a.val = l.val - c;
  a.dot = l.dot;
  return a;
}
  friend inline const doublet operator- (const doublet& l, const doublet& r) {
  doublet a;
  a.val = l.val - r.val;
  a.dot = l.dot - r.dot;
  return a;
}
  friend inline const doublet operator* (const double c,   const doublet& r) {
  doublet a;
  a.val = c*r.val;
  a.dot = c*r.dot;
  return a;
}
  friend inline const doublet operator* (const doublet& l, const double c) {
  doublet a;
  a.val = c*l.val;
  a.dot = c*l.dot;
  return a;
}
  friend inline const doublet operator* (const doublet& l, const doublet& r) {
  doublet a;
  a.val = l.val*r.val;
  a.dot = l.val*r.dot + l.dot*r.val;
  return a;
}
  friend inline const doublet operator/ (const doublet& l, const doublet& r) {
  doublet a;
  assert (r.val != 0.0);
  a.val = l.val / r.val;
  a.dot = 1/r.val*l.dot - l.val/(r.val*r.val)*r.dot;
  return a;
}
  friend inline const doublet operator/ (const doublet& l, const double& r) {
  doublet a;
  assert (r != 0.0);
  a.val = l.val / r;
  a.dot = 1/r*l.dot;
  return a;
}
  friend inline const doublet operator/ (const double& l, const doublet& r)  {
  doublet a;
  assert (r.val != 0.0);
  a.val = l / r.val; 
  a.dot = - l/(r.val*r.val)*r.dot;
  return a;
}

};
#endif
