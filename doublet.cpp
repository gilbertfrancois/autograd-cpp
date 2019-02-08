/*

File         : doublet.cpp

Author       : G.F. Duivesteijn
Organization : Delft University of Technology, 
               Faculty of Aerospace Technology
email        : info@gilbertfrancois.com

Copyright    : Copyright 2003 G.F. Duivesteijn, All rights reserved.

Description: 
  
  The doublet object is designed for evaluating algoritmic derivatives
  (AD) in forward mode. In the literature, this method is also called
  automatic or computational differentiation.  The tangent value for
  the dependent value will be automatically evaluated by replacing the
  "double" type with the "doublet" object.

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
  2003-04-28: Changed <iostream> header for compatibility with gcc 3.2.

*/

#include "doublet.h"
#include <iostream>
#include <assert.h>
#include <math.h>


// Default constructor

doublet::doublet() {
  val = 0.0;
  dot = 0.0;

// #ifdef DEBUG
//   printf("Doublet constructor called.\n");
// #endif
}


// Destructor

doublet::~doublet() {

// #ifdef DEBUG
//   printf("Doublet destructor called.\n");
// #endif
}


// Operator overloading


const doublet::doublet sin (const doublet& a) {
  doublet b;
  b.val = sin(a.val);
  b.dot = cos(a.val)*a.dot;
  return b;
}

const doublet::doublet cos (const doublet& a) {
  doublet b;
  b.val = cos(a.val);
  b.dot = -sin(a.val)*a.dot;
  return b;
}

const doublet::doublet exp (const doublet& a) {
  doublet b;
  b.val = exp(a.val);
  b.dot = exp(a.val)*a.dot;
  return b;
}

const doublet::doublet log (const doublet& a) {
  doublet b;
  b.val = log(a.val);
  b.dot = 1/a.val*a.dot;
  return b; 
}

const doublet::doublet pow (const doublet& l, double r) {
  doublet b;
  b.val = pow(l.val, r);
  b.dot = r*pow(l.val, r-1)*l.dot;
  return b;
}

const doublet::doublet sqrt (const doublet& a) {
  doublet b;
  b.val = sqrt(a.val);
  b.dot = 1/(2*sqrt(a.val))*a.dot;
  return b;
}
