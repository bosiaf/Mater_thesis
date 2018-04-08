//
// Created by iotn_ on 04.04.2018.
//

#ifndef HIVEVO_TRANSITIONMATRIX_H
#define HIVEVO_TRANSITIONMATRIX_H

#include <array>

namespace epi {
typedef std::array<const double, 4> transitionRowType;
typedef std::array<const transitionRowType, 4> transitionMatrixType;
struct TransitionMatrix {
  //TRANSITION MATRIX
  const transitionMatrixType transitionMatrix =
    {{
       {{0.545, 0.2475, 0.1137, 0.0937}},
       {{0.3836, 0.4357, 0.0891, 0.0915}},
       {{0.2192, 0.1108, 0.4274, 0.2425}},
       {{0.1823, 0.1149, 0.2448, 0.4581}}
     }}; //hardcoded transition matrix

  //Hardcoded cumulative sums
  // precalculated sums of the row of the
  // transition matrix.
  const transitionRowType cs_a =
    {
      {transitionMatrix[0][0], transitionMatrix[0][0] + transitionMatrix[0][1],
       transitionMatrix[0][0] + transitionMatrix[0][1] + transitionMatrix[0][2], 1.}
    }; // A nucleotide cumsum

  const transitionRowType cs_g =
    {
      {transitionMatrix[1][0], transitionMatrix[1][0] + transitionMatrix[1][1],
       transitionMatrix[1][0] + transitionMatrix[1][1] + transitionMatrix[1][2], 1.}
    }; // G nucleotide cumsum

  const transitionRowType cs_c =
    {
      {transitionMatrix[2][0], transitionMatrix[2][0] + transitionMatrix[2][1],
       transitionMatrix[2][0] + transitionMatrix[2][1] + transitionMatrix[2][2], 1.}
    }; // C nucleotide cumsum

  const transitionRowType cs_t =
    {
      {transitionMatrix[3][0], transitionMatrix[3][0] + transitionMatrix[3][1],
       transitionMatrix[3][0] + transitionMatrix[3][1] + transitionMatrix[3][2], 1.}
    }; // T nucleotide cumsum

  //Matrix containing the 4 cumulative sums
  const transitionMatrixType matrixOfCumulativeSums =
    {
      {cs_a, cs_g, cs_c, cs_t}
    };
  //END TRANSITION MATRIX

  const TransitionMatrix getTransitionMatrix() const {
    return *this;
  }
};
}

#endif //HIVEVO_TRANSITIONMATRIX_H
