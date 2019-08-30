pragma solidity ^0.5.0;

import "./EllipticCurve.sol";
import "./FastEcMul.sol";


/**
 * @title Secp256k1 Elliptic Curve
 * @dev Secp256k1 Elliptic Curve supporting point derivation function.
 * @author Witnet Foundation
 */
contract Secp256k1 is EllipticCurve {

  // Generator coordinate `x` of EC equation
  uint256 constant GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798;
  // Generator coordinate `y` of EC equation
  uint256 constant GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8;
  // Constant `a` of EC equation
  uint256 constant AA = 0;
  // Constant `B` of EC equation
  uint256 constant BB = 7;
  // Prime number of the curve
  uint256 constant PP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
  // Order of the curve
  uint256 constant NN = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
  // Square root of the characteristic polynomial of the endomorphism of the curve
  uint256 constant LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72;
  // Beta parameter of the endomorphism
  uint256 constant BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE;

  /// @dev Public key derivation from private key.
  /// @param _d The scalar
  /// @param _x The coordinate x
  /// @param _y The coordinate y
  /// @return (qx, qy) The derived point
  function derivePoint(uint256 _d, uint256 _x, uint256 _y) public pure returns(uint256 qx, uint256 qy) {
    (qx, qy) = EllipticCurve.ecMul(
      _d,
      _x,
      _y,
      AA,
      PP
    );
  }

  /// @dev Function to derive the `y` coordinate given the `x` coordinate and the parity byte (`0x03` for odd `y` and `0x04` for even `y`).
  /// @param _yByte The parity byte following the ec point compressed format
  /// @param _x The coordinate `x` of the point
  /// @return The coordinate `y` of the point
  function deriveY(uint8 _yByte, uint256 _x) public pure returns (uint256) {
    require(_yByte == 0x02 || _yByte == 0x03, "Invalid compressed EC point prefix");

    return EllipticCurve.deriveY(
      _yByte,
      _x,
      AA,
      BB,
      PP);
  }

  /// @dev Decomposition of the scalar k in two scalars k1 and k2 with half bit-length, such that k=k1+k2*LAMBDA (mod n)
  /// @param _k the scalar to be decompose
  /// @return k1 and k2  such that k=k1+k2*LAMBDA (mod n)
  function decomposeScalar (uint256 _k) public pure returns (int256, int256) {
    return FastEcMul.decomposeScalar(_k, LAMBDA, PP);
  }

  /// @dev Substracts two key derivation functionsas `s1*A - s2*B`.
  /// @param _scalar1 The scalar `s1`
  /// @param _a1 The `x` coordinate of point `A`
  /// @param _a2 The `y` coordinate of point `A`
  /// @param _scalar2 The scalar `s2`
  /// @param _b1 The `x` coordinate of point `B`
  /// @param _b2 The `y` coordinate of point `B`
  /// @return The derived point in affine cooridnates
  function ecMulSubMul(
    uint256 _scalar1,
    uint256 _a1,
    uint256 _a2,
    uint256 _scalar2,
    uint256 _b1,
    uint256 _b2)
  internal pure returns (uint256, uint256)
  {
    (uint256 m1, uint256 m2) = derivePoint(_scalar1, _a1, _a2);
    (uint256 n1, uint256 n2) = derivePoint(_scalar2, _b1, _b2);
    (uint256 r1, uint256 r2) = EllipticCurve.ecSub(
      m1,
      m2,
      n1,
      n2,
      AA,
      PP);

    return (r1, r2);
  }

  /// @dev Simultaneous multiplication of the form kP - lQ.
  /// Scalars k and l are expected to be decomposed such that k = k1 + k2 λ, and l = l1 + l2 λ,
  /// @param _scalars An array with the decomposition of k and l values, i.e., [k1, k2, l1, l2]
  /// @param _points An array with the affine coordinates of both P and Q, i.e., [P1, P2, Q1, Q2]
  /// @return (qx, qy) = P1-P2 in affine coordinates
  function fastEcMulSubMul(
    int256[4] memory _scalars,
    uint256[4] memory _points)
  internal pure returns(uint256, uint256)
  {
    uint256 minusY1;
    uint256 minusY2;
    (minusY1, minusY2) = EllipticCurve.ecInv(_points[2], _points[3], NN);
    int256[4] memory points = [
      _points[0],
      _points[1],
      minusY1,
      minusY2];

    return FastEcMul.ecSimMul(
      _scalars,
      points,
      AA,
      BETA,
      PP);
  }
}